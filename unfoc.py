#!/usr/bin/env python3
#
# unfoc.py: Unfocused Radar Processor
# Output file for pik1 (4-byte signed integer, network order)
#
# All output files are 4-byte network-order (big endian).
#
# see parse_channels.py for instructions on how to set channels
# TODO: make a graphviz diagram showing the pipeline

import argparse
from collections import namedtuple
import cProfile
import logging
import os
import gzip
import struct
import sys

import numpy as np
from scipy import signal

try:
    import typing
    from typing import Any, BinaryIO, Dict, Generator, List, Optional, Set, Tuple
except ImportError: #pragma: no cover
    pass  # Not installed on melt ...

from parse_channels import parse_channels, get_utig_channels, PIK1ChannelSpec

################################################
# Enable metadata index writing.
# This should normally be disabled only for legacy testing.
# It does not make things slower to output this data.
enable_meta_index=True
################################################

from radbxds import Trace, get_radar_stream, gen_ct


class IncoherentTrace:
    """
    Pairs channel and magnitude/phase data for an incoherent trace
    """
    def __init__(self, channel, magnitude, phase, ct):
        # type: (int, np.ndarray, np.ndarray, List[tuple]) -> None
        self.channel = channel
        self.magnitude = magnitude
        self.phase = phase
        self.ct = ct # type: tuple


class IncoherentStack:
    """
    Properly handles stacking incoherent magnitude and phase data.

    This needs to be separate from the utility for stacking complex traces
    because summing the phases would make no sense - instead, it just
    returns the phase of the central trace.
    """
    def __init__(self, depth, do_phase=False):
        # type: (int, bool) -> None
        self.depth = depth
        self.mag_stack = SingleStack(depth, False)
        self.phs_stack = SingleStack(depth, False, dtype=np.cdouble)
        self.do_phase = do_phase
        self.stack_center = self.depth // 2

    def add_trace(self, dechirped):
        # type: (Trace, int, List[tuple]) -> IncoherentTrace
        """
        Input trace is a complex, dechirped trace.
        Output is dechirped, broken into mag and phs
        """
        mag = self.mag_stack.add_trace(np.abs(dechirped.data), dechirped.ct)
        if self.do_phase:
            # Calculating angle later uses a bit more memory but decreases CPU usage.
            phs = self.phs_stack.add_trace(dechirped.data, dechirped.ct)

        if mag is not None:
            # Grab the seq and tim of the first raw trace
            # that was used to make this incoherent trace
            mmag, mct = mag
            mmag = np.mean(mmag, axis=0)

            if self.do_phase:
                pphs, pcts = phs
                pphs = np.angle(pphs[self.stack_center,...])
            else:
                pphs = None
            # Get seq and tim from the middle of the stack
            #seq, tim = mcts[self.stack_center]

            return IncoherentTrace(dechirped.channel, mmag, pphs, mct)

def cinterp(sweep_fft, index):
    # type: (np.ndarray, int) -> np.ndarray
    # sweep_fft is fft of sweep, and index is a bin affected by the LO noise.
    r = (np.abs(sweep_fft[index-1]) + np.abs(sweep_fft[index+1])) / 2
    t1 = np.angle(sweep_fft[index-1])
    t2 = np.angle(sweep_fft[index+1])
    if (np.abs(t1 - t2) > np.pi):
        t1 = t1 + 2 * np.pi
    theta = (t1 + t2) / 2
    sweep_fft[index] = r * (np.cos(theta) + 1j * np.sin(theta))
    return sweep_fft


# QUESTION: with numpy, is this modifying the input trace, or just the output?
def denoise_and_dechirp(trace, # type: np.ndarray
                        ref_chirp, # type: np.ndarray
                        blanking, # type: int
                        output_samples, # type: int
                        do_cinterp, # type: bool
                        detrend_type='linear'
                       ):
    # type: (...) -> np.ndarray

    # Input trace is a 1 x output_samples array.
    if (blanking >= 0):
        trace[:blanking] = np.zeros(blanking)
    else:
        trace[-blanking:] = np.zeros(len(trace)+blanking)

    #find peak energy below blanking samples
    ## [n,m]=sort(trace);
    ## shifter=abs((m(output_samples)));
    # GNG -- this line wasn't doing what you think it was doing.
    #shifter = int(np.median(np.argmax(trace)));
    # You have to do this:
    #shifter = np.median(np.argwhere(listy == np.amax(listy)))
    shifter = int(np.argmax(trace))
    trace = np.roll(trace, -shifter);

    DFT = np.fft.fft(signal.detrend(trace, type=detrend_type))

    if do_cinterp:
        # Remove five samples per cycle problem
        n1 = int(np.round(output_samples * (1.0/5)))
        n2 = output_samples - n1
        DFT = cinterp(DFT, n1)
        DFT = cinterp(DFT, n2)
        # Remove the first harmonic for five samples
        n1 = int(np.round(output_samples * (2.0/5)))
        n2 = output_samples - n1
        DFT = cinterp(DFT, n1)
        DFT = cinterp(DFT, n2)

    # Do the dechirp
    Product = np.multiply(ref_chirp, DFT)
    Dechirped = np.fft.ifft(Product)
    Dechirped = np.roll(Dechirped, shifter)
    return Dechirped


# TODO: LEL I don't like that this is a wrapper just to trim input data
# to the right length and then reattach the channel
# to the output of denoise_and_dechirp
def denoise_and_dechirp_gen(cohstacks, # type: Generator[Trace, None, None]
                            ref_chirp, # type: np.ndarray
                            blanking, # type: int
                            output_samples, # type: int
                            do_cinterp # type: bool
                           ):
    # type: (...) -> Generator[Trace, None, None]
    for trace in cohstacks:
        dechirped = denoise_and_dechirp(trace.data[0:output_samples], ref_chirp,
                                        blanking, output_samples, do_cinterp)

        yield Trace(trace.channel, dechirped, trace.ct)


class SingleStack:
    """
    Generates a stack for a single channel of complex (coherent) traces.
    """
    def __init__(self, depth, presum=True, dtype=np.float64):
        # type: (int, bool) -> None
        self.depth = depth # How many sweeps to add into a stack
        self.idx = 0 # Where to insert new sweep in the array
        self.count = 0 # Total number of sweeps processed
        self.presum = presum # Whether to sum the returned stack or not
        # these are autodetected from first call to add_trace
        self.samples = None # type: Optional[int]
        self.stacks = None # type: Optional[np.ndarray]
        # ct tim and sequence numbers of each stack processed
        self.cts = [] # type: List[tuple]
        self.dtype = dtype
        self.stack_center = self.depth // 2

    def add_trace(self, sweep, ct):
        # type: (np.ndarray, int, tuple) -> Optional([np.ndarray, tuple])
        # output will either be depth x samples or 1 x samples, depending
        # on presum
        if self.stacks is None:
            self.samples = sweep.size # How many samples long each sweep is
            self.stacks = np.zeros((self.depth, self.samples), dtype=self.dtype)

        assert sweep.size == self.samples

        self.stacks[self.idx, :] = sweep
        self.cts.append(ct)
        self.idx += 1
        self.count += 1

        if self.idx >= self.depth:
            self.idx = 0
            # This is coherent summing of a single trace
            # Select the center of the coherent stack as the ct
            logging.debug("Selecting CT {:d} of [0,{:d}, depth={:d})".format(self.stack_center, len(self.cts), self.depth))
            logging.debug(str(self.cts))
            ct = self.cts[self.stack_center]
            self.cts = []
            if self.presum:
                return (self.stacks.sum(axis=0), ct)
            else:
                return (self.stacks, ct)
        else:
            return None



class PairedStack:
    """
    Accepts multiple completed stacks and returns a stack if both
    input channel queues contain data
    """
    def __init__(self, channel_spec, depth):
        # type: (PIK1ChannelSpec, int) -> None
        self.channel_spec = channel_spec
        self.depth  = depth

        self.fullstacks0 = [] # type: List[np.ndarray]
        self.fullstacks1 = [] # type: List[np.ndarray]
        self.ct0 = []
        self.ct1 = []

        # Having the same input channels in the channel spec will cause
        # the stack queues to get confused.
        assert self.channel_spec.chan0in != self.channel_spec.chan1in

    def get_stack(self):
        # type: () -> Trace
        # Logic lifted from read_and_stack.cc
        bDo0 = (self.channel_spec.chan0in > 0)
        bDo1 = (self.channel_spec.chan1in > 0)
        bReady0 = bDo0 and len(self.fullstacks0) > 0
        bReady1 = bDo1 and len(self.fullstacks1) > 0
        scale0 = self.channel_spec.scalef0 / float(self.depth)
        scale1 = self.channel_spec.scalef1 / float(self.depth)
        array_out1 = None

        if not ((bDo0 and not bReady0) or (bDo1 and not bReady1)) :
            # assert self.fullstacks0[-1].shape[0] == self.depth?

            if bDo0 and bDo1:
                array_out1 = scale0 * self.fullstacks0[-1] + scale1 * self.fullstacks1[-1]
                ct = self.ct0[-1]
                #seq1, tim1 = self.ct1[-1]
                # Both stacks should be contemporaneous, but they will have
                # slightly different trace numbers.  So we choose the sequence
                # number of fullstacks0
                # we could try raising assertions that they are close to each other
            elif bDo0 and not bDo1:
                # Inconsistent with SingleStack, we choose the last ct.
                ct = self.ct0[-1]
                array_out1 = scale0 * self.fullstacks0[-1]
            elif not bDo0 and bDo1:
                # Inconsistent with SingleStack, we choose the last ct.
                ct = self.ct1[-1]
                array_out1 = scale1 * self.fullstacks1[-1]

            if bReady0:
                self.fullstacks0.pop()
                self.ct0.pop()
            if bReady1:
                self.fullstacks1.pop()
                self.ct1.pop()
            return Trace(self.channel_spec.chanout, array_out1, ct)


    def add_stack(self, trace):
        # type: (Trace) -> Optional[Trace]
        # Adds an existing stack to the state and returns the resulting
        # output stack if adding it enabled a new one to be computed.
        # Array dimensions are #channels.
        if len(self.fullstacks0) >= 1000 or len(self.fullstacks1) >= 1000:
            msg = ("overflow: p1cs={0:s} fullstacks0={0:d} fullstacks1={1:d}".format(
                       self.channel_spec, len(self.fullstacks0), len(self.fullstacks1)))
            raise RuntimeError(msg)
        if self.channel_spec.chan0in == trace.channel:
            self.fullstacks0.insert(0, trace.data)
            self.ct0.insert(0, trace.ct)
        elif self.channel_spec.chan1in == trace.channel:
            self.fullstacks1.insert(0, trace.data)
            self.ct1.insert(0, trace.ct)
        else:
            # no change in state, so no need for get_stack
            return None

        return self.get_stack()


# Yields a tuple of an output channel and a NDArray of values
# coherently stacked
# tracegen should be a sequence (generator) of tuples where traces[0]
# is the channel id, and traces[1] is an ndarray that is the trace
# TODO: sweep length is changing with RADnh5
def stacks_gen(traces, # type: List[Trace]
               channel_specs, # type: List[PIK1ChannelSpec]
               depth, # type: int
               ):
    # type: (...) -> Generator[Trace, None, None]
    stacks = [] # type: List[PairedStack]
    ss0 = {} # type: Dict[int, SingleStack]
    for (idx, p1cs) in enumerate(channel_specs):
        logging.debug("p1cs[%d]=%s"% (idx, str(p1cs)))
        stacks.append(PairedStack(p1cs, depth))
        for chan in (p1cs.chan0in, p1cs.chan1in):
            if chan > 0 and chan not in ss0:
                ss0[chan] = SingleStack(depth)

    ## TODO: make this read function selectable to accommodate for RADnh4
    for trace in traces:
        # Stack into appropriate bins
        # Iterate over all the channel specs and see which ones need to be stacked
        if trace.channel in ss0:
            cohstack = ss0[trace.channel].add_trace(np.float64(trace.data), trace.ct)

            if cohstack is not None:
                for (idx, p1cs) in enumerate(channel_specs):
                    # Test to prevent extraneous calls to add_stack()
                    if trace.channel in [p1cs.chan0in, p1cs.chan1in]:
                        stack, cts = cohstack
                        new_trace = Trace(trace.channel, stack, cts)
                        result = stacks[idx].add_stack(new_trace)
                        # If this stack state yielded a full stack, yield it
                        if result is not None:
                            yield result


def inco_stacks_gen(traces, # type: Generator[Trace, None, None]
                   channel_specs, # type: List[PIK1ChannelSpec]
                   depth, # type: int
                   output_samples=3200, # type: int
                   do_phase=False # type: bool
                  ):
    """

    """

    # type: (...) -> Generator[IncoherentTrace, None, None]
    ss0 = {} # type: Dict[int, IncoherentStack]
    for p1cs in channel_specs:
        chan = p1cs.chanout
        if chan > 0 and chan not in ss0:
            ss0[chan] = IncoherentStack(depth, do_phase)

    ## TODO: make this read function selectable to accommodate for RADnh4
    for trace in traces:
        # Stack into appropriate bins
        # Iterate over all the channel specs and see which ones need to be stacked
        if chan in ss0:
            stack = ss0[trace.channel].add_trace(trace)
            if stack is not None:
                yield stack


def get_radar_type(bxdsfile, nrecords=1000):
    """ Inspect a raw datafile and detect what type of radar it came from
    returns 'HiCARS2' if it is a 1-antenna radar, or 'MARFA' if it is a 2-antenna
    radar.
    """
    # Use a nominal channel spec that makes it return all records
    channel_specs = [PIK1ChannelSpec(chanout=1, chan0in=1, scalef0=1, chan1in=3, scalef1=1), # sum left and right low gain
                    PIK1ChannelSpec(chanout=2, chan0in=2, scalef0=1, chan1in=4, scalef1=1), # sum left and right high gain
                    ]

    choffs = {}
    gen = read_RADnhx_gen(bxdsfile, channel_specs)
    for ii, trace in enumerate(gen):
        choffs[trace.channel] = choffs.get(trace.channel, 0) + 1
        if ii >= nrecords or len(choffs) >= 3:
            break
    gen.close()

    if len(choffs) >= 3:
        return 'MARFA'
    else:
        assert len(choffs) == 2
        return 'HiCARS2'


# Read individual traces out of RADnh3 or RADnh5 file
# TODO: This doesn't yet filter on channels, which should be OK - the
# downstream users ALSO check channel.
def read_RADnhx_gen(input_filename, channel_specs):
    # type: (str, List[PIK1ChannelSpec]) -> Generator[Trace]


    stream = get_radar_stream(input_filename)
    if stream == "RADnh3":
        header_t = namedtuple('radnh3_header',
                              'nsamp nchan vr0 vr1 choff ver resvd2')
        fmtstr = '>HBBBBBB'
    elif stream == "RADnh5":
        header_t = namedtuple('radnh5_header',
                              'nsamp nchan vr0 vr1 choff ver resvd2 absix relix xinc rseq scount tscount')
        fmtstr = '>HBBBBBBddfLHL'
    else:
        raise ValueError("Invalid stream type for file: %s" % input_filename)

    fmtlen = struct.calcsize(fmtstr)
    # In theory, it's faster to do this here and only compile the format string once.
    header_struct = struct.Struct(fmtstr)

    # This way we don't waste time converting channel data that won't be used
    choffs = set() # type: Set[int]
    for p1cs in channel_specs:
        choffs.add(((p1cs.chan0in-1)//2)*2)
        choffs.add(((p1cs.chan1in-1)//2)*2)

    genct = gen_ct(input_filename)
    if not genct:
        logging.warning("Can't read corresponding ct file for " + input_filename)

    with open(input_filename, 'rb') as fd:
        while True:
            # read header
            buff = fd.read(fmtlen)
            if len(buff) < fmtlen:
                return
            ctdata = next(genct) if genct else (None, None)
            header = header_t._make(header_struct.unpack(buff))
            if stream == "RADnh5":
                # We don't do anything with the timestamps yet, but we need
                # to skip over 'em to get to the radar data.
                if header.tscount > 0:
                    # timestamp count should be on the order of the number of stacks.
                    assert(header.tscount < 100)
                    #tstamps = struct.unpack_from('>{:d}d'.format(tscount), fd.read(8*tscount))
                    fd.seek(8*header.tscount, os.SEEK_CUR)

            # both headers have choff
            # Legacy value of 0xff is equivalent to offset 0
            if header.choff == 0xff:
                choff = 0
            else:
                choff = 0x0f & header.choff

            # reclen = SweepLength
            input_samples = header.nsamp
            assert 0 < input_samples < 10000

            if choff in choffs:
                traces = np.fromfile(fd, dtype='>i2', count=input_samples*2)
                try:
                    traces.shape = (2, input_samples)
                except ValueError:
                    break # didn't get enough data.

                # no use yet for headers....
                yield Trace(choff+1, traces[0][...], ctdata) #, header, radnhx_header)
                yield Trace(choff+2, traces[1][...], ctdata) #, header, radnhx_header)
            else:
                # Skip the rest of this record without yielding any values
                fd.seek(4*input_samples, os.SEEK_CUR)






class PIK1OutputFile(object):
    """
    Outputs magnitude and phase information for NDArray that comes through.
    Optionally computes the mean before writing
    It can be used as a tee filter and tap output at any phase in the
    processing stream.
    """
    def __init__(self, outdir, channel, MagScale, StackDepth, IncoDepth):
        # type: (str, int, int, int, int) -> None
        self.MagFileName = "{0:s}/MagLoResInco{1:d}".format(outdir, channel)
        self.PhsFileName = "{0:s}/PhsLoResInco{1:d}".format(outdir, channel)
        self.MetaFileName = "{0:s}/MagLoResInco{1:d}.meta".format(outdir, channel)
        self.TracesFileName = "{0:s}/TraceNumbers{1:d}".format(outdir, channel)
        # File descriptors
        self.MagOutFD = None # type: Optional[BinaryIO]
        self.PhsOutFD = None # type: Optional[BinaryIO]
        self.MetaOutFD = None # type: Optional[BinaryIO]
        self.TraceNumbersFD = None # type: Optional[BinaryIO]

        self.ChannelNum = channel
        self.MagScale = MagScale
        self.StackDepth = StackDepth
        self.IncoDepth = IncoDepth

        self.record_increment = self.StackDepth*self.IncoDepth
        self.record_idx = self.record_increment/2

        # Cache result for whether we need to do byte swapping
        self.do_byteswap = sys.byteorder == 'little'

    def __del__(self):
        self.close()

    def open(self, input_filename, do_phase=True, do_index=True):
        # type: (str, bool, bool) -> None
        '''
        Opens all file descriptors and outputs the header for the metafile.
        Having args input is ugly, but it contains vars not used anywhere
        else in the class other than that header file.
        '''
        self.close()

        # TODO: remove references to args as possible
        self.enable_meta_idx = do_index

        ## Open Input and Output files
        self.MetaOutFD = open(self.MetaFileName, 'wt')
        self.MetaOutFD.write('#InputName = "' + input_filename + '"\n')
        logging.info("writing %s" % self.MagFileName)
        self.MagOutFD = open(self.MagFileName, 'wb')
        self.MetaOutFD.write('#MagName = "' + self.MagFileName + '"\n')

        if do_phase:
            self.PhsOutFD = open(self.PhsFileName, 'wt')
            self.MetaOutFD.write('#PhsName = "' + self.PhsFileName + '"\n')

        self.TraceNumbersFD = open(self.TracesFileName, 'wt')

        ###### FIXME? - request hdr file from xlob and read/forward
        # TODO: collect this information into a data structure rather than reading directly from 
        ###### Make this Meta similar
        self.MetaOutFD.write("#ChannelNum = " + str(self.ChannelNum) + "\n")
        self.MetaOutFD.write("#StackDepth = " + str(self.StackDepth) + "\n")
        self.MetaOutFD.write("#IncoDepth = "  + str(self.IncoDepth) + "\n")
        self.MetaOutFD.write("#Scale = "      + str(self.MagScale) + "\n")
        self.MetaOutFD.write("#Log = TRUE\n")

    def write_record(self, inco_trace):
        # type: (IncoherentTrace) -> None
        # Write component files if enabled
        if self.PhsOutFD is not None:
            Phase = np.int32(inco_trace.phase * 16777216)
            if self.do_byteswap:
                Phase.byteswap(True)
            Phase.tofile(self.PhsOutFD)

        if self.MagOutFD is not None:
            ScaledMag = np.int32(self.MagScale * np.log10(inco_trace.magnitude))
            if self.do_byteswap:
                ScaledMag.byteswap(True)
            ScaledMag.tofile(self.MagOutFD)

        if self.TraceNumbersFD is not None:
            self.TraceNumbersFD.write(str(inco_trace.ct.seq) + "\n")

        if self.enable_meta_idx:
            self.MetaOutFD.write(str(self.record_idx) + "\n")
        self.record_idx += self.record_increment

    def close(self):
        # type: () -> None
        # flush stacks,
        # TODO: warn about incomplete stacks
        # close file handles
        if self.MagOutFD is not None:
            self.MagOutFD.close()
            self.MagOutFD = None

        if self.PhsOutFD is not None:
            self.PhsOutFD.close()
            self.PhsOutFD = None

        if self.MetaOutFD is not None:
            self.MetaOutFD.close()
            self.MetaOutFD = None

        if self.TraceNumbersFD is not None:
            self.TraceNumbersFD.close()
            self.TraceNumbersFD = None


def get_ref_chirp(bandpass, trunc_sweep_length):
    # type: (bool, bool) -> np.ndarray
    I = np.array([-63, -92, -109, -75, -87, -50, -116, -154, -22, 68, 141,
                  -610, 1461, 3807, -6147, -5375, 10651, -4412, -9810, 15386,
                  -3070, -14499, 15130, 3677, -15935, 3743, 13362, -5884,
                  -13301, 8455, 12542, -8744, -11977, 5105, 13754, -961,
                  -14342, -5184, 10294, 12194, -2709, -14352, -8807, 5965,
                  15350, 8368, -6605, -14990, -11515, 196, 11276, 15490,
                  10300, -645, -10730, -15307, -13379, -7342, 377, 7264,
                  11662, 13435, 11530, 3243, -4865, -4427, -3233, -4000,
                  -2472, -2498, -2361, -1230, -1311, -618, -578, -569, -121,
                  -319, 206, 328, 436, 613, 318, 514, 353, 277, 221, 34, 250,
                  132, 199, 189, 75, 190, 65, 106, 19, -64, -14, -117])
    if not bandpass:
        rchirp = np.flipud(I)
    else:
        rchirp = I

    return np.fft.fft(rchirp, n=trunc_sweep_length)


def get_hfilter(trunc_sweep_length):
    # type: (int) -> np.ndarray
    # Convert MHz to samples
    min_freq = int(round(2.5 * trunc_sweep_length / 50))
    max_freq = int(round(17.5 * trunc_sweep_length / 50))
    dfreq = max_freq - min_freq + 1
    hamming = np.sin(np.linspace(0, 1, num=dfreq) * np.pi)
    hfilter = np.hstack((np.zeros(min_freq),
                         hamming*2,
                         np.zeros(trunc_sweep_length - 2*max_freq - 1),
                         np.zeros(hamming.size),
                         np.zeros(min_freq - 1)))
    return hfilter

def main():
    # type: () -> None


    parser = argparse.ArgumentParser(description='Pulse compress radar data')

    parser.add_argument('--outdir', required=True,
                        help='directory for output files')
    parser.add_argument('--infile', required=True,
                        help='filename of 2-byte radar file')

    cgroup = parser.add_mutually_exclusive_group(required=True)
    cgroup.add_argument('--channel_def',
                        help="Channel spec string. This option is deprecated. Use --channels")
    cgroup.add_argument('--channels',
                        help="comma-separated channels to produce, such as 'LoResInco1,LoResInco2'")

    parser.add_argument('--output_samples', required=True, type=int,
                        help='Length of each output sweep (in samples)')
    parser.add_argument('--StackDepth', required=True, type=int,
                        help='coherent stacking depth')
    parser.add_argument('--IncoDepth', required=True, type=int,
                        help='incoherent stacking depth')

    parser.add_argument('--Scale', type=int, default=20000,
                        help='Output scale default is 1000*dB')
    parser.add_argument('--blanking', type=int, default=50,
                        help='Samples at the top of the record to blank out. Negative number blanks bottom')

    parser.add_argument('--output_phases', action='store_true',
                        help="output phase, in addition to magnitude")
    parser.add_argument('--bandpass', action='store_true',
                        help='bandpass sampling, false is for legacy hicars. Disables cinterp and flips the chirp')
    parser.add_argument('--nmax', default=0, type=int,
                        help="Maximum number of stacks to output (usually used for testing")
    parser.add_argument('--debug', action='store_true',
                        help='Print debugging messages')

    args = parser.parse_args()


    LOGLEVEL = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=LOGLEVEL, stream=sys.stdout,
                    format='unfoc: [%(levelname)-5s] %(message)s',
                   )

    if args.channel_def: # legacy channel definitions
        channel_specs = parse_channels(args.channel_def)
        logging.debug("Channel spec: %r" % channel_specs)
    else:
        channel_specs = args.channels


    return unfoc(outdir=args.outdir,  infile=args.infile,
                 channels=channel_specs, output_samples=args.output_samples,
                stackdepth=args.StackDepth, incodepth=args.IncoDepth,
                blanking=args.blanking, bandpass=args.bandpass,
                scale=args.Scale, output_phases=args.output_phases, nmax=args.nmax)



def unfoc(outdir, infile, channels, output_samples, stackdepth, incodepth,
          blanking, bandpass, scale=20000, output_phases=False, nmax=0):


    # pass through if this is a legacy ChannelSpec object
    if type(channels) == 'str':
        radartype = get_radar_type(infile)
        logging.info("Radar type: " + radartype)
        channel_specs = get_utig_channels(channels, radar=radartype)
    else:
        channel_specs = channels

    # Obtain reference chirp
    ref_chirp = get_ref_chirp(bandpass, output_samples)

    # Compute Hamming Filter
    hfilter = get_hfilter(output_samples)

    ## Filter Chirp
    ref_chirp = np.multiply(ref_chirp, hfilter)


    """
    This shows how to construct a filter pipeline using generators.
    tracegen takes as input a bxds file (and ct file, implicitly), and generates raw Traces().

    """

    # Read traces from file
    tracegen = read_RADnhx_gen(infile, channel_specs)
    # Demultiplex stacks and generate coherently-stacked traces
    stackgen = stacks_gen(tracegen, channel_specs, stackdepth)
    # Dechirp coherent stacks
    dechirpgen = denoise_and_dechirp_gen(stackgen, ref_chirp, blanking,
                                         output_samples,
                                         not bandpass)
    # Incoherently stack
    istackgen = inco_stacks_gen(dechirpgen, channel_specs, incodepth,
                                output_samples, do_phase=output_phases)

    # Initialize output files
    os.makedirs(outdir, exist_ok=True)

    outfiles = {}
    for p1cs in channel_specs:
        outfiles[p1cs.chanout] = PIK1OutputFile(outdir,
                                                p1cs.chanout,
                                                scale,
                                                stackdepth,
                                                incodepth)
        outfiles[p1cs.chanout].open(infile, do_phase=output_phases)

    for ii, rec in enumerate(istackgen):
        if rec.channel in outfiles:
            outfiles[rec.channel].write_record(rec)
        else:
            raise ValueError("Invalid output channel %d" % rec.channel)
        if nmax > 0 and ii >= nmax:
            break

    for p1cs in channel_specs:
        outfiles[p1cs.chanout].close()



if __name__ == "__main__":
    sys.exit(main())
