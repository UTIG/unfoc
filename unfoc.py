#!/usr/bin/env python
#
# unfoc.py: Unfocused Radar Processor
# Output file for pik1 (4-byte signed integer, network order)
#
# All output files are 4-byte network-order.
#

import argparse
from collections import namedtuple
import cProfile
import logging
import numpy as np
import os
import Queue
from scipy import signal
import struct
import sys
try:
    import typing
    from typing import Any, BinaryIO, Dict, Generator, List, Optional, Set, Tuple
except:
    # Not installed on melt ...
    pass

from parse_channels import parse_channels, PIK1ChannelSpec

################################################
# Enable metadata index writing.
# This should normally be disabled only for legacy testing.
# It does not make things slower to output this data.
enable_meta_index=True
################################################

class Trace:
    """
    Pairs channel + data for a trace.
    Data is either raw amplitudes, or complex amplitudes
    """
    def __init__(self, channel, data):
        # type: (int, np.ndarray) -> None
        self.channel = channel
        self.data = data

class IncoherentTrace:
    """
    Pairs channel and magnitude/phase data for an incoherent trace
    """
    def __init__(self, channel, magnitude, phase):
        # type: (int, np.ndarray, Optional[np.ndarray]) -> None
        self.channel = channel
        self.magnitude = magnitude
        self.phase = phase


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
        self.phs_stack = SingleStack(depth, False)
        self.do_phase = do_phase
        self.StackCenter = int(np.fix((self.depth/2+0.5)))

    def add_trace(self, dechirped):
        # type: (Trace) -> IncoherentTrace
        """
        Input trace is a complex, dechirped trace.
        Output is dechirped, broken into mag and phs
        """
        mag = self.mag_stack.add_trace(np.abs(dechirped.data))
        if self.do_phase:
            phs = self.phs_stack.add_trace(np.angle(dechirped.data))

        if mag is not None:
            mag = np.mean(mag, axis=0)
            if self.do_phase:
                phs = phs[self.StackCenter,...]
            else:
                phs = None
            return IncoherentTrace(dechirped.channel, mag, phs)

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
                        do_cinterp # type: bool
                       ):
    # type: (...) -> np.ndarray

    # Input trace is a 1 x output_samples array.
    trace[0:blanking] = np.zeros(blanking)

    #find peak energy below blanking samples
    ## [n,m]=sort(trace);
    ## shifter=abs((m(output_samples)));
    shifter=int(np.median(np.argmax(trace)));
    trace=np.roll(trace,-shifter);

    DFT = np.fft.fft(signal.detrend(trace))

    if do_cinterp:
        # Remove five samples per cycle problem
        DFT = cinterp(DFT, output_samples * (1.0/5))
        DFT = cinterp(DFT, output_samples * (1 - 1.0/5))
        # Remove the first harmonic for five samples
        DFT = cinterp(DFT, output_samples * (2.0/5))
        DFT = cinterp(DFT, output_samples * (1 - 2.0/5))

    # Do the dechirp
    Product = np.multiply(ref_chirp, DFT)
    Dechirped = np.fft.ifft(Product)
    Dechirped = np.roll(Dechirped,shifter)
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
        yield Trace(trace.channel, dechirped)


class SingleStack:
    """
    Generates a stack for a single channel of complex (coherent) traces.
    """
    def __init__(self, depth, presum=True):
        # type: (int, bool) -> None
        self.depth = depth # How many sweeps to add into a stack
        self.idx = 0 # Where to insert new sweep in the array
        self.count = 0 # Total number of sweeps processed
        self.presum = presum # Whether to sum the returned stack or not
        # these are autodetected from first call to add_trace
        self.samples = None # type: Optional[int]
        self.stacks = None # type: Optional[np.ndarray]

    def add_trace(self, sweep):
        # type: (np.ndarray) -> Optional[np.ndarray]
        # output will either be depth x samples or 1 x samples, depending
        # on presum
        if self.stacks is None:
            self.samples = sweep.size # How many samples long each sweep is
            self.stacks = np.zeros((self.depth, self.samples), dtype=np.float64)

        assert sweep.size == self.samples

        self.stacks[self.idx,...] = sweep
        self.idx += 1
        self.count += 1

        if self.idx >= self.depth:
            self.idx = 0
            if self.presum:
                return self.stacks.sum(axis=0)
            else:
                return self.stacks
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
            elif bDo0 and not bDo1:
                array_out1 = scale0 * self.fullstacks0[-1]
            elif not bDo0 and bDo1:
                array_out1 = scale1 * self.fullstacks1[-1]

            if bReady0:
                self.fullstacks0.pop();
            if bReady1:
                self.fullstacks1.pop();
            return Trace(self.channel_spec.chanout, array_out1)


    def add_stack(self, trace):
        # type: (Trace) -> Optional[Trace]
        # Adds an existing stack to the state and returns the resulting
        # output stack if adding it enabled a new one to be computed.
        # Array dimensions are #channels.
        if len(self.fullstacks0) >= 1000 or len(self.fullstacks1) >= 1000:
            msg = ("overflow: p1cs={4:s} fullstacks0={0:d} fullstacks1={1:d}".format(
                       self.channel_spec, len(self.fullstacks0), len(self.fullstacks1)))
            raise Exception(msg)
        if self.channel_spec.chan0in == trace.channel:
            self.fullstacks0.insert(0, trace.data)
        elif self.channel_spec.chan1in == trace.channel:
            self.fullstacks1.insert(0, trace.data)
        else:
            # no change in state, so no need for get_stacks
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
            cohstack = ss0[trace.channel].add_trace(np.float64(trace.data))

            if cohstack is not None:
                for (idx, p1cs) in enumerate(channel_specs):
                    # Test to prevent extraneous calls to add_stack()
                    if trace.channel in [p1cs.chan0in, p1cs.chan1in]:
                        new_trace = Trace(trace.channel, cohstack)
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


def get_radar_stream(filename):
    # type: (str) -> str
    # both RADnh3 and RADnh5 are the same first elements in header format
    fmtstr='>HBBBBBB'
    fmtlen=struct.calcsize(fmtstr)
    with open(filename, 'rb') as fd:
        # read header
        buff = fd.read(fmtlen)
        if len(buff) < fmtlen:
            raise Exception("Unable to read %d bytes from %s."
                            % (fmtlen, filename))

        # TODO: This check should only  happen once, not every read!!
        # Check it then set a flag and move the file pointer back to the start!
        # Check the version byte and see if it is RADnh3 or RADnh5
        v = ord(buff[6])
        if v == 5:
            return "RADnh5"
        elif v == 0 or v == 255:
            return "RADnh3"
        else:
            raise Exception("Can't determine stream for file %s" % filename)

# Read individual traces out of RADnh3 or RADnh5 file
# TODO: This doesn't yet filter on channels, which should be OK - the
# downstream users ALSO check channel.
def read_RADnhx_gen(input_filename, channel_specs):
    # type: (str, List[PIK1ChannelSpec]) -> Generator[Trace, None, None]

    stream = get_radar_stream(input_filename)
    if stream == "RADnh3":
        header_t = namedtuple('radnh3_header',
                              'nsamp nchan vr0 vr1 choff ver resvd2')
        fmtstr='>HBBBBBB'
    elif stream == "RADnh5":
        header_t = namedtuple('radnh5_header',
                              'nsamp nchan vr0 vr1 choff ver resvd2 absix relix xinc rseq scount tscount')
        fmtstr='>HBBBBBBddfLHL'
    else:
        raise Exception("Invalid stream type for file: %s" % input_filename)

    fmtlen=struct.calcsize(fmtstr)
    # In theory, it's faster to do this here and only compile the format string once.
    header_struct = struct.Struct(fmtstr)

    # This way we don't waste time converting channel data that won't be used
    choffs = set() # type: Set[int]
    for p1cs in channel_specs:
        choffs.add(((p1cs.chan0in-1)//2)*2)
        choffs.add(((p1cs.chan1in-1)//2)*2)

    with open(input_filename, 'rb') as fd:
        while True:
            # read header
            buff = fd.read(fmtlen)
            if len(buff) < fmtlen:
                return

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
            assert(input_samples > 0 and input_samples < 10000)

            if choff in choffs:
                traces = np.fromfile(fd, dtype='>i2', count=input_samples*2)
                try:
                    traces.shape = (2, input_samples)
                except ValueError:
                    # didn't get enough data.
                    break

                # no use yet for headers....
                yield Trace(choff+1, traces[0][...])#, header, radnhx_header)
                yield Trace(choff+2, traces[1][...])#, header, radnhx_header)
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
        # File descriptors
        self.MagOutFD = None # type: Optional[BinaryIO]
        self.PhsOutFD = None # type: Optional[BinaryIO]
        self.MetaOutFD = None # type: Optional[BinaryIO]

        self.ChannelNum = channel
        self.MagScale = MagScale
        self.StackDepth = StackDepth
        self.IncoDepth = IncoDepth

        self.record_increment=self.StackDepth*self.IncoDepth
        self.record_idx=self.record_increment/2

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
        self.MetaOutFD = open(self.MetaFileName, 'w')
        self.MetaOutFD.write('#InputName = "' + input_filename + '"\n')
        logging.info("writing %s" % self.MagFileName)
        self.MagOutFD = open(self.MagFileName, 'wb')
        self.MetaOutFD.write('#MagName = "' + self.MagFileName + '"\n')

        if do_phase:
            self.PhsOutFD = open(self.PhsFileName, 'w')
            self.MetaOutFD.write('#PhsName = "' + self.PhsFileName + '"\n')

        ###### FIXME? - request hdr file from xlob and read/forward
        # TODO: collect this information into a data structure rather than reading directly from args.
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
            Phase.byteswap(True)
            Phase.tofile(self.PhsOutFD)

        if self.MagOutFD is not None:
            ScaledMag = np.int32(self.MagScale * np.log10(inco_trace.magnitude))
            ScaledMag.byteswap(True)
            ScaledMag.tofile(self.MagOutFD)

        if self.enable_meta_idx:
            self.MetaOutFD.write(str(self.record_idx) + "\n")
        self.record_idx += self.record_increment

    def close(self):
        # type: () -> None
        # flush stacks,
        # warn about incomplete stacks
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
    min_freq = round(2.5 * trunc_sweep_length / 50)
    max_freq = round(17.5 * trunc_sweep_length / 50)
    dfreq = max_freq - min_freq + 1
    hamming = np.sin(np.linspace(0, 1, num=dfreq) * np.pi)
    hfilter = np.hstack((np.zeros(min_freq),
                         hamming*2,
                         np.zeros(trunc_sweep_length - 2*max_freq - 1),
                         np.zeros(hamming.size),
                         np.zeros(min_freq - 1)))
    return hfilter

def main(args):
    # type: (Any) -> None
    if args.debug:
        LOGLEVEL=logging.DEBUG
    else:
        LOGLEVEL=logging.INFO
    logging.basicConfig(level=LOGLEVEL,
                    format='pik1: %(relativeCreated)8d [%(levelname)-5s] (%(process)d %(threadName)-10s) %(message)s',
                   )

    # Obtain reference chirp
    ref_chirp = get_ref_chirp(args.bandpass, args.output_samples)

    # Compute Hamming Filter
    hfilter = get_hfilter(args.output_samples)

    ## Filter Chirp
    # disable for production work
    # plt.plot(hfilter, range(0,filter.size))
    ref_chirp = np.multiply(ref_chirp, hfilter)

    channel_specs = parse_channels(args.channel_def)
    logging.debug("%r" % channel_specs)
    # Read traces from file
    # TODO: This is where we need to be able to read RADnh5 instead
    #tracegen = read_RADnh3_gen(args.infile, channel_specs, args.input_samples)
    tracegen = read_RADnhx_gen(args.infile, channel_specs)
    # Demultiplex stacks and generate coherently-stacked traces
    stackgen = stacks_gen(tracegen, channel_specs, args.StackDepth)
    # Dechirp coherent stacks
    dechirpgen = denoise_and_dechirp_gen(stackgen, ref_chirp, args.blanking,
                                         args.output_samples,
                                         not args.bandpass)
    # Incoherently stack
    istackgen = inco_stacks_gen(dechirpgen, channel_specs, args.IncoDepth,
                                args.output_samples, do_phase=args.output_phases)

    # Initialize output files
    outfiles = {}
    for p1cs in channel_specs:
        outfiles[p1cs.chanout] = PIK1OutputFile(args.outdir,
                                                p1cs.chanout,
                                                args.Scale,
                                                args.StackDepth,
                                                args.IncoDepth)
        outfiles[p1cs.chanout].open(args.infile, do_phase=args.output_phases)

    for rec in istackgen:
        if rec.channel in outfiles:
            outfiles[rec.channel].write_record(rec)
        else:
            raise Exception("Invalid output channel %d" % rec.channel)

    for p1cs in channel_specs:
        outfiles[p1cs.chanout].close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Pulse compress radar data')

    parser.add_argument('--outdir', required=True,
                        help='directory for output files')
    parser.add_argument('--infile', required=True,
                        help='filename of 2-byte radar file')
    # TODO(LEL): This argument is a PITA. Should be consistent channel spec,
    #    not optional int or string..
    parser.add_argument('--channel_def', required=True,
                        help="Channel Number to compress (counts from 1) OR string that's the channel spec")

    parser.add_argument('--output_samples', required=True, type=int,
                        help='Length of each output sweep (in samples)')
    parser.add_argument('--StackDepth', required=True, type=int,
                        help='coherent stacking depth for this output')
    parser.add_argument('--IncoDepth', required=True, type=int,
                        help='incoherent stacking depth for this output')

    parser.add_argument('--Scale', type=int, default=20000,
                        help='Output scale default is 1000*dB')
    parser.add_argument('--blanking', type=int, default=50,
                        help='Samples at the top of the record to blank out')

    parser.add_argument('--output_phases', action='store_true',
                        help="output phase, in addition to magnitude")
    parser.add_argument('--bandpass', action='store_true',
                        help='bandpass sampling, false is for legacy hicars. Disables cinterp and flips the chirp')
    parser.add_argument('--debug', action='store_true',
                        help='Print debugging messages')

    args = parser.parse_args()

    do_profile = False
    if do_profile:
        import os
        prof_file = "/tmp/pik1b.{0:d}.prof".format(os.getpid())
        cProfile.run('main(args)', prof_file)
        import pstats
        p = pstats.Stats(prof_file)
        p.sort_stats('cumulative').print_stats(50)
    else:
        main(args)
