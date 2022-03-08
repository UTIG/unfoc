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
import logging
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

# Pairs channel and magnitude/phase data for an incoherent trace
# type: (int, np.ndarray, np.ndarray, List[tuple]) -> None
IncoherentTrace = namedtuple('IncoherentTrace', 'channel magnitude phase ct')


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
