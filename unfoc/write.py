#!/usr/bin/env python3

"""
Handler for writing unfoc output files

The primary class in this module, PIK1Output, handles writing the
collection of unfocused data products to a specified directory in
the format and organization expected for PIK1 on the UTIG hierarchy.

"""

import argparse
from collections import namedtuple
import logging
import struct
import sys

import numpy as np

try:
    import typing
    from typing import Any, BinaryIO, Dict, Generator, List, Optional, Set, Tuple
except ImportError: #pragma: no cover
    pass  # Not installed on melt ...

from unfoc.read import Trace, get_radar_stream, gen_ct

# Pairs channel and magnitude/phase data for an incoherent trace
# type: (int, np.ndarray, np.ndarray, List[tuple]) -> None
# Channel is now superfluous but NBD.
IncoherentTrace = namedtuple('IncoherentTrace', 'channel magnitude phase ct')


class PIK1Output:
    """
    Outputs magnitude and phase information for NDArray that comes through.
    Optionally computes the mean before writing
    It can be used as a tee filter and tap output at any phase in the
    processing stream.
    TODO: make names PEP8 compliant
    TODO: use os.path.join
    TODO: open/close context manager
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
        # output data type 4-byte big endian (network order)
        # TODO: change from swapping to using built-in spec
        #self.dtype = '>i4'

    def __del__(self):
        self.close()

    def open(self, input_filename, do_phase=True, do_index=True):
        # type: (str, bool, bool) -> None
        '''
        Opens all file descriptors and outputs the header for the metafile.
        Having args input is ugly, but it contains vars not used anywhere
        else in the class other than that header file.
        # TODO: configure logging verbosity?
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
        metadata = """
#ChannelNum = {0.ChannelNum:d}
#StackDepth = {0.StackDepth:d}
#IncoDepth = {0.IncoDepth:d}
#Scale = {0.MagScale:d}
#Log = TRUE
""".lstrip().format(self)
        self.MetaOutFD.write(metadata)


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
        # TODO: change to print
        if self.TraceNumbersFD is not None:
            self.TraceNumbersFD.write("%d\n" % (inco_trace.ct.seq,))

        if self.enable_meta_idx:
            self.MetaOutFD.write("%d\n" % (self.record_idx,))
        self.record_idx += self.record_increment

    def close(self):
        # type: () -> None
        # flush stacks,
        # TODO: use a for loop
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


