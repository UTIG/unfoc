#!/usr/bin/env python3

"""
Handler for writing unfoc output files

The primary class in this module, PIK1Output, handles writing the
collection of unfocused data products to a specified directory in
the format and organization expected for PIK1 on the UTIG hierarchy.

"""

from collections import namedtuple
import logging
import sys
import os

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
    """
    def __init__(self, input_filename, outdir, channel, magscale,
                 stackdepth, incodepth, do_phase=True, do_index=True, loglevel=logging.INFO):
        # type: (str, int, int, int, int) -> None
        # File descriptors
        self.mag_fd = None # type: Optional[BinaryIO]
        self.phs_fd = None # type: Optional[BinaryIO]
        self.meta_fd = None # type: Optional[BinaryIO]
        self.tracenumbers_fd = None # type: Optional[BinaryIO]

        self.infile = input_filename
        self.channel = channel
        self.magscale = magscale
        self.stackdepth = stackdepth
        self.incodepth = incodepth

        self.record_increment = self.stackdepth*self.incodepth
        self.record_idx = self.record_increment/2

        self.do_phase = do_phase
        self.enable_meta_idx = do_index

        # type: (str, bool, bool) -> None
        '''
        Opens all file descriptors and outputs the header for the metafile.
        Having args input is ugly, but it contains vars not used anywhere
        else in the class other than that header file.
        '''
        self.close()

        channel = self.channel
        magfilename = os.path.join(outdir, "MagLoResInco{0:d}".format(channel))
        phsfilename = os.path.join(outdir, "PhsLoResInco{0:d}".format(channel))
        metafilename = os.path.join(outdir,"MagLoResInco{0:d}.meta".format(channel))
        tracesfilename = os.path.join(outdir, "TraceNumbers{0:d}".format(channel))


        ## Open Input and Output files
        self.meta_fd = open(metafilename, 'wt')
        self.meta_fd.write('#InputName = "' + self.infile + '"\n')
        logging.log(loglevel, "writing %s", magfilename)
        self.mag_fd = open(magfilename, 'wb')
        self.meta_fd.write('#MagName = "' + magfilename + '"\n')

        if self.do_phase:
            self.phs_fd = open(phsfilename, 'wt')
            self.meta_fd.write('#PhsName = "' + phsfilename + '"\n')

        self.tracenumbers_fd = open(tracesfilename, 'wt')

        ###### FIXME? - request hdr file from xlob and read/forward
        metadata = """
#ChannelNum = {0.channel:d}
#StackDepth = {0.stackdepth:d}
#IncoDepth = {0.incodepth:d}
#Scale = {0.magscale:d}
#Log = TRUE
""".lstrip().format(self)
        self.meta_fd.write(metadata)


    def write_record(self, inco_trace):
        # type: (IncoherentTrace) -> None
        # Write component files if enabled
        # output data type 4-byte big endian (network order)
        if self.phs_fd is not None:
            phase = (inco_trace.phase * 16777216).astype('>i4')
            phase.tofile(self.phs_fd)

        if self.mag_fd is not None:
            scaled_mag = (self.magscale * np.log10(inco_trace.magnitude)).astype('>i4')
            scaled_mag.tofile(self.mag_fd)

        if self.tracenumbers_fd is not None:
            self.tracenumbers_fd.write("%d\n" % inco_trace.ct.seq)

        if self.enable_meta_idx:
            self.meta_fd.write("%d\n" % self.record_idx)
        self.record_idx += self.record_increment

    def __del__(self):
        self.close()

    def close(self):
        # type: () -> None
        # close file handles
        for fd in ('mag_fd', 'phs_fd', 'meta_fd', 'tracenumbers_fd'):
            if getattr(self, fd) is not None:
                getattr(self, fd).close()
                setattr(self, fd, None)

    def __enter__(self):
        return self

    def __exit__(self, exc_type,exc_value, exc_traceback):
        self.close()
