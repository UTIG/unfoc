#!/usr/bin/env python3

"""
Handler for writing unfoc output files: unfocused radar data products
in the UTIG PIK1 format.

The primary class in this module, PIK1Output, handles writing the
collection of unfocused data products (magnitude, phase, metadata,
trace indices) to a specified directory in the format and organization
expected for PIK1 on the UTIG hierarchy.

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
    Writer for incoherent radar data products in UTIG-compatible PIK1 format.

    This class outputs unfocused magnitude and phase radar data along with metadata
    and trace indices to a directory. Output files include `.meta`, `.tracenumbers`,
    binary magnitude and optional phase files. Designed for use in radar data 
    processing pipelines, and can act as a tee or tap in streaming architectures.

    Parameters
    ----------
    input_filename : str
        Path to the input file (used for metadata reference).
    outdir : str
        Output directory for all files.
    channel : int
        Channel number to include in output filenames.
    magscale : int
        Scaling factor applied to log-magnitude values.
    stackdepth : int
        Number of stacks per output record.
    incodepth : int
        Number of incoherent traces per stack.
    do_phase : bool, optional
        If True, writes phase output. Default is True.
    do_index : bool, optional
        If True, writes trace indices to the metadata file. Default is True.
    tag : str, optional
        Custom tag for output filenames. Default is 'LoResInco'.
    loglevel : int, optional
        Logging level (e.g., logging.INFO). Default is INFO.

    Attributes
    ----------
    infile : str
        Input file name used in metadata headers.
    channel : int
        Channel number for trace assignment.
    magscale : int
        Logarithmic magnitude scaling factor.
    stackdepth : int
        Number of stacks per output record.
    incodepth : int
        Number of incoherent traces per stack.
    record_increment : int
        Trace index increment between records.
    record_idx : float
        Current trace index to be written into metadata.
    do_phase : bool
        Whether to write phase output.
    enable_meta_idx : bool
        Whether to write trace indices into the metadata.
    mag_fd, phs_fd, meta_fd, tracenumbers_fd : BinaryIO or None
        File descriptors for output binary and text files.

    Methods
    -------
    write_record(inco_trace)
        Write magnitude, phase, and trace number for one incoherent trace.

    close()
        Close all file descriptors.

    __enter__(), __exit__()
        Context manager methods for safely using this class with `with` blocks.

    """
    def __init__(self, input_filename, outdir, channel, magscale,
                 stackdepth, incodepth, do_phase=True, do_index=True,
                 tag='LoResInco', loglevel=logging.INFO):
        """
        Initializes the PIK1Output writer and opens output file descriptors.

        Opens all file descriptors and outputs the header for the metafile.
        Having args input is ugly, but it contains vars not used anywhere
        else in the class other than that header file.
            
        """
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

        self.close()

        channel = self.channel
        magfilename = os.path.join(outdir, "Mag{0:s}{1:d}".format(tag, channel))
        phsfilename = os.path.join(outdir, "Phs{0:s}{1:d}".format(tag, channel))
        metafilename = os.path.join(outdir,"Mag{0:s}{1:d}.meta".format(tag, channel))
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
        """
        Write one incoherent trace to magnitude, phase, trace number, and metadata files.

        Magnitude and phase arrays are scaled and saved as 4-byte big-endian integers.
        Trace numbers (`ct.seq`) and optional metadata indices are written as text.

        Parameters
        ----------
        inco_trace : IncoherentTrace
            Named tuple with fields: channel (int), magnitude (ndarray),
            phase (ndarray), and ct (object with `seq`).

        Notes
        -----
        Output binary format: 4-byte big-endian integers (network byte order).
        """

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
        """Close all output file descriptors, if open."""

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
