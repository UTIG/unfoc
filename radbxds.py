#!/usr/bin/env python3


"""
radbxds.py

Read a radar bxds file and return traces
Read the structure of a radar file
"""

import argparse
from collections import namedtuple
import logging
import os
#import gzip
import struct
import sys
import time
import itertools

import mmap

import numpy as np

import unfoc


HEADER_FORMATS = {
    'RADnh3': {
        'header_t': namedtuple('radnh3_header',
                              'nsamp nchan vr0 vr1 choff ver resvd2'),
        'fmtstr':  '>HBBBBBB',
    },
    'RADnh5': {
        'header_t': namedtuple('radnh5_header',
                              'nsamp nchan vr0 vr1 choff ver resvd2 absix relix xinc rseq scount tscount'),
        'fmtstr': '>HBBBBBBddfLHL',
    }
}


# Read individual traces out of RADnh3 or RADnh5 file
# TODO: This doesn't yet filter on channels, which should be OK - the
# downstream users ALSO check channel.
def index_RADnhx_bxds(input_filename, stream=None):
    # type: (str) -> Generator[tuple]
    """ Read the positions of packets within a RADnh3 and RADnh5 bxds file
    and return these as a generator
    Routines can then use this to seek to the correct location in a file.
    """

    if stream is None:
        stream = unfoc.get_radar_stream(input_filename)

    if stream not in HEADER_FORMATS: # pragma: no cover
        raise ValueError("Invalid stream type for file: %s" % input_filename)

    header_t = HEADER_FORMATS[stream]['header_t']
    fmtstr = HEADER_FORMATS[stream]['fmtstr']
    fmtlen = struct.calcsize(fmtstr)
    # In theory, it's faster to do this here and only compile the format string once.
    header_struct = struct.Struct(fmtstr)

    with open(input_filename, 'rb') as fd:
        fpos = 0
        while True:
            # read header
            fpos = fd.tell()
            buff = fd.read(fmtlen)
            if len(buff) < fmtlen:
                return
            #ctdata = next(genct) if genct else (None, None)
            header = header_t._make(header_struct.unpack(buff))
            headerlen = fmtlen
            if stream == "RADnh5":
                # We don't do anything with the timestamps yet, but we need
                # to skip over 'em to get to the radar data.
                if header.tscount > 0:
                    # timestamp count should be on the order of the number of stacks.
                    assert(header.tscount < 100)
                    #tstamps = struct.unpack_from('>{:d}d'.format(tscount), fd.read(8*tscount))
                    fd.seek(8*header.tscount, os.SEEK_CUR)
                    headerlen += 8*header.tscount

            input_samples = header.nsamp
            assert 0 < input_samples < 10000

            yield fpos, headerlen, header.choff, header.nsamp #, ctdata)
            # Skip the rest of this record without yielding any values
            fd.seek(4*input_samples, os.SEEK_CUR)



# Read individual traces out of RADnh3 or RADnh5 file
# TODO: This doesn't yet filter on channels, which should be OK - the
# downstream users ALSO check channel.
def index_RADnhx_bxds_mmap(input_filename, stream=None):
    # type: (str) -> Generator[tuple]
    """ Read the positions of packets within a RADnh3 and RADnh5 bxds file
    and return these as a generator
    Routines can then use this to seek to the correct location in a file.
    """

    if stream is None:
        stream = unfoc.get_radar_stream(input_filename)

    if stream not in HEADER_FORMATS: # pragma: no cover
        raise ValueError("Invalid stream type for file: %s" % input_filename)

    header_t = HEADER_FORMATS[stream]['header_t']
    fmtstr = HEADER_FORMATS[stream]['fmtstr']

    # In theory, it's faster to do this here and only compile the format string once.
    header_struct = struct.Struct(fmtstr)
    fmtlen = header_struct.size
    with open(input_filename, 'rb') as fd:
        mmbxds = mmap.mmap(fd.fileno(), 0, access=mmap.ACCESS_READ)
        filelen = mmbxds.size()
        fpos = 0
        while fpos + fmtlen < filelen:
            # read header
            headerdata = header_struct.unpack_from(mmbxds, fpos)
            fpos_next = fpos + header_struct.size
            #ctdata = next(genct) if genct else (None, None)
            header = header_t._make(headerdata)
            headerlen = header_struct.size
            if stream == "RADnh5":
                # We don't do anything with the timestamps yet, but we need
                # to skip over 'em to get to the radar data.
                if header.tscount > 0:
                    # timestamp count should be on the order of the number of stacks.
                    assert header.tscount < 100
                    fpos_next += 8*header.tscount
                    headerlen += 8*header.tscount

            assert 0 < header.nsamp < 10000

            yield fpos, headerlen, header.choff, header.nsamp #, ctdata)
            # Skip the rest of this record without yielding any values
            fpos_next += 4*header.nsamp
            fpos = fpos_next
        mmbxds.close()



def read_RADnhx_gen(bxds_filename, channel, stream=None):
    """ Return a sequence of traces from only one channel, from a bxds
    channel offset is a one-based index.

    Takes a bxds filename as an argument.
    Reads from this file and a ct file alongside it.

    """
    assert 1 <= channel # one-based channel number
    choff1 = (channel - 1) % 2
    # Filter for listed channel offset in radar record
    chfilter = (channel - 1) - choff1

    try:
        fd = open(bxds_filename, 'rb')
        mmbxds = mmap.mmap(fd.fileno(), 0, access=mmap.ACCESS_READ)

        ctgen = unfoc.gen_ct(bxds_filename)

        for radinfo, ctinfo in zip(index_RADnhx_bxds_mmap(bxds_filename, stream), ctgen):
            fpos, headerlen, rchoff, nsamp = radinfo

            if chfilter != rchoff:
                continue

            i0 = fpos + headerlen + choff1 * nsamp
            trace1 = np.frombuffer(mmbxds, dtype=np.int16, offset=i0, count=nsamp)

            yield unfoc.Trace(channel, trace1, ctinfo)
    finally:
        mmbxds.close()
        fd.close()


class BxdsMemmap:
    """ Reader for random access to traces in a RADnh3 or RADnh5 bxds as numpy arrays """
    def __init__(self, filename, channel, stream=None):
        """ Initialize the reader to access one channel's records """
        self.filename = filename
        self.index = None
        self.mmbxds = None
        self.fd = None
        self.channel = channel
        self.choff = 2*(self.channel & 0x01)
        self.fd = open(filename, 'rb')
        self.mmbxds = mmap.mmap(self.fd.fileno(), 0, access=mmap.ACCESS_READ)
        self.load_index(stream)

    def __del__(self):
        if self.mmbxds:
            self.mmbxds.close()
        if self.fd:
            self.fd.close()

    def load_index(self, stream=None):
        """ Load record index for bxds 
        This ignores sequence numbers and assumes that there are no
        out-of-order radar records within a channel.
        In theory we could sort the records by sequence number.
        """
        self.index = []
        choff = self.channel - self.channel & 1
        for item in index_RADnhx_bxds_mmap(self.filename, stream=stream):
            if item[2] == choff:
                self.index.append(item)

    def __getitem__(self, idx):
        """ Return a trace.
        Example:
        bxdsmm = BxdsMemmap(bxdsfile, channel=2)
        # Get the 6th trace for channel 5
        trace = bxdsmm[5]
         """
        if idx > len(self.index):
            return None
        fpos, headerlen, rchoff, nsamp = self.index[idx]
        # offset in bytes for even vs odd channels
        i0 = fpos + headerlen + self.choff * nsamp
        trace1 = np.frombuffer(self.mmbxds, dtype="<i2", offset=i0, count=nsamp)
        return trace1 #yield unfoc.Trace(self.channel, trace1, ctinfo)

def test_class(bxdsfile=None):
    if bxdsfile is None:
        bxdsfile = '/disk/kea/WAIS/orig/xlob/NIS4/IBH0e/X84b/RADnh5/bxds'
    for channel in (0, 1, 2, 3):
        print("Checking channel", channel)
        rread = BxdsMemmap(bxdsfile, channel=channel, stream='RADnh5')
        print(len(rread.index), " records")
        for n, idx in enumerate(rread.index):
            s = rread[n].shape
            if n == 0:
                print(s)

def main():
    # type: () -> None
    parser = argparse.ArgumentParser(description='Library for reading a radar BXDS file')

    parser.add_argument('-i', '--infile', required=True,
                        help='filename of 2-byte radar file (bxds file)')

    parser.add_argument('--debug', action='store_true',
                        help='Print debugging messages')

    parser.add_argument('--mmap', action='store_true',
                        help='Use mmap method')

    args = parser.parse_args()


    LOGLEVEL = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=LOGLEVEL, stream=sys.stdout,
                    format='unfoc: [%(levelname)-5s] %(message)s',
                   )

    t0 = time.time()
    records1 = []
    size1 = os.path.getsize(args.infile)

    if args.mmap:
        genfunc = index_RADnhx_bxds_mmap
    else:
        genfunc = index_RADnhx_bxds

    for ii, data in enumerate(genfunc(args.infile)):
        pass

    logging.info("%d records", ii+1)
    test_class(args.infile)

    elapsed_time = time.time() - t0
    logging.info("Elapsed time: %0.1f sec", elapsed_time)

    speed = size1 / elapsed_time / 1024 / 1024.
    logging.info("%0.3f MB/s", speed)




if __name__ == "__main__":
    sys.exit(main())
