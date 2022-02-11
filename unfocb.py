#!/usr/bin/env python3
#
# unfoc.py: Unfocused Radar Processor
# Output file for pik1 (4-byte signed integer, network order)

"""
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

try:
    import typing
    from typing import Any, BinaryIO, Dict, Generator, List, Optional, Set, Tuple
except ImportError: #pragma: no cover
    pass  # Not installed on melt ...


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

            yield (fpos, headerlen, header.choff) #, ctdata)
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



def dump_bxds_index(bxds_filename):
    """ Dump an index of the bxds to a file """
    #s_radinfo = 
    for radinfo in index_RADnhx_bxds_mmap(bxds_filename, stream):
        fpos, headerlen, rchoff, nsamp = radinfo
    

def read_RADnhx_gen2(bxds_filename, choff, stream=None):
    """ Return a sequence of traces from only one channel, from a bxds
    channel offset is a one-based index.

    Takes a bxds filename as an argument.
    Reads from this file and a ct file alongside it.

    """
    # Round to nearest 
    assert 1 <= choff
    choff1 = (choff - 1) % 2
    # Filter for listed channel offset in radar record
    chfilter = (choff - 1) - choff1 

    try:
        fd = open(bxds_filename, 'rb')
        mmbxds = mmap.mmap(fd.fileno(), 0, access=mmap.ACCESS_READ)

        ctgen = unfoc.gen_ct(bxds_filename)

        for radinfo, ctinfo in zip(index_RADnhx_bxds_mmap(bxds_filename, stream), ctgen):
            fpos, headerlen, rchoff, nsamp = radinfo

            if chfilter is not None and chfilter != rchoff:
                continue

            i0 = fpos + headerlen + choff1 * nsamp
            trace1 = np.frombuffer(mmbxds, dtype=np.int16, offset=i0, count=nsamp)

            yield unfoc.Trace(choff1, trace1, ctinfo)
    finally:
        mmbxds.close()
        fd.close()

def unfoc2(input_bxds, output_dir, channelspec):
    """ Generate one output channel of unfoc data and write it to the
    specified output directory with the standard file naming. """

def chunks(iterable, size=10):
    # https://stackoverflow.com/questions/24527006/\
    # split-a-generator-into-chunks-without-pre-walking-it/24527424
    iterator = iter(iterable)
    for first in iterator:
        yield itertools.chain([first], itertools.islice(iterator, size - 1))


def stack_coherent_chunk(coherent_chunk_gen, size=None, dtype=np.int32):
    """ Stack a chunk of traces coherently
    and return a Trace object and metadata from
    the input traces

    """
    stacked = None
    ctinfos = []
    for nrecs, rec in enumerate(coherent_chunk_gen):
        ctinfos.append(rec.ct)
        if stacked is None: # first stack
            stacked = rec.data.astype(dtype, copy=True)
        else:
            stacked += rec.data

    return stacked, nrecs+1, ctinfos

def stack_inco_chunk(inco_chunk_gen, dtype=np.int32):
    stacked = None

    # TODO: don't return the last stack if it's incomplete.
    for data, nrecs, ctinfos in inco_chunk_gen:
        #list_nrecs.append(nrecs)
        if stacked is None: # first stack
            stacked = data.astype(dtype, copy=True)
        else:
            stacked += data

    return stacked

def sum_traces(traces, dtype=np.int32, channel=-1):
    """
    TODO: make this part of the trace class?
    """
    data = traces[0].data.astype(np.int32, copy=False) + \
           traces[1].data.astype(np.int32, copy=False)
    return unfoc.Trace(channel=channel, data=data, ct=traces[0].ct)


def denoise_and_dechirp(data, args):
    # trace, nrecs+1, ctinfos
    return data


def main():
    # type: () -> None


    parser = argparse.ArgumentParser(description='Read ')

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
        genfunc = rad_traces
    else:
        genfunc = index_RADnhx_bxds

    #for ii, rec1 in enumerate(genfunc(args.infile, stream='RADnh5')):
    #    #records1.append(rec1)
    #    logging.debug("%r", rec1)

    # Read traces from channel 1
    channel_gen1 = read_RADnhx_gen2(args.infile, choff=1)
    # Read traces from channel 3
    channel_gen2 = read_RADnhx_gen2(args.infile, choff=3)

    # Combine two channels into one, and rewrite the output channel number
    f_sum_traces = lambda traces: sum_traces(traces, choff=1)
    sumchannel_gen = map(sum_traces, zip(channel_gen1, channel_gen2))

    # Stack traces coherently
    coherent_stack_size = 10
    coherent_chunker = chunks(sumchannel_gen, size=coherent_stack_size)
    coh_trace_gen = map(stack_coherent_chunk, coherent_chunker)

    # dechirp
    denoise_and_dechirp_lambda = lambda data: denoise_and_dechirp(data, args=None)
    denoise_gen = map(denoise_and_dechirp_lambda, coh_trace_gen)

    # incoherently stack
    incodepth = 5
    inco_chunker = chunks(denoise_gen, size=incodepth)
    igen1 = map(stack_inco_chunk, inco_chunker)

    # Run it all
    for ii, istack in enumerate(igen1):
        print(ii, istack.shape)


    logging.info("%d records", ii)
    elapsed_time = time.time() - t0
    logging.info("Elapsed time: %0.1f sec", elapsed_time)

    speed = size1 / elapsed_time / 1024 / 1024.
    logging.info("%0.3f MB/s", speed)




if __name__ == "__main__":
    sys.exit(main())
