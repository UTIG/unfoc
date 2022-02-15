#!/usr/bin/env python3
#
# unfoc.py: Unfocused Radar Processor
# Output file for pik1 (4-byte signed integer, network order)

"""
Perform unfocused processing on a radar file
"""

import argparse
from collections import namedtuple
import logging
import os
import sys
import time
import itertools

import mmap

import numpy as np

import unfoc as unfoc_mod
import radbxds


def unfoc(input_file, output_dir, channelspec):
    """ Generate one output channel of unfoc data and write it to the
    specified output directory with the standard file naming. """

    sumchannel_gen = setup_bxds_reader(input_file, p1cs)

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


    # TODO: should we write to a temporary directory first and then move?
    # Setup output files
    filename_mag = os.path.join(output_dir, 'MagLoresInco%d' % channelspec.chanout)
    filename_meta = os.path.join(output_dir, 'MagLoresInco%d.meta' % channelspec.chanout)
    filename_phs = os.path.join(output_dir, 'PhsLoresInco%d' % channelspec.chanout)
    filename_trc = os.path.join(output_dir, 'TraceNumbers%d' % channelspec.chanout)


    t0 = time.time()
    # Run it all
    for ii, istack in enumerate(igen1):
        print(ii, istack.shape)



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
    return unfoc_mod.Trace(channel=channel, data=data, ct=traces[0].ct)


def denoise_and_dechirp(data, args):
    # trace, nrecs+1, ctinfos
    return data


def setup_bxds_reader(bxdsfile, channel_specs):
    """
    Set up the generators for reading from a bxds file and producing
    pairwise-summed traces or individual unsummed traces
    """
    # PIK1ChannelSpec(chanout=1, chan0in=1, scalef0=1, chan1in=3, scalef1=1)
    assert channel_specs.scalef0 == 0 or channel_specs.scalef0 == 1
    assert channel_specs.scalef1 == 0 or channel_specs.scalef1 == 1

    if channel_specs.scalef0 == 1 and channel_specs.scalef1 == 1:
        # sum channels

        # Read traces from first channel of a bxds file
        channel_gen1 = radbxds.read_RADnhx_gen(bxdsfile, channel=channel_specs.chan0in)
        # Read traces from another channel of a bxds file
        channel_gen2 = radbxds.read_RADnhx_gen(bxdsfile, channel=channel_specs.chan1in)
        # Combine two channels into one, and rewrite the output channel number
        f_sum_traces = lambda traces: sum_traces(traces, choff=channel_specs.chanout)
        reader_gen = map(sum_traces, zip(channel_gen1, channel_gen2))
    else:
        # If we're only doing one channel, it better be chan0
        assert channel_specs.scalef0 == 1
        channel_gen1 = radbxds.read_RADnhx_gen2(bxdsfile, choff=channel_specs.chan0in)
        # Rewrite the channel number
        f_traceout = lambda trace: unfoc_mod.Trace(channel=channel_specs.chanout, data=trace.data, ct=trace.ct)
        reader_gen = map(f_traceout, channel_gen1)

    return reader_gen



def main():
    # type: () -> None


    parser = argparse.ArgumentParser(description='Unfocused processor ')

    parser.add_argument('-i', '--infile', required=True,
                        help='filename of 2-byte radar file (bxds file)')

    parser.add_argument('--debug', action='store_true',
                        help='Print debugging messages')

    #parser.add_argument('--mmap', action='store_true',
    #                    help='Use mmap method')

    args = parser.parse_args()


    LOGLEVEL = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=LOGLEVEL, stream=sys.stdout,
                    format='unfoc: [%(levelname)-5s] %(message)s',
                   )

    records1 = []
    size1 = os.path.getsize(args.infile)

    p1cs = unfoc_mod.PIK1ChannelSpec(chanout=1, chan0in=1, scalef0=1, chan1in=3, scalef1=1)
    sumchannel_gen = setup_bxds_reader(args.infile, p1cs)

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

    t0 = time.time()
    # Run it all
    for ii, istack in enumerate(igen1):
        print(ii, istack.shape)


    elapsed_time = time.time() - t0
    logging.info("Elapsed time: %0.1f sec", elapsed_time)

    speed = size1 / elapsed_time / 1024 / 1024.
    logging.info("%0.3f MB/s", speed)
    logging.info("%d records", ii)




if __name__ == "__main__":
    sys.exit(main())
