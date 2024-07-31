#!/usr/bin/env python3


"""
filtering operations for unfocused processing

Perform unfocused processing on a radar file

"""

import argparse
from collections import namedtuple
import logging
import os
import sys
import time
import itertools
import multiprocessing

import mmap

import numpy as np

from unfoc.parse_channels import get_utig_channels, PIK1ChannelSpec
from unfoc.write import PIK1Output, IncoherentTrace
import unfoc.dechirp as dechirp
import unfoc.read as read


def unfoc(outdir, infile, channels, output_samples, stackdepth, incodepth,
          blanking, bandpass, scale=20000, output_phases=False, nmax=0, processes=1):
    """ Generate all requested output channels of unfoc data and write it to the
    specified output directory with the standard file naming.

    The channels parameter is a string containing a comma-separated list
    of output channel specifications.  See parse_channels.get_utig_channels for
    definitions.

    The processes parameter specifies how many channels are processed in parallel.

    See unfoc_chan for other details on other input parameters.

    return value: None
    """

    # pass through if this is a legacy ChannelSpec object
    if isinstance(channels, str):
        radartype, data_channels = read.get_radar_type(infile)
        logging.info("Radar type: %s", radartype)
        channel_specs = get_utig_channels(channels, radar=radartype, input_channels=data_channels)
    else:
        channel_specs = channels
    delay = 0.01 if processes > 1 else 0.0 # cosmetically delay so channels output in same order
    unfoc_args = lambda p1cs: (outdir, infile, p1cs, output_samples, stackdepth, incodepth,
                       blanking, bandpass, scale, output_phases, nmax, (p1cs.chanout-1)*delay)
    gen_args = map(unfoc_args, channel_specs)
    if processes <= 1:
        for _ in map(unfoc_chan_, gen_args):
            pass
    else:
        with multiprocessing.Pool(processes) as pool:
            for _ in pool.map(unfoc_chan_, gen_args):
                pass

def unfoc_chan_(args):
    return unfoc_chan(*args)

def unfoc_chan(outdir, infile, p1cs, output_samples, stackdepth, incodepth,
          blanking, bandpass, scale=20000, output_phases=False, nmax=0, delay=0.):

    """ Generate one output channel of unfoc data
    outdir: output directory where data will be placed
    infile: input bxds file
    p1cs: channel specification for this output file (PIK1ChannelSpec object)
    output_samples: number of output samples to place into the output file
    stackdepth: coherent stacking depth
    incodepth: incoherent stacking depth
    blanking: samples to blank
    scale: output magnitude scaling factor
    output_phases: if true, also output phase data
    nmax: max number of output samples to prcess, then quit (usually for testing).


    """
    if delay > 0:
        time.sleep(delay)
    # Obtain reference chirp
    ref_chirp = dechirp.get_ref_chirp(bandpass, output_samples)

    # Compute Hamming Filter
    hfilter = dechirp.get_hfilter(output_samples)

    # Filter Chirp
    ref_chirp *= hfilter

    # Setup bxds file reader to get the correct files.
    sumchannel_gen, output_tag = setup_reader(infile, p1cs)

    # Stack traces coherently
    coherent_chunker = chunks(sumchannel_gen, size=stackdepth)
    # the last coherent trace is chosen when doing summed channels
    coh_trace_gen = map(stack_coherent_chunk, coherent_chunker)

    # dechirp
    f_denoise_and_dechirp = lambda data: denoise_and_dechirp(data,
                ref_chirp=ref_chirp, blanking=blanking, output_samples=output_samples, do_cinterp=not bandpass)

    denoise_gen = map(f_denoise_and_dechirp, coh_trace_gen)

    # incoherently stack
    inco_chunker = chunks(denoise_gen, size=incodepth, incomplete=False)
    f_stack_inco_chunk = lambda stack: stack_inco_chunk(stack, p1cs.chanout, output_phases=output_phases)
    igen1 = map(f_stack_inco_chunk, inco_chunker)


    os.makedirs(outdir, exist_ok=True)
    with PIK1Output(infile, outdir, channel=p1cs.chanout, magscale=scale,
                    stackdepth=stackdepth, incodepth=incodepth,
                    do_phase=output_phases, do_index=True, tag=output_tag) as p1out:

        for ii, istack in enumerate(igen1):
            p1out.write_record(istack)
            if nmax > 0 and ii >= nmax:
                break


def unfoc_1m_chan(outdir, infile, chanout, output_samples, stackdepth, incodepth,
          blanking, bandpass, scale=20000, output_phases=False, nmax=0):

    """ Generate one output channel of unfoc 1-meter-spaced data
    See unfoc_chan for parameters

    chanout: Output channel name -- used for setting output filename within outdir
    """

    p1cs = PIK1ChannelSpec(chanout, 0, 0, 0, 0)
    return unfoc_chan(outdir, infile, p1cs, output_samples, stackdepth, incodepth,
                      blanking, bandpass, scale, output_phases, nmax)


def chunks(iterable, size=10, incomplete=False):
    # https://stackoverflow.com/questions/24527006/\
    # split-a-generator-into-chunks-without-pre-walking-it/24527424
    iterator = iter(iterable)
    for first in iterator:
        #yield itertools.chain([first], itertools.islice(iterator, size - 1))
        data = list(itertools.chain([first], itertools.islice(iterator, size - 1)))
        if not incomplete and len(data) != size:
            # Don't return incomplete stacks
            continue
        yield data


def stack_coherent_chunk(coherent_chunk_gen, ct_type='mid', dtype=np.float64):
    """ Stack a chunk of traces coherently
    and return a Trace object, the number of records,
    and metadata from  the input traces

    ct_type can be any of 'mid' or 'all'

    """
    stacked = None
    ctinfos = []
    for nrecs, rec in enumerate(coherent_chunk_gen):
        ctinfos.append(rec.ct)
        if stacked is None: # first stack
            stacked = rec.data.astype(dtype, copy=True)
        else:
            stacked += rec.data

    stackdepth = len(ctinfos)
    if ct_type == 'mid':
        ctv = [ctinfos[stackdepth // 2],]
    else:
        assert ct_type == 'all'
        ctv = ctinfos

    stacked /= stackdepth
    return stacked, stackdepth, ctv


def stack_inco_chunk(inco_chunk_gen, channel, dtype=None, output_phases=False):
    """ Incoherently stack.
    Calculate magnitude and track phase.
    Return a summed magnitude trace,
    a ct record,
    and if phase is requested, the phase in the center of the stack
    """

    magnitude = None
    traces = []
    list_cts = []

    # TODO: don't return the last stack if it's incomplete.
    for data, _, ctinfos in inco_chunk_gen:
        traces.append(data)
        list_cts.extend(ctinfos)
        if magnitude is None: # first stack
            magnitude = np.abs(data)
        else:
            magnitude += np.abs(data)

    ntraces = len(traces)
    magnitude /= ntraces
    ct = list_cts[ntraces//2]
    if output_phases:
        phs = np.angle(traces[ntraces//2])
    else:
        phs = None

    itrace = IncoherentTrace(channel=channel, magnitude=magnitude, phase=phs, ct=ct)
    return itrace

def denoise_and_dechirp(coherent_data, *args, **kwargs):
    """ Assumes that the input coherent_data has the same output parameter
    structure as stack_coherent_chunk """
    stacked, nrecs, ctinfos = coherent_data
    output_samples = kwargs['output_samples']
    dechirped = dechirp.denoise_and_dechirp(stacked[0:output_samples].astype(np.double, copy=False), *args, **kwargs)

    return dechirped, nrecs, ctinfos


def setup_reader(bxdsfile, channel_specs, input_type=None):
    """ Return value: tuple(trace generator, output_tag) """
    if input_type is None:
        # auto-detect input type
        if bxdsfile.endswith('.i'):
            input_type = 'S2_FIL' # resampled data from S2_FIL
        else:
            input_type = 'orig'

    if input_type == 'S2_FIL':
        return read.read_1m_gen(bxdsfile, channel=channel_specs.chanout), 'HiResInco'
    else:
        assert input_type == 'orig'
        return setup_bxds_reader(bxdsfile, channel_specs), 'LoResInco'


def setup_bxds_reader(bxdsfile, channel_specs):
    """
    Set up the generators for reading from a bxds file and producing
    pairwise-summed traces or individual unsummed traces
    TODO: update PIK1ChannelSpec to have less redundant parameters since we
    make these assertions here.

    """
    # PIK1ChannelSpec(chanout=1, chan0in=1, scalef0=1, chan1in=3, scalef1=1)
    assert channel_specs.scalef0 == 0 or channel_specs.scalef0 == 1
    assert channel_specs.scalef1 == 0 or channel_specs.scalef1 == 1

    if channel_specs.scalef0 == 1 and channel_specs.scalef1 == 1:
        # sum channels

        # Read traces from first channel of a bxds file
        channel_gen1 = read.read_RADnhx_gen(bxdsfile, channel=channel_specs.chan0in)
        # Read traces from another channel of a bxds file
        channel_gen2 = read.read_RADnhx_gen(bxdsfile, channel=channel_specs.chan1in)
        # Combine two channels into one, and rewrite the output channel number
        reader_gen = map(sum_traces, channel_gen1, channel_gen2)
    else:
        # If we're only doing one channel, it better be chan0
        assert channel_specs.scalef0 == 1
        reader_gen = read.read_RADnhx_gen(bxdsfile, channel=channel_specs.chan0in)

    return reader_gen


def sum_traces(trace1, trace2, dtype=np.int32):
    data = trace1.data.astype(np.int32, copy=False) + \
           trace2.data.astype(np.int32, copy=False)
    return read.Trace(channel=-1, data=data, ct=trace1.ct)
