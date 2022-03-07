#!/usr/bin/env python3
#
# unfoc.py: Unfocused Radar Processor
# Output file for pik1 (4-byte signed integer, network order)

"""
Perform unfocused processing on a radar file

# TODO: make IncoherentTrace and Trace namedtuples instead of classes.
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

import unfoc as unfoc_mod
import radbxds
from radbxds import Trace
#Trace = namedtuple('Trace', 'data ct')


def unfoc(outdir, infile, channels, output_samples, stackdepth, incodepth,
          blanking, bandpass, scale=20000, output_phases=False, nmax=0, processes=1):
    """ Generate all requested output channels of unfoc data and write it to the
    specified output directory with the standard file naming. """

    # pass through if this is a legacy ChannelSpec object
    if isinstance(channels, str):
        radartype = unfoc_mod.get_radar_type(infile)
        logging.info("Radar type: " + radartype)
        channel_specs = unfoc_mod.get_utig_channels(channels, radar=radartype)
    else:
        channel_specs = channels

    unfoc_args = lambda p1cs: (outdir, infile, p1cs, output_samples, stackdepth, incodepth,
                       blanking, bandpass, scale, output_phases, nmax)
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
          blanking, bandpass, scale=20000, output_phases=False, nmax=0):

    """ Generate one output channel of unfoc data
    outdir: output directory where data will be placed
    infile: input bxds file
    p1cs: channel specification for this output file
    output_samples: number of output samples to place into the output file
    stackdepth: coherent stacking depth
    incodepth: incoherent stacking depth
    blanking: samples to blank
    scale: output magnitude scaling factor
    output_phases: if true, also output phase data
    nmax: max number of output samples to prcess, then quit (usually for testing).
    """

    # Obtain reference chirp
    ref_chirp = unfoc_mod.get_ref_chirp(bandpass, output_samples)

    # Compute Hamming Filter
    hfilter = unfoc_mod.get_hfilter(output_samples)

    # Filter Chirp
    ref_chirp *= hfilter

    sumchannel_gen = setup_bxds_reader(infile, p1cs)

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
    p1out = unfoc_mod.PIK1OutputFile(outdir, channel=p1cs.chanout,
                                     MagScale=scale, StackDepth=stackdepth, IncoDepth=incodepth)

    p1out.open(infile, do_phase=output_phases, do_index=True)
    for ii, istack in enumerate(igen1):
        p1out.write_record(istack)
        if nmax > 0 and ii >= nmax:
            break
    p1out.close()


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
        # convert to cdouble to match older version
        phs = np.angle(traces[ntraces//2].astype(np.cdouble))
    else:
        phs = None

    itrace = unfoc_mod.IncoherentTrace(channel=channel, magnitude=magnitude, phase=phs, ct=ct)
    return itrace

def sum_traces(trace1, trace2, dtype=np.int32):
    data = trace1.data.astype(np.int32, copy=False) + \
           trace2.data.astype(np.int32, copy=False)
    return Trace(channel=-1, data=data, ct=trace1.ct)


def denoise_and_dechirp(coherent_data, *args, **kwargs):
    """ Assumes that the input coherent_data has the same output parameter
    structure as stack_coherent_chunk """
    stacked, nrecs, ctinfos = coherent_data
    output_samples = kwargs['output_samples']
    dechirped = unfoc_mod.denoise_and_dechirp(stacked[0:output_samples].astype(np.double), *args, **kwargs)

    return dechirped, nrecs, ctinfos


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
        reader_gen = map(sum_traces, channel_gen1, channel_gen2)
    else:
        # If we're only doing one channel, it better be chan0
        assert channel_specs.scalef0 == 1
        reader_gen = radbxds.read_RADnhx_gen(bxdsfile, channel=channel_specs.chan0in)

    return reader_gen



def main():
    # type: () -> None



    parser = argparse.ArgumentParser(description='Pulse compress radar data with unfocused processor')

    parser.add_argument('-o', '--outdir', required=True,
                        help='directory for output files')
    parser.add_argument('-i', '--input', required=True,
                        help='filename of 2-byte radar file (bxds file)')

    cgroup = parser.add_mutually_exclusive_group(required=True)
    cgroup.add_argument('--channel_def',
                        help="Channel spec string. This option is deprecated. Use --channels")
    cgroup.add_argument('--channels',
                        help="comma-separated channels to produce, such as 'LoResInco1,LoResInco2'")

    parser.add_argument('--output_samples', default=3200, type=int,
                        help='Length of each output sweep (in samples)')
    parser.add_argument('--stackdepth', required=True, type=int,
                        help='coherent stacking depth')
    parser.add_argument('--incodepth', required=True, type=int,
                        help='incoherent stacking depth')

    parser.add_argument('--scale', type=int, default=20000,
                        help='Output scale default is 1000*dB')
    parser.add_argument('--blanking', type=int, default=50,
                        help='Samples at the top of the record to blank out. Negative number blanks bottom')

    parser.add_argument('--output_phases', action='store_true',
                        help="output phase, in addition to magnitude")
    parser.add_argument('--bandpass', action='store_true',
                        help='bandpass sampling, false is for legacy hicars. Disables cinterp and flips the chirp')
    parser.add_argument('-j', '--jobs', default=1, type=int,
                        help="Max number of CPUs to use for processing")
    parser.add_argument('--nmax', default=0, type=int,
                        help="Maximum number of stacks to output (usually used for testing)")
    parser.add_argument('--debug', action='store_true',
                        help='Print debugging messages')

    args = parser.parse_args()


    LOGLEVEL = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=LOGLEVEL, stream=sys.stdout,
                    format='unfoc: [%(levelname)-5s] %(message)s',
                   )

    unfoc(args.outdir, args.input, args.channels, args.output_samples, args.stackdepth, args.incodepth,
          args.blanking, args.bandpass, scale=args.scale, output_phases=False, nmax=args.nmax)


if __name__ == "__main__":
    sys.exit(main())
