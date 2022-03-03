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


def unfoc(outdir, infile, channels, output_samples, stackdepth, incodepth,
          blanking, bandpass, scale=20000, output_phases=False, nmax=0):
    """ Generate one output channel of unfoc data and write it to the
    specified output directory with the standard file naming. """

    # pass through if this is a legacy ChannelSpec object
    if isinstance(channels, str):
        radartype = unfoc_mod.get_radar_type(infile)
        logging.info("Radar type: " + radartype)
        channel_specs = unfoc_mod.get_utig_channels(channels, radar=radartype)
    else:
        channel_specs = channels

    # we can multiprocess here.
    for channel_spec in channel_specs:
        unfoc_chan(outdir, infile, channel_spec, output_samples, stackdepth, incodepth,
          blanking, bandpass, scale, output_phases, nmax)

def unfoc_chan(outdir, infile, p1cs, output_samples, stackdepth, incodepth,
          blanking, bandpass, scale=20000, output_phases=False, nmax=0):

    # Obtain reference chirp
    ref_chirp = unfoc_mod.get_ref_chirp(bandpass, output_samples)

    # Compute Hamming Filter
    hfilter = unfoc_mod.get_hfilter(output_samples)

    # Filter Chirp
    ref_chirp *= hfilter

    sumchannel_gen = setup_bxds_reader(infile, p1cs)

    # Stack traces coherently
    coherent_chunker = chunks(sumchannel_gen, size=stackdepth)
    coh_trace_gen = map(stack_coherent_chunk, coherent_chunker)

    # dechirp
    f_denoise_and_dechirp = lambda data: denoise_and_dechirp(data,
                ref_chirp=ref_chirp, blanking=blanking, output_samples=output_samples, do_cinterp=not bandpass)

    denoise_gen = map(f_denoise_and_dechirp, coh_trace_gen)

    # incoherently stack
    inco_chunker = chunks(denoise_gen, size=incodepth)
    igen1 = map(stack_inco_chunk, inco_chunker)


    p1out = unfoc_mod.PIK1OutputFile(outdir, channel=p1cs.chanout,
                                     MagScale=20000, StackDepth=stackdepth, IncoDepth=incodepth)

    p1out.open(infile, do_phase=True, do_index=True)

    t0 = time.time()
    # Run it all
    for ii, istack in enumerate(igen1):
        p1out.write_record(istack)
        if nmax > 0 and ii >= nmax:
            break
    p1out.close()


def chunks(iterable, size=10):
    # https://stackoverflow.com/questions/24527006/\
    # split-a-generator-into-chunks-without-pre-walking-it/24527424
    iterator = iter(iterable)
    for first in iterator:
        yield itertools.chain([first], itertools.islice(iterator, size - 1))


def stack_coherent_chunk(coherent_chunk_gen, size=None, dtype=np.float64):
    """ Stack a chunk of traces coherently
    and return a Trace object, the number of records,
    and metadata from  the input traces

    """
    stacked = None
    ctinfos = []
    for nrecs, rec in enumerate(coherent_chunk_gen):
        ctinfos.append(rec.ct)
        if stacked is None: # first stack
            stacked = rec.data.astype(dtype, copy=True)
        else:
            stacked += rec.data

    stacked /= (nrecs+1)
    return stacked, nrecs+1, ctinfos


def stack_inco_chunk(inco_chunk_gen, dtype=None):
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

    magnitude /= len(traces)
    ct = list_cts[0] # len(list_cts) // 2] apparently we just return the first one
    phs = np.angle(traces[len(traces)//2])

    itrace = unfoc_mod.IncoherentTrace(channel=0, magnitude=magnitude, phase=phs, ct=ct)
    return itrace

def sum_traces(traces, dtype=np.int32, channel=-1):
    """
    TODO: make this part of the trace class?
    """
    data = traces[0].data.astype(np.int32, copy=False) + \
           traces[1].data.astype(np.int32, copy=False)
    return unfoc_mod.Trace(channel=channel, data=data, ct=traces[0].ct)


def denoise_and_dechirp(coherent_data, *args, **kwargs):
    """ Assumes that the input coherent_data has the same output parameter
    structure as stack_coherent_chunk """
    stacked, nrecs, ctinfos = coherent_data
    dechirped = unfoc_mod.denoise_and_dechirp(stacked.astype(np.double), *args, **kwargs)

    # Should we use the first one?
    #return unfoc_mod.Trace(stacked.channel, dechirped, ctinfos[0])
    return dechirped, nrecs, ctinfos
    #return unfoc_mod.Trace(stacked.channel, dechirped, stacked.ct)

    # trace, nrecs+1, ctinfos
    #return data


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
        channel_gen1 = radbxds.read_RADnhx_gen(bxdsfile, channel=channel_specs.chan0in)
        # Rewrite the channel number
        f_traceout = lambda trace: unfoc_mod.Trace(channel=channel_specs.chanout, data=trace.data, ct=trace.ct)
        reader_gen = map(f_traceout, channel_gen1)

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
    parser.add_argument('--nmax', default=0, type=int,
                        help="Maximum number of stacks to output (usually used for testing")
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
