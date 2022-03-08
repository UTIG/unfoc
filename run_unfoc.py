#!/usr/bin/env python3

"""
Command line interface for unfocused processing
"""

import argparse
import logging
import os
import sys
import time

import mmap

import unfoc.filter as filter



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
                    format='unfoc: [%(levelname)-5s] %(message)s')

    filter.unfoc(args.outdir, args.input, args.channels, args.output_samples, args.stackdepth, args.incodepth,
          args.blanking, args.bandpass, scale=args.scale, output_phases=False, nmax=args.nmax)


if __name__ == "__main__":
    sys.exit(main())
