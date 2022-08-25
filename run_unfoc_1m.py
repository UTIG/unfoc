#!/usr/bin/env python3

"""
Command line interface for unfocused processing, to generate 1 meter, unfocused pik1
from S2_FIL 1-meter resampled radargrams.

Example usage:

./run_unfoc_1m.py -c 2 -i $WAIS/targ/xtra/ICP9/FOC/ngg_202205_f1/S2_FIL/TOT3/JKB2s/X15a/bxds2.i -o ./out1

# Defaults taken from $WAIS/code/xtra/ASE1/CMP/pik1.1m/Run_All_Hi_Res.sc:

  SweepLength=3200
  truncSweepLength=3200
  IName=""
  QName=""
  MagName="$radargram/MagHiResInco$OutChn"
  PhsName="$radargram/PhaseHiResInco$OutChn"
  MetaName="$radargram/MagHiResInco$OutChn.meta"
  StackDepth=1
  IncoDepth=1
  CenterMult=1
  MaxDepth=150
  StartSweep=1
  EndSweep=Inf
  StartSamp=1
  EndSamp=3437
  Scale=20000
  Log=TRUE

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
    parser = argparse.ArgumentParser(description='Pulse compress 1-m resampled radar data with unfocused processor')

    parser.add_argument('-o', '--outdir', required=True,
                        help='directory for output files')
    parser.add_argument('-i', '--input', required=True,
                        help='filename of 2-byte radar file (bxdsN.i)')

    parser.add_argument('-c', '--channel', required=True, type=int, help="Channel number")

    parser.add_argument('--output_samples', default=3200, type=int,
                        help='Length of each output sweep (in samples)')

    parser.add_argument('--stackdepth', default=1, type=int,
                        help='coherent stacking depth')
    parser.add_argument('--incodepth', default=1, type=int,
                        help='incoherent stacking depth')

    parser.add_argument('--scale', type=int, default=20000,
                        help='Output scale default is 1000*dB')
    parser.add_argument('--blanking', type=int, default=50,
                        help='Samples at the top of the record to blank out. Negative number blanks bottom')

    parser.add_argument('--output_phases', action='store_true',
                        help="output phase, in addition to magnitude")
    parser.add_argument('--bandpass', action='store_true',
                        help='Process bandpass-sampled data (for use with MARFA data, not for use with legacy HiCARS/HiCARS2 data). Disable cinterp and flips the chirp.')
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

    filter.unfoc_1m_chan(args.outdir, args.input, args.channel, args.output_samples, args.stackdepth, args.incodepth,
          args.blanking, args.bandpass, scale=args.scale, output_phases=args.output_phases, nmax=args.nmax)


if __name__ == "__main__":
    sys.exit(main())
