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
from pathlib import Path

import unfoc

p1 = Path(__file__).parent.absolute()
sys.path.insert(1, str(p1))
from run_unfoc import setup_common_args

def main():
    # type: () -> None
    parser = argparse.ArgumentParser(description='Pulse compress 1-m resampled radar data with unfocused processor')

    setup_common_args(parser)

    parser.add_argument('-c', '--channel', required=True, type=int, help="Channel number")

    args = parser.parse_args()


    LOGLEVEL = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=LOGLEVEL, stream=sys.stdout,
                    format='unfoc: [%(levelname)-5s] %(message)s')

    unfoc.unfoc_1m_chan(args.outdir, args.input, args.channel, args.output_samples, args.stackdepth, args.incodepth,
          args.blanking, args.bandpass, scale=args.scale, output_phases=args.output_phases, nmax=args.nmax,
          buffering=args.buffering)


if __name__ == "__main__":
    sys.exit(main())
