#!/usr/bin/env python
#
# test_unfoc.py: Run some test cases for unfocused processor
# Output file for pik1 (4-byte signed integer, network order)
#
# All output files are 4-byte network-order.
#


import argparse
#from collections import namedtuple
#import cProfile
import logging
import os
#import struct
import sys
import gzip

import numpy as np
from scipy import signal
from scipy import stats

try:
    import typing
    from typing import Any, BinaryIO, Dict, Generator, List, Optional, Set, Tuple
except ImportError: #pragma: no cover
    pass  # Not installed on melt ...

import unfoc

def main():
    # type: (Any) -> None
    parser = argparse.ArgumentParser(description='Pulse compress radar data')

    parser.add_argument('--output', default='./test_unfoc_data',
                        help='directory for output files')
    parser.add_argument('--nmax', default=0, type=int,
                        help="Maximum number of stacks to output (usually used for testing")
    parser.add_argument('--debug', action='store_true',
                        help='Print debugging messages')

    parser.add_argument('--WAIS', default=os.getenv('WAIS', '/disk/kea/WAIS'),
                        help='WAIS environment variable')

    args = parser.parse_args()

    LOGLEVEL = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=LOGLEVEL, stream=sys.stdout,
                    format='test_unfoc: %(relativeCreated)8d [%(levelname)-5s] %(message)s',
                   )

    test_err(args.WAIS)
    test_ct(args.WAIS)


def test_err(WAIS):
    """ Test with invalid data """
    # Non-radar data with a bxds
    infile = os.path.join(WAIS, "targ/xped/ICP8/breakout/ELSA/F05/KNX/JKB2q/X100a/GPSnc1/bxds")
    
    channel_specs = unfoc.parse_channels("1")
    try:
        tracegen = unfoc.read_RADnhx_gen(infile, channel_specs, b_parse_ct=True)
        for rec in tracegen:
            pass
        assert False
    except ValueError:
        pass

def test_ct(WAIS):
    """ Test reading ct and bxds simultaneously, and compare timestamps """
    channel_specs = unfoc.parse_channels("1")
    logging.info("channel_spec: %r" % channel_specs)
    # Read traces from file,  one from each season

    infiles = [
        os.path.join(WAIS, "targ/xped/ICP8/breakout/MARFA/F05/KNX/JKB2q/X100a/RADnh5/bxds"),
        os.path.join(WAIS, "targ/xped/SRH1/breakout/MARFA/F06/DEV/JKB2t/X80a/RADnh5/bxds"),
        os.path.join(WAIS, "targ/xped/ASE2/breakout/HERA/F01/GTZ/IBH0d/Y38a/RADnh5/bxds"),
        os.path.join(WAIS, "targ/xped/KPR2/breakout/HERA/J345/KPR2/IBH0a/J345c/RADnh5/bxds"),
        os.path.join(WAIS, "targ/xped/KRT2/breakout/HERA/F08/NIS4/IBH0e/X84b/RADnh5/bxds"),
        os.path.join(WAIS, "targ/xped/ICP10/breakout/MARFA/F04/DLY/JKB2u/Y61a/RADnh5/bxds"),
        os.path.join(WAIS, "targ/xped/ICP9/breakout/MARFA/F03/ICP9/JKB2s/F03T01a/RADnh5/bxds"),
    ]

    for infile in infiles:
        tracegen = unfoc.read_RADnhx_gen(infile, channel_specs, b_parse_ct=True)
        logging.info("Reading " + infile)
        t0 = None
        t1 = None
        tims = []
        tdigs = []
        assert 'RADnh5' in infile
        for i, rec in enumerate(tracegen):
            if i >= 100000:
                break
            #if (i % 10000 == 0):
            #    logging.info(i)
            # assume RADnh5
            if t0 is None:
                t0 = rec.header.absix
            if t1 is None:
                t1 = rec.ct.tim
            # Digitizer time
            tdig = (rec.header.absix - t0) + rec.header.relix
            # CT time
            tim = (rec.ct.tim - t1) / 100000.0
            tdigs.append(tdig)
            tims.append(tim)

        logging.info("Got {:d} records".format(len(tims)))

        slope, intercept, r_value, p_value, std_err = stats.linregress(tdigs, tims)
        logging.info("slope deviation from 1 = {:0.6g} err = {:0.6g}".format(1 - slope, std_err))


        assert abs(1.0 - slope) < 1e-7
        assert abs(1.0 - r_value) < 1e-7
        assert abs(std_err) < 1e-7

    # This should also work for RADnh3, but no secondary metadata.
    # Should we test that too?  Seems unnecessary


if __name__ == "__main__":
    main()
