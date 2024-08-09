#!/usr/bin/env python3


"""
Unit testing of unfoc

The main high-level test runs unfoc with --nmax 1000 (to limit size)
on a large number of transects of varying parameter combinations,
including undersampling/nonundersampling and checks whether

Checks to be implemented
- the output channels are each of the expected size (greater than zero, multiple of 4*3200bytes)
- output channels are approximately of the expected size (for expected number of channels)
- outputs of different channels are equal in size
- output contents of different channels differ?

# TODO: add multipol test examples
"""

import unittest
import sys
import os
import logging
import tempfile
from typing import Dict, Any
from pathlib import Path

from test_read import read_testlist

cwd = Path(__file__).parent
unfoc_path = cwd / '..' / 'src'
sys.path.insert(1, str(unfoc_path.absolute()))

import unfoc

#TESTLIST = os.path.join(cwd, 'test_lists/available_radnh_bxds.txt')
TESTLIST = cwd / 'test_lists' / 'tests_level0.txt'
OUTPUTDIR = os.path.abspath(os.path.join(cwd, 'covdata'))

# Normal output
DELETE_OUTPUT = True
LOGLEVEL = logging.WARNING
# Debugging
#DELETE_OUTPUT = False
#LOGLEVEL = logging.INFO

class UnfocBase(unittest.TestCase):
    def check_inputs_exist(self, bxds_input:str):
        """ Confirm that prerequisite inputs exist """
        self.assertTrue(os.path.exists(bxds_input))

    #def check_logging(self, cm, nwarnings, nerrors):
    #    # Count number of warnings
    #    cm_nwarnings = len([1 for m in cm.records if m.levelno == logging.WARNING])
    #    cm_nerrors = len([1 for m in cm.records if m.levelno == logging.ERROR])
    #    self.assertEqual(cm_nwarnings, nwarnings)
    #    self.assertEqual(cm_nerrors, nerrors, msg="Had %d errors" % nerrors)

    def check_outputs_exist(self, basepath:str, channels:str, has_phase:bool=True):
        """ Check that outputs were written for specified channels """
        for channel in channels.split(','):
            self.assertTrue(channel.startswith('LoResInco'))
            expected_files = [os.path.join(basepath, 'Mag' + channel)]
            if has_phase:
                expected_files.append(os.path.join(basepath, 'Phs' + channel))

            for file in expected_files:
                self.assertTrue(os.path.exists(file), msg=file)
                self.assertGreater(os.path.getsize(file), 0, msg=file)

    def run_unfoc(self, unfoc_params:Dict[str,Any], testlist:Path, outprefix:str):
        """ Run unfocused processor for tests and for all channels with given parameters """
        cdict = {
            'HiCARS2': 'LoResInco1,LoResInco2',
            'MARFA': 'LoResInco1,LoResInco2,LoResInco5,LoResInco6,LoResInco7,LoResInco8'
        }
        os.makedirs(outprefix, exist_ok=True)
        tempdir = tempfile.mkdtemp(prefix=outprefix)
        try:
            #with tempfile.TemporaryDirectory(prefix=outprefix, delete=False) as tempdir:
            # Use list so we make sure it fully reads the test list before starting.
            for pst, snm, bxds_input in list(read_testlist(testlist)):
                with self.subTest(pst=pst, snm=snm):
                    radartype, data_channels = unfoc.get_radar_type(bxds_input)
                    channels = cdict[radartype]
                    outdir = os.path.join(tempdir, pst)
                    os.makedirs(outdir, exist_ok=True)
                    self.check_inputs_exist(bxds_input)
                    unfoc.unfoc(infile=bxds_input, outdir=outdir, channels=channels,
                                 bandpass=is_bandpass(pst, snm), **unfoc_params)
                    self.check_outputs_exist(outdir, channels, unfoc_params['output_phases'])
        finally:
            if DELETE_OUTPUT:
                import shutil
                shutil.rmtree(tempdir)
        return tempdir


def is_bandpass(pst:str, snm:str):
    """ Return true if transect and stream is uses bandpass sampling """
    # this isn't exactly right but it's sometimes right.
    # and it'll at least let us try some different combinations
    return snm == 'RADnh5'

class TestUnfoc(UnfocBase):
    """ Test Unfocused processing functions """
    def test1(self):
        outprefix = os.path.join(OUTPUTDIR, 'unfoc_test1_')
        unfoc_params = {
            'output_samples': 3200,
            'stackdepth': 10,
            'incodepth': 5,
            'blanking': 50,
            'nmax': 0,
            'output_phases': True,
            'processes': 1,
        }
        return self.run_unfoc(unfoc_params, TESTLIST, outprefix)

    # TODO: test that single threaded comes out the same as multithreaded
    def test_multi(self):
        outprefix = os.path.join(OUTPUTDIR, 'unfoc_test1_')
        unfoc_params = {
            'output_samples': 3200,
            'stackdepth': 10,
            'incodepth': 5,
            'blanking': 50,
            'nmax': 0,
            'output_phases': True,
            'processes': 8,
        }
        return self.run_unfoc(unfoc_params, TESTLIST, outprefix)

    # TODO: test that single threaded comes out the same as multithreaded
    def test_reverseblank(self):
        outprefix = os.path.join(OUTPUTDIR, 'unfoc_test1_')
        unfoc_params = {
            'output_samples': 3200,
            'stackdepth': 10,
            'incodepth': 5,
            'blanking': -50,
            'nmax': 100,
            'output_phases': False,
            'processes': 1,
        }
        return self.run_unfoc(unfoc_params, TESTLIST, outprefix)

    # TODO: test that single threaded comes out the same as multithreaded
    def test_noblank(self):
        outprefix = os.path.join(OUTPUTDIR, 'unfoc_test1_')
        unfoc_params = {
            'output_samples': 3200,
            'stackdepth': 10,
            'incodepth': 5,
            'blanking': 0,
            'nmax': 100,
            'output_phases': False,
            'processes': 1,
        }
        return self.run_unfoc(unfoc_params, TESTLIST, outprefix)


    def test_denoise_burst(self):
        """ Enable denoising filter """
        outprefix = os.path.join(OUTPUTDIR, 'unfoc_test1_')
        unfoc_params = {
            'output_samples': 3200,
            'stackdepth': 10,
            'incodepth': 5,
            'blanking': 50,
            'nmax': 100,
            'output_phases': False,
            'processes': 1,
            'denoise': 'burst',
        }
        return self.run_unfoc(unfoc_params, TESTLIST, outprefix)

def main():
    logging.basicConfig(level=LOGLEVEL, stream=sys.stdout)
    unittest.main()


if __name__ == "__main__":
    main()
