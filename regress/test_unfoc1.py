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


"""

import unittest
import sys
import os
import logging
import tempfile

cwd = os.path.dirname(__file__)
p = os.path.abspath(os.path.join(cwd, ".."))
sys.path.insert(1, p)

import unfoc
import unfocb

#TESTLIST = os.path.join(cwd, 'test_lists/available_radnh_bxds.txt')
TESTLIST = os.path.join(cwd, 'test_lists/tests_level0.txt')
OUTPUTDIR = os.path.abspath(os.path.join(cwd, 'covdata'))
DELETE_OUTPUT = True


def read_testlist(testlist):
    with open(testlist, 'rt') as fin:
        for line in fin:
            filename = line.strip()
            if not filename or filename[0] == '#': # skip commented lines
                continue
            yield parse_fileinfo(filename)

def parse_fileinfo(filename):
    # line = /disk/kea/WAIS/orig/xlob/WSB/JKB2e/MAWG01c/RADnh3/bxds
    parts = filename.split('/')
    snm = parts[-2]
    pst = '/'.join(parts[-5:-2])
    return pst, snm, filename


class UnfocBase(unittest.TestCase):
    def check_inputs_exist(self, bxds_input):
        self.assertTrue(os.path.exists(bxds_input))

    def check_logging(self, cm, nwarnings, nerrors):
        # Count number of warnings
        cm_nwarnings = len([1 for m in cm.records if m.levelno == logging.WARNING])
        cm_nerrors = len([1 for m in cm.records if m.levelno == logging.ERROR])
        self.assertEqual(cm_nwarnings, nwarnings)
        self.assertEqual(cm_nerrors, nerrors, msg="Had %d errors" % nerrors)

    def check_outputs_exist(self, basepath, channels):
        for channel in channels.split(','):
            magfile = os.path.join(basepath, 'Mag' + channel)
            self.assertTrue(os.path.exists(magfile))
            self.assertGreater(os.path.getsize(magfile), 0)


def is_bandpass(pst, snm):
    # this isn't exactly right but it's sometimes right.
    # and it'll at least let us try some different combinations
    is_bandpass = (snm == 'RADnh5')
    return is_bandpass

class TestUnfoc(UnfocBase):
    def run_unfoc(self, unfoc_params, testlist, outprefix):
        channels_hc = 'LoResInco1,LoResInco2'
        channels_marfa = 'LoResInco1,LoResInco2,LoResInco5,LoResInco6,LoResInco7,LoResInco8'
        os.makedirs(outprefix, exist_ok=True)
        tempdir = tempfile.mkdtemp(prefix=outprefix)
        try:
            #with tempfile.TemporaryDirectory(prefix=outprefix, delete=False) as tempdir:
            # Use list so we make sure it fully reads the test list before starting.
            for pst, snm, bxds_input in list(read_testlist(testlist)):
                with self.subTest(pst=pst, snm=snm):
                    radartype = unfoc.get_radar_type(bxds_input)
                    channels = channels_marfa if radartype == 'MARFA' else channels_hc

                    outdir = os.path.join(tempdir, pst)
                    os.makedirs(outdir, exist_ok=True)
                    self.check_inputs_exist(bxds_input)
                    unfocb.unfoc(infile=bxds_input, outdir=outdir, channels=channels, bandpass=is_bandpass(pst, snm), **unfoc_params)
                    self.check_outputs_exist(outdir, channels)
        finally:
            if DELETE_OUTPUT:
                import shutil
                shutil.rmtree(tempdir)
        return tempdir

    def test1(self):
        testlist = os.path.join(cwd, 'test_lists/tests_level0.txt')
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
        return self.run_unfoc(unfoc_params, testlist, outprefix)

    def test_multi(self):
        testlist = os.path.join(cwd, 'test_lists/tests_level0.txt')
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
        return self.run_unfoc(unfoc_params, testlist, outprefix)


def main():
    logging.basicConfig(level=logging.INFO, stream=sys.stdout)
    unittest.main()


if __name__ == "__main__":
    main()
