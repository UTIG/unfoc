#!/usr/bin/env python3


"""
Unit testing of RADnh3 and RADnh5 bxds reader


"""

import unittest
import sys
import os
import logging
import itertools

import numpy as np

#from test_unfoc1 import read_testlist, parse_fileinfo, UnfocBase

cwd = os.path.dirname(__file__)
p = os.path.abspath(os.path.join(cwd, ".."))
sys.path.insert(1, p)

import radbxds

#TESTLIST = os.path.join(cwd, 'test_lists/available_radnh_bxds.txt')
TESTLIST = os.path.join(cwd, 'test_lists/tests_level1.txt')




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

class TestRadBxds(unittest.TestCase):
    def check_input_exists(self, bxds_input):
        self.assertTrue(os.path.exists(bxds_input))
        self.assertGreater(os.path.getsize(bxds_input), 0)

    def run_index_generator(self, genfunc):
        for pst, snm, bxds_input in list(read_testlist(TESTLIST)):
            with self.subTest(pst=pst, snm=snm):
                self.run_index_generator_bxds(genfunc, bxds_input)

    def run_index_generator_bxds(self, genfunc, bxds_input):
        self.check_input_exists(bxds_input)
        data_prev = None
        for ii, data in enumerate(genfunc(bxds_input)):
            if data_prev is not None:
                #fpos, headerlen, choff, nsamp = data
                self.assertLess(data_prev[0], data[0])
                self.assertGreater(data[1], 0)
                self.assertGreaterEqual(data[2], 0)
                self.assertGreater(data[3], 0)
            data_prev = data
        self.assertGreater(ii, 0)

    def test_index_generator1(self):
        return self.run_index_generator(radbxds.index_RADnhx_bxds_mmap)

    #@unittest.skip("Unneeded since we test equality later")
    #def test_index_generator2(self):
    #    return self.run_index_generator(radbxds.index_RADnhx_bxds)

    def test_read(self):
        for pst, snm, bxds_input in list(read_testlist(TESTLIST)):
            for channel in (1, 2):
                with self.subTest(pst=pst, snm=snm, channel=channel):
                    self.check_input_exists(bxds_input)
                    trace_p = None
                    for ii, trace in enumerate(radbxds.read_RADnhx_gen(bxds_input, channel=channel)):
                        self.assertEqual(channel, trace.channel)
                        self.assertEqual(len(trace.data.shape), 1)
                        self.assertGreaterEqual(trace.data.shape[0], 3200)
                        self.assertGreater(trace.ct.seq, 0)
                        if trace_p is not None:
                            self.assertGreater(trace.ct.seq, trace_p.ct.seq)
                            self.assertGreaterEqual(trace.ct.tim, trace_p.ct.tim)
                        trace_p = trace
                    self.assertGreater(ii, 0)


    def test_compare_gen_indexers(self):
        """ Check that indexers parse the bxds same """
        for pst, snm, bxds_input in list(read_testlist(TESTLIST)):
            with self.subTest(pst=pst, snm=snm):
                self.run_compare_gen_indexers(bxds_input)

    def run_compare_gen_indexers(self, bxds_input):
        self.check_input_exists(bxds_input)
        gen1 = radbxds.index_RADnhx_bxds_mmap(bxds_input)
        gen2 = radbxds.index_RADnhx_bxds(bxds_input)

        for ii, (data1, data2) in enumerate(itertools.zip_longest(gen1, gen2)):
            self.assertEqual(data1, data2, msg="parse mismatch at record %d" % (ii,))


    def test_memmap_class(self):
        for pst, snm, bxds_input in list(read_testlist(TESTLIST)):
            for channel in (1, 2):
                with self.subTest(pst=pst, snm=snm, channel=channel):
                    rread = radbxds.RadBxds(bxds_input, channel=channel)
                    self.assertGreater(len(rread.index), 1) # number of records

                    for n, idx in enumerate(rread.index):
                        s = rread[n].shape[0]
                        self.assertTrue(s == 3200 or s == 3437)
                    self.assertIsNone(rread[-1])

    def test_compare_reads(self):
        for pst, snm, bxds_input in list(read_testlist(TESTLIST)):
            for channel in (1, 2):
                with self.subTest(pst=pst, snm=snm, channel=channel):
                    rread = radbxds.RadBxds(bxds_input, channel=channel)
                    self.assertGreater(rread.size(), 1) # number of records
                    self.assertIsNone(rread[-1])
                    self.assertIsNone(rread[rread.size()])

                    for ii, trace1 in enumerate(radbxds.read_RADnhx_gen(bxds_input, channel=channel)):
                        msg = "mismatch at rread[%d]" % (ii,)
                        trace2 = rread[ii]
                        self.assertEqual(len(trace1.data), len(trace2.data), msg=msg)
                        self.assertTrue(np.array_equal(trace1.data, trace2.data))
                        # TODO: compare ct
                    self.assertEqual(len(rread.index), ii+1)


def main():
    logging.basicConfig(level=logging.INFO, stream=sys.stdout)
    unittest.main()


if __name__ == "__main__":
    main()
