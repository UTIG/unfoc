#!/usr/bin/env python3


"""
Unit testing of RADnh3 and RADnh5 bxds reader

TODO: cull some of the testing to only be done on one shorter transect rather than running through the really long ones.
(like indexing really only needs to be done on one transect)

"""

import unittest
import sys
import os
import logging
import itertools

import numpy as np

#from test_unfoc1 import read_testlist, parse_fileinfo, UnfocBase

cwd = os.path.dirname(__file__)
unfoc_path = os.path.join(cwd, '..')
sys.path.insert(1, os.path.abspath(unfoc_path))

import unfoc.read as read

#TESTLIST = os.path.join(cwd, 'test_lists/available_radnh_bxds.txt')
TESTLIST = os.path.join(cwd, 'test_lists/tests_level1.txt')

# List of tests for running just a single example of each type.
TESTLIST0 = os.path.join(cwd, 'test_lists/tests_level0.txt')




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


class RadBxdsBase(unittest.TestCase):
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

class TestParsers(RadBxdsBase):
    """ Test the functions and generators that parse the packet structure of the bxds file.
    Test these on a wide range of input stimulus """
    def run_compare_gen_parsers(self, bxds_input):
        self.check_input_exists(bxds_input)
        gen1 = read.index_RADnhx_bxds_mmap(bxds_input)
        gen2 = read.index_RADnhx_bxds(bxds_input)

        for ii, (data1, data2) in enumerate(itertools.zip_longest(gen1, gen2)):
            self.assertEqual(data1, data2, msg="parse mismatch at record %d" % (ii,))

    def test_compare_gen_parsers(self):
        """ Check that indexers parse the bxds same """
        for pst, snm, bxds_input in list(read_testlist(TESTLIST)):
            with self.subTest(pst=pst, snm=snm):
                self.run_compare_gen_parsers(bxds_input)


class TestClass1(RadBxdsBase):
    """ Test the RadBxds class. """
    def setUp(self):
        channel = 1
        pst, snm, bxds_input = list(read_testlist(TESTLIST0))[0]
        self.rread = read.RadBxds(bxds_input, channel=channel)

    def test_indexing(self):
        """ Check that slicing works consistently with how numpy does it. """
        rread = self.rread
        ntraces = rread.size()
        self.assertGreater(ntraces, 1) # number of records

        for n in range(ntraces):
            # Indexing a single element returns an ndarray of singleton dimension.
            trace1 = rread[n]
            self.assertTrue(trace1.shape == (3200,) or trace1.shape == (3437,))

            m = n - ntraces # Check that negative indices work.
            self.assertLess(m, 0)
            trace2 = rread[m]
            self.assertTrue(trace1.shape == (3200,) or trace1.shape == (3437,))

            self.assertTrue(np.array_equal(trace1, trace2))

        # Indexing beyond the end of the radargram raises an IndexError
        with self.assertRaises(IndexError):
            rread[ntraces]

    def test_slicing1(self):
        # Show that a slice comes out to the correct dimensions.
        traces1 = self.rread[0:10]
        self.assertEqual(traces1.shape[0], 10)
        self.assertTrue(traces1.shape[1] == 3200 or traces1.shape[1] == 3437)

        traces2 = self.rread[0:10]
        self.assertEqual(traces2.shape, (10, traces1.shape[1]))

        traces2 = self.rread[10:20]
        self.assertEqual(traces2.shape, (10, traces1.shape[1]))

        traces2 = self.rread[:10]
        self.assertEqual(traces2.shape, (10, traces1.shape[1]))

        traces2 = self.rread[-3:-1]
        self.assertEqual(traces2.shape, (2, traces1.shape[1]))



    def test_slicing_long(self):
        # Slicing beyond the end of the radargram doesn't raise an error
        # (This matches the behavior with lists and other sequences)

        ntraces = len(self.rread)
        expected_len = 10
        n1 = ntraces - expected_len
        traces1 = self.rread[n1:ntraces]

        self.assertEqual(traces1.shape[0], expected_len)
        # If you slice beyond the last element, the slice is automatically
        # limited to the last element.
        traces2 = self.rread[n1:ntraces + 15]
        self.assertEqual(traces2.shape[0], expected_len)
        self.assertTrue(np.array_equal(traces1, traces2))
        

    def test_slicing_backwards(self):
        expected_len = 10
        traces1 = self.rread[5:2]
        self.assertEqual(traces1.shape[0], 0)
        traces1 = self.rread[-1:-5]
        self.assertEqual(traces1.shape[0], 0)

    def test_slicing_stride(self):
        # Slicing with a stride/step
        traces1 = self.rread[0:10]
        s = traces1.shape
        self.assertEqual(s[0], 10)
        self.assertTrue(s[1] == 3200 or s[1] == 3437)
        
        traces2 = self.rread[0:10:2]
        s = traces2.shape
        self.assertEqual(s[0], 5)
        self.assertEqual(traces1.shape[1], traces2.shape[1])
        self.assertTrue(np.array_equal(traces1[0, :], traces2[0, :]))
        self.assertTrue(np.array_equal(traces1[2, :], traces2[1, :]))

    def test_ct(self):
        for ii in range(len(self.rread)):
            ct1 = self.rread.ct(ii)
            self.assertEqual(len(ct1), 2)
        self.assertGreaterEqual(len(self.rread.cts_), len(self.rread))

class TestRADjh1Class(TestClass1):
    def setUp(self):
        channel = 1
        testlist1 = os.path.join(cwd, 'test_lists/tests_radjh1.txt')
        pst, snm, bxds_input = list(read_testlist(testlist1))[0]
        self.rread = read.RADjh1Bxds(bxds_input, channel=channel)

class TestRadBxds(RadBxdsBase):
    """ Run tests on many different bxdses """
    def test_index_generator1(self):
        return self.run_index_generator(read.index_RADnhx_bxds_mmap)

    #@unittest.skip("Unneeded since we test equality later")
    #def test_index_generator2(self):
    #    return self.run_index_generator(read.index_RADnhx_bxds)

    def test_read(self):
        for pst, snm, bxds_input in list(read_testlist(TESTLIST)):
            for channel in (1, 2):
                with self.subTest(pst=pst, snm=snm, channel=channel):
                    self.check_input_exists(bxds_input)
                    trace_p = None
                    for ii, trace in enumerate(read.read_RADnhx_gen(bxds_input, channel=channel)):
                        self.assertEqual(channel, trace.channel)
                        self.assertEqual(len(trace.data.shape), 1)
                        self.assertGreaterEqual(trace.data.shape[0], 3200)
                        self.assertGreater(trace.ct.seq, 0)
                        if trace_p is not None:
                            self.assertGreater(trace.ct.seq, trace_p.ct.seq)
                            self.assertGreaterEqual(trace.ct.tim, trace_p.ct.tim)
                        trace_p = trace
                    self.assertGreater(ii, 0)


    def test_compare_reads(self):
        """ Check that we get the same data out of the function as from the class interface """
        for pst, snm, bxds_input in list(read_testlist(TESTLIST0)):
            for channel in (1, 2):
                with self.subTest(pst=pst, snm=snm, channel=channel):
                    rread = read.RadBxds(bxds_input, channel=channel)
                    self.assertGreater(len(rread), 1) # number of records

                    for ii, trace1 in enumerate(read.read_RADnhx_gen(bxds_input, channel=channel)):
                        msg = "mismatch at rread[%d]" % (ii,)
                        trace2 = rread[ii]
                        self.assertEqual(trace1.data.shape, trace2.shape, msg=msg)
                        self.assertTrue(np.array_equal(trace1.data, trace2), msg=msg)
                        self.assertEqual(trace1.ct, rread.ct(ii), msg=msg)

                    self.assertEqual(len(rread), ii+1)

    def test_get_rad_stream(self):
        # Try detecting bxds.
        pass



def main():
    logging.basicConfig(level=logging.INFO, stream=sys.stdout)
    unittest.main()


if __name__ == "__main__":
    main()
