#!/usr/bin/env python3


"""
Unit testing of RADnh3 and RADnh5 bxds reader

"""

import unittest
import sys
import os
import logging
import itertools
from pathlib import Path
from typing import Tuple
import tempfile
import shutil
import gzip

import numpy as np

#from test_unfoc1 import read_testlist, parse_fileinfo, UnfocBase

cwd = Path(__file__).parent
unfoc_path = cwd / '..' / 'src'
sys.path.insert(1, str(unfoc_path.absolute()))

import unfoc

#TESTLIST = os.path.join(cwd, 'test_lists/available_radnh_bxds.txt')
TESTLIST = cwd / 'test_lists' / 'tests_level1.txt'

# List of tests for running just a single example of each type.
TESTLIST0 = cwd / 'test_lists' / 'tests_level0.txt'

WAIS = Path(os.getenv('WAIS', '/disk/kea/WAIS'))


def read_testlist(testlist: Path):
    """ Yield file information from a test list """
    with testlist.open('rt') as fin:
        for line in fin:
            filename = line.strip()
            if not filename or filename[0] == '#': # skip commented lines
                continue
            yield parse_fileinfo(filename)

def parse_fileinfo(filename: str) -> Tuple[str,str,str]:
    """ Split the pst, streamname, and whole filename """
    # line = /disk/kea/WAIS/orig/xlob/WSB/JKB2e/MAWG01c/RADnh3/bxds
    parts = filename.split('/')
    snm = parts[-2]
    pst = '/'.join(parts[-5:-2])
    return pst, snm, filename


class RadBxdsBase(unittest.TestCase):
    def check_input_exists(self, bxds_input):
        """ Make sure inputs exist as precondition for running test """
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

class TestFileFunctions(unittest.TestCase):
    """ Check that the file information reading is as-expected """
    def check_type_and_stream(self, pst:str, snm:str, radartype:str):
        filename = WAIS / 'orig/xlob' / pst / snm / 'bxds'
        restype = unfoc.get_radar_type(filename)
        resstream = unfoc.get_radar_stream(filename)
        assert restype[0] == radartype
        assert resstream == snm
    def test_type_mpol(self):
        """ Check a multipol transect """
        pst, snm, radartype = 'DOT/HVU0a/Y05a', 'RADnh5', 'MPOL'
        self.check_type_and_stream(pst, snm, radartype)

    def test_type_marfa(self):
        """ Check a multipol transect """
        pst, snm, radartype = 'DEV2/JKB2t/Y91a', 'RADnh5', 'MARFA'
        self.check_type_and_stream(pst, snm, radartype)
    def test_type_hicars(self):
        """ Check a multipol transect """
        pst, snm, radartype = 'AGAW/JKB2k/JEVANSb', 'RADnh3', 'HiCARS2'
        self.check_type_and_stream(pst, snm, radartype)

class TestParsers(RadBxdsBase):
    """ Test the functions and generators that parse the packet structure of the bxds file.
    Test these on a wide range of input stimulus """
    def run_compare_gen_parsers(self, bxds_input, **kwargs):
        self.check_input_exists(bxds_input)
        gen1 = unfoc.index_RADnhx_bxds_mmap_(bxds_input, **kwargs)
        gen2 = unfoc.index_RADnhx_bxds(bxds_input, **kwargs)

        for ii, (data1, data2) in enumerate(itertools.zip_longest(gen1, gen2)):
            self.assertEqual(data1, data2, msg="parse mismatch at record %d" % (ii,))

    def test_compare_gen_parsers(self):
        """ Check that indexers parse the bxds same """
        for ii, (pst, snm, bxds_input) in enumerate(list(read_testlist(TESTLIST))):
            with self.subTest(pst=pst, snm=snm):
                logging.debug("pst=%s snm=%s", pst, snm)
                self.run_compare_gen_parsers(bxds_input)
                if ii == 0:
                    self.run_compare_gen_parsers(bxds_input, full_header=True)


class TestGenCT(unittest.TestCase):
    """ Hit specific input cases for genct """
    def test_gzip(self):
        """ Make sure we can read a gzip file alongside a gzip file """
        for _, _, bxdsfile0 in read_testlist(TESTLIST):
            with self.subTest(file=bxdsfile0):
                bxdsfile = Path(bxdsfile0)
                ctfile = bxdsfile.with_name('ct.gz')
                if ctfile.exists():
                    for _ in unfoc.gen_ct(bxdsfile, full=False):
                        pass
                    for _ in unfoc.gen_ct(bxdsfile, full=True):
                        pass
                else:
                    # make a gzip and try it
                    ctfile = ctfile.with_name('ct')
                    with tempfile.TemporaryDirectory() as tempdir:
                        ptemp = Path(tempdir)
                        bxdsfile2 = ptemp / bxdsfile.name
                        bxdsfile2.touch()
                        gzip_text(ctfile, ptemp / 'ct.gz')
                        for _ in unfoc.gen_ct(bxdsfil2):
                            pass


    def test_nongzip(self):
        """ Make sure we can read a non-gzip file alongside a gzip file """
        for _, _, bxdsfile0 in read_testlist(TESTLIST):
            with self.subTest(file=bxdsfile0):
                bxdsfile = Path(bxdsfile0)
                ctfile = bxdsfile.with_name('ct')
                if ctfile.exists():
                    for _ in unfoc.gen_ct(bxdsfile):
                        pass
                else:
                    # make a gzip and try it
                    ctfile = ctfile.with_name('ct.gz')
                    with tempfile.TemporaryDirectory() as tempdir:
                        ptemp = Path(tempdir)
                        bxdsfile2 = ptemp / bxdsfile.name
                        bxdsfile2.touch()
                        self.zcat_text(ctfile, ptemp / 'ct')
                        for _ in unfoc.gen_ct(bxdsfile2, full=False):
                            pass
                        for _ in unfoc.gen_ct(bxdsfile2, full=True):
                            pass

    def zcat_text(self, file_compressed: Path, file_plain: Path):
        """ Read a gzipped file and write it to a non-gzipped file """
        with gzip.open(file_compressed, 'rt') as fin, \
             open(file_plain, 'wt') as fout:
            for line in fin:
                print(line, file=fout, end='')

    def gzip_text(self, file_plain: Path, file_compressed: Path):
        """ Compress text to gzip format """
        with open(file_plain, 'rt') as fin, \
             gzip.open(file_compressed, 'wt') as fout:
            for line in fin:
                print(line, file=fout)

class RadBxdsTestLoader:
    @classmethod
    def setUpClass(cls):
        """ we can reuse the rread object for all tests since it is readonly """
        channel = 1
        pst, snm, bxds_input = list(read_testlist(TESTLIST0))[0]
        cls.rread = unfoc.RadBxds(bxds_input, channel=channel)
        cls.expected_shape = (134676, 3437)

    @classmethod
    def tearDownClass(cls):
        del cls.rread

class IteratorTestComponent:
    """ Tests for classes that support iteration """
    def test_iteration(self):
        """  make sure we can run two iterators on the same array concurrently,
        independently """
        for _ in range(2):
            counter = 0
            for trace1, trace2 in zip(self.rread, self.rread):
                counter += 1
                self.assertTrue(isinstance(trace1, np.ndarray))
                self.assertTrue(isinstance(trace2, np.ndarray))
                np.testing.assert_array_equal(trace1, trace2)

        self.assertEqual(counter, len(self.rread))




class TestClass1(RadBxdsTestLoader, RadBxdsBase, IteratorTestComponent):
    """ Test the RadBxds class. """

    def test_indexing(self):
        """ Check that slicing works consistently with how numpy does it. """
        rread = self.rread
        ntraces = len(rread)
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
        """ Show that a slice comes out to the correct dimensions. """
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
        """ Slicing beyond the end of the radargram doesn't raise an error
         (This matches the behavior with lists and other sequences) """

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
        """ Slicing with a stride/step """
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



    def test_slicing_2d_a(self):
        """
        2d slicing with slices for both axes
        Test for issue #4, Support 2D ndarray slicing syntax in RadBxds.getitem
        Test multiple input files but no need to read the full file
        """
        # Check that behavior of this class is consistent with numpy
        arr = np.zeros((100, 3200), dtype=np.int16)
        traces0 = arr[0:10][:, 20:120]
        expected_shape = (10, 100)
        self.assertEqual(traces0.shape, expected_shape)

        traces0 = arr[0:10, 20:120]
        self.assertEqual(traces0.shape, expected_shape)

        # 1d way
        traces1 = self.rread[0:10][:, 20:120]
        self.assertIsNotNone(traces1)
        self.assertEqual(traces1.shape, expected_shape)

        # 2d way
        traces2 = self.rread[0:10, 20:120]
        self.assertIsNotNone(traces2)
        self.assertEqual(traces2.shape, expected_shape)

        # Check that they come out to the same value
        np.testing.assert_equal(traces1, traces2)


    def test_slicing_2d_b(self):
        """
        2d slicing with singleton first axis, slice second axis
        Test for issue #4, Support 2D ndarray slicing syntax in RadBxds.getitem
        Test multiple input files but no need to read the full file
        """
        # Check that behavior of this class is consistent with numpy
        arr = np.zeros((100, 3200), dtype=np.int16)
        traces0 = arr[5][20:120]
        expected_shape = (100,)
        self.assertEqual(traces0.shape, expected_shape)

        traces0 = arr[5, 20:120]
        self.assertEqual(traces0.shape, expected_shape)


        # 1d way
        traces1 = self.rread[5][20:120]
        self.assertIsNotNone(traces1)
        self.assertEqual(traces1.shape, expected_shape)

        # 2d way
        traces2 = self.rread[5, 20:120]
        self.assertIsNotNone(traces2)
        self.assertEqual(traces2.shape, expected_shape)
        #self.assertEqual(s[1], 100)

        # Check that they come out to the same value
        np.testing.assert_equal(traces1, traces2)



    def test_slicing_2d_c(self):
        """ 2d slicing with slice first axis, singleton 2nd axis
        Test for issue #4, Support 2D ndarray slicing syntax in RadBxds.getitem
        Test multiple input files but no need to read the full file
        """

        # Check that behavior of this class is consistent with numpy
        arr = np.zeros((100, 3200), dtype=np.int16)
        traces0 = arr[5:20][:, 5]
        expected_shape = (15,)
        self.assertEqual(traces0.shape, expected_shape)

        traces0 = arr[5:20, 5]
        self.assertEqual(traces0.shape, expected_shape)


        # 1d way
        traces1 = self.rread[5:20][:, 5]
        self.assertIsNotNone(traces1)
        self.assertEqual(traces1.shape, expected_shape)

        # 2d way
        traces2 = self.rread[5:20, 5]
        self.assertIsNotNone(traces2)
        self.assertEqual(traces2.shape, expected_shape)

        # Check that they come out to the same value
        np.testing.assert_equal(traces1, traces2)

    def test_slicing_2d_d(self):
        """ 2d slicing with singleton dimensions on both axes
        Test for issue #4, Support 2D ndarray slicing syntax in RadBxds.getitem
        Test multiple input files but no need to read the full file
        """


        # Check that behavior of this class is consistent with numpy
        arr = np.zeros((100, 3200), dtype=np.int16)
        traces0 = arr[29][5]
        expected_shape = (15, 1)
        self.assertEqual(type(traces0), np.int16)

        traces0 = arr[29, 5]
        self.assertEqual(type(traces0), np.int16)

        # 1d way
        traces1 = self.rread[29][5]
        self.assertEqual(type(traces1), np.int16)

        # 2d way
        traces2 = self.rread[29, 5]
        self.assertEqual(type(traces2), np.int16)

        # Check that they come out to the same value
        np.testing.assert_equal(traces1, traces2)


    def test_slicing_2d_stride(self):
        """ Test for issue #4, Support 2D ndarray slicing syntax in RadBxds.getitem
        Test multiple input files but no need to read the full file
        """
        # 1d way
        traces1 = self.rread[0:10][:, 20:120:2]
        s = traces1.shape
        self.assertEqual(s[0], 10)
        self.assertTrue(s[1] == 50)

        # 2d way
        traces2 = self.rread[0:10, 20:120:2]
        s = traces2.shape
        self.assertEqual(s[0], 10)
        self.assertTrue(s[1] == 50)

        # Check that they come out to the same value
        np.testing.assert_equal(traces1, traces2)




    def test_ct_index(self):
        for ii in range(len(self.rread)):
            ct1 = self.rread.ct(ii)
            self.assertEqual(len(ct1), 2)
        self.assertGreaterEqual(len(self.rread.cts_), len(self.rread))

    def test_ct_slice(self):
        for ii in range(0, len(self.rread), 10):
            ct1 = self.rread.ct(slice(ii, ii+10))
            self.assertLessEqual(len(ct1), 10)


class TestClassAttributes(RadBxdsTestLoader, unittest.TestCase):
    """ Test dimensional attributes on a dataset of known size """

    def test_size(self):
        expected_size = self.expected_shape[0]*self.expected_shape[1]
        self.assertEqual(self.rread.size, expected_size)

    def test_shape(self):
        self.assertEqual(self.rread.shape, self.expected_shape)

    def test_ndim(self):
        self.assertEqual(self.rread.ndim, 2)

    def test_dtype(self):
        dtype1 = self.rread.dtype
        if isinstance(dtype1, str):
            self.assertEqual(dtype1[1:], 'i2')

    def test_nbytes(self):
        itemsize = getattr(self, 'itemsize', 2)
        self.assertEqual(self.rread.nbytes, self.expected_shape[0]*self.expected_shape[1]*itemsize)


    def test_unimplemented(self):
        with self.assertRaises(AttributeError):
            self.rread.doesnotexist
        with self.assertRaises(NotImplementedError):
            self.rread.T

class Radjh1TestLoader:
    """ Component to load a RADjh1 test """
    @classmethod
    def setUpClass(cls):
        channel = 1
        testlist1 = cwd / 'test_lists' / 'tests_radjh1.txt'
        pst, snm, bxds_input = list(read_testlist(testlist1))[0]
        cls.rread = unfoc.RADjh1Bxds(bxds_input, channel=channel)
        cls.expected_shape = (509136, 3200)

    @classmethod
    def tearDownClass(cls):
        del cls.rread

class TestRADjh1Class(Radjh1TestLoader, TestClass1, IteratorTestComponent):
    """ Test general class methods on RADjh1 dataset """
    pass


class TestRADjh1ClassAttributes(Radjh1TestLoader, TestClassAttributes):
    """ Test dimensional attributes on a RADjh1 dataset """
    pass


class TestRadBxds(RadBxdsBase):
    """ Run tests on many different bxdses """
    def test_index_generator1(self):
        return self.run_index_generator(unfoc.index_RADnhx_bxds_mmap_)

    #@unittest.skip("Unneeded since we test equality later")
    #def test_index_generator2(self):
    #    return self.run_index_generator(read.index_RADnhx_bxds)

    def test_read(self):
        for pst, snm, bxds_input in list(read_testlist(TESTLIST)):
            for channel in (1, 2):
                with self.subTest(pst=pst, snm=snm, channel=channel):
                    self.check_input_exists(bxds_input)
                    trace_p = None
                    for ii, trace in enumerate(unfoc.read_RADnhx_gen(bxds_input, channel=channel)):
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
                    rread = unfoc.RadBxds(bxds_input, channel=channel)
                    self.assertGreater(len(rread), 1) # number of records

                    for ii, trace1 in enumerate(unfoc.read_RADnhx_gen(bxds_input, channel=channel)):
                        msg = "mismatch at rread[%d]" % (ii,)
                        trace2 = rread[ii]
                        self.assertEqual(trace1.data.shape, trace2.shape, msg=msg)
                        self.assertTrue(np.array_equal(trace1.data, trace2), msg=msg)
                        self.assertEqual(trace1.ct, rread.ct(ii), msg=msg)

                    self.assertEqual(len(rread), ii+1)

    def test_get_rad_stream(self):
        # Try detecting bxds.
        pass


class RadBxdsEx_2ch_Loader:
    """ Load a RadBxdsEx object with a two-channel transect. """
    @classmethod
    def setUpClass(cls):
        p1cs = unfoc.get_utig_channels('LoResInco1', radar='MARFA')[0]
        pst, snm = 'DEV2/JKB2t/Y91a', 'RADnh5'
        bxds_input = WAIS / 'orig/xlob' / pst / snm / 'bxds'
        cls.rread = unfoc.RadBxdsEx(bxds_input, channels=p1cs)

class RadBxdsExBase:
    """ Additional tests for the RadBxdsEx class """
    def test_ex_index(self):
        # single index
        for ii in range(len(self.rread)):
            arr1 = self.rread[ii]
            self.assertEqual(len(arr1.shape), 1)
            ct1 = self.rread.ct(ii)
            self.assertTrue(isinstance(ct1, tuple))

    def test_ex_slice(self):
        slicelen = 10
        len1 = len(self.rread)
        for ii in range(0, len(self.rread), slicelen):
            msg = "mismatch at %d of %d" % (ii, len1)
            arr1 = self.rread[ii:(ii+slicelen)]
            self.assertEqual(len(arr1.shape), 2)
            self.assertLessEqual(arr1.shape[0], slicelen)
            ct1 = self.rread.ct(slice(ii,ii+slicelen) )
            self.assertEqual(len(ct1), arr1.shape[0], msg=msg)

class TestClassEx2ch(RadBxdsEx_2ch_Loader, RadBxdsExBase, IteratorTestComponent, RadBxdsBase):
    """ Test RadBxdsBase with two channels """
    pass


class RadBxdsEx_1ch_Loader:
    """ Load a RadBxdsEx object with a one-channel transect. """
    @classmethod
    def setUpClass(cls):
        p1cs = unfoc.get_utig_channels('LoResInco5', radar='MARFA')[0]
        pst, snm = 'DEV2/JKB2t/Y91a', 'RADnh5'
        bxds_input = WAIS / 'orig/xlob' / pst / snm / 'bxds'
        cls.rread = unfoc.RadBxdsEx(bxds_input, channels=p1cs, dtype=np.double)


class TestClassEx1ch(RadBxdsEx_1ch_Loader, RadBxdsExBase, IteratorTestComponent, RadBxdsBase):
    """ Test RadBxdsBase with one channel.  Same tests but different channel spec """
    pass

class RadBxdsExTestLoader:
    @classmethod
    def setUpClass(cls):
        p1cs = unfoc.get_utig_channels('LoResInco5', radar='MARFA')[0]
        pst, snm = 'DEV2/JKB2t/Y91a', 'RADnh5'
        bxds_input = WAIS / 'orig/xlob' / pst / snm / 'bxds'
        cls.rread = unfoc.RadBxdsEx(bxds_input, channels=p1cs, dtype=np.double)
        cls.expected_shape = (55916, 3200)
        cls.itemsize = 8

    @classmethod
    def tearDownClass(cls):
        del cls.rread

class TestRadBxdsExClassAttributes(RadBxdsExTestLoader, TestClassAttributes):
    """ Test dimensional attributes on a RadBxdsEx object """
    pass


class NoiseTestLoader:
    """ Load a RadBxdsEx object configured for noise suppresion """
    @classmethod
    def setUpClass(cls):
        p1cs = unfoc.get_utig_channels('LoResInco2', radar='MARFA')[0]
        p1cs = unfoc.enable_burstnoise(p1cs)
        pst, snm = 'DEV2/JKB2t/Y91a', 'RADnh5'
        bxds_input = WAIS / 'orig/xlob' / pst / snm / 'bxds'
        cls.rread = unfoc.RadBxdsEx(bxds_input, channels=p1cs)

class TestClassExNoise(NoiseTestLoader, RadBxdsExBase, RadBxdsBase):
    """ Test RadBxdsBase with the high sum gain channel and enable denoising """
    pass

class TestOneMeter(unittest.TestCase):
    def test_1m_a(self):
        # pst = 'NIS4/IBH0e/X84b'
        pst = 'D2DG/IBH0e/X30a'
        bxds_path = WAIS / 'targ/xtra/KRT2/FOC/Best_Versions/S2_FIL' / pst
        for chan in (1, 2, 5, 6, 7, 8):
            with self.subTest(chan=chan):
                bxds_input = os.path.join(bxds_path, 'bxds{:d}.i'.format(chan))
                prevct = None
                for t in unfoc.read_1m_gen(bxds_input, chan):
                    if prevct is not None:
                        self.assertLess(prevct.seq, t.ct.seq)
                    prevct = t.ct


def main():
    #logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)
    unittest.main()


if __name__ == "__main__":
    main()
