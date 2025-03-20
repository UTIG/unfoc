#!/usr/bin/env python3

from pathlib import Path
import json
import os
import unittest
import sys
import logging
from pathlib import Path
from collections import Counter
import tempfile

cwd = Path(__file__).parent
unfoc_path = cwd / '..' / 'src'
sys.path.insert(1, str(unfoc_path.absolute()))

import unfoc

WAIS = Path(os.getenv('WAIS', '/disk/kea/WAIS'))

class TestMultichannelRect(unittest.TestCase):

    def psts(self):
        """ List of nonrectangular PSTs """
        nrpsts = ['TOT3/JKB2s/X15a', # ICP9, orphan at beginning
                  'ASB/JKB2s/GL0360b', # ICP9, orphan at beginning
                  'THW2/UBH0c/X245a'] # ASE4

    def test1(self):
        bxds = WAIS / 'orig/xlob/THW2/UBH0c/X245a/RADnh5/bxds'

        s = unfoc.get_radar_stream(str(bxds))
        assert s == 'RADnh5'

        x = unfoc.sync_radar_start(str(bxds), stream=s)
        rec0, fpos0, rseq0, nchan = x

        # for this one, the first sample is 
        assert rec0 == 0
        assert fpos0 == 0
        assert nchan == 2



    def test_nonrect_tot3_jkb2s_x15a(self):
        """ Check with a nonrectangular PST """
        pst = 'TOT3/JKB2s/X15a'        # ICP9, orphan at beginning
        snm = 'RADnh5'
        bxds = WAIS / 'orig/xlob' / pst / snm / 'bxds'

        s = unfoc.get_radar_stream(str(bxds))
        assert s == snm


        x = unfoc.sync_radar_start(str(bxds), stream=s)
        rec0, fpos0, rseq0, nchan = x
        # for this one, the first sample is not the first
        assert rec0 == 1
        assert fpos0 == 12870

        info = unfoc.radar_index_summary(bxds)
        assert info['incomplete_records']

        for item in unfoc.index_RADnhx_bxds(bxds, stream=s, filepos=fpos0):
            pass


    def test_nonrect_sba_jkb2s_x38a(self):
        """ This is a weird one with a lot of mismatched data
        SBA/JKB2s/X38a  {
        75254016: 1,
        75268288: 1,
        75268352: 1,
        75268512: 1,
        75284160: 1,
        75284192: 1,
        75284224: 1,
        75284256: 1,
        75284288: 1,
        75284320: 1,
        75284352: 1,
        75284512: 1,
        75284544: 1,
        75284576: 1,
        75284608: 1,
        75386880: 1})

        """

        pst = 'SBA/JKB2s/X38a'
        snm = 'RADnh5'

        bxds = WAIS / 'orig/xlob' / pst / snm / 'bxds'

        s = unfoc.get_radar_stream(str(bxds))
        assert s == snm

        x = unfoc.sync_radar_start(str(bxds), stream=s)
        rec0, fpos0, rseq0, nchan = x
        # for this one, the first sample is not the first
        assert rec0 == 0, x
        assert fpos0 == 0

        info = unfoc.radar_index_summary(bxds, s)
        assert info['incomplete_records'], "These shouldn't have the same number of records"


    #@unittest.skip("Takes a long time but this works")
    def test_index_all_mismatched_radnh5(self):
        psts = []
        with open('mismatched.txt', 'rt') as fin:
            for line in fin:
                pst, rest = line.rstrip().split(';', 2)
                psts.append(pst)

        p_testout = Path(__file__).parent / 'mismatchout'
        p_testout.mkdir(parents=True, exist_ok=True)
        for pst in psts[0:10]:
            logging.debug("pst=%s", pst)

            snm = 'RADnh5'
            bxds = WAIS / 'orig/xlob' / pst / snm / 'bxds'
            with self.subTest(pst=pst):
                info = unfoc.radar_index_summary(bxds)
                outfile = p_testout / (pst.replace('/', '_') + '.json')
                with outfile.open('wt') as fhjson:
                    json.dump(info, fhjson, indent="\t")
                assert info['incomplete_records'], "Should have mismatched records"


class TestSync(unittest.TestCase):
    """ Test the sync_radar_start function """

    def run_sync(self, pst:str, snm:str):
        """ Exercise the sync_radar_start function """
        bxds = WAIS / 'orig/xlob' / pst / snm / 'bxds'
        s = unfoc.get_radar_stream(str(bxds))
        assert s == snm

        x1 = unfoc.sync_radar_start(str(bxds), stream=s)
        x2 = unfoc.sync_radar_start(str(bxds))
        assert x1 == x2
        return x1

    def test_radnh3_1channel(self):
        """ Check that sync_radar_start works with RADnh3, but only one set of channels"""
        rec0, fpos0, rseq0, nchan = self.run_sync(pst='CLH/JKB2h/Y32a', snm='RADnh3')
        assert nchan == 1
        assert rec0 == 0
        assert fpos0 == 0

    def test_radnh3_2channel_mismatched(self):
        """ Test on a 2-channel version with a mismatch """
        rec0, fpos0, rseq0, nchan = self.run_sync(pst='GOG3/JKB2j/BWN01a', snm='RADnh3')
        assert nchan == 2
        assert rec0 > 0
        assert fpos0 > 0

    # Can't seem to find any occurrences that don't have a mismatch?
    def test_radnh3_2channel_matched(self):
        """ Test on a 2-channel version without a mismatch """
        rec0, fpos0, rseq0, nchan = self.run_sync(pst='NAQLK/JKB2j/X24a', snm='RADnh3')
        assert nchan == 2
        assert rec0 == 0
        assert fpos0 == 0




class TestBxds(unittest.TestCase):
    """ Make sure behavior of bxds class handles valid/invalid properly """

    def test_radnh5(self):
        """ First make sure it has different lengths if valid_only is false, then same length (the lesser) for all channels if valid_only is True """

        pst = 'TOT3/JKB2s/X15a'        # ICP9, orphan at beginning
        snm = 'RADnh5'
        bxds = WAIS / 'orig/xlob' / pst / snm / 'bxds'

        # In this part of the test case, channel 1 should have 83704 records,
        # and channel 3 should have 83703 records
        validonly = False
        rada0 = unfoc.RadBxds(str(bxds), channel=1, validonly=validonly)
        rada1 = unfoc.RadBxds(str(bxds), channel=3, validonly=validonly)

        assert len(rada0) == len(rada1) + 1 , \
               "validonly=%r filename=%s channel_1_len=%d channel_3_len=%d" % (validonly, str(bxds), len(rada0), len(rada1))

        # But then after running specifying valid only, it should work ok and be min length
        validonly = True
        radb0 = unfoc.RadBxds(str(bxds), channel=1, validonly=validonly)
        radb1 = unfoc.RadBxds(str(bxds), channel=3, validonly=validonly)

        assert len(radb0) == len(rada1), \
           "validonly=%r filename=%s channel_1_len=%d channel_3_len=%d" % (validonly, str(bxds), len(radb0), len(rada1))
        assert len(radb1) == len(rada1), \
           "validonly=%r filename=%s channel_1_len=%d channel_3_len=%d" % (validonly, str(bxds), len(radb1), len(rada1))

    def test_radnh3_mismatched1(self):
        """ Mismatched before and after """
        pst = 'GOG3/JKB2j/BWN01a'
        snm = 'RADnh3'
        bxds = WAIS / 'orig/xlob' / pst / snm / 'bxds'


        validonly = False
        rada0 = unfoc.RadBxds(str(bxds), channel=1, validonly=validonly)
        rada1 = unfoc.RadBxds(str(bxds), channel=3, validonly=validonly)

        # This is in valid_records
        assert len(rada0) == 77827
        assert len(rada1) == 77827

        del rada0
        del rada1

        # But then after running specifying valid only, one at the beginning
        # and at the end should get truncated off
        validonly = True
        radb0 = unfoc.RadBxds(str(bxds), channel=1, validonly=validonly)
        radb1 = unfoc.RadBxds(str(bxds), channel=3, validonly=validonly)

        assert len(radb0) == 77826
        assert len(radb1) == 77826



    def test_radnh3_a(self):
        """ First make sure it has different lengths if valid_only is false, then same length (the lesser) for all channels if valid_only is True """

        pst = 'CLH/JKB2h/Y32a'        # One channel
        snm = 'RADnh3'
        bxds = WAIS / 'orig/xlob' / pst / snm / 'bxds'

        # This should do the same thing with validonly == 1 because it doesn't have a rada1

        # In this part of the test case, channel 1 should have 83704 records,
        # and channel 3 should have 83703 records
        validonly = False
        rada0 = unfoc.RadBxds(str(bxds), channel=1, validonly=validonly)
        rada1 = unfoc.RadBxds(str(bxds), channel=3, validonly=validonly)

        assert len(rada1) == 0, "Shouldn't be a second set of channels for this data"

        # But then after running specifying valid only, it should work ok and be min length
        validonly = True
        radb0 = unfoc.RadBxds(str(bxds), channel=1, validonly=validonly)
        radb1 = unfoc.RadBxds(str(bxds), channel=3, validonly=validonly)

        assert len(radb0) == len(rada0)
        assert len(radb1) == 0

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, stream=sys.stdout)
    unittest.main()
