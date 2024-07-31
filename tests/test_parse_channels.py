#!/usr/bin/env python3


import unittest
import sys
import logging
from pathlib import Path

cwd = Path(__file__).parent
unfoc_path = cwd / '..' / 'src'
sys.path.insert(1, str(unfoc_path.absolute()))

import unfoc

class TestParseChannels(unittest.TestCase):
    def test_invalid_values(self):
        with self.assertRaises(ValueError):
            unfoc.parse_channels("1.2")
        with self.assertRaises(ValueError):
            unfoc.parse_channels("[1,2,3,4]")

    def test_multi1(self):
        x = unfoc.parse_channels("[1,2,3,4,5;6,7,8.0,9,10]")
        self.assertEqual(len(x), 2)
        self.assertEqual(x[1].chanout, 6)

    def test_str1(self):
        x = unfoc.parse_channels("[1,2,3,4,5]")
        self.assertEqual(len(x), 1)
        self.assertEqual(x[0].chanout, 1)

    def test_many(self):
        tests = "[1,1,1,3,1] [2,2,1,4,1] [5,1,1,0,0;7,3,1,0,0] [6,2,1,0,0;8,4,1,0,0]".split(" ")

        for s in tests:
            x = unfoc.parse_channels(s)
            self.assertGreaterEqual(len(x), 1)

    def test_deprecated(self):
        for s in '1 2 3 4 5 6 7 8'.split():
            with self.assertLogs(level=logging.WARNING) as cm:
                x = unfoc.parse_channels(s)
                self.assertEqual(len(x), 1)
                self.assertEqual(x[0].chanout, int(s))

    def test_utig_channels(self):
        p1cs = unfoc.get_utig_channels('LoResInco5')
        self.assertEqual(len(p1cs), 1)

        p1cs = unfoc.get_utig_channels('LoResInco1,LoResInco2,LoResInco3')
        self.assertEqual(len(p1cs), 3)

    def test_utig_channels_invalid(self):
        p1cs = unfoc.get_utig_channels('LoResInco5,LoResInco88')
        self.assertEqual(len(p1cs), 1)

    def test_utig_channels_multipol(self):
        # In the case of multipol, we might ask for some extended channels that we don't have
        p1cs = unfoc.get_utig_channels('LoResInco1,LoResInco3,LoResInco6,LoResInco9', input_channels=[1,2,3,4])
        self.assertEqual(len(p1cs), 2, msg=str(p1cs))



if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, stream=sys.stdout)
    unittest.main()
