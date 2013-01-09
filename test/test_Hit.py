'''
Created on Nov 30, 2012

@author: bneron
'''

import sys
import os

TXSSCAN_HOME = os.path.abspath('..')
if not TXSSCAN_HOME in sys.path: 
    sys.path.append(os.path.abspath('..') )

import unittest

from txsscanlib.report import Hit
from txsscanlib.config import Config
from txsscanlib.gene import Gene
from txsscanlib.system import System


class Test(unittest.TestCase):


    def setUp(self):
        self.cfg = Config( hmmer_exe = "",
                           e_value_res = 1,
                           i_evalue_sel = 0.5,
                           def_dir = "../data/DEF",
                           res_search_dir = ".",
                           res_search_suffix = "",
                           profile_dir = "../data/profiles",
                           profile_suffix = ".fasta-aln_edit.hmm",
                           res_extract_suffix = "",
                           log_level = 30
                           )

    def test_cmp(self):
        gene_name = "gspD"
        gene = Gene("gspD", self.cfg)
        system = System("T2SS", self.cfg)
        gene.system = system
        h0 = Hit( gene, system, "PSAE001c01_006940", "PSAE001c01", 694 , float(1.2e-234), float(779.2), float(1.000000))
        h1 = Hit( gene, system, "PSAE001c01_013980", "PSAE001c01", 1398, float(3.7e-76), float(255.8), float(1.000000))
        self.assertGreater(h1, h0)
        self.assertLess(h0, h1)

    def test_eq(self):
        gene_name = "gspD"
        gene = Gene("gspD", self.cfg)
        system = System("T2SS", self.cfg)
        gene.system = system
        h0 = Hit( gene, system, "PSAE001c01_006940", "PSAE001c01", 694 , float(1.2e-234), float(779.2), float(1.000000))
        h1 = Hit( gene, system, "PSAE001c01_006940", "PSAE001c01", 694 , float(1.2e-234), float(779.2), float(1.000000))
        self.assertEqual(h0, h1)
        
        
if __name__ == "__main__":
    unittest.main()