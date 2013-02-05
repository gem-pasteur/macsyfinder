# -*- coding: utf-8 -*-

#===============================================================================
# Created on Nov 30, 2012
# 
# @author: bneron
# @contact: user_email
# @organization: organization_name
# @license: license
#===============================================================================


import sys
import os

TXSSCAN_HOME = os.path.abspath('..')
if not TXSSCAN_HOME in sys.path: 
    sys.path.append(os.path.abspath('..') )

import unittest
import shutil
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
                           log_level = 30,
                           log_file = '/dev/null'
                           )
  
    def tearDown(self):
        shutil.rmtree(self.cfg.working_dir)

    def test_cmp(self):
        system = System("T2SS", self.cfg)
        gene_name = "gspD"
        gene = Gene(self.cfg, "gspD",system)
        h0 = Hit( gene, system, "PSAE001c01_006940", "PSAE001c01", 694 , float(1.2e-234), float(779.2), float(1.000000))
        h1 = Hit( gene, system, "PSAE001c01_013980", "PSAE001c01", 1398, float(3.7e-76), float(255.8), float(1.000000))
        self.assertGreater(h1, h0)
        self.assertLess(h0, h1)

    def test_eq(self):
        system = System("T2SS", self.cfg)
        gene_name = "gspD"
        gene = Gene(self.cfg, "gspD", system)
        h0 = Hit( gene, system, "PSAE001c01_006940", "PSAE001c01", 694 , float(1.2e-234), float(779.2), float(1.000000))
        h1 = Hit( gene, system, "PSAE001c01_006940", "PSAE001c01", 694 , float(1.2e-234), float(779.2), float(1.000000))
        self.assertEqual(h0, h1)
        
    def test_str(self):
        system = System("T2SS", self.cfg)
        gene_name = "gspD"
        gene = Gene(self.cfg, "gspD", system)
        hit_prop={'id' : "PSAE001c01_006940",
                  'replicon_name' : "PSAE001c01",
                  'position' : 694,
                  'i_eval' : float(1.2e-234),
                  'score' : float(779.2),
                  'coverage' : float(1.0),
                  'gene_name' : gene.name,
                  'system_name' : system.name 
                  }
        
        hit = Hit( gene, system, hit_prop['id'], hit_prop['replicon_name'], hit_prop['position'] , hit_prop['i_eval'], hit_prop['score'], hit_prop['coverage'])
        s = "%(id)s\t%(replicon_name)s\t%(position)d\t%(gene_name)s\t%(system_name)s\t%(i_eval)s\t%(score)s\t%(coverage)f\n" % hit_prop
        self.assertEqual(s,str(hit))
        
if __name__ == "__main__":
    unittest.main()