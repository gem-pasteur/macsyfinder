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
import shutil
from txsscanlib.report import HMMReport, UnOrderedHMMReport, OrderedHMMReport, Hit
from txsscanlib.gene import Gene
from txsscanlib.system import System
from txsscanlib.config import Config

class Test(unittest.TestCase):

    _working_dir = "./working_dir"
    _data_dir = "./datatest/res_search" 
    
    def setUp(self):
        if os.path.exists(self._working_dir):
            shutil.rmtree(self._working_dir)
        os.mkdir(self._working_dir)
        
        self.cfg = Config( hmmer_exe = "hmmsearch",
                           sequence_db = "./datatest/prru_psae.001.c01.fasta",
                           ordered_db = True,
                           e_value_res = 1,
                           i_evalue_sel = 0.5,
                           def_dir = "../data/DEF",
                           res_search_dir = self._working_dir,
                           res_search_suffix = ".search_hmm.out",
                           profile_dir = "../data/profiles",
                           profile_suffix = ".fasta-aln_edit.hmm",
                           res_extract_suffix = "",
                           log_level = 30
                           )

    def tear_down(self):
        shutil.rmtree(self._working_dir)

    def test_HMMReport(self):
        gene_name = "gspD"
        gene = Gene("gspD", self.cfg)
        shutil.copy( os.path.join(self._data_dir, gene_name + self.cfg.res_search_suffix), self._working_dir)
        report_path = os.path.join(self._working_dir, gene_name + self.cfg.res_search_suffix)
        self.assertRaises(TypeError, HMMReport, gene, report_path, self.cfg)

    def test_OrderedHMMReport_extract(self):
        gene_name = "gspD"
        gene = Gene("gspD", self.cfg)
        system = System("T2SS", self.cfg)
        gene.system = system
        shutil.copy(os.path.join(self._data_dir, gene_name + self.cfg.res_search_suffix), self._working_dir)
        report_path = os.path.join(self._working_dir, gene_name + self.cfg.res_search_suffix)
        report = OrderedHMMReport(gene, report_path, self.cfg)
        report.extract()
        self.assertEqual(len(report.hits), 5)
        
        hits=[ Hit( gene, system, "PSAE001c01_006940", "PSAE001c01", 694 , float(1.2e-234), float(779.2), float(1.000000)),
               Hit( gene, system, "PSAE001c01_013980", "PSAE001c01", 1398, float(3.7e-76), float(255.8), float(1.000000)),
               Hit( gene, system, "PSAE001c01_017350", "PSAE001c01", 1735, float(3.2e-27), float(94.2), float(0.500000)),
               Hit( gene, system, "PSAE001c01_018920", "PSAE001c01", 1892, float(6.1e-183), float(608.4), float(1.000000)),
               Hit( gene, system, "PSAE001c01_031420", "PSAE001c01", 3142, float(1.8e-210), float(699.3), float(1.000000))
        ]
        self.assertListEqual(hits, report.hits)
            
        
        
if __name__ == "__main__":
    unittest.main()