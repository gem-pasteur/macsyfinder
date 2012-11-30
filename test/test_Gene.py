'''
Created on Nov 30, 2012

@author: bneron
'''

import sys
import os
TXSSCAN_HOME = os.path.abspath( '..')
if not TXSSCAN_HOME in sys.path: 
    sys.path.append( os.path.abspath( '..') )

import unittest

from txsscanlib.secretion import Gene
from txsscanlib.config import Config

class Test(unittest.TestCase):


    def setUp(self):
        self.cfg = Config( hmmer_exe = "",
                           e_value_res = 1,
                           e_value_sel = 0.5,
                           def_dir = "../data/DEF",
                           res_search_dir = ".",
                           res_search_suffix = "",
                           profile_dir = "../data/profiles",
                           profile_suffix = ".fasta-aln_edit.hmm",
                           res_extract_suffix = "",
                           log_level = 30
                           )


    def test_add_homolog(self):
        gene = Gene('sctJ_FLG' , self.cfg)
        homolog = Gene( 'sctJ', self.cfg)
        gene.add_homolog( homolog )
        self.assertEqual(len( gene.homologs), 1)
        self.assertEqual(gene.homologs[0], homolog)
    
    def test_get_homologs(self):
        gene = Gene('sctJ_FLG' , self.cfg)
        homolog_1 = Gene( 'sctJ', self.cfg)
        gene.add_homolog( homolog_1 )
        homolog_2 = Gene( 'sctN_FLG', self.cfg)
        gene.add_homolog( homolog_2 )
        self.assertEqual( gene.get_homologs(), [homolog_1, homolog_2] )
                         
                         
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()