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

from txsscanlib.gene import Homolog
from txsscanlib.gene import Gene
from txsscanlib.config import Config


class Test(unittest.TestCase):


    def setUp(self):
        self.cfg = Config( sequence_db = ".",
                           hmmer_exe = "",
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



    def test_gene_ref(self):
        gene = Gene('sctJ_FLG' , self.cfg)
        homolog_1 = Homolog( 'sctJ', self.cfg, gene)
        self.assertEqual( homolog_1.gene_ref , gene)
 
    def test_is_aligned(self):
        gene = Gene('sctJ_FLG' , self.cfg)
        homolog = Homolog( 'sctJ', self.cfg, gene )
        self.assertFalse( homolog.is_aligned() )
        homolog = Homolog( 'sctJ', self.cfg, gene, aligned = True  )
        self.assertTrue( homolog.is_aligned() )


if __name__ == "__main__":
    unittest.main()