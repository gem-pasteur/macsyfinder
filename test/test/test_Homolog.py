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
from txsscanlib.gene import Homolog
from txsscanlib.gene import Gene
from txsscanlib.system import System
from txsscanlib.config import Config


class Test(unittest.TestCase):


    def setUp(self):
        self.cfg = Config( sequence_db = ".",
                           db_type = "gembase",
                           hmmer_exe = "",
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


    def test_gene_ref(self):
        system = System("T2SS", self.cfg)
        gene_ref = Gene(self.cfg, 'sctJ_FLG', system)
        gene = Gene(self.cfg, 'sctJ', system)
        homolog_1 = Homolog(gene, gene_ref)
        self.assertEqual( homolog_1.gene_ref , gene_ref)
 
    def test_is_aligned(self):
        system = System("T2SS", self.cfg)
        gene_ref = Gene(self.cfg, 'sctJ_FLG', system)
        gene = Gene(self.cfg, 'sctJ', system)
        homolog = Homolog( gene, gene_ref)
        self.assertFalse( homolog.is_aligned() )
        homolog = Homolog(gene, gene_ref, aligned = True  )
        self.assertTrue( homolog.is_aligned() )

    def test_delegation(self):
        system = System("T2SS", self.cfg)
        gene_ref = Gene(self.cfg, 'sctJ_FLG', system)
        gene = Gene(self.cfg, 'sctJ', system)
        homolog = Homolog( gene, gene_ref)
        self.assertEqual( homolog.system , system )

if __name__ == "__main__":
    unittest.main()