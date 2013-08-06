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
from txsscanlib.gene import Gene
from txsscanlib.gene import Homolog
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
                           res_search_dir = "/tmp",
                           res_search_suffix = "",
                           profile_dir = "../data/profiles",
                           profile_suffix = ".fasta-aln_edit.hmm",
                           res_extract_suffix = "",
                           log_level = 30,
                           log_file = '/dev/null'
                           )
    def tearDown(self):
        shutil.rmtree(self.cfg.working_dir)

    def test_add_homolog(self):
        system_foo = System(self.cfg, "foo", 10)
        system_bar = System(self.cfg, "bar", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo)
        gene_ref = Gene(self.cfg, 'sctJ', system_bar)
        homolog = Homolog(self.cfg, gene, gene_ref)
        gene.add_homolog( homolog )
        self.assertEqual(len( gene.homologs), 1)
        self.assertEqual(gene.homologs[0], homolog)
    
    
    def test_get_homologs(self):
        system_foo = System(self.cfg, "foo", 10)
        system_bar = System(self.cfg, "bar", 10)
        gene = Gene(self.cfg, 'sctN', system_foo)
        sctJ_FLG = Gene(self.cfg, 'sctJ_FLG', system_foo)
        sctJ = Gene(self.cfg, 'sctJ', system_bar)
        homolog_1 = Homolog(sctJ_FLG, gene)
        gene.add_homolog(homolog_1)
        homolog_2 = Homolog(sctJ, gene)
        gene.add_homolog(homolog_2)
        self.assertEqual(gene.get_homologs(), [homolog_1, homolog_2] )
    
        
    def test_system(self):
        """
        test getter/setter for system property
        """
        system_foo = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo)
        self.assertEqual(gene.system, system_foo)
    
    
    def test_loner(self):
        """
        test getter for loner property
        """
        system_foo = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo)
        self.assertFalse(gene.loner)
        gene = Gene(self.cfg, 'sctJ', system_foo, loner = True)
        self.assertTrue(gene.loner)


    def test_exchangeable(self):
        """
        test getter for exchangeable property
        """
        system_foo = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo)
        self.assertFalse(gene.exchangeable)
        gene = Gene(self.cfg, 'sctJ', system_foo, exchangeable = True)
        self.assertTrue(gene.exchangeable)
 
    def test_multi_system(self):
        """
        test getter for multi_system property
        """
        system_foo = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo)
        self.assertFalse(gene.multi_system)
        gene = Gene(self.cfg, 'sctJ', system_foo, multi_system = True)
        self.assertTrue(gene.multi_system)


    def test_inter_gene_max_space(self):
        """
        test getter for inter_gene_max_space property
        """
        system_inter_gene_max_space = 40
        gene_inter_gene_max_space = 50
        system_foo = System(self.cfg, "foo", system_inter_gene_max_space)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo)
        self.assertEqual(gene.inter_gene_max_space, system_inter_gene_max_space)
        gene = Gene(self.cfg, 'sctJ', system_foo, inter_gene_max_space = gene_inter_gene_max_space)
        self.assertEqual(gene.inter_gene_max_space, gene_inter_gene_max_space)


    def test_str(self):
        """
        """
        system_foo = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo)
        system_bar = System(self.cfg, "bar", 20)
        gene_homolog = Gene(self.cfg, 'sctJ', system_bar)
        homolog = Homolog( gene_homolog, gene, self.cfg)
        gene.add_homolog( homolog )
        s = """name : sctJ_FLG
    homologs: sctJ"""
        self.assertEqual( str(gene) , s )


if __name__ == "__main__":
    unittest.main()
