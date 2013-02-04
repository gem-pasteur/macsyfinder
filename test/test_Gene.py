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

    def test_add_homolog(self):
        system_foo = System( "foo", self.cfg)
        system_bar = System( "bar", self.cfg)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo)
        gene_ref = Gene(self.cfg, 'sctJ', system_bar)
        homolog = Homolog(self.cfg, gene, gene_ref)
        gene.add_homolog( homolog )
        self.assertEqual(len( gene.homologs), 1)
        self.assertEqual(gene.homologs[0], homolog)
    
    
    def test_get_homologs(self):
        system_foo = System( "foo", self.cfg)
        system_bar = System( "bar", self.cfg)
        gene = Gene(self.cfg, 'sctN', system_foo)
        sctJ_FLG = Gene(self.cfg, 'sctJ_FLG', system_foo)
        sctJ = Gene(self.cfg, 'sctJ', system_bar)
        homolog_1 = Homolog( sctJ_FLG, gene)
        gene.add_homolog( homolog_1 )
        homolog_2 = Homolog(sctJ, gene)
        gene.add_homolog( homolog_2 )
        self.assertEqual( gene.get_homologs(), [homolog_1, homolog_2] )
    
        
    def test_system(self):
        """
        test getter/setter for system property
        """
        system_foo = System( "foo", self.cfg)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo)
        self.assertEqual(gene.system, system_foo)
    
    
    def test_loner(self):
        """
        test getter for loner property
        """
        system_foo = System( "foo", self.cfg)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo)
        self.assertFalse(gene.loner)
        gene = Gene(self.cfg, 'sctJ', system_foo, loner = True)
        self.assertTrue(gene.loner)

   
    def test_exchangeable(self):
        """
        test getter for exchangeable property
        """
        system_foo = System( "foo", self.cfg)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo)
        self.assertFalse(gene.exchangeable)
        gene = Gene(self.cfg, 'sctJ', system_foo, exchangeable = True)
        self.assertTrue(gene.exchangeable)
        
        
    def test_str(self):
        """
        """
        system_foo = System( "foo", self.cfg)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo)
        
        system_bar = System( "bar", self.cfg)
        gene_homolog = Gene(self.cfg, 'sctJ', system_bar)
        homolog = Homolog( gene_homolog, gene, self.cfg)
        gene.add_homolog( homolog )
        s = """name : sctJ_FLG
    homologs: sctJ"""
        self.assertEqual( str(gene) , s )


if __name__ == "__main__":
    unittest.main()