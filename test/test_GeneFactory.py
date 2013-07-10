# -*- coding: utf-8 -*-

#===============================================================================
# Created on Jan 14, 2012
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
from txsscanlib.gene import gene_factory
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
        gene_factory._genes_bank = {}
        shutil.rmtree(self.cfg.working_dir)

    def test_get_gene(self):
        system_foo = System("foo", self.cfg)
        gene_name = 'sctJ_FLG'
        gene = gene_factory.get_gene(self.cfg, gene_name, system_foo)
        self.assertTrue(isinstance(gene, Gene))
        self.assertEqual(gene.name, gene_name)

    def test_get_loner_gene(self):
        system_foo = System("foo", self.cfg)
        gene = gene_factory.get_gene(self.cfg, 'sctJ_FLG', system_foo)
        self.assertFalse(gene.loner)
        gene = gene_factory.get_gene(self.cfg, 'sctJ', system_foo, loner = True)
        self.assertTrue(gene.loner)

    def test_get_exchangeable_gene(self):
        system_foo = System("foo", self.cfg)
        gene = gene_factory.get_gene(self.cfg, 'sctJ_FLG', system_foo)
        self.assertFalse(gene.exchangeable)
        gene = gene_factory.get_gene(self.cfg, 'sctJ', system_foo, exchangeable = True)
        self.assertTrue(gene.exchangeable)

    def test_get_multi_system_gene(self):
        system_foo = System("foo", self.cfg)
        gene = gene_factory.get_gene(self.cfg, 'sctJ_FLG', system_foo)
        self.assertFalse(gene.multi_system)
        gene = gene_factory.get_gene(self.cfg, 'sctJ', system_foo, multi_system = True)
        self.assertTrue(gene.multi_system)

    def test_inter_gene_max_space(self):
        """
        test getter for nter_gene_max_space property
        """
        system_foo = System("foo", self.cfg)
        system_inter_gene_max_space = 40
        gene_inter_gene_max_space = 50
        system_foo.inter_gene_max_space = system_inter_gene_max_space
        gene = gene_factory.get_gene(self.cfg, 'sctJ_FLG', system_foo)
        self.assertEqual(gene.inter_gene_max_space, system_inter_gene_max_space)
        gene = gene_factory.get_gene(self.cfg, 'sctJ', system_foo, gene_inter_gene_max_space)
        self.assertEqual(gene.inter_gene_max_space, system_inter_gene_max_space)


    def test_get_uniq_object(self):
        system_foo = System( "foo", self.cfg)
        gene1 = gene_factory.get_gene(self.cfg, 'sctJ_FLG', system_foo)
        gene2 = gene_factory.get_gene(self.cfg, 'sctJ_FLG', system_foo)
        self.assertEqual( gene1, gene2 )



if __name__ == "__main__":
    unittest.main()