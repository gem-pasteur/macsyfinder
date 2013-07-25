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
from txsscanlib.gene import gene_bank
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
        gene_bank._genes_bank = {}
        shutil.rmtree(self.cfg.working_dir)
    

    def test_add_get_gene(self):
        gene_name = 'sctJ_FLG'
        self.assertRaises(KeyError, gene_bank.__getitem__, gene_name)
        system_foo = System( "foo", self.cfg, 10)
        gene = Gene(self.cfg, gene_name, system_foo)
        gene_bank.add_gene(gene)
        gene_from_bank = gene_bank[gene_name]
        self.assertTrue(isinstance(gene_from_bank, Gene))
        self.assertEqual(gene_from_bank, gene)

    def test_contains(self):
        system_foo = System(self.cfg, "foo", 10)
        gene_in = Gene(self.cfg, 'sctJ_FLG', system_foo)
        gene_bank.add_gene(gene_in)
        self.assertIn(gene_in, gene_bank)
        gene_out = Gene(self.cfg, 'abc', system_foo)
        self.assertNotIn( gene_out, gene_bank)

    def test_iter(self):
        system_foo = System(self.cfg, "foo", 10)
        genes = [Gene(self.cfg, 'sctJ_FLG', system_foo), Gene(self.cfg, 'abc', system_foo)]
        for g in genes:
            gene_bank.add_gene(g)
        i = 0
        for g in gene_bank:
            self.assertIn(g, genes)
            i = i + 1
        self.assertEqual(i, len(genes))

    def test_get_uniq_object(self):
        system_foo = System(self.cfg, "foo", 10)
        gene_in = Gene(self.cfg, 'sctJ_FLG', system_foo)
        gene_bank.add_gene(gene_in)
        gene1 = gene_bank['sctJ_FLG']
        gene2 = gene_bank['sctJ_FLG']
        self.assertEqual(gene1, gene2)
        self.assertIs(gene1, gene2)


if __name__ == "__main__":
    unittest.main()