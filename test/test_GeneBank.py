# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Néron                                         #
# Copyright © 2014  Institut Pasteur, Paris.                                   #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################




import os
import unittest
import shutil
import tempfile
import platform
import logging
from macsypy.gene import gene_bank
from macsypy.gene import Gene
from macsypy.system import System
from macsypy.config import Config
from macsypy.registries import ProfilesRegistry


class Test(unittest.TestCase):

    _data_dir = os.path.join(os.path.dirname(__file__), "datatest")


    def setUp(self):
        self.cfg = Config( sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type = "gembase",
                           hmmer_exe = "",
                           e_value_res = 1,
                           i_evalue_sel = 0.5,
                           def_dir = os.path.join(self._data_dir, 'DEF'),
                           res_search_dir = tempfile.gettempdir(),
                           res_search_suffix = "",
                           profile_dir = os.path.join(self._data_dir, 'profiles'),
                           profile_suffix = ".hmm",
                           res_extract_suffix = "",
                           log_level = 30,
                           log_file = 'NUL' if platform.system() == 'Windows' else '/dev/null'
                           )
        self.profile_registry = ProfilesRegistry(self.cfg)
        

    def tearDown(self):
        # close loggers filehandles, so they don't block file deletion
        # in shutil.rmtree calls in Windows
        logging.shutdown()
        gene_bank._genes_bank = {}
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass

    def test_add_get_gene(self):
        gene_name = 'sctJ_FLG'
        self.assertRaises(KeyError, gene_bank.__getitem__, gene_name)
        system_foo = System( "foo", self.cfg, 10)
        gene = Gene(self.cfg, gene_name, system_foo, self.profile_registry)
        gene_bank.add_gene(gene)
        gene_from_bank = gene_bank[gene_name]
        self.assertTrue(isinstance(gene_from_bank, Gene))
        self.assertEqual(gene_from_bank, gene)

    def test_contains(self):
        system_foo = System(self.cfg, "foo", 10)
        gene_in = Gene(self.cfg, 'sctJ_FLG', system_foo, self.profile_registry)
        gene_bank.add_gene(gene_in)
        self.assertIn(gene_in, gene_bank)
        gene_out = Gene(self.cfg, 'abc', system_foo, self.profile_registry)
        self.assertNotIn( gene_out, gene_bank)

    def test_iter(self):
        system_foo = System(self.cfg, "foo", 10)
        genes = [Gene(self.cfg, 'sctJ_FLG', system_foo, self.profile_registry), Gene(self.cfg, 'abc', system_foo, self.profile_registry)]
        for g in genes:
            gene_bank.add_gene(g)
        i = 0
        for g in gene_bank:
            self.assertIn(g, genes)
            i = i + 1
        self.assertEqual(i, len(genes))

    def test_get_uniq_object(self):
        system_foo = System(self.cfg, "foo", 10)
        gene_in = Gene(self.cfg, 'sctJ_FLG', system_foo, self.profile_registry)
        gene_bank.add_gene(gene_in)
        gene1 = gene_bank['sctJ_FLG']
        gene2 = gene_bank['sctJ_FLG']
        self.assertEqual(gene1, gene2)
        self.assertIs(gene1, gene2)
