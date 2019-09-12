# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Neron                                         #
# Copyright (c) 2014-2019  Institut Pasteur (Paris) and CNRS.                  #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################


import os
import shutil
import tempfile
import argparse

from macsypy.gene import GeneBank
from macsypy.gene import Gene
from macsypy.profile import ProfileFactory
from macsypy.model import Model
from macsypy.config import Config, MacsyDefaults
from macsypy.registries import ModelLocation
from tests import MacsyTest


class Test(MacsyTest):

    def setUp(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.res_search_dir = tempfile.gettempdir()
        args.log_level = 30
        self.cfg = Config(MacsyDefaults(), args)

        self.model_name = 'foo'
        self.models_location = ModelLocation(path=os.path.join(args.models_dir, self.model_name))
        self.gene_bank = GeneBank()
        self.profile_factory = ProfileFactory(self.cfg)

    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass

    def test_add_get_gene(self):
        gene_name = 'sctJ_FLG'
        self.assertRaises(KeyError, self.gene_bank.__getitem__, gene_name)
        system_foo = Model("foo/bar", 10)
        gene = Gene(self.profile_factory, gene_name, system_foo, self.models_location)
        self.gene_bank.add_gene(gene)
        gene_from_bank = self.gene_bank[(self.model_name, gene_name)]
        self.assertTrue(isinstance(gene_from_bank, Gene))
        self.assertEqual(gene_from_bank, gene)
        with self.assertRaises(KeyError) as ctx:
            self.gene_bank.add_gene(gene)
        self.assertEqual(str(ctx.exception),
                         '"a gene named \'{0}/{1}\' is already registered"'.format("foo", gene.name))

    def test_contains(self):
        system_foo = Model("foo/bar", 10)
        gene_in = Gene(self.profile_factory, 'sctJ_FLG', system_foo, self.models_location)
        self.gene_bank.add_gene(gene_in)
        self.assertIn(gene_in, self.gene_bank)
        gene_out = Gene(self.profile_factory, 'abc', system_foo, self.models_location)
        self.assertNotIn(gene_out, self.gene_bank)

    def test_iter(self):
        system_foo = Model("foo/bar", 10)
        genes = [Gene(self.profile_factory, 'sctJ_FLG', system_foo, self.models_location),
                 Gene(self.profile_factory, 'abc', system_foo, self.models_location)]
        for g in genes:
            self.gene_bank.add_gene(g)
        i = 0
        for g in self.gene_bank:
            self.assertIn(g, genes)
            i += 1
        self.assertEqual(i, len(genes))


    def test_get_uniq_object(self):
        system_foo = Model("foo", 10)
        gene_in = Gene(self.profile_factory, 'sctJ_FLG', system_foo, self.models_location)
        self.gene_bank.add_gene(gene_in)
        gene1 = self.gene_bank[(self.model_name, 'sctJ_FLG')]
        gene2 = self.gene_bank[(self.model_name, 'sctJ_FLG')]
        self.assertEqual(gene1, gene2)
        self.assertIs(gene1, gene2)
