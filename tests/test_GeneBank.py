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
import shutil
import tempfile
import logging
import argparse

from macsypy.gene import gene_bank
from macsypy.gene import Gene
from macsypy.system import System
from macsypy.config import Config, MacsyDefaults
from macsypy.registries import ModelRegistry
from tests import MacsyTest


class Test(MacsyTest):

    def setUp(self):
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        
        # add only one handler to the macsypy logger
        from macsypy.gene import _log
        macsy_log = _log.parent
        log_file = os.devnull
        log_handler = logging.FileHandler(log_file)
        macsy_log.addHandler(log_handler)

        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.res_search_dir = tempfile.gettempdir()
        args.log_level = 30
        args.log_file = log_file
        self.cfg = Config(MacsyDefaults(), args)

        models_registry = ModelRegistry(self.cfg)
        self.model_name = 'foo'
        self.models_location = models_registry[self.model_name]


    def tearDown(self):
        # close loggers filehandles, so they don't block file deletion
        # in shutil.rmtree calls in Windows
        logging.shutdown()
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        gene_bank._genes_bank = {}
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass

    def test_add_get_gene(self):
        gene_name = 'sctJ_FLG'
        self.assertRaises(KeyError, gene_bank.__getitem__, gene_name)
        system_foo = System(self.cfg, "foo/bar", 10)
        gene = Gene(self.cfg, gene_name, system_foo, self.models_location)
        gene_bank.add_gene(gene)
        gene_from_bank = gene_bank[(self.model_name, gene_name)]
        self.assertTrue(isinstance(gene_from_bank, Gene))
        self.assertEqual(gene_from_bank, gene)
        with self.assertRaises(KeyError) as ctx:
            gene_bank.add_gene(gene)
        self.assertEqual(str(ctx.exception),
                         '"a gene named \'{0}/{1}\' is already registered"'.format("foo", gene.name))

    def test_contains(self):
        system_foo = System(self.cfg, "foo/bar", 10)
        gene_in = Gene(self.cfg, 'sctJ_FLG', system_foo, self.models_location)
        gene_bank.add_gene(gene_in)
        self.assertIn(gene_in, gene_bank)
        gene_out = Gene(self.cfg, 'abc', system_foo, self.models_location)
        self.assertNotIn(gene_out, gene_bank)

    def test_iter(self):
        system_foo = System(self.cfg, "foo/bar", 10)
        genes = [Gene(self.cfg, 'sctJ_FLG', system_foo, self.models_location),
                 Gene(self.cfg, 'abc', system_foo, self.models_location)]
        for g in genes:
            gene_bank.add_gene(g)
        i = 0
        for g in gene_bank:
            self.assertIn(g, genes)
            i += 1
        self.assertEqual(i, len(genes))


    def test_get_uniq_object(self):
        system_foo = System(self.cfg, "foo", 10)
        gene_in = Gene(self.cfg, 'sctJ_FLG', system_foo, self.models_location)
        gene_bank.add_gene(gene_in)
        gene1 = gene_bank[(self.model_name, 'sctJ_FLG')]
        gene2 = gene_bank[(self.model_name, 'sctJ_FLG')]
        self.assertEqual(gene1, gene2)
        self.assertIs(gene1, gene2)
