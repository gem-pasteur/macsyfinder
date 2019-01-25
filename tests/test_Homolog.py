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


import shutil
import tempfile
import argparse

from macsypy.gene import Homolog
from macsypy.gene import Gene
from macsypy.model import Model
from macsypy.config import Config, MacsyDefaults
from macsypy.registries import ModelRegistry
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

        models_registry = ModelRegistry(self.cfg)
        self.model_name = 'foo'
        self.models_location = models_registry[self.model_name]


    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass

    def test_gene_ref(self):
        system = Model(self.cfg, "T2SS", 10)
        gene_ref = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        gene = Gene(self.cfg, 'sctJ', system, self.models_location)
        homolog_1 = Homolog(gene, gene_ref)
        self.assertEqual(homolog_1.gene_ref, gene_ref)
 
    def test_is_aligned(self):
        system = Model(self.cfg, "T2SS", 10)
        gene_ref = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        gene = Gene(self.cfg, 'sctJ', system, self.models_location)
        homolog = Homolog(gene, gene_ref)
        self.assertFalse(homolog.is_aligned())
        homolog = Homolog(gene, gene_ref, aligned=True)
        self.assertTrue(homolog.is_aligned())

    def test_delegation(self):
        system = Model(self.cfg, "T2SS", 10)
        gene_ref = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        gene = Gene(self.cfg, 'sctJ', system, self.models_location)
        homolog = Homolog(gene, gene_ref)
        self.assertEqual(homolog.system, system)

