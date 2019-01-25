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

from macsypy.gene import Gene
from macsypy.gene import Homolog, Analog
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
            shutil.rmtree(self.cfg.working_dir())
        except:
            pass


    def test_add_homolog(self):
        system_foo = Model(self.cfg, "foo", 10)
        system_bar = Model(self.cfg, "bar", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.models_location)
        gene_ref = Gene(self.cfg, 'sctJ', system_bar, self.models_location)
        homolog = Homolog(self.cfg, gene, gene_ref)
        gene.add_homolog(homolog)
        self.assertEqual(len(gene.homologs), 1)
        self.assertEqual(gene.homologs[0], homolog)


    def test_get_homologs(self):
        system_foo = Model(self.cfg, "foo", 10)
        system_bar = Model(self.cfg, "bar", 10)
        gene = Gene(self.cfg, 'sctN', system_foo, self.models_location)
        sctJ_FLG = Gene(self.cfg, 'sctJ_FLG', system_foo, self.models_location)
        sctJ = Gene(self.cfg, 'sctJ', system_bar, self.models_location)
        homolog_1 = Homolog(sctJ_FLG, gene)
        gene.add_homolog(homolog_1)
        homolog_2 = Homolog(sctJ, gene)
        gene.add_homolog(homolog_2)
        self.assertEqual(gene.get_homologs(), [homolog_1, homolog_2])


    def test_is_homolog(self):
        system_foo = Model(self.cfg, "foo", 10)
        system_bar = Model(self.cfg, "bar", 10)
        gene = Gene(self.cfg, 'sctN', system_foo, self.models_location)
        sctJ_FLG = Gene(self.cfg, 'sctJ_FLG', system_foo, self.models_location)
        sctJ = Gene(self.cfg, 'sctJ', system_bar, self.models_location)
        homolog = Homolog(sctJ_FLG, gene)
        gene.add_homolog(homolog)
        self.assertTrue(gene.is_homolog(gene))
        self.assertTrue(gene.is_homolog(homolog))
        self.assertFalse(gene.is_homolog(sctJ))


    def test_add_analog(self):
        system_foo = Model(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.models_location)
        analog = Analog(self.cfg, gene)
        gene.add_analog(analog)
        self.assertEqual(len(gene.analogs), 1)
        self.assertEqual(gene.analogs[0], analog)


    def test_get_analogs(self):
        system_foo = Model(self.cfg, "foo", 10)
        system_bar = Model(self.cfg, "bar", 10)
        gene = Gene(self.cfg, 'sctN', system_foo, self.models_location)
        sctJ_FLG = Gene(self.cfg, 'sctJ_FLG', system_foo, self.models_location)
        sctJ = Gene(self.cfg, 'sctJ', system_bar, self.models_location)
        analog_1 = Analog(sctJ_FLG, gene)
        gene.add_analog(analog_1)
        analog_2 = Analog(sctJ, gene)
        gene.add_analog(analog_2)
        self.assertEqual(gene.get_analogs(), [analog_1, analog_2])


    def test_is_analog(self):
        system_foo = Model(self.cfg, "foo", 10)
        system_bar = Model(self.cfg, "bar", 10)
        gene = Gene(self.cfg, 'sctN', system_foo, self.models_location)
        sctJ_FLG = Gene(self.cfg, 'sctJ_FLG', system_foo, self.models_location)
        sctJ = Gene(self.cfg, 'sctJ', system_bar, self.models_location)
        analog = Analog(sctJ_FLG, gene)
        gene.add_analog(analog)
        self.assertTrue(gene.is_analog(gene))
        self.assertTrue(gene.is_analog(analog))
        self.assertFalse(gene.is_analog(sctJ))


    def test_system(self):
        """
        test getter/setter for system property
        """
        system_foo = Model(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.models_location)
        self.assertEqual(gene.system, system_foo)


    def test_loner(self):
        """
        test getter for loner property
        """
        system_foo = Model(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.models_location)
        self.assertFalse(gene.loner)
        gene = Gene(self.cfg, 'sctJ', system_foo, self.models_location, loner=True)
        self.assertTrue(gene.loner)


    def test_exchangeable(self):
        """
        test getter for exchangeable property
        """
        system_foo = Model(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.models_location)
        self.assertFalse(gene.exchangeable)
        gene = Gene(self.cfg, 'sctJ', system_foo, self.models_location, exchangeable=True)
        self.assertTrue(gene.exchangeable)


    def test_multi_system(self):
        """
        test getter for multi_system property
        """
        system_foo = Model(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.models_location)
        self.assertFalse(gene.multi_system)
        gene = Gene(self.cfg, 'sctJ', system_foo, self.models_location, multi_system=True)
        self.assertTrue(gene.multi_system)


    def test_inter_gene_max_space(self):
        """
        test getter for inter_gene_max_space property
        """
        system_inter_gene_max_space = 40
        gene_inter_gene_max_space = 50
        system_foo = Model(self.cfg, "foo", system_inter_gene_max_space)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.models_location)
        self.assertEqual(gene.inter_gene_max_space, system_inter_gene_max_space)
        gene = Gene(self.cfg, 'sctJ', system_foo, self.models_location, inter_gene_max_space=gene_inter_gene_max_space)
        self.assertEqual(gene.inter_gene_max_space, gene_inter_gene_max_space)


    def test_str(self):
        """
        """
        system_foo = Model(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.models_location)
        system_bar = Model(self.cfg, "bar", 20)
        gene_homolog = Gene(self.cfg, 'sctJ', system_bar, self.models_location)
        homolog = Homolog(gene_homolog, gene, self.cfg)
        gene.add_homolog(homolog)
        analog = Gene(self.cfg, 'sctN', system_foo, self.models_location)
        gene.add_analog(analog)
        s = """name : sctJ_FLG
inter_gene_max_space: 10
    homologs: sctJ
    analogs: sctN"""
        self.assertEqual(str(gene), s)

        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.models_location,
                    loner=True, exchangeable=True, multi_system=True, inter_gene_max_space=10)
        s = """name : sctJ_FLG
inter_gene_max_space: 10
loner
multi_system
exchangeable"""
        self.assertEqual(str(gene), s)
