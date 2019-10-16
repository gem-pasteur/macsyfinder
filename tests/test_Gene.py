#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2020  Institut Pasteur (Paris) and CNRS.           #
# See the COPYRIGHT file for details                                    #
#                                                                       #
# This file is part of MacSyFinder package.                             #
#                                                                       #
# MacSyFinder is free software: you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or     #
# (at your option) any later version.                                   #
#                                                                       #
# MacSyFinder is distributed in the hope that it will be useful,        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
# GNU General Public License for more details .                         #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with MacSyFinder (COPYING).                                     #
# If not, see <https://www.gnu.org/licenses/>.                          #
#########################################################################

import os
import shutil
import tempfile
import argparse

from macsypy.gene import CoreGene, ModelGene, Homolog, Analog, GeneStatus
from macsypy.model import Model
from macsypy.config import Config, MacsyDefaults
from macsypy.registries import ModelLocation
from tests import MacsyTest

from unittest import skipIf


class TestCoreGene(MacsyTest):

    def setUp(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.res_search_dir = tempfile.gettempdir()
        args.log_level = 30
        self.cfg = Config(MacsyDefaults(), args)
        self.model_name = 'foo'
        self.model_location = ModelLocation(path=os.path.join(args.models_dir, self.model_name))

    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir())
        except:
            pass

    def test_core_gene(self):
        model_fqn = "foo/bar"
        model = Model(model_fqn, 10)
        gene_name = 'toto'
        profile = self.model_location.get_profile(gene_name)
        cg = CoreGene(gene_name, model.family_name, profile)
        self.assertEqual(cg.name, gene_name)
        self.assertEqual(cg.model_family_name, model.family_name)
        self.assertEqual(cg.profile, profile)
        cg2 = CoreGene(gene_name, model.family_name, profile)
        self.assertTrue(isinstance(hash(cg), int))
        self.assertEqual(hash(cg), hash(cg2))
        gene_name = 'totote'
        profile = self.model_location.get_profile(gene_name)
        cg3 = CoreGene(gene_name, model.family_name, profile)
        self.assertNotEqual(hash(cg), hash(cg3))


class TestModelGene(MacsyTest):

    def setUp(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.res_search_dir = tempfile.gettempdir()
        args.log_level = 30
        self.cfg = Config(MacsyDefaults(), args)
        self.model_name = 'foo'
        self.model_location = ModelLocation(path=os.path.join(args.models_dir, self.model_name))


    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir())
        except:
            pass

    def test_hash(self):
        model_foo = Model("foo", 10)
        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        gene_1 = ModelGene(c_gene, model_foo)
        gene_2 = ModelGene(c_gene, model_foo)

        self.assertTrue(isinstance(hash(gene_1), int))
        self.assertEqual(hash(gene_1), hash(gene_1))
        self.assertNotEqual(hash(gene_1), hash(gene_2))

    def test_unknown_attribute(self):
        model_foo = Model("foo", 10)
        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        gene = ModelGene(c_gene, model_foo)
        with self.assertRaises(AttributeError) as ctx:
            gene.foo
        self.assertEqual(str(ctx.exception), "'ModelGene' object has no attribute 'foo'")


    def test_add_homolog(self):
        model_foo = Model("foo", 10)
        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        gene = ModelGene(c_gene, model_foo)
        gene_name = 'sctJ'
        profile = self.model_location.get_profile(gene_name)
        c_gene_ref = CoreGene(gene_name, model_foo.family_name, profile)
        gene_ref = ModelGene(c_gene_ref,  model_foo)
        homolog = Homolog(gene, gene_ref)
        gene.add_homolog(homolog)
        self.assertEqual(len(gene.homologs), 1)
        self.assertEqual(gene.homologs[0], homolog)


    def test_homologs(self):
        model_foo = Model("foo", 10)
        gene_name = 'sctN'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        gene = ModelGene(c_gene, model_foo)

        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ_FLG = ModelGene(c_gene, model_foo)
        gene_name = 'sctJ'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ = ModelGene(c_gene, model_foo)

        homolog_1 = Homolog(sctJ_FLG, gene)
        gene.add_homolog(homolog_1)
        homolog_2 = Homolog(sctJ, gene)
        gene.add_homolog(homolog_2)
        self.assertEqual(gene.homologs, [homolog_1, homolog_2])


    def test_is_homolog(self):
        model_foo = Model("foo", 10)

        gene_name = 'sctN'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctN = ModelGene(c_gene, model_foo)

        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ_FLG = ModelGene(c_gene, model_foo)

        gene_name = 'sctJ'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ = ModelGene(c_gene, model_foo)
        homolog = Homolog(sctJ_FLG, sctJ)
        sctJ.add_homolog(homolog)

        self.assertTrue(sctN.is_homolog(sctN))
        self.assertFalse(sctJ_FLG.is_homolog(sctJ))
        self.assertTrue(sctJ.is_homolog(sctJ_FLG))
        self.assertFalse(sctN.is_homolog(sctJ))


    def test_add_analog(self):
        model_foo = Model("foo", 10)

        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ_FLG = ModelGene(c_gene, model_foo)

        gene_name = 'sctJ'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ = ModelGene(c_gene, model_foo)

        analog = Analog(sctJ, sctJ_FLG)
        sctJ_FLG.add_analog(analog)
        self.assertEqual(len(sctJ_FLG.analogs), 1)
        self.assertEqual(sctJ_FLG.analogs[0], analog)


    def test_analogs(self):
        model_foo = Model("foo", 10)
        gene_name = 'sctN'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        gene = ModelGene(c_gene, model_foo)

        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ_FLG = ModelGene(c_gene, model_foo)
        gene_name = 'sctJ'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ = ModelGene(c_gene, model_foo)

        analog_1 = Analog(sctJ_FLG, gene)
        gene.add_analog(analog_1)
        analog_2 = Analog(sctJ, gene)
        gene.add_analog(analog_2)
        self.assertEqual(gene.analogs, [analog_1, analog_2])


    def test_is_analog(self):
        model_foo = Model("foo", 10)

        gene_name = 'sctN'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctN = ModelGene(c_gene, model_foo)

        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ_FLG = ModelGene(c_gene, model_foo)

        gene_name = 'sctJ'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ = ModelGene(c_gene, model_foo)
        analog = Analog(sctJ_FLG, sctJ)
        sctJ.add_analog(analog)

        self.assertTrue(sctN.is_analog(sctN))
        self.assertFalse(sctJ_FLG.is_analog(sctJ))
        self.assertTrue(sctJ.is_analog(sctJ_FLG))
        self.assertFalse(sctN.is_analog(sctJ))


    def test_gene_ref(self):
        model_foo = Model("foo", 10)

        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ_FLG = ModelGene(c_gene, model_foo)

        gene_name = 'sctJ'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ = ModelGene(c_gene, model_foo)
        analog = Analog(sctJ_FLG, sctJ)
        sctJ.add_analog(analog)
        self.assertEqual(analog.gene_ref, sctJ)
        self.assertIsNone(sctJ_FLG.gene_ref)


    def test_model(self):
        """
        test getter/setter for model property
        """
        model_foo = Model("foo", 10)
        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ_FLG = ModelGene(c_gene, model_foo)
        self.assertEqual(sctJ_FLG.model, model_foo)


    def test_loner(self):
        """
        test getter for loner property
        """
        model_foo = Model("foo", 10)
        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ_FLG = ModelGene(c_gene, model_foo)
        self.assertFalse(sctJ_FLG.loner)

        gene_name = 'sctJ'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ = ModelGene(c_gene, model_foo, loner=True)
        self.assertTrue(sctJ.loner)


    def test_is_mandatory(self):
        """
        test if gene belong to model mandatory genes
        """
        model_foo = Model("foo", 10)

        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ_FLG = ModelGene(c_gene, model_foo)
        model_foo.add_mandatory_gene(sctJ_FLG)
        self.assertTrue(sctJ_FLG.is_mandatory(model_foo))

        gene_name = 'sctJ'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ = ModelGene(c_gene, model_foo)
        model_foo.add_accessory_gene(sctJ)
        self.assertFalse(sctJ.is_mandatory(model_foo))


    def test_is_accessory(self):
        """
        test if gene belong to model mandatory genes
        """
        model_foo = Model("foo", 10)

        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ_FLG = ModelGene(c_gene, model_foo)
        model_foo.add_mandatory_gene(sctJ_FLG)
        self.assertFalse(sctJ_FLG.is_accessory(model_foo))

        gene_name = 'sctJ'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ = ModelGene(c_gene, model_foo)
        model_foo.add_accessory_gene(sctJ)
        self.assertTrue(sctJ.is_accessory(model_foo))


    def test_is_Forbidden(self):
        """
        test if gene belong to model mandatory genes
        """
        model_foo = Model("foo", 10)
        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ_FLG = ModelGene(c_gene, model_foo)
        model_foo.add_mandatory_gene(sctJ_FLG)
        self.assertFalse(sctJ_FLG.is_forbidden(model_foo))

        gene_name = 'sctJ'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ = ModelGene(c_gene, model_foo)
        model_foo.add_forbidden_gene(sctJ)
        self.assertTrue(sctJ.is_forbidden(model_foo))


    def test_exchangeable(self):
        """
        test getter for exchangeable property
        """
        model_foo = Model("foo", 10)
        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ_FLG = ModelGene(c_gene, model_foo)
        self.assertFalse(sctJ_FLG.exchangeable)

        gene_name = 'sctJ'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ = ModelGene(c_gene, model_foo, exchangeable=True)
        self.assertTrue(sctJ.exchangeable)


    def test_multi_system(self):
        """
        test getter for multi_system property
        """
        model_foo = Model("foo", 10)

        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ_FLG = ModelGene(c_gene, model_foo)
        self.assertFalse(sctJ_FLG.multi_system)

        gene_name = 'sctJ'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ = ModelGene(c_gene, model_foo, multi_system=True)
        self.assertTrue(sctJ.multi_system)


    def test_inter_gene_max_space(self):
        """
        test getter for inter_gene_max_space property
        """
        system_inter_gene_max_space = 40
        gene_inter_gene_max_space = 50
        model_foo = Model("foo", system_inter_gene_max_space)

        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ_FLG = ModelGene(c_gene, model_foo)
        self.assertEqual(sctJ_FLG.inter_gene_max_space, system_inter_gene_max_space)

        gene_name = 'sctJ'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ = ModelGene(c_gene, model_foo, inter_gene_max_space=gene_inter_gene_max_space)
        self.assertEqual(sctJ.inter_gene_max_space, gene_inter_gene_max_space)


    def test_str(self):
        """
        """
        model_foo = Model("foo", 10)

        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ_FLG = ModelGene(c_gene, model_foo)

        gene_name = 'sctJ'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ = ModelGene(c_gene, model_foo)

        homolog = Homolog(sctJ, sctJ_FLG)
        sctJ_FLG.add_homolog(homolog)

        gene_name = 'sctN'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctN = ModelGene(c_gene, model_foo)

        analog = Analog(sctN, sctJ_FLG)
        sctJ_FLG.add_analog(analog)
        s = """name : sctJ_FLG
inter_gene_max_space: 10
    homologs: sctJ
    analogs: sctN"""
        self.assertEqual(str(sctJ_FLG), s)

        gene_name = 'sctJ_FLG'
        profile = self.model_location.get_profile(gene_name)
        c_gene = CoreGene(gene_name, model_foo.family_name, profile)
        sctJ_FLG = ModelGene(c_gene, model_foo, loner=True, exchangeable=True, multi_system=True, inter_gene_max_space=10)
        s = """name : sctJ_FLG
inter_gene_max_space: 10
loner
multi_system
exchangeable"""
        self.assertEqual(str(sctJ_FLG), s)


class TestGeneStatus(MacsyTest):

    def test_str(self):
        self.assertEqual(str(GeneStatus.MANDATORY), 'mandatory')
