#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2021  Institut Pasteur (Paris) and CNRS.           #
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

from macsypy.gene import Exchangeable
from macsypy.gene import CoreGene, ModelGene
from macsypy.model import Model
from macsypy.profile import ProfileFactory
from macsypy.config import Config, MacsyDefaults
from macsypy.registries import ModelLocation
from macsypy.error import MacsypyError
from tests import MacsyTest


class TestExchangeable(MacsyTest):

    def setUp(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.res_search_dir = tempfile.gettempdir()
        args.log_level = 30
        self.cfg = Config(MacsyDefaults(), args)

        self.model_name = 'foo'
        self.model_location = ModelLocation(path=os.path.join(args.models_dir, self.model_name))
        self.profile_factory = ProfileFactory(self.cfg)


    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass


    def test_alternate_of(self):
        model = Model("T2SS", 10)

        gene_name = 'sctJ_FLG'
        c_gene_ref = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene_ref = ModelGene(c_gene_ref, model)

        gene_name = 'sctJ'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        homolog_1 = Exchangeable(c_gene, gene_ref)
        gene_ref.add_exchangeable(homolog_1)

        self.assertEqual(homolog_1.alternate_of(), gene_ref)

    def test_is_exchangeable(self):
        model = Model("T2SS", 10)
        gene_name = 'sctJ_FLG'
        c_gene_ref = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene_ref = ModelGene(c_gene_ref, model)

        gene_name = 'sctJ'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        homolog_1 = Exchangeable(c_gene, gene_ref)

        self.assertTrue(homolog_1.is_exchangeable)

    def test_add_exchangeable(self):
        model = Model("T2SS", 10)
        gene_name = 'sctJ_FLG'
        c_gene_ref = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene_ref = ModelGene(c_gene_ref, model)

        gene_name = 'sctJ'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        homolog_1 = Exchangeable(c_gene, gene_ref)
        homolog_2 = Exchangeable(c_gene, gene_ref)

        with self.assertRaises(MacsypyError) as ctx:
            homolog_1.add_exchangeable(homolog_2)
        self.assertEqual(str(ctx.exception),
                         "Cannot add 'Exchangeable' to an Exchangeable")

    def test_model(self):
        model = Model("T2SS", 10)
        gene_name = 'sctJ_FLG'
        c_gene_ref = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene_ref = ModelGene(c_gene_ref, model)

        gene_name = 'sctJ'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        homolog_1 = Exchangeable(c_gene, gene_ref)

        self.assertEqual(homolog_1.model, gene_ref.model)


    def test_loner(self):
        model = Model("T2SS", 10)
        gene_name = 'sctJ_FLG'
        c_gene_ref = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene_ref = ModelGene(c_gene_ref, model)
        gene_ref_loner = ModelGene(c_gene_ref, model, loner=True)

        gene_name = 'sctJ'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        homolog_1 = Exchangeable(c_gene, gene_ref)
        homolog_2 = Exchangeable(c_gene, gene_ref_loner)

        self.assertFalse(homolog_1.loner)
        self.assertTrue(homolog_2.loner)

    def test_multi_system(self):
        model = Model("T2SS", 10)
        gene_name = 'sctJ_FLG'
        c_gene_ref = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene_ref = ModelGene(c_gene_ref, model)
        gene_ref_multi_system = ModelGene(c_gene_ref, model, multi_system=True)

        gene_name = 'sctJ'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        homolog_1 = Exchangeable(c_gene, gene_ref)
        homolog_2 = Exchangeable(c_gene, gene_ref_multi_system)

        self.assertFalse(homolog_1.multi_system)
        self.assertTrue(homolog_2.multi_system)

