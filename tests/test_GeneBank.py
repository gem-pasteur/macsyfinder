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

from macsypy.gene import GeneBank
from macsypy.gene import CoreGene, ModelGene
from macsypy.model import Model
from macsypy.config import Config, MacsyDefaults
from macsypy.registries import ModelLocation
from macsypy.error import MacsypyError
from macsypy.profile import ProfileFactory

from tests import MacsyTest


class Test(MacsyTest):

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
        self.gene_bank = GeneBank()
        self.profile_factory = ProfileFactory(self.cfg)

    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass


    def test_add_get_gene(self):
        gene_name = 'sctJ_FLG'
        with self.assertRaises(KeyError) as ctx:
            self.gene_bank[f"foo/{gene_name}"]
        self.assertEqual(str(ctx.exception),
                         f"\"No such gene 'foo/{gene_name}' in this bank\"")
        model_foo = Model(self.model_name, 10)

        self.gene_bank.add_new_gene(self.model_location, gene_name, self.profile_factory)

        gene_from_bank = self.gene_bank[(model_foo.family_name, gene_name)]
        self.assertTrue(isinstance(gene_from_bank, CoreGene))
        self.assertEqual(gene_from_bank.name, gene_name)
        gbk_contains_before = list(self.gene_bank)
        self.gene_bank.add_new_gene(self.model_location, gene_name, self.profile_factory)
        gbk_contains_after = list(self.gene_bank)
        self.assertEqual(gbk_contains_before, gbk_contains_after)

        gene_name = "bar"
        with self.assertRaises(MacsypyError) as ctx:
            self.gene_bank.add_new_gene(self.model_location, gene_name, self.profile_factory)
        self.assertEqual(str(ctx.exception),
                         f"'{self.model_name}/{gene_name}': No such profile")


    def test_contains(self):
        model_foo = Model("foo/bar", 10)
        gene_name = 'sctJ_FLG'

        self.gene_bank.add_new_gene(self.model_location, gene_name, self.profile_factory)
        gene_in = self.gene_bank[(model_foo.family_name, gene_name)]
        self.assertIn(gene_in, self.gene_bank)

        gene_name = 'abc'
        c_gene_out = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene_out = ModelGene(c_gene_out, model_foo)
        self.assertNotIn(gene_out, self.gene_bank)


    def test_iter(self):
        genes_names = ['sctJ_FLG', 'abc']
        for g in genes_names:
            self.gene_bank.add_new_gene(self.model_location, g, self.profile_factory)
        self.assertListEqual([g.name for g in self.gene_bank],
                             genes_names)


    def test_get_uniq_object(self):
        gene_name = 'sctJ_FLG'
        self.gene_bank.add_new_gene(self.model_location, gene_name, self.profile_factory)
        self.gene_bank.add_new_gene(self.model_location, gene_name, self.profile_factory)
        self.assertEqual(len(self.gene_bank), 1)
