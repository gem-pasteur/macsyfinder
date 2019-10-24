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

from macsypy.gene import Homolog
from macsypy.gene import CoreGene, ModelGene
from macsypy.model import Model
from macsypy.profile import ProfileFactory
from macsypy.config import Config, MacsyDefaults
from macsypy.registries import ModelLocation
from tests import MacsyTest


class TestHomolog(MacsyTest):

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
        homolog_1 = Homolog(c_gene, gene_ref)
        gene_ref.add_homolog(homolog_1)

        self.assertEqual(homolog_1.alternate_of(), gene_ref)


    def test_is_aligned(self):
        model = Model("T2SS", 10)

        gene_name = 'sctJ_FLG'
        c_gene_ref = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene_ref = ModelGene(c_gene_ref, model)

        gene_name = 'sctJ'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene = ModelGene(c_gene, model)

        homolog = Homolog(gene, gene_ref)
        self.assertFalse(homolog.is_aligned())

        homolog = Homolog(gene, gene_ref, aligned=True)
        self.assertTrue(homolog.is_aligned())
