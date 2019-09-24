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
import unittest
import shutil
import tempfile
import argparse

from macsypy.config import Config, MacsyDefaults
from macsypy.gene import Gene
import macsypy.gene
from macsypy.model import Model
from macsypy.hit import Hit
from macsypy.registries import ModelLocation
from macsypy.database import Indexes
from macsypy.search_genes import search_genes
import macsypy.search_genes
from tests import MacsyTest, which


class TestSearchGenes(MacsyTest):

    def setUp(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.res_search_dir = tempfile.gettempdir()
        args.log_level = 30
        args.out_dir = os.path.join(args.res_search_dir,
                                    'test_macsyfinder_search_genes')
        if os.path.exists(args.out_dir):
            shutil.rmtree(args.out_dir)
        os.mkdir(args.out_dir)

        self.cfg = Config(MacsyDefaults(), args)

        self.model_name = 'foo'
        self.models_location = ModelLocation(path=os.path.join(args.models_dir, self.model_name))

        idx = Indexes(self.cfg)
        idx._build_my_indexes()
        self.profile_factory = macsypy.profile.ProfileFactory(self.cfg)

    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir())
            pass
        except:
            pass

    @unittest.skipIf(not which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_search(self):
        model = Model("foo/T2SS", 10)
        gene_abc = Gene(self.profile_factory, "abc", model, self.models_location)
        report = search_genes([gene_abc], self.cfg)
        expected_hit = [Hit(gene_abc, model, "ESCO030p01_000260", 706, "ESCO030p01",
                            26, float(1.000e-200), float(660.800), float(1.000), float(0.714), 160, 663
                            )]
        self.assertEqual(len(report), 1)
        self.assertEqual(expected_hit[0], report[0].hits[0])
