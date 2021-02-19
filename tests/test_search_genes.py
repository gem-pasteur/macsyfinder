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
import unittest
import shutil
import tempfile
import argparse

import macsypy
from macsypy.config import Config, MacsyDefaults
from macsypy.gene import CoreGene
from macsypy.hit import Hit
from macsypy.registries import ModelLocation
from macsypy.profile import ProfileFactory
from macsypy.database import Indexes
from macsypy.search_genes import search_genes
from tests import MacsyTest, which


class TestSearchGenes(MacsyTest):

    def setUp(self):
        self.tmp_dir = os.path.join(tempfile.gettempdir(),
                                    'test_macsyfinder_search_genes')
        if os.path.exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.mkdir(self.tmp_dir)

        macsypy.init_logger()
        macsypy.logger_set_level(30)

        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.log_level = 30
        args.out_dir = os.path.join(self.tmp_dir, 'job_1')
        args.res_search_dir = args.out_dir
        args.no_cut_ga = True
        os.mkdir(args.out_dir)

        self.cfg = Config(MacsyDefaults(), args)

        self.model_name = 'foo'
        self.model_location = ModelLocation(path=os.path.join(args.models_dir, self.model_name))

        idx = Indexes(self.cfg)
        idx._build_my_indexes()
        self.profile_factory = ProfileFactory(self.cfg)

    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir())
            #pass
        except:
            pass


    @unittest.skipIf(not which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_search(self):
        gene_name = "abc"
        c_gene_abc = CoreGene(self.model_location, gene_name, self.profile_factory)
        report = search_genes([c_gene_abc], self.cfg)
        expected_hit = [Hit(c_gene_abc, "ESCO030p01_000260", 706, "ESCO030p01",
                            26, float(1.000e-200), float(660.800), float(1.000), float(0.714), 160, 663
                            )]
        self.assertEqual(len(report), 1)
        self.assertEqual(expected_hit[0], report[0].hits[0])


    @unittest.skipIf(not which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_search_recover(self):
        # first job searching using hmmsearch
        gene_name = "abc"
        c_gene_abc = CoreGene(self.model_location, gene_name, self.profile_factory)
        report = search_genes([c_gene_abc], self.cfg)
        expected_hit = [Hit(c_gene_abc, "ESCO030p01_000260", 706, "ESCO030p01",
                            26, float(1.000e-200), float(660.800), float(1.000), float(0.714), 160, 663
                            )]

        # second job using recover
        # disable hmmer to be sure that test use the recover inner function
        self.cfg.hmmer = lambda: "hmmer_disable"
        # and create a new dir for the second job
        previous_job_path = self.cfg.working_dir()
        self.cfg.previous_run = lambda: previous_job_path
        self.cfg.out_dir = lambda: os.path.join(self.tmp_dir, 'job_2')
        os.mkdir(self.cfg.out_dir())

        # rerun with previous run
        # but we have to reset the profile attached to the gene gene._profile._report
        self.profile_factory = ProfileFactory(self.cfg)
        c_gene_abc = CoreGene(self.model_location, gene_name, self.profile_factory)
        report = search_genes([c_gene_abc], self.cfg)
        self.assertEqual(len(report), 1)
        self.assertEqual(expected_hit[0], report[0].hits[0])
