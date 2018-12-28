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
import unittest
import shutil
import tempfile
import logging
import argparse

from macsypy.config import Config, MacsyDefaults
from macsypy.gene import Gene
import macsypy.gene
from macsypy.system import System
from macsypy.report import Hit
from macsypy.registries import ModelRegistry
from macsypy.database import Indexes
from macsypy.utils import which
from macsypy.search_genes import search_genes
import macsypy.search_genes
from tests import MacsyTest


class TestSearchGenes(MacsyTest):

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
        args.out_dir = os.path.join(args.res_search_dir,
                                    'test_macsyfinder_search_genes')
        if os.path.exists(args.out_dir):
            shutil.rmtree(args.out_dir)
        os.mkdir(args.out_dir)

        self.cfg = Config(MacsyDefaults(), args)

        models_registry = ModelRegistry(self.cfg)
        self.model_name = 'foo'
        self.models_location = models_registry[self.model_name]
        idx = Indexes(self.cfg)
        idx._build_my_indexes()
        macsypy.gene.profile_factory = macsypy.gene.ProfileFactory()
        macsypy.gene.profile_factory._profiles = {}

    def tearDown(self):
        # close loggers filehandles, so they don't block file deletion
        # in shutil.rmtree calls in Windows
        logging.shutdown()
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        try:
            shutil.rmtree(self.cfg.working_dir())
            pass
        except:
            pass
        macsypy.gene.profile_factory = macsypy.gene.ProfileFactory()
        macsypy.gene.profile_factory._profiles = {}

    @unittest.skipIf(not which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_search(self):
        system = System(self.cfg, "foo/T2SS", 10)
        gene_abc = Gene(self.cfg, "abc", system, self.models_location)
        report = search_genes([gene_abc], self.cfg)
        expected_hit = [Hit(gene_abc, system, "ESCO030p01_000260", 706, "ESCO030p01",
                            26, float(1.000e-200), float(660.800), float(1.000), float(0.714), 160, 663
                            )]
        self.assertEqual(len(report), 1)
        self.assertEqual(expected_hit[0], report[0].hits[0])
