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
import argparse

from macsypy.config import MacsyDefaults, Config
from macsypy.model import Model, ModelBank
from macsypy.search_systems import SystemOccurence, systemDetectionReport
from macsypy.registries import ModelRegistry
from macsypy.database import RepliconDB, Indexes
from tests import MacsyTest


class TestModelDetectionReport(MacsyTest):

    def setUp(self):
        args = argparse.Namespace()
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.res_search_dir = tempfile.gettempdir()
        args.log_level = 30
        args.out_dir = os.path.join(args.res_search_dir,
                                    'test_macsyfinder_systemDetectionReport')
        if os.path.exists(args.out_dir):
            shutil.rmtree(args.out_dir)
        os.mkdir(args.out_dir)

        seq_db = self.find_data("base", "test_base.fa")
        shutil.copy(seq_db, args.out_dir)
        args.sequence_db = os.path.join(args.out_dir, os.path.basename(seq_db))
        self.cfg = Config(MacsyDefaults(), args)

        idx = Indexes(self.cfg)
        idx._build_my_indexes()

        models_registry = ModelRegistry(self.cfg)
        self.model_name = 'foo'
        self.models_location = models_registry[self.model_name]
        self.model_bank = ModelBank()

        # hack to test abstract methods (more info =>
        # https://stackoverflow.com/questions/36413844/writing-unittests-for-abstract-classes)
        self.abstractmethods = systemDetectionReport.__abstractmethods__
        systemDetectionReport.__abstractmethods__ = frozenset()


    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir())
        except:
            pass

        systemDetectionReport.__abstractmethods__ = self.abstractmethods


    def test_init(self):
        system = Model(self.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        os.environ['MACSY_DEBUG'] = '1'
        sdr = systemDetectionReport([system_occurence], self.cfg, self.model_bank)
        del os.environ['MACSY_DEBUG']
        self.assertEqual(sdr._indent, 2)


    def test_report_output(self):
        system = Model(self.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        sdr = systemDetectionReport([system_occurence], self.cfg, self.model_bank)
        sdr.report_output('foo')


    def test_json_output(self):
        system = Model(self.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        sdr = systemDetectionReport([system_occurence], self.cfg, self.model_bank)
        db = RepliconDB(self.cfg)
        sdr.json_output('foo', db)


    def test_summary_output(self):
        system = Model(self.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        sdr = systemDetectionReport([system_occurence], self.cfg, self.model_bank)
        db = RepliconDB(self.cfg)
        rep_info = db['NC_xxxxx_xx']
        sdr.summary_output('foo', rep_info)
