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
import logging
from macsypy.config import Config
from macsypy.system import System
from macsypy.search_systems import SystemOccurence, systemDetectionReport
from macsypy.registries import ModelRegistry
from macsypy.database import RepliconDB, Indexes
from tests import MacsyTest


class TestSystemDetectionReport(MacsyTest):

    def setUp(self):
        l = logging.getLogger()
        l.manager.loggerDict.clear()

        # add only one handler to the macsypy logger
        from macsypy.report import _log
        macsy_log = _log.parent
        log_file = os.devnull
        log_handler = logging.FileHandler(log_file)
        macsy_log.addHandler(log_handler)

        self.cfg = Config(hmmer_exe="hmmsearch",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type="gembase",
                          e_value_res=1,
                          i_evalue_sel=0.5,
                          res_search_dir=tempfile.gettempdir(),
                          res_search_suffix=".search_hmm.out",
                          profile_suffix=".hmm",
                          res_extract_suffix="",
                          log_level=30,
                          models_dir=self.find_data('models'),
                          log_file=log_file)

        shutil.copy(self.cfg.sequence_db, self.cfg.working_dir)
        self.cfg.options['sequence_db'] = os.path.join(self.cfg.working_dir, os.path.basename(self.cfg.sequence_db))

        idx = Indexes(self.cfg)
        idx._build_my_indexes()

        models_registry = ModelRegistry(self.cfg)
        self.model_name = 'foo'
        self.models_location = models_registry[self.model_name]

        # hack to test abstract methods (more info =>
        # https://stackoverflow.com/questions/36413844/writing-unittests-for-abstract-classes)
        self.abstractmethods = systemDetectionReport.__abstractmethods__
        systemDetectionReport.__abstractmethods__ = frozenset()

    def tearDown(self):
        # close loggers filehandles, so they don't block file deletion
        # in shutil.rmtree calls in Windows
        MacsyTest.close_loggers_filehandles()
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass

        systemDetectionReport.__abstractmethods__ = self.abstractmethods

    def test_init(self):
        system = System(self.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        os.environ['MACSY_DEBUG'] = '1'
        sdr = systemDetectionReport([system_occurence], self.cfg)
        del os.environ['MACSY_DEBUG']
        self.assertEqual(sdr._indent, 2)

    def test_report_output(self):
        system = System(self.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        sdr = systemDetectionReport([system_occurence], self.cfg)
        sdr.report_output('foo')

    def test_json_output(self):
        system = System(self.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        sdr = systemDetectionReport([system_occurence], self.cfg)
        db = RepliconDB(self.cfg)
        sdr.json_output('foo', db)

    def test_summary_output(self):
        system = System(self.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        sdr = systemDetectionReport([system_occurence], self.cfg)
        db = RepliconDB(self.cfg)
        rep_info = db['NC_xxxxx_xx']
        sdr.summary_output('foo', rep_info)
