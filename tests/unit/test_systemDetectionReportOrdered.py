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
from macsypy.search_systems import SystemOccurence, systemDetectionReportOrdered
from macsypy.registries import ModelRegistry
from macsypy.database import Indexes, RepliconDB
from tests import MacsyTest, md5sum


class Test(MacsyTest):

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

        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        # close loggers filehandles, so they don't block file deletion
        # in shutil.rmtree calls in Windows
        logging.shutdown()
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        try:
            shutil.rmtree(self.cfg.working_dir)
            shutil.rmtree(self.test_dir)
        except:
            pass

    def test_counter_output(self):
        system = System(self.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        sdro = systemDetectionReportOrdered('bar', [system_occurence], self.cfg)
        c = sdro.counter_output()
        self.assertEqual(c['foo_empty'], 1)

    def test_tabulated_output_header(self):
        system_occurences_states = ['single_locus', 'multi_loci']
        system = System(self.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        sdro = systemDetectionReportOrdered('bar', [system_occurence], self.cfg)
        out = sdro.tabulated_output_header(system_occurences_states, [system.name])
        expected_output = '#Replicon	foo_single_locus	foo_multi_loci\n'
        self.assertEqual(out, expected_output)

    def test_summary_output(self):
        test_file = os.path.join(self.test_dir, 'test.txt')
        system = System(self.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        sdro = systemDetectionReportOrdered('bar', [system_occurence], self.cfg)
        db = RepliconDB(self.cfg)
        rep_info = db['NC_xxxxx_xx']
        sdro.summary_output(test_file, rep_info)
        self.assertEqual(md5sum(test_file), '81d35e845603dfc1124c348231f46547')

    def test_tabulated_output(self):
        test_file = os.path.join(self.test_dir, 'test.txt')
        system_occurences_states = ['single_locus', 'multi_loci']
        system = System(self.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        sdro = systemDetectionReportOrdered('bar', [system_occurence], self.cfg)
        sdro.tabulated_output(system_occurences_states, [system.name], test_file)
        self.assertEqual(md5sum(test_file), '44feb14f77227e9a6dcb77905723d866')
