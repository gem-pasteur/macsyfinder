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
from macsypy.system import System
from macsypy.search_systems import SystemOccurence, systemDetectionReportOrdered, system_name_generator
from macsypy.database import RepliconDB
from tests import MacsyTest, md5sum
from tests.unit import MacsyTestEnv


class Test(MacsyTest):

    def setUp(self):
        self.macsy_test_env = MacsyTestEnv()
        self.macsy_test_env.load("env_001")

        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        self.macsy_test_env.unload("env_001")

        shutil.rmtree(self.test_dir)

    def test_counter_output(self):
        system = System(self.macsy_test_env.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        sdro = systemDetectionReportOrdered('bar', [system_occurence], self.macsy_test_env.cfg)
        c = sdro.counter_output()
        self.assertEqual(c['foo_empty'], 1)

    def test_tabulated_output_header(self):
        system_occurences_states = ['single_locus', 'multi_loci']
        system = System(self.macsy_test_env.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        sdro = systemDetectionReportOrdered('bar', [system_occurence], self.macsy_test_env.cfg)
        out = sdro.tabulated_output_header(system_occurences_states, [system.name])
        expected_output = '#Replicon	foo_single_locus	foo_multi_loci\n'
        self.assertEqual(out, expected_output)

    def test_summary_output(self):
        test_file = os.path.join(self.test_dir, 'test.txt')
        system = System(self.macsy_test_env.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        sdro = systemDetectionReportOrdered('bar', [system_occurence], self.macsy_test_env.cfg)
        db = RepliconDB(self.macsy_test_env.cfg)
        rep_info = db['NC_xxxxx_xx']
        sdro.summary_output(test_file, rep_info)
        self.assertEqual(md5sum(test_file), '81d35e845603dfc1124c348231f46547')

    def test_tabulated_output(self):
        test_file = os.path.join(self.test_dir, 'test.txt')
        system_occurences_states = ['single_locus', 'multi_loci']
        system = System(self.macsy_test_env.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        sdro = systemDetectionReportOrdered('bar', [system_occurence], self.macsy_test_env.cfg)
        sdro.tabulated_output(system_occurences_states, [system.name], test_file)
        self.assertEqual(md5sum(test_file), '44feb14f77227e9a6dcb77905723d866')
