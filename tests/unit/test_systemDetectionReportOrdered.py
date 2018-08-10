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
from macsypy.search_systems import systemDetectionReportOrdered
from macsypy.database import RepliconDB
from tests import MacsyTest, md5sum
from tests.unit import MacsyTestEnv


class Test(MacsyTest):

    def setUp(self):
        self.macsy_test_env = MacsyTestEnv()
        self.macsy_test_env.load("env_002")

        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        self.macsy_test_env.unload("env_002")

        shutil.rmtree(self.test_dir)

    def test_counter_output(self):
        so = self.macsy_test_env.system_occurence
        sdro = systemDetectionReportOrdered('bar', [so], self.macsy_test_env.cfg)
        c = sdro.counter_output()
        self.assertEqual(c['T9SS_multi_loci'], 1)

    def test_tabulated_output_header(self):
        system_occurences_states = ['single_locus', 'multi_loci']
        so = self.macsy_test_env.system_occurence
        sdro = systemDetectionReportOrdered('bar', [so], self.macsy_test_env.cfg)
        out = sdro.tabulated_output_header(system_occurences_states, [self.macsy_test_env.system.name])
        expected_output = '#Replicon\tT9SS_single_locus\tT9SS_multi_loci\n'
        self.assertEqual(out, expected_output)

    def test_summary_output(self):
        test_file = os.path.join(self.test_dir, 'test.txt')
        so = self.macsy_test_env.system_occurence
        sdro = systemDetectionReportOrdered('bar', [so], self.macsy_test_env.cfg)
        db = RepliconDB(self.macsy_test_env.cfg)
        rep_info = db['AESU001c01a']
        sdro.summary_output(test_file, rep_info)
        self.assertEqual(md5sum(test_file), 'ff58f3dbff33338db2f272ad0b3923b8')

    def test_tabulated_output(self):
        test_file = os.path.join(self.test_dir, 'test.txt')
        system_occurences_states = ['single_locus', 'multi_loci']
        so = self.macsy_test_env.system_occurence
        sdro = systemDetectionReportOrdered('bar', [so], self.macsy_test_env.cfg)
        sdro.tabulated_output(system_occurences_states, [self.macsy_test_env.system.name], test_file)
        self.assertEqual(md5sum(test_file), 'aa175729fe8a3b787cf84b9c09cf5e07')

    def test_report_output(self):
        test_file = os.path.join(self.test_dir, 'test_foo.txt')
        so = self.macsy_test_env.system_occurence
        sdro = systemDetectionReportOrdered('bar', [so], self.macsy_test_env.cfg)
        sdro.report_output(test_file)
        self.assertEqual(md5sum(test_file), '10791ef7c77ec198863d645e00bc4d6e')

    def test_system_2_json(self):
        so = self.macsy_test_env.system_occurence
        so.get_system_unique_name('mew')
        sdro = systemDetectionReportOrdered('bar', [so], self.macsy_test_env.cfg)
        rep_db = RepliconDB(self.macsy_test_env.cfg)
        out = sdro.system_2_json(rep_db)
        self.assertEqual(md5sum(str_=str(out)), '0a21f19165da3de2241f7602c61377cc')

    def test_match2json(self):
        so = self.macsy_test_env.system_occurence
        valid_hit = so.valid_hits[0]
        sdro = systemDetectionReportOrdered('bar', [so], self.macsy_test_env.cfg)
        gene = sdro._match2json(valid_hit, so)
        self.assertEqual(md5sum(str_=str(gene)), '8c00bb9dfb43f4595f2f33a1a643c202')

    def test_json_output(self):
        test_file = os.path.join(self.test_dir, 'test_foo.txt')
        so = self.macsy_test_env.system_occurence
        so.get_system_unique_name('mew')
        sdro = systemDetectionReportOrdered('bar', [so], self.macsy_test_env.cfg)
        rep_db = RepliconDB(self.macsy_test_env.cfg)
        json_all_systems = sdro.system_2_json(rep_db)
        sdro.json_output(test_file, json_all_systems)
        self.assertEqual(md5sum(test_file), '662683cc95dfddbbbd7807027b7c3327')

    def test_gene2json(self):
        so = self.macsy_test_env.system_occurence
        sdro = systemDetectionReportOrdered('bar', [so], self.macsy_test_env.cfg)
        gene = sdro._gene2json('foobar', 44, 72)
        self.assertEqual(md5sum(str_=str(gene)), '6fce3f2b977ed8ddf267efcd313acba9')
