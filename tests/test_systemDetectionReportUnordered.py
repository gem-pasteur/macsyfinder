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
from macsypy.search_systems import systemDetectionReportUnordered
from tests import MacsyTest
from tests.macsy_test_env import MacsyEnvManager


class Test(MacsyTest, MacsyEnvManager):

    def setUp(self):
        self.load_env("env_002")
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        self.unload_env("env_002")
        shutil.rmtree(self.test_dir)

    def test_json_output(self):
        test_file = os.path.join(self.test_dir, 'test_foo.txt')
        so = self.macsy_test_env.system_occurence
        sdru = systemDetectionReportUnordered([so], self.macsy_test_env.cfg)
        sdru.json_output(test_file)
        self.assertFileEqual(test_file, self.output_control_file('001'))

    def test_summary_output(self):
        test_file = os.path.join(self.test_dir, 'test_bar.txt')
        so = self.macsy_test_env.system_occurence
        sdru = systemDetectionReportUnordered([so], self.macsy_test_env.cfg)
        sdru.summary_output(test_file, print_header=True)
        self.assertFileEqual(test_file, self.output_control_file('001'))

    def test_report_output(self):
        test_file = os.path.join(self.test_dir, 'test_foo.txt')
        so = self.macsy_test_env.system_occurence
        sdru = systemDetectionReportUnordered([so], self.macsy_test_env.cfg)
        sdru.report_output(test_file, print_header=True)
        self.assertFileEqual(test_file, self.output_control_file('001'))
