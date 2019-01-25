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


from macsypy.search_systems import validSystemHit
from tests import MacsyTest
from tests.macsy_test_env import MacsyEnvManager


class Test(MacsyTest, MacsyEnvManager):

    def setUp(self):
        self.load_env("env_003", log_out=False,)

    def tearDown(self):
        self.unload_env("env_003")
        pass

    def test_str(self):
        hit = self.macsy_test_env.all_hits[0]
        valid_system_hit = validSystemHit(hit, self.macsy_test_env.system, "mandatory")
        self.assertEqual(str(valid_system_hit), self.output_control_str('001'))

    def test_output_system(self):
        hit = self.macsy_test_env.all_hits[0]
        valid_system_hit = validSystemHit(hit, self.macsy_test_env.system, "mandatory")
        report_str = valid_system_hit.output_system("foo132", "all_clear")
        self.assertEqual(str(report_str), self.output_control_str('001'))

    def test_getattr(self):
        hit = self.macsy_test_env.all_hits[0]
        valid_system_hit = validSystemHit(hit, self.macsy_test_env.system, "mandatory")
        i_eval = valid_system_hit.i_eval
        self.assertEqual(i_eval, 7.5e-103)

    def test_output_system_header(self):
        hit = self.macsy_test_env.all_hits[0]
        valid_system_hit = validSystemHit(hit, self.macsy_test_env.system, "mandatory")
        header_str = valid_system_hit.output_system_header()
        self.assertEqual(header_str, self.output_control_str('001'))
