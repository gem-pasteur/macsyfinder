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
        self.model = [m for m in self.macsy_test_env.models if m.name == 'T9SS'][0]
        self.hit = [h for h in self.macsy_test_env.all_hits if h.gene.name == 'T9SS_sprT'][0]

    def tearDown(self):
        self.unload_env("env_003")

    def test_str(self):
        valid_system_hit = validSystemHit(self.hit, self.model, "mandatory")
        self.assertEqual(str(valid_system_hit), self.output_control_str('001'))

    def test_output_system(self):
        valid_system_hit = validSystemHit(self.hit, self.model, "mandatory")
        report_str = valid_system_hit.output_system("foo132", "all_clear")
        self.assertEqual(str(report_str), self.output_control_str('001'))

    def test_getattr(self):
        valid_system_hit = validSystemHit(self.hit, self.model, "mandatory")
        i_eval = valid_system_hit.i_eval
        self.assertEqual(i_eval, 7.5e-103)

    def test_output_system_header(self):
        valid_system_hit = validSystemHit(self.hit, self.model, "mandatory")
        header_str = valid_system_hit.output_system_header()
        self.assertEqual(header_str, self.output_control_str('001'))
