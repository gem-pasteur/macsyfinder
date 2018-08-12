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
from tests import MacsyTest, md5sum
from tests.unit import MacsyTestEnv


class Test(MacsyTest):

    def setUp(self):
        self.macsy_test_env = MacsyTestEnv()
        self.macsy_test_env.load("env_003")

    def tearDown(self):
        self.macsy_test_env.unload("env_003")

    def test_str(self):
        hit = self.macsy_test_env.all_hits[0]
        valid_system_hit = validSystemHit(hit, self.macsy_test_env.system, "mandatory")
        self.assertEqual(md5sum(str_=str(valid_system_hit)), '55efaa255d201e98cc104cfefe7ea8c0')
