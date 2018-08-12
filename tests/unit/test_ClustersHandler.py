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


from macsypy.search_systems import ClustersHandler
from macsypy.macsypy_error import SystemDetectionError
from tests import MacsyTest, md5sum
from tests.unit import MacsyTestEnv


class Test(MacsyTest):

    def setUp(self):
        self.macsy_test_env = MacsyTestEnv()
        self.macsy_test_env.load("env_002")

    def tearDown(self):
        self.macsy_test_env.unload("env_002")

    def test_add(self):
        ch = ClustersHandler()
        cluster = self.macsy_test_env.cluster

        ch.replicon_name = None
        cluster.replicon_name = 'foo'
        ch.add(cluster)
        self.assertEqual(ch.replicon_name, 'foo')
        self.assertEqual(len(ch.clusters), 1)

        ch.replicon_name = 'foo'
        cluster.replicon_name = 'foo'
        ch.add(cluster)
        self.assertEqual(len(ch.clusters), 2)

        ch.replicon_name = 'foo'
        cluster.replicon_name = 'bla'
        with self.assertRaises(SystemDetectionError):
            ch.add(cluster)

    def test_str(self):
        ch = ClustersHandler()
        cluster = self.macsy_test_env.cluster
        ch.add(cluster)
        ch.add(cluster)
        str_ = str(ch)
        self.assertEqual(md5sum(str_=str_), 'f08677489ff1dd2ed89f27995427233a')
