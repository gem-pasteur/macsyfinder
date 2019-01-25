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


from macsypy.database import RepliconInfo
from macsypy.report import Hit
from macsypy.gene import Gene
from macsypy.search_systems import ClustersHandler, Cluster
from macsypy.error import SystemDetectionError
from tests import MacsyTest
from tests.macsy_test_env import MacsyEnvManager


class TestCircularizeData(object):
    """Helper class used in test_circularize() method."""

    def __init__(self, cfg, models_location, system):
        self.cfg = cfg
        self.system = system
        self.models_location = models_location

        self._hits_group = (None, None, None, None, None, None)

    def build_clusters_group(self):
        c1 = self.build_cluster(10, 15)
        c2 = self.build_cluster(11, 16)
        c3 = self.build_cluster(12, 17)

        return [c1, c2, c3]

    def build_hit(self, pos, max_space):
        gene_name = 'T9SS_gldJ_TIGR03524'
        gene = Gene(self.cfg, gene_name, self.system, self.models_location, inter_gene_max_space=max_space)
        hit = Hit(gene, self.system, None, None, '', pos, None, None, None, None, None, None)

        return hit

    def set_hits_group(self, p1, s1, p2, s2, p3, s3):
        self._hits_group = (p1, s1, p2, s2, p3, s3)

    def get_hits_group(self):
        (p1, s1, p2, s2, p3, s3) = self._hits_group

        h1 = self.build_hit(p1, s1)
        h2 = self.build_hit(p2, s2)
        h3 = self.build_hit(p3, s3)

        return [h1, h2, h3]

    def build_cluster(self, begin, end):
        c = Cluster([self.system])
        c.end = end
        c.begin = begin
        c.hits = self.get_hits_group()
        return c

    def get_end_hits(self, p1, s1, p2, s2):
        h1 = self.build_hit(p1, s1)
        h2 = self.build_hit(p2, s2)
        return [h1, h2]


class Test(MacsyTest, MacsyEnvManager):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_add(self):
        self.load_env("env_002", log_out=False)

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

        self.unload_env("env_002")


    def test_str(self):
        self.load_env("env_002", log_out=False)

        ch = ClustersHandler()
        cluster = self.macsy_test_env.cluster
        ch.add(cluster)
        ch.add(cluster)
        str_ = str(ch)
        self.assertEqual(str(str_), self.output_control_str('001'))

        self.unload_env("env_002")


    def test_circularize(self):
        self.load_env("env_004", log_out=False)

        system = self.macsy_test_env.model
        cfg = self.macsy_test_env.cfg
        models_location = self.macsy_test_env.models_location
        ch = ClustersHandler()
        tcd = TestCircularizeData(cfg, models_location, system)

        # case 1
        tcd.set_hits_group(1, 3, 6, 3, 15, 3)
        ch.clusters = tcd.build_clusters_group()
        end_hits = tcd.get_end_hits(11, 3, 6, 3)
        rep_info = RepliconInfo('circular', 200, 400, [])
        ch.circularize(rep_info, end_hits, [system])
        str_ = str(ch)
        self.assertEqual(str(str_), self.output_control_str('001'))

        # case 2
        tcd.set_hits_group(1, 3, 6, 3, 15, 3)
        ch.clusters = tcd.build_clusters_group()
        end_hits = tcd.get_end_hits(11, 50, 6, 50)
        rep_info = RepliconInfo('circular', 200, 220, [])
        ch.circularize(rep_info, end_hits, [system])
        str_ = str(ch)
        self.assertEqual(str(str_), self.output_control_str('002'))

        # case 3
        tcd.set_hits_group(1, 350, 6, 30, 15, 30)
        ch.clusters = tcd.build_clusters_group()
        end_hits = tcd.get_end_hits(11, 30, 6, 30)
        rep_info = RepliconInfo('circular', 200, 400, [])
        ch.circularize(rep_info, end_hits, [system])
        str_ = str(ch)
        self.assertEqual(str(str_), self.output_control_str('003'))

        # case 4
        tcd.set_hits_group(30, 3, 6, 3, 15, 3)
        ch.clusters = tcd.build_clusters_group()
        end_hits = tcd.get_end_hits(11, 195, 6, 128)
        rep_info = RepliconInfo('circular', 200, 400, [])
        ch.circularize(rep_info, end_hits, [system])
        str_ = str(ch)
        self.assertEqual(str(str_), self.output_control_str('004'))

        # case 5
        tcd.set_hits_group(1, 3, 6, 3, 15, 400)
        ch.clusters = tcd.build_clusters_group()
        end_hits = tcd.get_end_hits(11, 3, 6, 3)
        rep_info = RepliconInfo('circular', 200, 400, [])
        ch.circularize(rep_info, end_hits, [system])
        str_ = str(ch)
        self.assertEqual(str(str_), self.output_control_str('005'))

        self.unload_env("env_004")
