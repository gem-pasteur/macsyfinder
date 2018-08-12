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


import shutil
import tempfile
from macsypy.search_systems import build_clusters, get_compatible_systems, get_best_hits, disambiguate_cluster, analyze_clusters_replicon
from macsypy.database import RepliconDB
from tests import MacsyTest, md5sum
from tests.unit import MacsyTestEnv


class Test(MacsyTest):

    def setUp(self):
        self.macsy_test_env = MacsyTestEnv()
        self.macsy_test_env.load("env_003")

        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        self.macsy_test_env.unload("env_003")

        shutil.rmtree(self.test_dir)

    def test_build_clusters(self):

        # case 1

        (clusters, multi_syst_genes) = build_clusters(self.macsy_test_env.all_hits,
                                                      [self.macsy_test_env.system],
                                                      self.macsy_test_env.rep_info)
        self.assertEqual(md5sum(str_=str(clusters)), '5e8f78636b7d480ce6f6b6db949755ce')
        self.assertEqual(len(multi_syst_genes), 0)

        # case 2

        for h in self.macsy_test_env.all_hits:
            h.gene._multi_system = True
        (clusters, multi_syst_genes) = build_clusters(self.macsy_test_env.all_hits,
                                                      [self.macsy_test_env.system],
                                                      self.macsy_test_env.rep_info)
        self.assertEqual(md5sum(str_=str(clusters)), '5e8f78636b7d480ce6f6b6db949755ce')
        self.assertEqual(len(multi_syst_genes), 1)

        # FIXME: parts below are not tested
        """
        if positions.count(first.position) == 0:
            hitstoconsider.append(first)

        if len(cur_cluster) > 1 or (len(cur_cluster) == 1 and prev.gene.loner):
            #print "Recap clusters"
        """

    def test_get_compatible_systems(self):
        inter = get_compatible_systems([1, 2, 3], [3, 4])
        self.assertEqual(inter, [3])

    def test_get_best_hits(self):
        hits = self.macsy_test_env.all_hits[0:2]

        hits[0].position = hits[1].position

        # debug
        """
        for h in hits:
            print h.position
            print h.i_eval
            print h.score
            print h.profile_coverage
            print h.gene.name
        """

        best_hits = get_best_hits(hits, criterion="score")
        self.assertEqual(best_hits[0].gene.name, "T9SS_sprT")

        hits[1].i_eval = 7.5e-106
        best_hits = get_best_hits(hits, criterion="i_eval")
        self.assertEqual(best_hits[0].gene.name, "T9SS_gldN_TIGR03523")

        hits[0].profile_coverage = 0.92
        best_hits = get_best_hits(hits, criterion="profile_coverage")
        self.assertEqual(best_hits[0].gene.name, "T9SS_gldN_TIGR03523")

    def test_disambiguate_cluster(self):
        rep_db = RepliconDB(self.macsy_test_env.cfg)
        rep_info = rep_db['AESU001c01a']
        (clusters, multi_syst_genes) = build_clusters(self.macsy_test_env.all_hits, [self.macsy_test_env.system], rep_info)
        cluster = clusters.clusters[1]

        with self.catch_io(out=True, err=False) as stdxxx:
            clusters = disambiguate_cluster(cluster)

        stdout = stdxxx[0].getvalue()

        self.assertEqual(md5sum(str_=stdout), '2567ef1661050233b190073a57400d42')

        # FIXME
        # in disambiguate_cluster func, block starting with comment below is
        # not tested
        """
        Deal with "accessory foreign genes",
        """

    def test_analyze_clusters_replicon(self):
        rep_db = RepliconDB(self.macsy_test_env.cfg)
        rep_info = rep_db['AESU001c01a']
        (clusters, multi_syst_genes) = build_clusters(self.macsy_test_env.all_hits, [self.macsy_test_env.system], rep_info)
        systems_occurences_list = analyze_clusters_replicon(clusters, [self.macsy_test_env.system], multi_syst_genes)
        self.assertEqual(len(systems_occurences_list), 1)

        # FIXME
        # many part of analyze_clusters_replicon func not tested
