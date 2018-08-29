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
from macsypy.search_systems import build_clusters, get_compatible_systems, get_best_hits, disambiguate_cluster, analyze_clusters_replicon, search_systems
from macsypy.database import RepliconDB
from tests import MacsyTest, md5sum
from tests.unit import MacsyTestEnv


class Test(MacsyTest):

    def setUp(self):
        self.macsy_test_env = MacsyTestEnv()

    def tearDown(self):
        self.macsy_test_env = None

        # reset static members (hacked in test_search_systems func)
        RepliconDB.ordered_replicon_name = 'UserReplicon'

    def test_build_clusters(self):
        self.macsy_test_env.load("env_003")

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

        self.macsy_test_env.unload("env_003")

    def test_get_compatible_systems(self):
        self.macsy_test_env.load("env_003")

        inter = get_compatible_systems([1, 2, 3], [3, 4])
        self.assertEqual(inter, [3])

        self.macsy_test_env.unload("env_003")

    def test_get_best_hits(self):
        self.macsy_test_env.load("env_003")

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

        self.macsy_test_env.unload("env_003")

    def test_disambiguate_cluster(self):
        self.macsy_test_env.load("env_003")

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

        self.macsy_test_env.unload("env_003")

    def test_analyze_clusters_replicon(self):
        self.macsy_test_env.load("env_003")

        rep_db = RepliconDB(self.macsy_test_env.cfg)
        rep_info = rep_db['AESU001c01a']
        (clusters, multi_syst_genes) = build_clusters(self.macsy_test_env.all_hits, [self.macsy_test_env.system], rep_info)
        systems_occurences_list = analyze_clusters_replicon(clusters, [self.macsy_test_env.system], multi_syst_genes)
        self.assertEqual(len(systems_occurences_list), 1)

        # FIXME
        # many part of analyze_clusters_replicon func not tested

        self.macsy_test_env.unload("env_003")

    def test_search_systems(self):

        # case 1

        self.macsy_test_env.load("env_003")

        tabfilename = os.path.join(self.macsy_test_env.cfg.working_dir, 'macsyfinder.tab')
        reportfilename = os.path.join(self.macsy_test_env.cfg.working_dir, 'macsyfinder.report')
        summaryfilename = os.path.join(self.macsy_test_env.cfg.working_dir, 'macsyfinder.summary')

        search_systems(self.macsy_test_env.all_hits, [self.macsy_test_env.system], self.macsy_test_env.cfg)

        # debug
        """
        with open(tabfilename) as f:
            print f.read()
        """

        self.assertEqual(md5sum(tabfilename), 'ecfb8fd9705884fff9773adb16bb22e0')
        self.assertEqual(md5sum(reportfilename), 'f6fc34319ef2e97c6f6b837fd7093709')
        self.assertEqual(md5sum(summaryfilename), '678a182c63ba693d08ebe5f263096f42')

        self.macsy_test_env.unload("env_003")

        # case 2

        self.macsy_test_env.load("env_005")

        tabfilename = os.path.join(self.macsy_test_env.cfg.working_dir, 'macsyfinder.tab')
        reportfilename = os.path.join(self.macsy_test_env.cfg.working_dir, 'macsyfinder.report')
        summaryfilename = os.path.join(self.macsy_test_env.cfg.working_dir, 'macsyfinder.summary')

        search_systems(self.macsy_test_env.all_hits, [self.macsy_test_env.system], self.macsy_test_env.cfg)

        self.assertEqual(md5sum(tabfilename), '9641c5b07036c25187d4e81300572771')
        self.assertEqual(md5sum(reportfilename), 'f6fc34319ef2e97c6f6b837fd7093709')
        self.assertEqual(md5sum(summaryfilename), '678a182c63ba693d08ebe5f263096f42')

        # case 3

        RepliconDB.ordered_replicon_name = 'AESU001c01a'

        self.macsy_test_env.cfg.options['db_type'] = "ordered_replicon"
        search_systems(self.macsy_test_env.all_hits, [self.macsy_test_env.system], self.macsy_test_env.cfg)

        self.assertEqual(md5sum(tabfilename), '650f42eca24eb2357fec9f63b5c41c05')
        self.assertEqual(md5sum(reportfilename), 'dbe75a1bab3be68b56c4cb95b18c10fc')
        self.assertEqual(md5sum(summaryfilename), '0c0f3a682d393e1408f2e00d0cacb624')

        # case 4

        self.macsy_test_env.cfg.options['db_type'] = "unordered_replicon"
        search_systems(self.macsy_test_env.all_hits, [self.macsy_test_env.system], self.macsy_test_env.cfg)

        self.assertEqual(md5sum(tabfilename), '650f42eca24eb2357fec9f63b5c41c05')
        self.assertEqual(md5sum(reportfilename), 'bba614dc69362195d5a264e2339dfe09')
        self.assertEqual(md5sum(summaryfilename), '4af8657e1877ba123783db7088c7f132')

        # case 5

        self.macsy_test_env.cfg.options['db_type'] = "unordered_replicon"
        forbidden_gene = self.macsy_test_env.all_hits[0].gene
        self.macsy_test_env.system._forbidden_genes.append(forbidden_gene)
        search_systems(self.macsy_test_env.all_hits, [self.macsy_test_env.system], self.macsy_test_env.cfg)

        self.assertEqual(md5sum(tabfilename), '650f42eca24eb2357fec9f63b5c41c05')
        self.assertEqual(md5sum(reportfilename), 'a7d305b8194e1b0b58d0bef6f127d47f')
        self.assertEqual(md5sum(summaryfilename), '4d0453c5a03d7e56fa5231a88ba2b7f6')

        # case 6

        self.macsy_test_env.cfg.options['db_type'] = "foobar"
        with self.assertRaises(ValueError) as context:
            search_systems(self.macsy_test_env.all_hits, [self.macsy_test_env.system], self.macsy_test_env.cfg)
        self.assertEqual(context.exception.message, 'Invalid database type. ')

        self.macsy_test_env.unload("env_005")
