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
import copy
from macsypy.search_systems import build_clusters, get_compatible_systems, get_best_hits, disambiguate_cluster, analyze_clusters_replicon, search_systems
from macsypy.database import RepliconDB
from macsypy.registries import ModelRegistry
from macsypy.gene import Gene
from tests import MacsyTest
from tests.macsy_test_env import MacsyEnvManager


class TestSearchSystem(MacsyTest, MacsyEnvManager):

    def setUp(self):
        pass


    def tearDown(self):

        # reset static members (hacked in test_search_systems func)
        RepliconDB.ordered_replicon_name = 'UserReplicon'


    def test_build_clusters(self):
        # case 1
        self.load_env("env_003")
        try:
            clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits,
                                                          [self.macsy_test_env.system],
                                                          self.macsy_test_env.rep_info)
            self.assertEqual(str(clusters), self.output_control_str('001'))
            self.assertEqual(len(multi_syst_genes), 0)
        finally:
            self.unload_env("env_003")

        # case 2
        self.load_env("env_003")
        try:
            for h in self.macsy_test_env.all_hits:
                h.gene._multi_system = True
            (clusters, multi_syst_genes) = build_clusters(self.macsy_test_env.all_hits,
                                                          [self.macsy_test_env.system],
                                                          self.macsy_test_env.rep_info)
            self.assertEqual(str(clusters), self.output_control_str('002'))
            self.assertEqual(len(multi_syst_genes), 1)
        finally:
            self.unload_env("env_003")


    def test_get_compatible_systems(self):
        self.load_env("env_003")
        try:
            inter = get_compatible_systems([1, 2, 3], [3, 4])
            self.assertEqual(inter, [3])
        finally:
            self.unload_env("env_003")


    def test_get_best_hits(self):
        self.load_env("env_003")
        try:
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
        finally:
            self.unload_env("env_003")

    def dc_helper(self, cluster): # 'dc' stands for Disambiguate Cluster
        """
        This method
            - calls disambiguate_cluster()
            - returns stdxxx
        """
        with self.catch_io(out=True, err=True) as stdxxx:
            clusters = disambiguate_cluster(cluster)
        stdout = stdxxx[0].getvalue()
        stderr = stdxxx[1].getvalue()
        buffer_ = os.linesep.join([stdout, stderr])
        return buffer_

    def test_disambiguate_cluster_case01(self):

        self.load_env("env_003")
        try:
            rep_db = RepliconDB(self.macsy_test_env.cfg)
            rep_info = rep_db['AESU001c01a']
            clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits,
                                                        [self.macsy_test_env.system], rep_info)
            cluster = clusters.clusters[1]
            str_ = self.dc_helper(cluster)
            self.assertEqual(str_, self.output_control_str('001'))
        finally:
            self.unload_env("env_003")

    def test_disambiguate_cluster_case02(self):

        self.load_env("env_003")
        try:
            rep_db = RepliconDB(self.macsy_test_env.cfg)
            rep_info = rep_db['AESU001c01a']
            (clusters, multi_syst_genes) = build_clusters(self.macsy_test_env.all_hits,
                                                          [self.macsy_test_env.system], rep_info)
            cluster = clusters.clusters[1]
            cluster.systems_to_detect = []
            str_ = self.dc_helper(cluster)
            self.assertEqual(str_, self.output_control_str('002'))
        finally:
            self.unload_env("env_003")


    def test_disambiguate_cluster_case03(self):

        self.load_env("env_009")
        try:
            rep_db = RepliconDB(self.macsy_test_env.cfg)
            rep_info = rep_db['VICH001.B.00001.C001']
            (clusters, multi_syst_genes) = build_clusters(self.macsy_test_env.all_hits,
                                                          [self.macsy_test_env.system], rep_info)
            cluster = clusters.clusters[6]

            str_ = self.dc_helper(cluster)
            self.assertEqual(str_, self.output_control_str('003'))
        finally:
            self.unload_env("env_009")


    def test_disambiguate_cluster_case04(self):

        self.load_env("env_009")
        try:
            rep_db = RepliconDB(self.macsy_test_env.cfg)
            rep_info = rep_db['VICH001.B.00001.C001']
            (clusters, multi_syst_genes) = build_clusters(self.macsy_test_env.all_hits,
                                                          [self.macsy_test_env.system], rep_info)
            cluster = clusters.clusters[6]

            for h in cluster.hits:
                h.gene._loner = True

            li = copy.copy(cluster.hits[1:3])
            li[1].gene = Gene(self.macsy_test_env.cfg, 'T4SS_MOBV', self.macsy_test_env.system,
                              self.macsy_test_env.models_location)
            cluster.hits.extend(li)

            str_ = self.dc_helper(cluster)
            self.assertEqual(str_, self.output_control_str('004'))
        finally:
            self.unload_env("env_009")


    def test_disambiguate_cluster_case05(self):

        self.load_env("env_009")
        try:
            rep_db = RepliconDB(self.macsy_test_env.cfg)
            rep_info = rep_db['VICH001.B.00001.C001']
            (clusters, multi_syst_genes) = build_clusters(self.macsy_test_env.all_hits,
                                                          [self.macsy_test_env.system], rep_info)
            cluster = clusters.clusters[6]

            for h in cluster.hits:
                h.gene._loner = True

            cluster.hits[1].gene = Gene(self.macsy_test_env.cfg, 'T4SS_MOBV',
                                        self.macsy_test_env.system, self.macsy_test_env.models_location)

            str_ = self.dc_helper(cluster)
            self.assertEqual(str_, self.output_control_str('005'))
        finally:
            self.unload_env("env_009")


    def test_disambiguate_cluster_case06(self):

        self.load_env("env_009")
        try:
            rep_db = RepliconDB(self.macsy_test_env.cfg)
            rep_info = rep_db['VICH001.B.00001.C001']
            clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits,
                                                        [self.macsy_test_env.system], rep_info)
            cluster = clusters.clusters[6]

            for h in cluster.hits:
                h.gene._loner = True

            li = copy.copy(cluster.hits[0:2])
            cluster.hits.extend(li)
            cluster.hits[1].gene = Gene(self.macsy_test_env.cfg, 'T4SS_MOBV',
                                        self.macsy_test_env.system, self.macsy_test_env.models_location)

            str_ = self.dc_helper(cluster)
            self.assertEqual(str_, self.output_control_str('006'))
        finally:
            self.unload_env("env_009")


    def test_disambiguate_cluster_case07(self):

        self.load_env("env_003")
        try:
            models_registry = ModelRegistry(self.macsy_test_env.cfg)
            model_name = 'set_1'
            models_location = models_registry[model_name]

            rep_db = RepliconDB(self.macsy_test_env.cfg)
            rep_info = rep_db['AESU001c01a']
            (clusters, multi_syst_genes) = build_clusters(self.macsy_test_env.all_hits,
                                                          [self.macsy_test_env.system], rep_info)
            cluster = clusters.clusters[1]

            li = copy.copy(cluster.hits[0:3])
            cluster.hits.extend(li)
            new_gene = Gene(self.macsy_test_env.cfg, 'T4SS_MOBV', self.macsy_test_env.system, models_location)
            cluster.hits[4].gene = new_gene

            str_ = self.dc_helper(cluster)
            self.assertEqual(str_, self.output_control_str('007'))
        finally:
            self.unload_env("env_003")


    def test_disambiguate_cluster_case08(self):

        self.load_env("env_009")
        try:
            rep_db = RepliconDB(self.macsy_test_env.cfg)
            rep_info = rep_db['VICH001.B.00001.C001']
            (clusters, multi_syst_genes) = build_clusters(self.macsy_test_env.all_hits,
                                                          [self.macsy_test_env.system], rep_info)
            cluster = clusters.clusters[3]

            str_ = self.dc_helper(cluster)
            self.assertEqual(str_, self.output_control_str('008'))
        finally:
            self.unload_env("env_009")


    def test_analyze_clusters_replicon(self):
        def stringify(so_list):
            buffer_ = os.linesep.join([str(so) for so in so_list])
            return buffer_

        # case 1
        self.load_env("env_003")
        try:
            (clusters, multi_syst_genes) = build_clusters(self.macsy_test_env.all_hits, [self.macsy_test_env.system], self.macsy_test_env.rep_info)
            multi_syst_genes['set_1/T9SS'] = self.macsy_test_env.all_hits[:4]
            systems_occurences_list = analyze_clusters_replicon(clusters, [self.macsy_test_env.system], multi_syst_genes)
            str_ = stringify(systems_occurences_list)
            self.assertEqual(str_, self.output_control_str('001'))
        finally:
            self.unload_env("env_003")

        # case 2
        self.load_env("env_004")
        try:
            cfg = self.macsy_test_env.cfg
            models_location = self.macsy_test_env.models_location
            system = self.macsy_test_env.system

            system._min_mandatory_genes_required = 1
            system._min_genes_required = 1

            rep_db = RepliconDB(cfg)
            rep_info = rep_db['AESU001c01a']

            gene = Gene(cfg, 'T4SS_t4cp2', system, models_location)
            system.add_mandatory_gene(gene)
            gene = Gene(cfg, 'T4SS_MOBH', system, models_location)
            system.add_accessory_gene(gene)

            (clusters, multi_syst_genes) = build_clusters(self.macsy_test_env.all_hits, [system], rep_info)
            clusters.clusters = clusters.clusters[:1] # keep only one cluster

            multi_syst_genes['set_1/T9SS'] = self.macsy_test_env.all_hits[:4]

            systems_occurences_list = analyze_clusters_replicon(clusters, [self.macsy_test_env.system], multi_syst_genes)

            str_ = stringify(systems_occurences_list)
            self.assertEqual(str_, self.output_control_str('002'))
        finally:
            self.unload_env("env_004")

        # case 3
        self.load_env("env_004")
        try:
            cfg = self.macsy_test_env.cfg
            models_location = self.macsy_test_env.models_location
            system = self.macsy_test_env.system

            system._min_mandatory_genes_required = 1
            system._min_genes_required = 1

            rep_db = RepliconDB(cfg)
            rep_info = rep_db['AESU001c01a']

            gene = Gene(cfg, 'T4SS_t4cp2', system, models_location)
            system.add_mandatory_gene(gene)
            gene = Gene(cfg, 'T4SS_MOBH', system, models_location)
            system.add_accessory_gene(gene)

            (clusters, multi_syst_genes) = build_clusters(self.macsy_test_env.all_hits, [system], rep_info)
            clusters.clusters = clusters.clusters[:1] # keep only one cluster
            cluster = clusters.clusters[0]
            cluster._state = "ambiguous"

            multi_syst_genes['set_1/T9SS'] = self.macsy_test_env.all_hits[:4]

            with self.catch_io(out=True, err=True) as stdxxx:
                systems_occurences_list = analyze_clusters_replicon(clusters, [self.macsy_test_env.system], multi_syst_genes)

            str_ = stringify(systems_occurences_list)
            self.assertEqual(str_, self.output_control_str('003'))
        finally:
            self.unload_env("env_004")

        # case 4
        self.load_env("env_004")
        try:
            cfg = self.macsy_test_env.cfg
            models_location = self.macsy_test_env.models_location
            system = self.macsy_test_env.system

            system._min_mandatory_genes_required = 1
            system._min_genes_required = 1

            rep_db = RepliconDB(cfg)
            rep_info = rep_db['AESU001c01a']

            gene = Gene(cfg, 'T4SS_t4cp2', system, models_location)
            system.add_mandatory_gene(gene)
            gene = Gene(cfg, 'T4SS_MOBH', system, models_location)
            system.add_accessory_gene(gene)

            (clusters, multi_syst_genes) = build_clusters(self.macsy_test_env.all_hits, [system], rep_info)
            clusters.clusters = clusters.clusters[:1] # keep only one cluster
            cluster = clusters.clusters[0]
            cluster._state = "ineligible"

            multi_syst_genes['set_1/T9SS'] = self.macsy_test_env.all_hits[:4]

            systems_occurences_list = analyze_clusters_replicon(clusters, [self.macsy_test_env.system], multi_syst_genes)

            str_ = stringify(systems_occurences_list)
            self.assertEqual(str_, self.output_control_str('004'))
        finally:
            self.unload_env("env_004")


    def test_search_systems_case01(self):
        self.load_env("env_003")
        try:
            tabfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.tab')
            reportfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.report')
            summaryfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.summary')

            search_systems(self.macsy_test_env.all_hits, [self.macsy_test_env.system], self.macsy_test_env.cfg)


            self.assertFileEqual(tabfilename, self.output_control_file('tabfilename_001'))
            self.assertFileEqual(reportfilename, self.output_control_file('reportfilename_001'))
            self.assertFileEqual(summaryfilename, self.output_control_file('summaryfilename_001'))
        finally:
            self.unload_env("env_003")


    def test_search_systems_case02(self):
        self.load_env("env_005")
        try:
            tabfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.tab')
            reportfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.report')
            summaryfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.summary')

            search_systems(self.macsy_test_env.all_hits, [self.macsy_test_env.system], self.macsy_test_env.cfg)

            self.assertFileEqual(tabfilename, self.output_control_file('tabfilename_002'))
            self.assertFileEqual(reportfilename, self.output_control_file('reportfilename_002'))
            self.assertFileEqual(summaryfilename, self.output_control_file('summaryfilename_002'))
        finally:
            self.unload_env("env_005")


    def test_search_systems_case03(self):
        self.load_env("env_005", db_type="ordered_replicon")
        try:
            tabfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.tab')
            reportfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.report')
            summaryfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.summary')

            RepliconDB.ordered_replicon_name = 'AESU001c01a'
            search_systems(self.macsy_test_env.all_hits, [self.macsy_test_env.system], self.macsy_test_env.cfg)

            self.assertFileEqual(tabfilename, self.output_control_file('tabfilename_003'))
            self.assertFileEqual(reportfilename, self.output_control_file('reportfilename_003'))
            self.assertFileEqual(summaryfilename, self.output_control_file('summaryfilename_003'))
        finally:
            self.unload_env("env_005")


    def test_search_systems_case04(self):
        self.load_env("env_005", db_type="unordered_replicon")
        try:
            tabfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.tab')
            reportfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.report')
            summaryfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.summary')

            RepliconDB.ordered_replicon_name = 'AESU001c01a'
            search_systems(self.macsy_test_env.all_hits, [self.macsy_test_env.system], self.macsy_test_env.cfg)

            self.assertFileEqual(tabfilename, self.output_control_file('tabfilename_004'))
            self.assertFileEqual(reportfilename, self.output_control_file('reportfilename_004'))
            self.assertFileEqual(summaryfilename, self.output_control_file('summaryfilename_004'))
        finally:
            self.unload_env("env_005")


    def test_search_systems_case05(self):
        self.load_env("env_005", db_type="unordered_replicon")
        try:
            tabfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.tab')
            reportfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.report')
            summaryfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.summary')

            RepliconDB.ordered_replicon_name = 'AESU001c01a'

            forbidden_gene = self.macsy_test_env.all_hits[0].gene
            self.macsy_test_env.system._forbidden_genes.append(forbidden_gene)
            search_systems(self.macsy_test_env.all_hits, [self.macsy_test_env.system], self.macsy_test_env.cfg)

            self.assertFileEqual(tabfilename, self.output_control_file('tabfilename_005'))
            self.assertFileEqual(reportfilename, self.output_control_file('reportfilename_005'))
            self.assertFileEqual(summaryfilename, self.output_control_file('summaryfilename_005'))
        finally:
            self.unload_env("env_005")


    def test_search_systems_case06(self):
        self.load_env("env_005", db_type="foobar")
        try:
            RepliconDB.ordered_replicon_name = 'AESU001c01a'
            with self.assertRaises(ValueError) as context:
                search_systems(self.macsy_test_env.all_hits, [self.macsy_test_env.system], self.macsy_test_env.cfg)
            self.assertEqual(str(context.exception), 'Invalid database type. ')
        finally:
            self.unload_env("env_005")
