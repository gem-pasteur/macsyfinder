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

    def _stringify(self, so_list):
        buffer_ = os.linesep.join([str(so) for so in so_list])
        return buffer_

    # def test_build_clusters_case01(self):
    #     self.load_env("env_003", log_out=False)
    #     try:
    #         clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits,
    #                                                     self.macsy_test_env.models,
    #                                                      self.macsy_test_env.rep_info)
    #         self.maxDiff = None
    #         self.assertEqual(str(clusters), self.output_control_str('001'))
    #         self.assertEqual(len(multi_syst_genes), 0)
    #     finally:
    #         self.unload_env("env_003")
    #
    #
    # def test_build_clusters_case02(self):
    #     self.load_env("env_003", log_out=False)
    #     try:
    #         for h in self.macsy_test_env.all_hits:
    #             h.gene._multi_system = True
    #         clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits,
    #                                                     self.macsy_test_env.models,
    #                                                     self.macsy_test_env.rep_info)
    #         self.assertEqual(str(clusters), self.output_control_str('002'))
    #         self.assertEqual(len(multi_syst_genes), 1)
    #     finally:
    #         self.unload_env("env_003")
    #
    #
    # def test_get_compatible_systems(self):
    #     self.load_env("env_003", log_out=False)
    #     try:
    #         inter = get_compatible_systems([1, 2, 3], [3, 4])
    #         self.assertEqual(inter, [3])
    #     finally:
    #         self.unload_env("env_003")
    #
    #
    # def test_get_best_hits(self):
    #     self.load_env("env_003", log_out=False)
    #     try:
    #         hits = self.macsy_test_env.all_hits[0:2]
    #         hits[0].position = hits[1].position
    #         best_hits = get_best_hits(hits, criterion="score")
    #         self.assertEqual(best_hits[0].gene.name, "T9SS_sprT")
    #
    #         hits[1].i_eval = 7.5e-106
    #         best_hits = get_best_hits(hits, criterion="i_eval")
    #         self.assertEqual(best_hits[0].gene.name, "T9SS_sprT")
    #
    #         hits[0].profile_coverage = 1.92
    #         best_hits = get_best_hits(hits, criterion="profile_coverage")
    #         self.assertEqual(best_hits[0].gene.name, "T3SS_sctN")
    #     finally:
    #         self.unload_env("env_003")
    #
    #
    # def test_disambiguate_cluster_case01(self):
    #     self.load_env("env_003", log_out=False)
    #     try:
    #         rep_db = RepliconDB(self.macsy_test_env.cfg)
    #         rep_info = rep_db['AESU001c01a']
    #         clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits,
    #                                                     self.macsy_test_env.models, rep_info)
    #         cluster = clusters.clusters[1]
    #         with self.catch_log() as log:
    #             dc_clusters = disambiguate_cluster(cluster)
    #             catch_msg = log.get_value().strip()
    #
    #         self.maxDiff = None
    #         self.assertEqual(catch_msg, self.output_control_str('001'))
    #         self.assertEqual(len(dc_clusters), 1)
    #         self.assertEqual(str(dc_clusters[0]), str(cluster))
    #     finally:
    #         self.unload_env("env_003")
    #
    # def test_disambiguate_cluster_case02(self):
    #     self.load_env("env_003", log_out=False)
    #     try:
    #         rep_db = RepliconDB(self.macsy_test_env.cfg)
    #         rep_info = rep_db['AESU001c01a']
    #         clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits,
    #                                                     self.macsy_test_env.models, rep_info)
    #         cluster = clusters.clusters[1]
    #         cluster.models_to_detect = []
    #         with self.catch_log() as log:
    #             dc_clusters = disambiguate_cluster(cluster)
    #             catch_msg = log.get_value().strip()
    #
    #         self.assertEqual(catch_msg, self.output_control_str('002').strip())
    #         self.assertListEqual(dc_clusters, [])
    #     finally:
    #         self.unload_env("env_003")
    #
    #
    # def test_disambiguate_cluster_case03(self):
    #     self.load_env("env_009", log_out=False)
    #     try:
    #         rep_db = RepliconDB(self.macsy_test_env.cfg)
    #         rep_info = rep_db['VICH001.B.00001.C001']
    #         clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits,
    #                                                     self.macsy_test_env.models, rep_info)
    #         cluster = clusters.clusters[6]
    #         with self.catch_log() as log:
    #             dc_clusters = disambiguate_cluster(cluster)
    #             catch_msg = log.get_value().strip()
    #         del (cluster.hits[-1])
    #         self.assertEqual(catch_msg, self.output_control_str('003').strip())
    #         self.assertEqual(len(dc_clusters), 1)
    #         self.assertEqual(str(dc_clusters[0]), str(cluster))
    #     finally:
    #         self.unload_env("env_009")
    #
    #
    # def test_disambiguate_cluster_case04(self):
    #     self.load_env("env_009", log_out=False)
    #     try:
    #         rep_db = RepliconDB(self.macsy_test_env.cfg)
    #         rep_info = rep_db['VICH001.B.00001.C001']
    #         clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits,
    #                                                     self.macsy_test_env.models, rep_info)
    #         cluster = clusters.clusters[6]
    #         for h in cluster.hits:
    #             h.gene._loner = True
    #         li = copy.copy(cluster.hits[1:3])
    #         t4p_model = self.macsy_test_env.model_bank['set_1/T4P']
    #         li[1].gene = Gene(self.macsy_test_env.cfg, 'T4SS_MOBV', t4p_model,
    #                           self.macsy_test_env.models_location)
    #         cluster.hits.extend(li)
    #
    #         with self.catch_log() as log:
    #             dc_clusters = disambiguate_cluster(cluster)
    #             catch_msg = log.get_value().strip()
    #
    #         self.assertEqual(catch_msg, self.output_control_str('004').strip())
    #         self.assertListEqual(dc_clusters, [])
    #     finally:
    #         self.unload_env("env_009")
    #
    #
    # def test_disambiguate_cluster_case05(self):
    #     self.load_env("env_009", log_out=False)
    #     try:
    #         rep_db = RepliconDB(self.macsy_test_env.cfg)
    #         rep_info = rep_db['VICH001.B.00001.C001']
    #         clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits,
    #                                                     self.macsy_test_env.models, rep_info)
    #         cluster = clusters.clusters[6]
    #         for h in cluster.hits:
    #             h.gene._loner = True
    #         t4p_model = self.macsy_test_env.model_bank['set_1/T4P']
    #         cluster.hits[1].gene = Gene(self.macsy_test_env.cfg, 'T4SS_MOBV',
    #                                     t4p_model, self.macsy_test_env.models_location)
    #         with self.catch_log() as log:
    #             dc_clusters = disambiguate_cluster(cluster)
    #             catch_msg = log.get_value().strip()
    #         del(cluster.hits[-2:])
    #         self.assertEqual(catch_msg, self.output_control_str('005').strip())
    #         self.assertEqual(len(dc_clusters), 1)
    #         self.assertEqual(str(dc_clusters[0]), str(cluster))
    #     finally:
    #         self.unload_env("env_009")
    #
    #
    # def test_disambiguate_cluster_case06(self):
    #     self.load_env("env_009", log_out=False)
    #     try:
    #         rep_db = RepliconDB(self.macsy_test_env.cfg)
    #         rep_info = rep_db['VICH001.B.00001.C001']
    #         clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits,
    #                                                     self.macsy_test_env.models, rep_info)
    #         cluster = clusters.clusters[6]
    #         for h in cluster.hits:
    #             h.gene._loner = True
    #         li = copy.copy(cluster.hits[0:2])
    #         cluster.hits.extend(li)
    #         t4p_model = self.macsy_test_env.model_bank['set_1/T4P']
    #         cluster.hits[1].gene = Gene(self.macsy_test_env.cfg, 'T4SS_MOBV',
    #                                     t4p_model,
    #                                     self.macsy_test_env.models_location)
    #         with self.catch_log() as log:
    #             dc_clusters = disambiguate_cluster(cluster)
    #             catch_msg = log.get_value().strip()
    #
    #         self.assertEqual(catch_msg, self.output_control_str('006').strip())
    #         self.assertListEqual(dc_clusters, [])
    #     finally:
    #         self.unload_env("env_009")
    #
    #
    # def test_disambiguate_cluster_case07(self):
    #
    #     self.load_env("env_003", log_out=False)
    #     try:
    #         models_registry = ModelRegistry(self.macsy_test_env.cfg)
    #         model_name = 'set_1'
    #         models_location = models_registry[model_name]
    #
    #         rep_db = RepliconDB(self.macsy_test_env.cfg)
    #         rep_info = rep_db['AESU001c01a']
    #         clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits,
    #                                                     self.macsy_test_env.models,
    #                                                     rep_info)
    #         cluster = clusters.clusters[1]
    #         li = copy.copy(cluster.hits[0:3])
    #         cluster.hits.extend(li)
    #         t4ss_model = self.macsy_test_env.model_bank['set_1/T4SS_typeI']
    #         new_gene = Gene(self.macsy_test_env.cfg, 'T4SS_MOBV', t4ss_model, models_location)
    #         cluster.hits[4].gene = new_gene
    #         with self.catch_log() as log:
    #             dc_clusters = disambiguate_cluster(cluster)
    #             catch_msg = log.get_value().strip()
    #
    #         self.assertEqual(catch_msg, self.output_control_str('007').strip())
    #         self.assertListEqual(dc_clusters, [])
    #     finally:
    #         self.unload_env("env_003")
    #
    #
    # def test_disambiguate_cluster_case08(self):
    #
    #     self.load_env("env_009", log_out=False)
    #     try:
    #         rep_db = RepliconDB(self.macsy_test_env.cfg)
    #         rep_info = rep_db['VICH001.B.00001.C001']
    #         clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits,
    #                                                     self.macsy_test_env.models, rep_info)
    #         cluster = clusters.clusters[3]
    #         with self.catch_log() as log:
    #             dc_clusters = disambiguate_cluster(cluster)
    #             catch_msg = log.get_value().strip()
    #
    #         self.assertEqual(catch_msg, self.output_control_str('008').strip())
    #         self.assertListEqual(dc_clusters, [])
    #     finally:
    #         self.unload_env("env_009")
    #
    #
    # def test_analyze_clusters_replicon_case01(self):
    #     self.load_env("env_003", log_out=False)
    #     try:
    #         clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits,
    #                                                     self.macsy_test_env.models,
    #                                                     self.macsy_test_env.rep_info)
    #         multi_syst_genes['set_1/T9SS'] = self.macsy_test_env.all_hits[:4]
    #         systems_occurences_list = analyze_clusters_replicon(clusters,
    #                                                             self.macsy_test_env.models,
    #                                                             multi_syst_genes)
    #         str_ = self._stringify(systems_occurences_list)
    #         self.assertEqual(str_, self.output_control_str('001'))
    #     finally:
    #         self.unload_env("env_003")
    #
    #
    # def test_analyze_clusters_replicon_case02(self):
    #     self.load_env("env_004", log_out=False)
    #     try:
    #         cfg = self.macsy_test_env.cfg
    #         models_location = self.macsy_test_env.models_location
    #         t9ss_model = self.macsy_test_env.models[0]
    #
    #         t9ss_model._min_mandatory_genes_required = 1
    #         t9ss_model._min_genes_required = 1
    #
    #         rep_db = RepliconDB(cfg)
    #         rep_info = rep_db['AESU001c01a']
    #
    #         gene = Gene(cfg, 'T4SS_t4cp2', t9ss_model, models_location)
    #         t9ss_model.add_mandatory_gene(gene)
    #         gene = Gene(cfg, 'T4SS_MOBH', t9ss_model, models_location)
    #         t9ss_model.add_accessory_gene(gene)
    #
    #         clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits, [t9ss_model], rep_info)
    #         clusters.clusters = clusters.clusters[:1]  # keep only one cluster
    #
    #         multi_syst_genes['set_1/T9SS'] = self.macsy_test_env.all_hits[:4]
    #         systems_occurences_list = analyze_clusters_replicon(clusters,
    #                                                             self.macsy_test_env.models,
    #                                                             multi_syst_genes)
    #
    #         str_ = self._stringify(systems_occurences_list)
    #         self.assertEqual(str_, self.output_control_str('002'))
    #     finally:
    #         self.unload_env("env_004")
    #
    #
    # def test_analyze_clusters_replicon_case03(self):
    #     self.load_env("env_004", log_out=False)
    #     try:
    #         cfg = self.macsy_test_env.cfg
    #         models_location = self.macsy_test_env.models_location
    #         t9ss_model = self.macsy_test_env.models[0]
    #
    #         t9ss_model._min_mandatory_genes_required = 1
    #         t9ss_model._min_genes_required = 1
    #
    #         rep_db = RepliconDB(cfg)
    #         rep_info = rep_db['AESU001c01a']
    #
    #         gene = Gene(cfg, 'T4SS_t4cp2', t9ss_model, models_location)
    #         t9ss_model.add_mandatory_gene(gene)
    #         gene = Gene(cfg, 'T4SS_MOBH', t9ss_model, models_location)
    #         t9ss_model.add_accessory_gene(gene)
    #
    #         clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits,
    #                                                     [t9ss_model], rep_info)
    #         clusters.clusters = clusters.clusters[:1] # keep only one cluster
    #         cluster = clusters.clusters[0]
    #         cluster._state = "ambiguous"
    #
    #         multi_syst_genes['set_1/T9SS'] = self.macsy_test_env.all_hits[:4]
    #
    #         systems_occurences_list = analyze_clusters_replicon(clusters,
    #                                                             self.macsy_test_env.models,
    #                                                             multi_syst_genes)
    #
    #         str_ = self._stringify(systems_occurences_list)
    #         self.assertEqual(str_, self.output_control_str('003'))
    #     finally:
    #         self.unload_env("env_004")
    #
    #
    # def test_analyze_clusters_replicon_case04(self):
    #     self.load_env("env_004", log_out=False)
    #     try:
    #         cfg = self.macsy_test_env.cfg
    #         models_location = self.macsy_test_env.models_location
    #         t9ss_model = self.macsy_test_env.models[0]
    #
    #         t9ss_model._min_mandatory_genes_required = 1
    #         t9ss_model._min_genes_required = 1
    #
    #         rep_db = RepliconDB(cfg)
    #         rep_info = rep_db['AESU001c01a']
    #
    #         gene = Gene(cfg, 'T4SS_t4cp2', t9ss_model, models_location)
    #         t9ss_model.add_mandatory_gene(gene)
    #         gene = Gene(cfg, 'T4SS_MOBH', t9ss_model, models_location)
    #         t9ss_model.add_accessory_gene(gene)
    #
    #         clusters, multi_syst_genes = build_clusters(self.macsy_test_env.all_hits, [t9ss_model], rep_info)
    #         clusters.clusters = clusters.clusters[:1]  # keep only one cluster
    #         cluster = clusters.clusters[0]
    #         cluster._state = "ineligible"
    #
    #         multi_syst_genes['set_1/T9SS'] = self.macsy_test_env.all_hits[:4]
    #         systems_occurences_list = analyze_clusters_replicon(clusters, self.macsy_test_env.models, multi_syst_genes)
    #
    #         str_ = self._stringify(systems_occurences_list)
    #         self.assertEqual(str_, self.output_control_str('004'))
    #     finally:
    #         self.unload_env("env_004")
    #
    #
    # def test_search_systems_case01(self):
    #     self.load_env("env_003", log_out=False)
    #     try:
    #         tabfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.tab')
    #         reportfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.report')
    #         json_expected = self.output_control_file('results.macsyfinder.json')
    #
    #         search_systems(self.macsy_test_env.all_hits,
    #                        self.macsy_test_env.models,
    #                        self.macsy_test_env.cfg)
    #         self.assertFileEqual(tabfilename, self.output_control_file('tabfilename_001'))
    #         self.assertFileEqual(reportfilename, self.output_control_file('reportfilename_001'))
    #         json_result = os.path.join(self.macsy_test_env.cfg.working_dir(), "results.macsyfinder.json")
    #         self.assertJsonEqual(json_expected, json_result)
    #     finally:
    #         self.unload_env("env_003")
    #
    #
    # def test_search_systems_case02(self):
    #     self.load_env("env_005", log_out=False)
    #     try:
    #         tabfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.tab')
    #         reportfilename = os.path.join(self.macsy_test_env.cfg.working_dir(), 'macsyfinder.report')
    #         json_expected = self.output_control_file('results.macsyfinder.json')
    #         search_systems(self.macsy_test_env.all_hits,
    #                        self.macsy_test_env.models,
    #                        self.macsy_test_env.cfg)
    #         self.assertFileEqual(tabfilename, self.output_control_file('tabfilename_002'))
    #         self.assertFileEqual(reportfilename, self.output_control_file('reportfilename_002'))
    #         json_result = os.path.join(self.macsy_test_env.cfg.working_dir(), "results.macsyfinder.json")
    #         self.assertJsonEqual(json_expected, json_result, max_diff=None)
    #     finally:
    #         self.unload_env("env_005")
    #
    # def test_search_systems_case03(self):
    #     self.load_env("env_013", log_out=False, db_type="ordered_replicon")
    #     try:
    #         json_expected = self.output_control_file('results.macsyfinder.json')
    #         search_systems(self.macsy_test_env.all_hits,
    #                        self.macsy_test_env.models,
    #                        self.macsy_test_env.cfg)
    #         json_result = os.path.join(self.macsy_test_env.cfg.working_dir(), "results.macsyfinder.json")
    #         self.assertJsonEqual(json_expected, json_result)
    #     finally:
    #         self.unload_env("env_013")
    #
    #
    # def test_search_systems_case04(self):
    #     self.load_env("env_011", log_out=False, db_type="unordered_replicon")
    #     try:
    #         json_expected = self.output_control_file('results.macsyfinder.json')
    #         search_systems(self.macsy_test_env.all_hits,
    #                        self.macsy_test_env.models,
    #                        self.macsy_test_env.cfg)
    #         json_result = os.path.join(self.macsy_test_env.cfg.working_dir(), "results.macsyfinder.json")
    #         self.assertJsonEqual(json_expected, json_result)
    #     finally:
    #         self.unload_env("env_011")


    # def test_search_systems_case05(self):
    #     self.load_env("env_012", log_out=False, db_type="unordered_replicon")
    #     try:
    #         json_expected = self.output_control_file('results.macsyfinder.json')
    #         search_systems(self.macsy_test_env.all_hits,
    #                        self.macsy_test_env.models,
    #                        self.macsy_test_env.cfg,
    #                        self.macsy_test_env.model_bank)
    #         json_result = os.path.join(self.macsy_test_env.cfg.working_dir(), "results.macsyfinder.json")
    #         print(json_expected)
    #         self.assertJsonEqual(json_expected, json_result, max_diff=None)
    #     finally:
    #         # self.unload_env("env_012")
    #         pass


class TestSearchSystem4x4(MacsyTest, MacsyEnvManager):
    """SearchSystem systematic test (4x4 matrix based).

    4x4 test matrix is made of tests environment from env_017 to env_032.
    """

    def exec_wrapper(self):
        """Helper func."""
        json_result = os.path.join(self.macsy_test_env.cfg.working_dir(), "results.macsyfinder.json")
        search_systems(self.macsy_test_env.all_hits,
                       self.macsy_test_env.models,
                       self.macsy_test_env.cfg,
                       self.macsy_test_env.model_bank)
        return json_result

    def test_env_017(self):
        self.load_env("env_017", log_out=False)
        try:
            json_expected = self.output_control_file('results.macsyfinder.json')
            json_result = self.exec_wrapper()
            self.assertJsonEqual(json_expected, json_result)
        finally:
            self.unload_env("env_017")

    def test_env_018(self):
        self.load_env("env_018", log_out=False)
        try:
            json_expected = self.output_control_file('results.macsyfinder.json')
            json_result = self.exec_wrapper()
            self.assertJsonEqual(json_expected, json_result)
        finally:
            self.unload_env("env_018")

    def test_env_019(self):
        self.load_env("env_019", log_out=False)
        try:
            json_expected = self.output_control_file('results.macsyfinder.json')
            json_result = self.exec_wrapper()
            self.assertJsonEqual(json_expected, json_result)
        finally:
            self.unload_env("env_019")

    def test_env_020(self):
        self.load_env("env_020", log_out=False)
        try:
            json_expected = self.output_control_file('results.macsyfinder.json')
            json_result = self.exec_wrapper()
            self.assertJsonEqual(json_expected, json_result)
        finally:
            self.unload_env("env_020")

    def test_env_021(self):
        self.load_env("env_021", log_out=False)
        try:
            json_expected = self.output_control_file('results.macsyfinder.json')
            json_result = self.exec_wrapper()
            self.assertJsonEqual(json_expected, json_result)
        finally:
            self.unload_env("env_021")

    def test_env_022(self):
        self.load_env("env_022", log_out=False)
        try:
            json_expected = self.output_control_file('results.macsyfinder.json')
            json_result = self.exec_wrapper()
            self.assertJsonEqual(json_expected, json_result)
        finally:
            self.unload_env("env_022")

    def test_env_023(self):
        self.load_env("env_023", log_out=False)
        try:
            json_expected = self.output_control_file('results.macsyfinder.json')
            json_result = self.exec_wrapper()
            self.assertJsonEqual(json_expected, json_result)
        finally:
            self.unload_env("env_023")

    def test_env_024(self):
        self.load_env("env_024", log_out=False)
        try:
            json_expected = self.output_control_file('results.macsyfinder.json')
            json_result = self.exec_wrapper()
            self.assertJsonEqual(json_expected, json_result)
        finally:
            self.unload_env("env_024")

    def test_env_025(self):
        self.load_env("env_025", log_out=False)
        try:
            json_expected = self.output_control_file('results.macsyfinder.json')
            json_result = self.exec_wrapper()
            self.assertJsonEqual(json_expected, json_result)
        finally:
            self.unload_env("env_025")

    def test_env_026(self):
        self.load_env("env_026", log_out=False)
        try:
            json_expected = self.output_control_file('results.macsyfinder.json')
            json_result = self.exec_wrapper()
            self.assertJsonEqual(json_expected, json_result)
        finally:
            self.unload_env("env_026")

    def test_env_027(self):
        self.load_env("env_027", log_out=False)
        try:
            json_expected = self.output_control_file('results.macsyfinder.json')
            json_result = self.exec_wrapper()
            self.assertJsonEqual(json_expected, json_result)
        finally:
            self.unload_env("env_027")

    def test_env_028(self):
        self.load_env("env_028", log_out=False)
        try:
            json_expected = self.output_control_file('results.macsyfinder.json')
            json_result = self.exec_wrapper()
            self.assertJsonEqual(json_expected, json_result)
        finally:
            self.unload_env("env_028")

    def test_env_029(self):
        self.load_env("env_029", log_out=False)
        try:
            json_expected = self.output_control_file('results.macsyfinder.json')
            json_result = self.exec_wrapper()
            self.assertJsonEqual(json_expected, json_result)
        finally:
            self.unload_env("env_029")

    def test_env_030(self):
        self.load_env("env_030", log_out=False)
        try:
            json_expected = self.output_control_file('results.macsyfinder.json')
            json_result = self.exec_wrapper()
            self.assertJsonEqual(json_expected, json_result)
        finally:
            self.unload_env("env_030")

    def test_env_031(self):
        self.load_env("env_031", log_out=False)
        try:
            json_expected = self.output_control_file('results.macsyfinder.json')
            json_result = self.exec_wrapper()
            self.assertJsonEqual(json_expected, json_result)
        finally:
            self.unload_env("env_031")

    def test_env_032(self):
        self.load_env("env_032", log_out=False)
        try:
            json_expected = self.output_control_file('results.macsyfinder.json')
            json_result = self.exec_wrapper()
            self.assertJsonEqual(json_expected, json_result)
        finally:
            self.unload_env("env_032")

    """
    def test_all_4x4_env(self):
        for i in xrange(17,33):
            env_id = 'env_{0:03d}'.format(i)

            self.load_env(env_id, log_out=False)
            try:
                json_expected = self.output_control_file('results.macsyfinder.json')
                json_result = self.exec_wrapper()
                self.assertJsonEqual(json_expected, json_result)
            finally:
                self.unload_env(env_id)
    """
