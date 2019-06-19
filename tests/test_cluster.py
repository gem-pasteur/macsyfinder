# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Neron                                         #
# Copyright (c) 2014-2019  Institut Pasteur (Paris) and CNRS.                  #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################

import argparse
import random

from macsypy.error import MacsypyError
from macsypy.config import Config, MacsyDefaults
from macsypy.registries import ModelRegistry
from macsypy.gene import Gene, ProfileFactory
from macsypy.hit import Hit
from macsypy.model import Model
from macsypy.database import RepliconInfo
from macsypy.cluster import Cluster, build_clusters, RejectedClusters, get_loners, filter_loners
from tests import MacsyTest


class TestBuildCluster(MacsyTest):

    def setUp(self) -> None:
        self.args = argparse.Namespace()
        self.args.sequence_db = self.find_data("base", "test_base.fa")
        self.args.db_type = 'gembase'
        self.args.models_dir = self.find_data('models')
        self.args.res_search_dir = "blabla"

        self.cfg = Config(MacsyDefaults(), self.args)
        models_registry = ModelRegistry(self.cfg)
        self.model_name = 'foo'
        self.models_location = models_registry[self.model_name]
        self.profile_factory = ProfileFactory()

    def test_build_clusters(self):
        model = Model(self.cfg, "foo/T2SS", 11)
        # handle name, topology type, and min/max positions in the sequence dataset for a replicon and list of genes.
        # each genes is representing by a tuple (seq_id, length)"""
        rep_info = RepliconInfo('linear', 1, 60, [("g_{}".format(i), i * 10) for i in range(1, 7)])

        gene_1 = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location)
        gene_2 = Gene(self.cfg, self.profile_factory, "sctC", model, self.models_location)
        gene_3 = Gene(self.cfg, self.profile_factory, "sctJ", model, self.models_location)
        gene_4 = Gene(self.cfg, self.profile_factory, "sctN", model, self.models_location)
        gene_5 = Gene(self.cfg, self.profile_factory, "abc", model, self.models_location, loner=True)

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(gene_1, model, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        h11 = Hit(gene_1, model, "h11", 10, "replicon_1", 10, 1.0, 11.0, 1.0, 1.0, 10, 20)
        h20 = Hit(gene_2, model, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        h21 = Hit(gene_3, model, "h21", 10, "replicon_1", 20, 1.0, 21.0, 1.0, 1.0, 10, 20)
        h30 = Hit(gene_3, model, "h30", 10, "replicon_1", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        h31 = Hit(gene_2, model, "h31", 10, "replicon_1", 30, 1.0, 31.0, 1.0, 1.0, 10, 20)
        h50 = Hit(gene_3, model, "h50", 10, "replicon_1", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        h51 = Hit(gene_3, model, "h51", 10, "replicon_1", 50, 1.0, 51.0, 1.0, 1.0, 10, 20)
        h60 = Hit(gene_3, model, "h60", 10, "replicon_1", 60, 1.0, 60.0, 1.0, 1.0, 10, 20)
        h61 = Hit(gene_4, model, "h61", 10, "replicon_1", 60, 1.0, 61.0, 1.0, 1.0, 10, 20)

        # case replicon is linear, 2 clusters
        hits = [h10, h11, h20, h21, h30, h31, h50, h51, h60, h61]
        random.shuffle(hits)
        clusters = build_clusters(hits, rep_info, model)
        self.assertEqual(len(clusters), 2)
        self.assertListEqual(clusters[0].hits, [h11, h21, h31])
        self.assertListEqual(clusters[1].hits, [h51, h61])

        # case replicon is linear with a single hit (not loner) between 2 clusters
        h70 = Hit(gene_4, model, "h70", 10, "replicon_1", 70, 1.0, 80.0, 1.0, 1.0, 10, 20)
        h80 = Hit(gene_5, model, "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        hits = [h10, h11, h20, h21, h50, h51, h70, h80]
        random.shuffle(hits)
        clusters = build_clusters(hits, rep_info, model)
        self.assertEqual(len(clusters), 2)
        self.assertListEqual(clusters[0].hits, [h11, h21])
        self.assertListEqual(clusters[1].hits, [h70, h80])

        # replicon is linear, 3 clusters, the last one contains only one hit (loner)
        rep_info = RepliconInfo('linear', 1, 100, [("g_{}".format(i), i*10) for i in range(1, 101)])
        h80 = Hit(gene_5, model, "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        hits = [h10, h11, h20, h21, h30, h31, h50, h51, h60, h61, h80]
        random.shuffle(hits)
        clusters = build_clusters(hits, rep_info, model)
        self.assertEqual(len(clusters), 3)
        self.assertListEqual(clusters[0].hits, [h11, h21, h31])
        self.assertListEqual(clusters[1].hits, [h51, h61])
        self.assertListEqual(clusters[2].hits, [h80])

        # replicon is circular contains only one cluster
        rep_info = RepliconInfo('circular', 1, 60, [("g_{}".format(i), i*10) for i in range(1, 7)])
        hits = [h10, h20, h30]
        clusters = build_clusters(hits, rep_info, model)
        self.assertEqual(len(clusters), 1)
        self.assertListEqual(clusters[0].hits, [h10, h20, h30])

        # replicon is circular the last cluster is merge  with the first So we have only one cluster
        rep_info = RepliconInfo('circular', 1, 60, [("g_{}".format(i), i*10) for i in range(1, 7)])
        hits = [h10, h11, h20, h21, h30, h31, h50, h51, h60, h61]
        clusters = build_clusters(hits, rep_info, model)
        self.assertEqual(len(clusters), 1)
        self.assertListEqual(clusters[0].hits, [h51, h61, h11, h21, h31])

        # relicon is circular the last hit is incorporate to the first cluster
        rep_info = RepliconInfo('circular', 1, 80, [("g_{}".format(i), i*10) for i in range(1, 9)])
        h80 = Hit(gene_4, model, "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        hits = [h10, h11, h20, h21, h30, h31, h50, h51, h60, h61, h80]
        random.shuffle(hits)
        clusters = build_clusters(hits, rep_info, model)
        self.assertEqual(len(clusters), 2)
        self.assertListEqual(clusters[0].hits, [h80, h11, h21, h31])
        self.assertListEqual(clusters[1].hits, [h51, h61])

        # replicon is circular the last hit is not merged with the first cluster
        rep_info = RepliconInfo('linear', 1, 80, [("g_{}".format(i), i*10) for i in range(1, 9)])
        hits = [h10, h11, h20, h21, h30, h31, h50, h51, h60, h61, h80]
        random.shuffle(hits)
        clusters = build_clusters(hits, rep_info, model)
        self.assertEqual(len(clusters), 2)
        self.assertEqual(len(clusters), 2)
        self.assertListEqual(clusters[0].hits, [h11, h21, h31])
        self.assertListEqual(clusters[1].hits, [h51, h61])

        # case replicon is linear, 2 clusters, the hits 11,21,31 and 51,61 are contiguous
        h10 = Hit(gene_1, model, "h10", 10, "replicon_1", 10, 1.0, 11.0, 1.0, 1.0, 10, 20)
        h11 = Hit(gene_3, model, "h11", 10, "replicon_1", 11, 1.0, 21.0, 1.0, 1.0, 10, 20)
        h12 = Hit(gene_2, model, "h12", 10, "replicon_1", 12, 1.0, 31.0, 1.0, 1.0, 10, 20)
        h50 = Hit(gene_3, model, "h50", 10, "replicon_1", 50, 1.0, 51.0, 1.0, 1.0, 10, 20)
        h51 = Hit(gene_4, model, "h51", 10, "replicon_1", 51, 1.0, 61.0, 1.0, 1.0, 10, 20)
        hits = [h10, h11, h12, h50, h51]
        random.shuffle(hits)
        clusters = build_clusters(hits, rep_info, model)
        self.assertEqual(len(clusters), 2)
        self.assertListEqual(clusters[0].hits, [h10, h11, h12])
        self.assertListEqual(clusters[1].hits, [h50, h51])


class TestGetFilterLoners(MacsyTest):

    def setUp(self) -> None:
        self.args = argparse.Namespace()
        self.args.sequence_db = self.find_data("base", "test_base.fa")
        self.args.db_type = 'gembase'
        self.args.models_dir = self.find_data('models')
        self.args.res_search_dir = "blabla"

        self.cfg = Config(MacsyDefaults(), self.args)
        models_registry = ModelRegistry(self.cfg)
        self.model_name = 'foo'
        self.models_location = models_registry[self.model_name]
        self.profile_factory = ProfileFactory()

    def test_get_loners(self):
        model = Model(self.cfg, "foo/T2SS", 11)
        # handle name, topology type, and min/max positions in the sequence dataset for a replicon and list of genes.
        # each genes is representing by a tuple (seq_id, length)"""
        rep_info = RepliconInfo('linear', 1, 60, [("g_{}".format(i), i * 10) for i in range(1, 7)])

        gene_1 = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location)
        gene_2 = Gene(self.cfg, self.profile_factory, "sctC", model, self.models_location)
        gene_3 = Gene(self.cfg, self.profile_factory, "sctJ", model, self.models_location)
        gene_4 = Gene(self.cfg, self.profile_factory, "sctN", model, self.models_location, loner=True)
        gene_5 = Gene(self.cfg, self.profile_factory, "abc", model, self.models_location, loner=True)

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(gene_1, model, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        h20 = Hit(gene_2, model, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        h30 = Hit(gene_3, model, "h30", 10, "replicon_1", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        h61 = Hit(gene_4, model, "h61", 10, "replicon_1", 60, 1.0, 61.0, 1.0, 1.0, 10, 20)
        h80 = Hit(gene_5, model, "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)

        # loners are clusters of one hit
        loners = get_loners([h10, h20, h30, h61, h80], model)
        hit_from_clusters = [h.hits[0] for h in loners]
        self.assertListEqual(hit_from_clusters, [h61, h80])

    def test_filter_loners(self):
        model = Model(self.cfg, "foo/T2SS", 11)

        gene_1 = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location, loner=True)
        gene_2 = Gene(self.cfg, self.profile_factory, "sctC", model, self.models_location)
        gene_3 = Gene(self.cfg, self.profile_factory, "sctJ", model, self.models_location, loner=True)

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(gene_1, model, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        h20 = Hit(gene_2, model, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        h30 = Hit(gene_3, model, "h30", 10, "replicon_1", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        h50 = Hit(gene_3, model, "h50", 10, "replicon_1", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        c1 = Cluster([h10, h20], model)

        l10 = Cluster([h10], model)
        l30 = Cluster([h30], model)
        l50 = Cluster([h50], model)
        loners = [l10, l30, l50]
        filtered_loners = filter_loners(c1, loners)
        self.assertListEqual(filtered_loners, [l30, l50])


class TestCluster(MacsyTest):

    def setUp(self) -> None:
        self.args = argparse.Namespace()
        self.args.sequence_db = self.find_data("base", "test_base.fa")
        self.args.db_type = 'gembase'
        self.args.models_dir = self.find_data('models')
        self.args.res_search_dir = "blabla"

        self.cfg = Config(MacsyDefaults(), self.args)
        models_registry = ModelRegistry(self.cfg)
        self.model_name = 'foo'
        self.models_location = models_registry[self.model_name]
        self.profile_factory = ProfileFactory()


    def test_init(self):
        model_1 = Model(self.cfg, "foo/T2SS", 11)
        model_2 = Model(self.cfg, "foo/T3SS", 11)

        gene_1 = Gene(self.cfg, self.profile_factory, "gspD", model_1, self.models_location)
        gene_2 = Gene(self.cfg, self.profile_factory, "sctC", model_1, self.models_location)
        gene_3 = Gene(self.cfg, self.profile_factory, "sctJ", model_2, self.models_location)

        h10 = Hit(gene_1, model_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        h20 = Hit(gene_2, model_1, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        h30 = Hit(gene_3, model_2, "h30", 10, "replicon_2", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        h50 = Hit(gene_3, model_2, "h50", 10, "replicon_2", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)

        with self.assertRaises(MacsypyError) as ctx:
            with self.catch_log():
                Cluster([h10, h20, h30, h50], model_1)
        msg = "Cannot build a cluster from hits coming from different replicons"
        self.assertEqual(str(ctx.exception), msg)

    def test_contains(self):
        model = Model(self.cfg, "foo/T2SS", 11)

        gene_1 = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location)
        gene_2 = Gene(self.cfg, self.profile_factory, "sctC", model, self.models_location)
        gene_3 = Gene(self.cfg, self.profile_factory, "sctJ", model, self.models_location)

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(gene_1, model, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        h20 = Hit(gene_2, model, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        h30 = Hit(gene_3, model, "h30", 10, "replicon_1", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        h50 = Hit(gene_3, model, "h50", 10, "replicon_1", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        c1 = Cluster([h10, h20], model)

        self.assertTrue(h10 in c1)
        self.assertFalse(h30 in c1)


    def test_merge(self):
        model = Model(self.cfg, "foo/T2SS", 11)

        gene_1 = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location)
        gene_2 = Gene(self.cfg, self.profile_factory, "sctC", model, self.models_location)
        gene_3 = Gene(self.cfg, self.profile_factory, "sctJ", model, self.models_location)

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(gene_1, model, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        h20 = Hit(gene_2, model, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        h30 = Hit(gene_3, model, "h30", 10, "replicon_1", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        h50 = Hit(gene_3, model, "h50", 10, "replicon_1", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)

        c1 = Cluster([h10, h20], model)
        c2 = Cluster([h30, h50], model)
        c1.merge(c2)
        self.assertListEqual(c1.hits, [h10, h20, h30, h50])

        c1 = Cluster([h10, h20], model)
        c2 = Cluster([h30, h50], model)
        c2.merge(c1)
        self.assertListEqual(c2.hits, [h30, h50, h10, h20])

        c1 = Cluster([h10, h20], model)
        c2 = Cluster([h30, h50], model)
        c1.merge(c2, before=True)
        self.assertListEqual(c1.hits, [h30, h50, h10, h20])

        model_2 = Model(self.cfg, "foo/T3SS", 11)
        gene_3 = Gene(self.cfg, self.profile_factory, "sctJ", model_2, self.models_location)
        h30 = Hit(gene_3, model_2, "h30", 10, "replicon_2", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        h50 = Hit(gene_3, model_2, "h50", 10, "replicon_2", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        c3 = Cluster([h30, h50], model_2)
        with self.assertRaises(MacsypyError) as ctx:
            c1.merge(c3)
        self.assertEqual(str(ctx.exception), "Try to merge Clusters from different model")


    def test_str(self):
        model = Model(self.cfg, "foo/T2SS", 11)

        gene_1 = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location)
        gene_2 = Gene(self.cfg, self.profile_factory, "sctC", model, self.models_location)

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(gene_1, model, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        h20 = Hit(gene_2, model, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        c1 = Cluster([h10, h20], model)
        s ="""Cluster:
    - model: T2SS
    - hits: (h10, gspD, 10), (h20, sctC, 20)"""
        self.assertEqual(str(c1), s)


class TestRejectedCluster(MacsyTest):

    def setUp(self) -> None:
        self.args = argparse.Namespace()
        self.args.sequence_db = self.find_data("base", "test_base.fa")
        self.args.db_type = 'gembase'
        self.args.models_dir = self.find_data('models')
        self.args.res_search_dir = "blabla"

        self.cfg = Config(MacsyDefaults(), self.args)
        models_registry = ModelRegistry(self.cfg)
        self.model_name = 'foo'
        self.models_location = models_registry[self.model_name]
        self.profile_factory = ProfileFactory()


    def test_str(self):
        model = Model(self.cfg, "foo/T2SS", 11)

        gene_1 = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location)
        gene_2 = Gene(self.cfg, self.profile_factory, "sctC", model, self.models_location)

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(gene_1, model, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        h20 = Hit(gene_2, model, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        h40 = Hit(gene_1, model, "h10", 10, "replicon_1", 40, 1.0, 10.0, 1.0, 1.0, 10, 20)
        h50 = Hit(gene_2, model, "h20", 10, "replicon_1", 50, 1.0, 20.0, 1.0, 1.0, 10, 20)
        c1 = Cluster([h10, h20], model)
        c2 = Cluster([h40, h50], model)
        r_c = RejectedClusters(model, [c1, c2], "bla")

        expected_str = """Cluster:
    - model: T2SS
    - hits: (h10, gspD, 10), (h20, sctC, 20)
Cluster:
    - model: T2SS
    - hits: (h10, gspD, 40), (h20, sctC, 50)
These clusters has been rejected because:
bla"""
        self.assertEqual(expected_str, str(r_c))
