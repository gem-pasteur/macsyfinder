#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2021  Institut Pasteur (Paris) and CNRS.           #
# See the COPYRIGHT file for details                                    #
#                                                                       #
# This file is part of MacSyFinder package.                             #
#                                                                       #
# MacSyFinder is free software: you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or     #
# (at your option) any later version.                                   #
#                                                                       #
# MacSyFinder is distributed in the hope that it will be useful,        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
# GNU General Public License for more details .                         #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with MacSyFinder (COPYING).                                     #
# If not, see <https://www.gnu.org/licenses/>.                          #
#########################################################################

import os
import argparse
import random

from macsypy.error import MacsypyError
from macsypy.config import Config, MacsyDefaults
from macsypy.registries import ModelLocation
from macsypy.gene import CoreGene, ModelGene, Exchangeable, GeneStatus
from macsypy.profile import ProfileFactory
from macsypy.hit import Hit, ValidHit, HitWeight
from macsypy.model import Model
from macsypy.database import RepliconInfo
from macsypy.cluster import Cluster, build_clusters, get_loners, filter_loners
from tests import MacsyTest


class TestBuildCluster(MacsyTest):

    def setUp(self) -> None:
        self.args = argparse.Namespace()
        self.args.sequence_db = self.find_data("base", "test_1.fasta")
        self.args.db_type = 'gembase'
        self.args.models_dir = self.find_data('models')
        self.args.res_search_dir = "blabla"

        self.cfg = Config(MacsyDefaults(), self.args)
        self.model_name = 'foo'
        self.model_location = ModelLocation(path=os.path.join(self.args.models_dir, self.model_name))
        self.profile_factory = ProfileFactory(self.cfg)
        self.hit_weights = HitWeight(**self.cfg.hit_weights())


    def test_build_clusters(self):
        # handle name, topology type, and min/max positions in the sequence dataset for a replicon and list of genes.
        # each genes is representing by a tuple (seq_id, length)"""
        rep_info = RepliconInfo('linear', 1, 60, [(f"g_{i}", i * 10) for i in range(1, 7)])

        model = Model("foo/T2SS", 11)

        core_genes = []
        model_genes = []
        for g_name in ('gspD', 'sctC', 'sctJ', 'sctN', 'abc'):
            core_gene = CoreGene(self.model_location, g_name, self.profile_factory)
            core_genes.append(core_gene)
            model_genes.append(ModelGene(core_gene, model))
        model_genes[4]._loner = True

        model.add_mandatory_gene(model_genes[0])
        model.add_mandatory_gene(model_genes[1])
        model.add_accessory_gene(model_genes[2])
        model.add_accessory_gene(model_genes[3])
        model.add_neutral_gene(model_genes[4])

        #     Hit(gene, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(core_genes[0], "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        h11 = Hit(core_genes[0], "h11", 10, "replicon_1", 10, 1.0, 11.0, 1.0, 1.0, 10, 20)
        h20 = Hit(core_genes[1], "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        h21 = Hit(core_genes[2], "h21", 10, "replicon_1", 20, 1.0, 21.0, 1.0, 1.0, 10, 20)
        h30 = Hit(core_genes[2], "h30", 10, "replicon_1", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        h31 = Hit(core_genes[1], "h31", 10, "replicon_1", 30, 1.0, 31.0, 1.0, 1.0, 10, 20)
        h50 = Hit(core_genes[2], "h50", 10, "replicon_1", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        h51 = Hit(core_genes[2], "h51", 10, "replicon_1", 50, 1.0, 51.0, 1.0, 1.0, 10, 20)
        h60 = Hit(core_genes[2], "h60", 10, "replicon_1", 60, 1.0, 60.0, 1.0, 1.0, 10, 20)
        h61 = Hit(core_genes[3], "h61", 10, "replicon_1", 60, 1.0, 61.0, 1.0, 1.0, 10, 20)

        # case replicon is linear, 2 clusters
        hits = [h10, h11, h20, h21, h30, h31, h50, h51, h60, h61]
        random.shuffle(hits)
        clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(clusters), 2)
        self.assertListEqual(clusters[0].hits, [h11, h21, h31])
        self.assertListEqual(clusters[1].hits, [h51, h61])

        # case replicon is linear with a single hit (not loner) between 2 clusters
        h70 = Hit(core_genes[3], "h70", 10, "replicon_1", 70, 1.0, 80.0, 1.0, 1.0, 10, 20)
        h80 = Hit(core_genes[4], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        hits = [h10, h11, h20, h21, h50, h51, h70, h80]
        random.shuffle(hits)
        clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(clusters), 2)
        self.assertListEqual(clusters[0].hits, [h11, h21])
        self.assertListEqual(clusters[1].hits, [h70, h80])

        # replicon is linear, 3 clusters, the last one contains only one hit (loner)
        rep_info = RepliconInfo('linear', 1, 100, [(f"g_{i}", i*10) for i in range(1, 101)])
        h80 = Hit(core_genes[4], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        hits = [h10, h11, h20, h21, h30, h31, h50, h51, h60, h61, h80]
        random.shuffle(hits)
        clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(clusters), 3)
        self.assertListEqual(clusters[0].hits, [h11, h21, h31])
        self.assertListEqual(clusters[1].hits, [h51, h61])
        self.assertListEqual(clusters[2].hits, [h80])

        # replicon is circular contains only one cluster
        rep_info = RepliconInfo('circular', 1, 60, [(f"g_{i}", i*10) for i in range(1, 7)])
        hits = [h10, h20, h30]
        clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(clusters), 1)
        self.assertListEqual(clusters[0].hits, [h10, h20, h30])

        # replicon is circular the last cluster is merge  with the first So we have only one cluster
        rep_info = RepliconInfo('circular', 1, 60, [(f"g_{i}", i*10) for i in range(1, 7)])
        hits = [h10, h11, h20, h21, h30, h31, h50, h51, h60, h61]
        clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(clusters), 1)
        self.assertListEqual(clusters[0].hits, [h51, h61, h11, h21, h31])

        # replicon is circular the last hit is incorporate to the first cluster
        rep_info = RepliconInfo('circular', 1, 80, [(f"g_{i}", i*10) for i in range(1, 9)])
        h80 = Hit(core_genes[3], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        hits = [h10, h11, h20, h21, h30, h31, h50, h51, h60, h61, h80]
        random.shuffle(hits)
        clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(clusters), 2)
        self.assertListEqual(clusters[0].hits, [h80, h11, h21, h31])
        self.assertListEqual(clusters[1].hits, [h51, h61])

        # replicon is circular the last hit is not merged with the first cluster
        rep_info = RepliconInfo('linear', 1, 80, [(f"g_{i}", i*10) for i in range(1, 9)])
        hits = [h10, h11, h20, h21, h30, h31, h50, h51, h60, h61, h80]
        random.shuffle(hits)
        clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(clusters), 2)
        self.assertEqual(len(clusters), 2)
        self.assertListEqual(clusters[0].hits, [h11, h21, h31])
        self.assertListEqual(clusters[1].hits, [h51, h61])

        # case replicon is linear, 2 clusters, the hits 11,21,31 and 51,61 are contiguous
        h10 = Hit(core_genes[0], "h10", 10, "replicon_1", 10, 1.0, 11.0, 1.0, 1.0, 10, 20)
        h11 = Hit(core_genes[2], "h11", 10, "replicon_1", 11, 1.0, 21.0, 1.0, 1.0, 10, 20)
        h12 = Hit(core_genes[1], "h12", 10, "replicon_1", 12, 1.0, 31.0, 1.0, 1.0, 10, 20)
        h50 = Hit(core_genes[2], "h50", 10, "replicon_1", 50, 1.0, 51.0, 1.0, 1.0, 10, 20)
        h51 = Hit(core_genes[3], "h51", 10, "replicon_1", 51, 1.0, 61.0, 1.0, 1.0, 10, 20)
        hits = [h10, h11, h12, h50, h51]
        random.shuffle(hits)
        clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(clusters), 2)
        self.assertListEqual(clusters[0].hits, [h10, h11, h12])
        self.assertListEqual(clusters[1].hits, [h50, h51])

        # case replicon is linear
        # one cluster with one hit loner
        h80 = Hit(core_genes[4], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        hits = [h80]
        clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(clusters), 1)
        self.assertListEqual(clusters[0].hits, [h80])

        # case replicon is linear
        # one cluster with one hit min_gene_required == 1
        # last hit is a cluster
        model = Model("foo/T2SS", 11, min_mandatory_genes_required=1, min_genes_required=1)
        model.add_mandatory_gene(model_genes[0])
        model.add_accessory_gene(model_genes[1])

        h80 = Hit(core_genes[1], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        hits = [h80]
        clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(clusters), 1)
        self.assertListEqual(clusters[0].hits, [h80])

        # case replicon is linear
        # one cluster with one hit min_gene_required == 1
        # first hit is a cluster
        model = Model("foo/T2SS", 11, min_mandatory_genes_required=1, min_genes_required=1)
        model.add_mandatory_gene(model_genes[0])
        model.add_accessory_gene(model_genes[1])
        h10 = Hit(core_genes[0], "h10", 10, "replicon_1", 10, 1.0, 11.0, 1.0, 1.0, 10, 20)
        h80 = Hit(core_genes[1], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        hits = [h10, h80]
        random.shuffle(hits)
        clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(clusters), 2)
        self.assertListEqual(clusters[0].hits, [h10])
        self.assertListEqual(clusters[1].hits, [h80])


        # case replicon is linear
        # one cluster with one hit min_gene_required != 1
        model = Model("foo/T2SS", 11, min_mandatory_genes_required=1, min_genes_required=2)
        model.add_mandatory_gene(model_genes[0])
        model.add_accessory_gene(model_genes[1])
        h80 = Hit(core_genes[1], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        hits = [h80]
        clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertListEqual(clusters, [])

        # case replicon is linear, no hits
        clusters = build_clusters([], rep_info, model, self.hit_weights)
        self.assertListEqual(clusters, [])



class TestHitFunc(MacsyTest):

    def setUp(self) -> None:
        self.args = argparse.Namespace()
        self.args.sequence_db = self.find_data("base", "test_1.fasta")
        self.args.db_type = 'gembase'
        self.args.models_dir = self.find_data('models')
        self.args.res_search_dir = "blabla"

        self.cfg = Config(MacsyDefaults(), self.args)
        self.model_name = 'foo'
        self.model_location = ModelLocation(path=os.path.join(self.args.models_dir, self.model_name))
        self.profile_factory = ProfileFactory(self.cfg)
        self.hit_weights = HitWeight(**self.cfg.hit_weights())

    def test_get_loners(self):
        model = Model("foo/T2SS", 11)
        # handle name, topology type, and min/max positions in the sequence dataset for a replicon and list of genes.
        # each genes is representing by a tuple (seq_id, length)"""
        rep_info = RepliconInfo('linear', 1, 60, [(f"g_{i}", i * 10) for i in range(1, 7)])

        core_genes = []
        model_genes = []
        for g_name in ('gspD', 'sctC', 'sctJ', 'sctN', 'abc'):
            core_gene = CoreGene(self.model_location, g_name, self.profile_factory)
            core_genes.append(core_gene)
            model_genes.append(ModelGene(core_gene, model))
        model_genes[3]._loner = True
        model_genes[4]._loner = True

        model.add_mandatory_gene(model_genes[0])
        model.add_mandatory_gene(model_genes[1])
        model.add_accessory_gene(model_genes[2])
        model.add_accessory_gene(model_genes[3])
        model.add_neutral_gene(model_genes[4])

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(core_genes[0], "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        h20 = Hit(core_genes[1], "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        h30 = Hit(core_genes[2], "h30", 10, "replicon_1", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        h61 = Hit(core_genes[3], "h61", 10, "replicon_1", 60, 1.0, 61.0, 1.0, 1.0, 10, 20)
        h80 = Hit(core_genes[4], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)

        # loners are clusters of one hit
        loners = get_loners([h10, h20, h30, h61, h80], model, self.hit_weights)
        hit_from_clusters = [h.hits[0] for h in loners]
        self.assertListEqual(hit_from_clusters, [h61, h80])


    def test_filter_loners(self):
        model = Model("foo/T2SS", 11)

        core_genes = []
        model_genes = []
        for g_name in ('gspD', 'sctC', 'sctJ', 'sctN', 'abc'):
            core_gene = CoreGene(self.model_location, g_name, self.profile_factory)
            core_genes.append(core_gene)
            model_genes.append(ModelGene(core_gene, model))
        model_genes[2]._loner = True
        model_genes[3]._loner = True
        model_genes[4]._loner = True

        model.add_mandatory_gene(model_genes[0])
        model.add_mandatory_gene(model_genes[1])
        model.add_accessory_gene(model_genes[2])
        model.add_accessory_gene(model_genes[3])
        model.add_neutral_gene(model_genes[4])

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(core_genes[0], "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        h20 = Hit(core_genes[1], "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        h30 = Hit(core_genes[2], "h30", 10, "replicon_1", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        h40 = Hit(core_genes[3], "h40", 10, "replicon_1", 40, 1.0, 61.0, 1.0, 1.0, 10, 20)
        h50 = Hit(core_genes[4], "h50", 10, "replicon_1", 50, 1.0, 80.0, 1.0, 1.0, 10, 20)

        c1 = Cluster([h10, h20, h30, h40, h50], model, self.hit_weights)
        filtered_loners = filter_loners(c1, [Cluster([h30], model, self.hit_weights),
                                             Cluster([h40], model, self.hit_weights),
                                             Cluster([h50], model, self.hit_weights)]
                                        )
        self.assertListEqual(filtered_loners, [])
        c1 = Cluster([h10, h20, h40], model, self.hit_weights)
        c30 = Cluster([h30], model, self.hit_weights)
        c40 = Cluster([h40], model, self.hit_weights)
        c50 = Cluster([h50], model, self.hit_weights)
        filtered_loners = filter_loners(c1, [c30, c40, c50])
        self.assertListEqual(filtered_loners, [c30, c50])


class TestCluster(MacsyTest):

    def setUp(self) -> None:
        self.args = argparse.Namespace()
        self.args.sequence_db = self.find_data("base", "test_1.fasta")
        self.args.db_type = 'gembase'
        self.args.models_dir = self.find_data('models')
        self.args.res_search_dir = "blabla"

        self.cfg = Config(MacsyDefaults(), self.args)
        self.model_name = 'foo'
        self.model_location = ModelLocation(path=os.path.join(self.args.models_dir, self.model_name))
        self.profile_factory = ProfileFactory(self.cfg)
        self.hit_weights = HitWeight(**self.cfg.hit_weights())


    def test_init(self):
        model_1 = Model("foo/T2SS", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)
        c_gene_3 = CoreGene(self.model_location, "sctJ", self.profile_factory)

        gene_1 = ModelGene(c_gene_1, model_1)

        h10 = Hit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        v_h10 = ValidHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = Hit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        v_h20 = ValidHit(h20, gene_1, GeneStatus.MANDATORY)
        h30 = Hit(c_gene_3, "h30", 10, "replicon_2", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        v_h30 = ValidHit(h30, gene_1, GeneStatus.ACCESSORY)
        h50 = Hit(c_gene_3, "h50", 10, "replicon_2", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        v_h50 = ValidHit(h50, gene_1, GeneStatus.ACCESSORY)

        with self.assertRaises(MacsypyError) as ctx:
            with self.catch_log():
                Cluster([v_h10, v_h20, v_h30, v_h50], model_1, self.hit_weights)
        msg = "Cannot build a cluster from hits coming from different replicons"
        self.assertEqual(str(ctx.exception), msg)


    def test_replicon_name(self):
        model = Model("foo/T2SS", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)

        gene_1 = ModelGene(c_gene_1, model)
        gene_2 = ModelGene(c_gene_2, model)

        replicon_name = "replicon_1"
        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(c_gene_1, "h10", 10, replicon_name, 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        v_h10 = ValidHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = Hit(c_gene_2, "h20", 10, replicon_name, 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        v_h20 = ValidHit(h20, gene_2, GeneStatus.MANDATORY)

        c1 = Cluster([v_h10, v_h20], model, self.hit_weights)
        self.assertEqual(c1.replicon_name, replicon_name)


    def test_len(self):
        model = Model("foo/T2SS", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)

        gene_1 = ModelGene(c_gene_1, model)
        gene_2 = ModelGene(c_gene_2, model)

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        v_h10 = ValidHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = Hit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        v_h20 = ValidHit(h20, gene_2, GeneStatus.MANDATORY)

        c1 = Cluster([v_h10, v_h20], model, self.hit_weights)
        self.assertEqual(len(c1), 2)

    def test_loner(self):
        model = Model("foo/bar", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)

        gene_1 = ModelGene(c_gene_1, model)
        gene_2 = ModelGene(c_gene_2, model)

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        v_h10 = ValidHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = Hit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        v_h20 = ValidHit(h20, gene_2, GeneStatus.MANDATORY)

        c1 = Cluster([v_h10], model, self.hit_weights)
        c2 = Cluster([v_h10, v_h20], model, self.hit_weights)
        self.assertTrue(c1.loner())
        self.assertFalse(c2.loner())

    def test_contains(self):
        model = Model("foo/T2SS", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)
        c_gene_3 = CoreGene(self.model_location, "sctJ", self.profile_factory)

        gene_1 = ModelGene(c_gene_1, model)
        gene_2 = ModelGene(c_gene_2, model)
        gene_3 = ModelGene(c_gene_3, model)

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        v_h10 = ValidHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = Hit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        v_h20 = ValidHit(h20, gene_2, GeneStatus.MANDATORY)
        h30 = Hit(c_gene_3, "h30", 10, "replicon_1", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        v_h30 = ValidHit(h30, gene_3, GeneStatus.ACCESSORY)
        h50 = Hit(c_gene_3, "h50", 10, "replicon_1", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        v_h50 = ValidHit(h50, gene_3, GeneStatus.ACCESSORY)
        c1 = Cluster([v_h10, v_h20, v_h50], model, self.hit_weights)

        self.assertTrue(v_h10 in c1)
        self.assertFalse(v_h30 in c1)


    def test_fulfilled_function(self):
        model = Model("foo/T2SS", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)
        c_gene_3 = CoreGene(self.model_location, "sctJ", self.profile_factory)
        c_gene_4 = CoreGene(self.model_location, "sctJ_FLG", self.profile_factory)

        gene_1 = ModelGene(c_gene_1, model)
        gene_2 = ModelGene(c_gene_2, model)
        gene_3 = ModelGene(c_gene_3, model)
        gene_4 = Exchangeable(c_gene_4, gene_3)
        gene_3.add_exchangeable(gene_4)

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        v_h10 = ValidHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = Hit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        v_h20 = ValidHit(h20, gene_2, GeneStatus.MANDATORY)

        c = Cluster([v_h10, v_h20], model, self.hit_weights)

        self.assertTrue(c.fulfilled_function(gene_1))
        self.assertFalse(c.fulfilled_function(gene_3))

        h50 = Hit(c_gene_4, "h50", 10, "replicon_1", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        v_h50 = ValidHit(h50, gene_4, GeneStatus.ACCESSORY)
        c = Cluster([v_h10, v_h50], model, self.hit_weights)
        self.assertTrue(c.fulfilled_function(gene_3))

    def test_score(self):
        model = Model("foo/T2SS", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)

        c_gene_tadZ = CoreGene(self.model_location, "tadZ", self.profile_factory)
        gene_tadZ = ModelGene(c_gene_tadZ, model)
        model.add_mandatory_gene(gene_tadZ)

        c_gene_sctj = CoreGene(self.model_location, "sctC", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)

        c_gene_sctJ_FLG = CoreGene(self.model_location, "sctJ_FLG", self.profile_factory)

        analog_sctJ_FLG = Exchangeable(c_gene_sctJ_FLG, gene_sctj)
        gene_sctj.add_exchangeable(analog_sctJ_FLG)
        model.add_accessory_gene(gene_sctj)

        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model, loner=True)
        c_gene_sctn_FLG = CoreGene(self.model_location, "sctN_FLG", self.profile_factory)
        homolog_sctn_FLG = Exchangeable(c_gene_sctn_FLG, gene_sctn)
        gene_sctn.add_exchangeable(homolog_sctn_FLG)
        model.add_accessory_gene(gene_sctn)

        c_gene_toto = CoreGene(self.model_location, "toto", self.profile_factory)
        gene_toto = ModelGene(c_gene_toto, model)
        model.add_neutral_gene(gene_toto)

        c_gene_flie = CoreGene(self.model_location, "fliE", self.profile_factory)
        gene_flie = ModelGene(c_gene_flie, model, loner=True, multi_system=True)
        model.add_mandatory_gene(gene_flie)

        h_gspd = Hit(c_gene_gspd, "h_gspd", 10, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_h_gspd = ValidHit(h_gspd, gene_gspd, GeneStatus.MANDATORY)
        h_tadz = Hit(c_gene_tadZ, "h_tadz", 20, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_h_tadz = ValidHit(h_tadz, gene_tadZ, GeneStatus.MANDATORY)

        h_sctj = Hit(c_gene_sctj, "h_sctj", 30, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_h_sctj = ValidHit(h_sctj, gene_sctj, GeneStatus.ACCESSORY)
        h_sctj_an = Hit(c_gene_sctJ_FLG, "h_sctj_an", 30, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_h_sctj_an = ValidHit(h_sctj_an, analog_sctJ_FLG, GeneStatus.ACCESSORY)

        h_sctn = Hit(c_gene_sctn, "sctn", 40, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_h_sctn = ValidHit(h_sctn, gene_sctn, GeneStatus.ACCESSORY)
        h_sctn_hom = Hit(c_gene_sctn_FLG, "h_scth_hom", 30, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_h_sctn_hom = ValidHit(h_sctn_hom, homolog_sctn_FLG, GeneStatus.ACCESSORY)

        h_toto = Hit(c_gene_sctn, "toto", 50, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_h_toto = ValidHit(h_toto, gene_toto, GeneStatus.NEUTRAL)

        h_flie = Hit(c_gene_flie, "h_flie", 100, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_h_flie = ValidHit(h_flie, gene_flie, GeneStatus.MANDATORY)

        # 2 mandatory, 2 accessory no analog/homolog
        c1 = Cluster([v_h_gspd, v_h_tadz, v_h_sctj, v_h_sctn], model, self.hit_weights)
        self.assertEqual(c1.score, 3.0)

        # 2 mandatory, 2 accessory 1 neutral, no analog/homolog
        c1 = Cluster([v_h_gspd, v_h_tadz, v_h_sctj, v_h_sctn, v_h_toto], model, self.hit_weights)
        self.assertEqual(c1.score, 3.0)

        # 1 mandatory + 1 mandatory duplicated 1 time
        # 1 accessory + 1 accessory duplicated 1 times
        # no analog/homolog
        c1 = Cluster([v_h_gspd, v_h_tadz, v_h_sctj, v_h_sctn, v_h_gspd, v_h_sctn], model, self.hit_weights)
        self.assertEqual(c1.score, 3.0)

        # 2 mandatory
        # 1 accessory + 1 accessory homolog
        c1 = Cluster([v_h_gspd, v_h_tadz, v_h_sctj, v_h_sctn_hom], model, self.hit_weights)
        self.assertEqual(c1.score, 2.9)

        # # 2 mandatory
        # # 1 accessory + 1 accessory analog of the 1rst accessory
        # c1 = Cluster([v_h_gspd, v_h_tadz, v_h_sctj, v_h_sctj_an], model, self.hit_weights)
        # self.assertEqual(c1.score, 2.5)

        # test loners multi system
        c1 = Cluster([v_h_flie], model, self.hit_weights)
        self.assertEqual(c1.score, 0.7)

        # test the cache score
        self.assertEqual(c1.score, 0.7)

        non_valid_hit = ValidHit(h_sctn, gene_sctn, GeneStatus.FORBIDDEN)
        c1 = Cluster([v_h_gspd, non_valid_hit, v_h_tadz], model, self.hit_weights)
        with self.assertRaises(MacsypyError) as ctx:
            c1.score
        self.assertEqual(str(ctx.exception),
                         "a Cluster contains hit which is neither mandatory nor accessory")


    def test_merge(self):
        model = Model("foo/T2SS", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)
        c_gene_3 = CoreGene(self.model_location, "sctJ", self.profile_factory)

        gene_1 = ModelGene(c_gene_1, model)
        gene_2 = ModelGene(c_gene_2, model)
        gene_3 = ModelGene(c_gene_3, model)

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        v_h10 = ValidHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = Hit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        v_h20 = ValidHit(h20, gene_2, GeneStatus.MANDATORY)
        h30 = Hit(c_gene_3, "h30", 10, "replicon_1", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        v_h30 = ValidHit(h30, gene_3, GeneStatus.ACCESSORY)
        h50 = Hit(c_gene_3, "h50", 10, "replicon_1", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        v_h50 = ValidHit(h50, gene_3, GeneStatus.ACCESSORY)

        c1 = Cluster([v_h10, v_h20], model, self.hit_weights)
        c2 = Cluster([v_h30, v_h50], model, self.hit_weights)
        c1.merge(c2)
        self.assertListEqual(c1.hits, [v_h10, v_h20, v_h30, v_h50])

        c1 = Cluster([v_h10, v_h20], model, self.hit_weights)
        c2 = Cluster([v_h30, v_h50], model, self.hit_weights)
        c2.merge(c1)
        self.assertListEqual(c2.hits, [v_h30, v_h50, v_h10, v_h20])

        c1 = Cluster([v_h10, v_h20], model, self.hit_weights)
        c2 = Cluster([v_h30, v_h50], model, self.hit_weights)
        c1.merge(c2, before=True)
        self.assertListEqual(c1.hits, [v_h30, v_h50, v_h10, v_h20])

        model_2 = Model("foo/T3SS", 11)
        c_gene_3 = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_3 = ModelGene(c_gene_3, model)

        h30 = Hit(c_gene_3, "h30", 10, "replicon_2", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        v_h30 = ValidHit(h30, gene_3, GeneStatus.ACCESSORY)
        h50 = Hit(c_gene_3, "h50", 10, "replicon_2", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        v_h50 = ValidHit(h50, gene_3, GeneStatus.ACCESSORY)
        c3 = Cluster([v_h30, v_h50], model_2, self.hit_weights)
        with self.assertRaises(MacsypyError) as ctx:
            c1.merge(c3)
        self.assertEqual(str(ctx.exception), "Try to merge Clusters from different model")


    def test_str(self):
        model = Model("foo/T2SS", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)

        gene_1 = ModelGene(c_gene_1, model)
        gene_2 = ModelGene(c_gene_2, model)

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        v_h10 = ValidHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = Hit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        v_h20 = ValidHit(h20, gene_2, GeneStatus.MANDATORY)
        c1 = Cluster([v_h10, v_h20], model, self.hit_weights)
        s ="""Cluster:
- model = T2SS
- replicon = replicon_1
- hits = (h10, gspD, 10), (h20, sctC, 20)"""
        self.assertEqual(str(c1), s)
