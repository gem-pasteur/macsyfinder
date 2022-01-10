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
from macsypy.hit import CoreHit, ModelHit, Loner, MultiSystem, LonerMultiSystem, HitWeight
from macsypy.model import Model
from macsypy.database import RepliconInfo
from macsypy.cluster import Cluster, build_clusters, _colocates, _clusterize, _get_multi_system_hits, _get_true_loners
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

        #              fqn      , inter_gene_max_sapce
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

        #     CoreHit(gene, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        #                                                     pos     score
        h10 = CoreHit(core_genes[0], "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        m_h10 = ModelHit(h10, gene_ref=model_genes[0], gene_status=GeneStatus.MANDATORY)
        h11 = CoreHit(core_genes[0], "h11", 10, "replicon_1", 10, 1.0, 11.0, 1.0, 1.0, 10, 20)
        m_h11 = ModelHit(h11, gene_ref=model_genes[0], gene_status=GeneStatus.MANDATORY)

        h20 = CoreHit(core_genes[1], "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        m_h20 = ModelHit(h20, gene_ref=model_genes[1], gene_status=GeneStatus.MANDATORY)
        h21 = CoreHit(core_genes[2], "h21", 10, "replicon_1", 20, 1.0, 21.0, 1.0, 1.0, 10, 20)
        m_h21 = ModelHit(h21, gene_ref=model_genes[1], gene_status=GeneStatus.MANDATORY)

        h30 = CoreHit(core_genes[2], "h30", 10, "replicon_1", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        m_h30 = ModelHit(h30, gene_ref=model_genes[2], gene_status=GeneStatus.ACCESSORY)
        h31 = CoreHit(core_genes[1], "h31", 10, "replicon_1", 30, 1.0, 31.0, 1.0, 1.0, 10, 20)
        m_h31 = ModelHit(h31, gene_ref=model_genes[1], gene_status=GeneStatus.MANDATORY)

        h50 = CoreHit(core_genes[2], "h50", 10, "replicon_1", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        m_h50 = ModelHit(h50, gene_ref=model_genes[2], gene_status=GeneStatus.ACCESSORY)
        h51 = CoreHit(core_genes[2], "h51", 10, "replicon_1", 50, 1.0, 51.0, 1.0, 1.0, 10, 20)
        m_h51 = ModelHit(h51, gene_ref=model_genes[2], gene_status=GeneStatus.ACCESSORY)

        h60 = CoreHit(core_genes[2], "h60", 10, "replicon_1", 60, 1.0, 60.0, 1.0, 1.0, 10, 20)
        m_h60 = ModelHit(h60, gene_ref=model_genes[2], gene_status=GeneStatus.ACCESSORY)
        h61 = CoreHit(core_genes[3], "h61", 10, "replicon_1", 60, 1.0, 61.0, 1.0, 1.0, 10, 20)
        m_h61 = ModelHit(h61, gene_ref=model_genes[3], gene_status=GeneStatus.ACCESSORY)

        # case replicon is linear, 2 clusters
        hits = [m_h10, m_h11, m_h20, m_h21, m_h30, m_h31, m_h50, m_h51, m_h60, m_h61]
        random.shuffle(hits)
        true_clusters, special_clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(true_clusters), 2)
        self.assertListEqual(true_clusters[0].hits, [m_h11, m_h21, m_h31])
        self.assertListEqual(true_clusters[1].hits, [m_h51, m_h61])
        self.assertEqual(special_clusters, {})

        # case replicon is linear with a single hit (not loner) between 2 clusters
        h70 = CoreHit(core_genes[3], "h70", 10, "replicon_1", 70, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h70 = ModelHit(h70, gene_ref=model_genes[3], gene_status=GeneStatus.ACCESSORY)
        h80 = CoreHit(core_genes[4], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[4], gene_status=GeneStatus.NEUTRAL)
        hits = [m_h10, m_h11, m_h20, m_h21, m_h50, m_h51, m_h70, m_h80]
        random.shuffle(hits)
        true_clusters, special_clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(true_clusters), 2)
        self.assertListEqual(true_clusters[0].hits, [m_h11, m_h21])
        self.assertListEqual(true_clusters[1].hits, [m_h70, m_h80])
        self.assertEqual(special_clusters, {})

        # replicon is linear, 3 clusters, the last one contains only one hit (loner h80)
        rep_info = RepliconInfo('linear', 1, 100, [(f"g_{i}", i*10) for i in range(1, 101)])
        h80 = CoreHit(core_genes[4], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[4], gene_status=GeneStatus.NEUTRAL)
        hits = [m_h10, m_h11, m_h20, m_h21, m_h30, m_h31, m_h50, m_h51, m_h60, m_h61, m_h80]
        random.shuffle(hits)
        true_clusters, special_clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(true_clusters), 2)
        self.assertListEqual(true_clusters[0].hits, [m_h11, m_h21, m_h31])
        self.assertListEqual(true_clusters[1].hits, [m_h51, m_h61])
        self.assertEqual(len(special_clusters), 1)
        self.assertListEqual(special_clusters['abc'].hits, [m_h80])

        # replicon is circular contains only one cluster
        rep_info = RepliconInfo('circular', 1, 60, [(f"g_{i}", i*10) for i in range(1, 7)])
        hits = [m_h10, m_h20, m_h30]
        true_clusters, special_clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(true_clusters), 1)
        self.assertListEqual(true_clusters[0].hits, [m_h10, m_h20, m_h30])
        self.assertEqual(special_clusters, {})

        # replicon is circular the last cluster is merge  with the first So we have only one cluster
        rep_info = RepliconInfo('circular', 1, 60, [(f"g_{i}", i*10) for i in range(1, 7)])
        hits = [m_h10, m_h11, m_h20, m_h21, m_h30, m_h31, m_h50, m_h51, m_h60, m_h61]
        true_clusters, special_clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(true_clusters), 1)
        self.assertListEqual(true_clusters[0].hits, [m_h51, m_h61, m_h11, m_h21, m_h31])
        self.assertEqual(special_clusters, {})

        # replicon is circular the last hit is incorporate to the first cluster
        # m_h80 is not considered as loner it is included in a cluster
        rep_info = RepliconInfo('circular', 1, 80, [(f"g_{i}", i*10) for i in range(1, 9)])
        h80 = CoreHit(core_genes[3], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[4], gene_status=GeneStatus.NEUTRAL)
        hits = [m_h10, m_h11, m_h20, m_h21, m_h30, m_h31, m_h50, m_h51, m_h60, m_h61, m_h80]
        random.shuffle(hits)
        true_clusters, special_clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(true_clusters), 2)
        self.assertListEqual(true_clusters[0].hits, [m_h80, m_h11, m_h21, m_h31])
        self.assertListEqual(true_clusters[1].hits, [m_h51, m_h61])
        self.assertEqual(special_clusters, {})

        # replicon is linear the last hit is not merged with the first cluster
        # m_h80 is link to gene_4 'abc'. So, in this test, it's not a loner.
        rep_info = RepliconInfo('linear', 1, 62, [(f"g_{i}", i*10) for i in range(1, 7)])
        hits = [m_h10, m_h11, m_h20, m_h21, m_h30, m_h31, m_h50, m_h51, m_h60, m_h61]
        random.shuffle(hits)
        true_clusters, special_clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(true_clusters), 2)
        self.assertListEqual(true_clusters[0].hits, [m_h11, m_h21, m_h31])
        self.assertListEqual(true_clusters[1].hits, [m_h51, m_h61])
        self.assertEqual(special_clusters, {})

        # case replicon is linear, 2 clusters, the hits 11,21,31 and 51,61 are contiguous
        #                                                              pos
        h10 = CoreHit(core_genes[0], "h10", 10, "replicon_1", 10, 1.0, 11.0, 1.0, 1.0, 10, 20)
        m_h10 = ModelHit(h10, gene_ref=model_genes[0], gene_status=GeneStatus.MANDATORY)
        h11 = CoreHit(core_genes[2], "h11", 10, "replicon_1", 11, 1.0, 21.0, 1.0, 1.0, 10, 20)
        m_h11 = ModelHit(h11, gene_ref=model_genes[2], gene_status=GeneStatus.ACCESSORY)
        h12 = CoreHit(core_genes[1], "h12", 10, "replicon_1", 12, 1.0, 31.0, 1.0, 1.0, 10, 20)
        m_h12 = ModelHit(h12, gene_ref=model_genes[1], gene_status=GeneStatus.MANDATORY)
        h50 = CoreHit(core_genes[2], "h50", 10, "replicon_1", 50, 1.0, 51.0, 1.0, 1.0, 10, 20)
        m_h50 = ModelHit(h50, gene_ref=model_genes[2], gene_status=GeneStatus.ACCESSORY)
        h51 = CoreHit(core_genes[3], "h51", 10, "replicon_1", 51, 1.0, 61.0, 1.0, 1.0, 10, 20)
        m_h51 = ModelHit(h51, gene_ref=model_genes[3], gene_status=GeneStatus.ACCESSORY)
        hits = [m_h10, m_h11, m_h12, m_h50, m_h51]
        random.shuffle(hits)
        true_clusters, special_clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(true_clusters), 2)
        self.assertListEqual(true_clusters[0].hits, [m_h10, m_h11, m_h12])
        self.assertListEqual(true_clusters[1].hits, [m_h50, m_h51])
        self.assertEqual(special_clusters, {})

        # case replicon is linear
        # one cluster with one hit loner
        h80 = CoreHit(core_genes[4], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[4], gene_status=GeneStatus.NEUTRAL)
        hits = [m_h80]
        true_clusters, special_clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(true_clusters), 0)
        self.assertListEqual(special_clusters['abc'].hits, [m_h80])

        # case replicon is linear
        # one cluster with one hit min_gene_required == 1
        # last hit is a cluster
        model = Model("foo/T2SS", 11, min_mandatory_genes_required=1, min_genes_required=1)
        core_genes = []
        model_genes = []
        for g_name in ('gspD', 'sctC'):
            core_gene = CoreGene(self.model_location, g_name, self.profile_factory)
            core_genes.append(core_gene)
            model_genes.append(ModelGene(core_gene, model))
        model.add_mandatory_gene(model_genes[0])
        model.add_accessory_gene(model_genes[1])

        h80 = CoreHit(core_genes[1], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[1], gene_status=GeneStatus.ACCESSORY)
        hits = [m_h80]
        true_clusters, special_clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(true_clusters), 1)
        self.assertListEqual(true_clusters[0].hits, [m_h80])
        self.assertEqual(special_clusters, {})

        # case replicon is linear
        # one cluster with one hit min_gene_required == 1
        # first hit is a cluster
        model = Model("foo/T2SS", 11, min_mandatory_genes_required=1, min_genes_required=1)
        core_genes = []
        model_genes = []
        for g_name in ('gspD', 'sctC'):
            core_gene = CoreGene(self.model_location, g_name, self.profile_factory)
            core_genes.append(core_gene)
            model_genes.append(ModelGene(core_gene, model))
        model.add_mandatory_gene(model_genes[0])
        model.add_accessory_gene(model_genes[1])

        h10 = CoreHit(core_genes[0], "h10", 10, "replicon_1", 10, 1.0, 11.0, 1.0, 1.0, 10, 20)
        m_h10 = ModelHit(h10, gene_ref=model_genes[0], gene_status=GeneStatus.MANDATORY)
        h80 = CoreHit(core_genes[1], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[1], gene_status=GeneStatus.ACCESSORY)
        hits = [m_h10, m_h80]
        random.shuffle(hits)
        true_clusters, special_clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertEqual(len(true_clusters), 2)
        self.assertListEqual(true_clusters[0].hits, [m_h10])
        self.assertListEqual(true_clusters[1].hits, [m_h80])
        self.assertEqual(special_clusters, {})

        # case replicon is linear
        # one cluster with one hit min_gene_required != 1
        model = Model("foo/T2SS", 11, min_mandatory_genes_required=1, min_genes_required=2)
        core_genes = []
        model_genes = []
        for g_name in ('gspD', 'sctC'):
            core_gene = CoreGene(self.model_location, g_name, self.profile_factory)
            core_genes.append(core_gene)
            model_genes.append(ModelGene(core_gene, model))
        model.add_mandatory_gene(model_genes[0])
        model.add_accessory_gene(model_genes[1])
        h80 = CoreHit(core_genes[1], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[1], gene_status=GeneStatus.ACCESSORY)
        hits = [m_h80]
        true_clusters, special_clusters = build_clusters(hits, rep_info, model, self.hit_weights)
        self.assertListEqual(true_clusters, [])
        self.assertEqual(special_clusters, {})

        # case replicon is linear, no hits
        true_clusters, special_clusters = build_clusters([], rep_info, model, self.hit_weights)
        self.assertListEqual(true_clusters, [])
        self.assertEqual(special_clusters, {})


    def test_colocates(self):
        rep_info = RepliconInfo('linear', 1, 60, [(f"g_{i}", i * 10) for i in range(1, 7)])

        #              fqn      , inter_gene_max_sapce
        model = Model("foo/T2SS", 7)

        core_genes = {}
        model_genes = {}
        for g_name in ('gspD', 'sctC', 'sctJ', 'sctN'):
            core_gene = CoreGene(self.model_location, g_name, self.profile_factory)
            core_genes[g_name] = core_gene
            model_genes[g_name] = ModelGene(core_gene, model)

        model.add_mandatory_gene(model_genes['gspD'])
        model.add_mandatory_gene(model_genes['sctC'])
        model.add_accessory_gene(model_genes['sctJ'])
        model.add_accessory_gene(model_genes['sctN'])

        h10 = CoreHit(core_genes['gspD'], "h10", 10, "replicon_1", 10, 1.0, 11.0, 1.0, 1.0, 10, 20)
        m_h10 = ModelHit(h10, gene_ref=model_genes['gspD'], gene_status=GeneStatus.MANDATORY)

        h15 = CoreHit(core_genes['sctC'], "h15", 10, "replicon_1", 15, 1.0, 21.0, 1.0, 1.0, 10, 20)
        m_h15 = ModelHit(h15, gene_ref=model_genes['sctC'], gene_status=GeneStatus.ACCESSORY)
        self.assertTrue(_colocates(m_h10, m_h15, rep_info))

        h20 = CoreHit(core_genes['sctJ'], "h20", 10, "replicon_1", 20, 1.0, 21.0, 1.0, 1.0, 10, 20)
        m_h20 = ModelHit(h20, gene_ref=model_genes['sctJ'], gene_status=GeneStatus.ACCESSORY)
        self.assertFalse(_colocates(m_h10, m_h20, rep_info))

        model_genes['sctJ']._inter_gene_max_space = 10
        self.assertTrue(_colocates(m_h10, m_h20, rep_info))

        # case inter_gene_max_space is define at gene level sctJ and > than model intergene_max_space
        h30 = CoreHit(core_genes['sctJ'], "h30", 10, "replicon_1", 30, 1.0, 21.0, 1.0, 1.0, 10, 20)
        m_h30 = ModelHit(h30, gene_ref=model_genes['sctJ'], gene_status=GeneStatus.ACCESSORY)
        model_genes['sctJ']._inter_gene_max_space = 30
        self.assertTrue(_colocates(m_h10, m_h30, rep_info))

        # same case but inter_gene_max_space is define in first gene
        h30 = CoreHit(core_genes['sctJ'], "h30", 10, "replicon_1", 30, 1.0, 21.0, 1.0, 1.0, 10, 20)
        m_h30 = ModelHit(h30, gene_ref=model_genes['sctJ'], gene_status=GeneStatus.ACCESSORY)
        model_genes['sctJ']._inter_gene_max_space = 5
        h35 = CoreHit(core_genes['sctN'], "h35", 10, "replicon_1", 35, 1.0, 21.0, 1.0, 1.0, 10, 20)
        m_h35 = ModelHit(h35, gene_ref=model_genes['sctN'], gene_status=GeneStatus.ACCESSORY)
        self.assertTrue(_colocates(m_h30, m_h35, rep_info))

        # case inter_gene_max_sapce is define at gene level sctJ and > than model intergene_max_space
        h30 = CoreHit(core_genes['sctJ'], "h30", 10, "replicon_1", 15, 1.0, 21.0, 1.0, 1.0, 10, 20)
        m_h30 = ModelHit(h30, gene_ref=model_genes['sctJ'], gene_status=GeneStatus.ACCESSORY)
        model_genes['sctJ']._inter_gene_max_space = 5
        self.assertTrue(_colocates(m_h10, m_h30, rep_info))

        h30 = CoreHit(core_genes['sctJ'], "h30", 10, "replicon_1", 17, 1.0, 21.0, 1.0, 1.0, 10, 20)
        m_h30 = ModelHit(h30, gene_ref=model_genes['sctJ'], gene_status=GeneStatus.ACCESSORY)
        model_genes['sctJ']._inter_gene_max_space = 5
        self.assertFalse(_colocates(m_h10, m_h30, rep_info))

        # case inter_gene_max_sapce is define at gene level sctJ and gspD use the smallest one
        model_genes['gspD']._inter_gene_max_space = 7
        h30 = CoreHit(core_genes['sctJ'], "h30", 10, "replicon_1", 17, 1.0, 21.0, 1.0, 1.0, 10, 20)
        m_h30 = ModelHit(h30, gene_ref=model_genes['sctJ'], gene_status=GeneStatus.ACCESSORY)
        self.assertFalse(_colocates(m_h10, m_h30, rep_info))

        h30 = CoreHit(core_genes['sctJ'], "h30", 10, "replicon_1", 15, 1.0, 21.0, 1.0, 1.0, 10, 20)
        m_h30 = ModelHit(h30, gene_ref=model_genes['sctJ'], gene_status=GeneStatus.ACCESSORY)
        model_genes['sctJ']._inter_gene_max_space = 5
        self.assertTrue(_colocates(m_h10, m_h30, rep_info))


    def test_clusterize(self):
        # handle name, topology type, and min/max positions in the sequence dataset for a replicon and list of genes.
        # each genes is representing by a tuple (seq_id, length)"""
        rep_info = RepliconInfo('linear', 1, 60, [(f"g_{i}", i * 10) for i in range(1, 7)])

        #              fqn      , inter_gene_max_sapce
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

        #     CoreHit(gene, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        #                                                     pos     score
        h10 = CoreHit(core_genes[0], "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        m_h10 = ModelHit(h10, gene_ref=model_genes[0], gene_status=GeneStatus.MANDATORY)
        h11 = CoreHit(core_genes[0], "h11", 10, "replicon_1", 10, 1.0, 11.0, 1.0, 1.0, 10, 20)
        m_h11 = ModelHit(h11, gene_ref=model_genes[0], gene_status=GeneStatus.MANDATORY)

        h20 = CoreHit(core_genes[1], "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        m_h20 = ModelHit(h20, gene_ref=model_genes[1], gene_status=GeneStatus.MANDATORY)
        h21 = CoreHit(core_genes[2], "h21", 10, "replicon_1", 20, 1.0, 21.0, 1.0, 1.0, 10, 20)
        m_h21 = ModelHit(h21, gene_ref=model_genes[1], gene_status=GeneStatus.MANDATORY)

        h30 = CoreHit(core_genes[2], "h30", 10, "replicon_1", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        m_h30 = ModelHit(h30, gene_ref=model_genes[2], gene_status=GeneStatus.ACCESSORY)
        h31 = CoreHit(core_genes[1], "h31", 10, "replicon_1", 30, 1.0, 31.0, 1.0, 1.0, 10, 20)
        m_h31 = ModelHit(h31, gene_ref=model_genes[1], gene_status=GeneStatus.MANDATORY)

        h50 = CoreHit(core_genes[2], "h50", 10, "replicon_1", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        m_h50 = ModelHit(h50, gene_ref=model_genes[2], gene_status=GeneStatus.ACCESSORY)
        h51 = CoreHit(core_genes[2], "h51", 10, "replicon_1", 50, 1.0, 51.0, 1.0, 1.0, 10, 20)
        m_h51 = ModelHit(h51, gene_ref=model_genes[2], gene_status=GeneStatus.ACCESSORY)

        h60 = CoreHit(core_genes[2], "h60", 10, "replicon_1", 60, 1.0, 60.0, 1.0, 1.0, 10, 20)
        m_h60 = ModelHit(h60, gene_ref=model_genes[2], gene_status=GeneStatus.ACCESSORY)
        h61 = CoreHit(core_genes[3], "h61", 10, "replicon_1", 60, 1.0, 61.0, 1.0, 1.0, 10, 20)
        m_h61 = ModelHit(h61, gene_ref=model_genes[3], gene_status=GeneStatus.ACCESSORY)

        # case replicon is linear, 2 clusters
        hits = [m_h10, m_h11, m_h20, m_h21, m_h30, m_h31, m_h50, m_h51, m_h60, m_h61]
        random.shuffle(hits)
        got_clusters = _clusterize(hits, model, self.hit_weights, rep_info)
        self.assertEqual(len(got_clusters), 2)
        self.assertListEqual(got_clusters[0].hits, [m_h11, m_h21, m_h31])
        self.assertListEqual(got_clusters[1].hits, [m_h51, m_h61])

        # case replicon is linear with a single hit (not loner) between 2 clusters
        h70 = CoreHit(core_genes[3], "h70", 10, "replicon_1", 70, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h70 = ModelHit(h70, gene_ref=model_genes[3], gene_status=GeneStatus.ACCESSORY)
        h80 = CoreHit(core_genes[4], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[4], gene_status=GeneStatus.NEUTRAL)
        hits = [m_h10, m_h11, m_h20, m_h21, m_h50, m_h51, m_h70, m_h80]
        random.shuffle(hits)
        got_clusters = _clusterize(hits, model, self.hit_weights, rep_info)
        self.assertEqual(len(got_clusters), 2)
        self.assertListEqual(got_clusters[0].hits, [m_h11, m_h21])
        self.assertListEqual(got_clusters[1].hits, [m_h70, m_h80])

        # replicon is linear, 3 clusters, the last one contains only one hit (loner h80)
        rep_info = RepliconInfo('linear', 1, 100, [(f"g_{i}", i * 10) for i in range(1, 101)])
        h80 = CoreHit(core_genes[4], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[4], gene_status=GeneStatus.NEUTRAL)
        hits = [m_h10, m_h11, m_h20, m_h21, m_h30, m_h31, m_h50, m_h51, m_h60, m_h61, m_h80]
        random.shuffle(hits)
        got_clusters = _clusterize(hits, model, self.hit_weights, rep_info)
        self.assertEqual(len(got_clusters), 3)
        self.assertListEqual(got_clusters[0].hits, [m_h11, m_h21, m_h31])
        self.assertListEqual(got_clusters[1].hits, [m_h51, m_h61])
        self.assertListEqual(got_clusters[2].hits, [m_h80])

        # replicon is circular contains only one cluster
        rep_info = RepliconInfo('circular', 1, 60, [(f"g_{i}", i * 10) for i in range(1, 7)])
        hits = [m_h10, m_h20, m_h30]
        got_clusters = _clusterize(hits, model, self.hit_weights, rep_info)
        self.assertEqual(len(got_clusters), 1)
        self.assertListEqual(got_clusters[0].hits, [m_h10, m_h20, m_h30])

        # replicon is circular the last cluster is merge  with the first So we have only one cluster
        rep_info = RepliconInfo('circular', 1, 60, [(f"g_{i}", i * 10) for i in range(1, 7)])
        hits = [m_h10, m_h11, m_h20, m_h21, m_h30, m_h31, m_h50, m_h51, m_h60, m_h61]
        got_clusters = _clusterize(hits, model, self.hit_weights, rep_info)
        self.assertEqual(len(got_clusters), 1)
        self.assertListEqual(got_clusters[0].hits, [m_h51, m_h61, m_h11, m_h21, m_h31])

        # replicon is circular the last hit is incorporate to the first cluster
        # m_h80 is not considered as loner it is included in a cluster
        rep_info = RepliconInfo('circular', 1, 80, [(f"g_{i}", i * 10) for i in range(1, 9)])
        h80 = CoreHit(core_genes[3], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[4], gene_status=GeneStatus.NEUTRAL)
        hits = [m_h10, m_h11, m_h20, m_h21, m_h30, m_h31, m_h50, m_h51, m_h60, m_h61, m_h80]
        random.shuffle(hits)
        got_clusters = _clusterize(hits, model, self.hit_weights, rep_info)
        self.assertEqual(len(got_clusters), 2)
        self.assertListEqual(got_clusters[0].hits, [m_h80, m_h11, m_h21, m_h31])
        self.assertListEqual(got_clusters[1].hits, [m_h51, m_h61])

        # replicon is linear the last hit is not merged with the first cluster
        # m_h80 is link to gene_4 'abc'. So, in this test, it's not a loner.
        rep_info = RepliconInfo('linear', 1, 62, [(f"g_{i}", i * 10) for i in range(1, 7)])
        hits = [m_h10, m_h11, m_h20, m_h21, m_h30, m_h31, m_h50, m_h51, m_h60, m_h61]
        random.shuffle(hits)
        got_clusters = _clusterize(hits, model, self.hit_weights, rep_info)
        self.assertEqual(len(got_clusters), 2)
        self.assertListEqual(got_clusters[0].hits, [m_h11, m_h21, m_h31])
        self.assertListEqual(got_clusters[1].hits, [m_h51, m_h61])

        # case replicon is linear, 2 clusters, the hits 11,21,31 and 51,61 are contiguous
        #                                                              pos
        h10 = CoreHit(core_genes[0], "h10", 10, "replicon_1", 10, 1.0, 11.0, 1.0, 1.0, 10, 20)
        m_h10 = ModelHit(h10, gene_ref=model_genes[0], gene_status=GeneStatus.MANDATORY)
        h11 = CoreHit(core_genes[2], "h11", 10, "replicon_1", 11, 1.0, 21.0, 1.0, 1.0, 10, 20)
        m_h11 = ModelHit(h11, gene_ref=model_genes[2], gene_status=GeneStatus.ACCESSORY)
        h12 = CoreHit(core_genes[1], "h12", 10, "replicon_1", 12, 1.0, 31.0, 1.0, 1.0, 10, 20)
        m_h12 = ModelHit(h12, gene_ref=model_genes[1], gene_status=GeneStatus.MANDATORY)
        h50 = CoreHit(core_genes[2], "h50", 10, "replicon_1", 50, 1.0, 51.0, 1.0, 1.0, 10, 20)
        m_h50 = ModelHit(h50, gene_ref=model_genes[2], gene_status=GeneStatus.ACCESSORY)
        h51 = CoreHit(core_genes[3], "h51", 10, "replicon_1", 51, 1.0, 61.0, 1.0, 1.0, 10, 20)
        m_h51 = ModelHit(h51, gene_ref=model_genes[3], gene_status=GeneStatus.ACCESSORY)
        hits = [m_h10, m_h11, m_h12, m_h50, m_h51]
        random.shuffle(hits)
        got_clusters = _clusterize(hits, model, self.hit_weights, rep_info)
        self.assertEqual(len(got_clusters), 2)
        self.assertListEqual(got_clusters[0].hits, [m_h10, m_h11, m_h12])
        self.assertListEqual(got_clusters[1].hits, [m_h50, m_h51])

        # case replicon is linear
        # one cluster with one hit loner
        h80 = CoreHit(core_genes[4], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[4], gene_status=GeneStatus.NEUTRAL)
        hits = [m_h80]
        got_clusters = _clusterize(hits, model, self.hit_weights, rep_info)
        self.assertEqual(len(got_clusters), 1)
        self.assertListEqual(got_clusters[0].hits, [m_h80])

        # case replicon is linear
        # one cluster with one hit min_gene_required == 1
        # last hit is a cluster
        model = Model("foo/T2SS", 11, min_mandatory_genes_required=1, min_genes_required=1)
        core_genes = []
        model_genes = []
        for g_name in ('gspD', 'sctC'):
            core_gene = CoreGene(self.model_location, g_name, self.profile_factory)
            core_genes.append(core_gene)
            model_genes.append(ModelGene(core_gene, model))
        model.add_mandatory_gene(model_genes[0])
        model.add_accessory_gene(model_genes[1])

        h80 = CoreHit(core_genes[1], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[1], gene_status=GeneStatus.ACCESSORY)
        hits = [m_h80]
        got_clusters = _clusterize(hits, model, self.hit_weights, rep_info)
        self.assertEqual(len(got_clusters), 1)
        self.assertListEqual(got_clusters[0].hits, [m_h80])

        # case replicon is linear
        # one cluster with one hit min_gene_required == 1
        # first hit is a cluster
        model = Model("foo/T2SS", 11, min_mandatory_genes_required=1, min_genes_required=1)
        core_genes = []
        model_genes = []
        for g_name in ('gspD', 'sctC'):
            core_gene = CoreGene(self.model_location, g_name, self.profile_factory)
            core_genes.append(core_gene)
            model_genes.append(ModelGene(core_gene, model))
        model.add_mandatory_gene(model_genes[0])
        model.add_accessory_gene(model_genes[1])

        h10 = CoreHit(core_genes[0], "h10", 10, "replicon_1", 10, 1.0, 11.0, 1.0, 1.0, 10, 20)
        m_h10 = ModelHit(h10, gene_ref=model_genes[0], gene_status=GeneStatus.MANDATORY)
        h80 = CoreHit(core_genes[1], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[1], gene_status=GeneStatus.ACCESSORY)
        hits = [m_h10, m_h80]
        random.shuffle(hits)
        got_clusters = _clusterize(hits, model, self.hit_weights, rep_info)
        self.assertEqual(len(got_clusters), 2)
        self.assertListEqual(got_clusters[0].hits, [m_h10])
        self.assertListEqual(got_clusters[1].hits, [m_h80])

        # case replicon is linear
        # one cluster with one hit min_gene_required != 1
        model = Model("foo/T2SS", 11, min_mandatory_genes_required=1, min_genes_required=2)
        core_genes = []
        model_genes = []
        for g_name in ('gspD', 'sctC'):
            core_gene = CoreGene(self.model_location, g_name, self.profile_factory)
            core_genes.append(core_gene)
            model_genes.append(ModelGene(core_gene, model))
        model.add_mandatory_gene(model_genes[0])
        model.add_accessory_gene(model_genes[1])
        h80 = CoreHit(core_genes[1], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[1], gene_status=GeneStatus.ACCESSORY)
        hits = [m_h80]
        got_clusters = _clusterize(hits, model, self.hit_weights, rep_info)
        self.assertListEqual(got_clusters, [])

        # case replicon is linear, no hits
        got_clusters = _clusterize([], model, self.hit_weights, rep_info)
        self.assertListEqual(got_clusters, [])


    def test_get_multi_system_hits(self):
        #              fqn      , inter_gene_max_sapce
        model = Model("foo/T2SS", 11)

        core_genes = []
        model_genes = []
        for g_name in ('gspD', 'sctC', 'sctJ', 'sctN', 'abc', 'tadZ'):
            core_gene = CoreGene(self.model_location, g_name, self.profile_factory)
            core_genes.append(core_gene)
            model_genes.append(ModelGene(core_gene, model))
        model_genes[4]._loner = True
        model_genes[5]._multi_system = True

        model.add_mandatory_gene(model_genes[0])
        model.add_mandatory_gene(model_genes[1])
        model.add_accessory_gene(model_genes[2])
        model.add_accessory_gene(model_genes[3])
        model.add_neutral_gene(model_genes[4])
        model.add_accessory_gene(model_genes[5])

        #     CoreHit(gene, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        #                                                     pos     score
        h11 = CoreHit(core_genes[0], "h11", 10, "replicon_1", 10, 1.0, 11.0, 1.0, 1.0, 10, 20)
        m_h11 = ModelHit(h11, gene_ref=model_genes[0], gene_status=GeneStatus.MANDATORY)

        h21 = CoreHit(core_genes[2], "h21", 10, "replicon_1", 20, 1.0, 21.0, 1.0, 1.0, 10, 20)
        m_h21 = ModelHit(h21, gene_ref=model_genes[1], gene_status=GeneStatus.MANDATORY)

        h31 = CoreHit(core_genes[1], "h31", 10, "replicon_1", 30, 1.0, 31.0, 1.0, 1.0, 10, 20)
        m_h31 = ModelHit(h31, gene_ref=model_genes[1], gene_status=GeneStatus.MANDATORY)

        h51 = CoreHit(core_genes[2], "h51", 10, "replicon_1", 50, 1.0, 51.0, 1.0, 1.0, 10, 20)
        m_h51 = ModelHit(h51, gene_ref=model_genes[2], gene_status=GeneStatus.ACCESSORY)

        h61 = CoreHit(core_genes[3], "h61", 10, "replicon_1", 60, 1.0, 61.0, 1.0, 1.0, 10, 20)
        m_h61 = ModelHit(h61, gene_ref=model_genes[3], gene_status=GeneStatus.ACCESSORY)

        # case replicon is linear with a single hit (not loner) between 2 clusters
        # no multi system hit
        h70 = CoreHit(core_genes[3], "h70", 10, "replicon_1", 70, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h70 = ModelHit(h70, gene_ref=model_genes[3], gene_status=GeneStatus.ACCESSORY)
        h80 = CoreHit(core_genes[4], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[4], gene_status=GeneStatus.NEUTRAL)
        c0 = Cluster([m_h11, m_h21], model, self.hit_weights)
        c1 = Cluster([m_h70, m_h80], model, self.hit_weights)
        clusters, ms_hits = _get_multi_system_hits([c0, c1])
        self.assertEqual(len(clusters), 2)
        self.assertListEqual(clusters[0].hits, [m_h11, m_h21])
        self.assertListEqual(clusters[1].hits, [m_h70, m_h80])
        self.assertEqual(ms_hits, {})

        # replicon is linear, 1 cluster, containing one multi_system hit
        h32 = CoreHit(core_genes[5], "h32", 10, "replicon_1", 32, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h32 = ModelHit(h32, gene_ref=model_genes[5], gene_status=GeneStatus.ACCESSORY)

        c0 = Cluster([m_h11, m_h21, m_h31, m_h32], model, self.hit_weights)
        clusters, ms_hits = _get_multi_system_hits([c0])
        self.assertEqual(len(clusters), 1)
        self.assertListEqual(clusters[0].hits, [m_h11, m_h21, m_h31, m_h32])
        self.assertEqual(len(ms_hits), 1)
        self.assertTrue(isinstance(clusters[0].hits[-1], MultiSystem))

        # replicon is linear,
        # 2 clusters, containing one multi_system hit each (same gene)
        h32 = CoreHit(core_genes[5], "h32", 10, "replicon_1", 32, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h32 = ModelHit(h32, gene_ref=model_genes[5], gene_status=GeneStatus.ACCESSORY)
        h62 = CoreHit(core_genes[5], "h62", 10, "replicon_1", 62, 1.0, 80.0, 1.0, 1.0, 10, 30)
        m_h62 = ModelHit(h62, gene_ref=model_genes[5], gene_status=GeneStatus.ACCESSORY)

        c0 = Cluster([m_h11, m_h21, m_h32], model, self.hit_weights)
        c1 = Cluster([m_h51, m_h61, m_h62], model, self.hit_weights)
        clusters, ms_hits = _get_multi_system_hits([c0, c1])
        self.assertEqual(len(clusters), 2)
        # ModelHit implement __gt_, __lt__, __eq__ so MultiSystem inherits
        # and we can compare ModelHit and MultiSystem
        self.assertListEqual(clusters[0].hits, [m_h11, m_h21, m_h32])
        self.assertTrue(isinstance(clusters[0].hits[-1], MultiSystem))
        self.assertListEqual(clusters[0].hits[-1].counterpart,
                             [m_h62])

        self.assertListEqual(clusters[1].hits, [m_h51, m_h61, m_h62])
        self.assertTrue(isinstance(clusters[1].hits[-1], MultiSystem))
        self.assertListEqual(clusters[1].hits[-1].counterpart,
                             [m_h32])

        self.assertEqual(len(ms_hits), 1)
        # ms_hits values are Clusters of one MultiSystem hit
        # check that the ms_hit in ms_hits is the best (higest score) multisystem hit
        self.assertEqual(ms_hits['tadZ'].hits[0],
                         clusters[1].hits[-1].counterpart[-1])


    def test_get_true_loners(self):
        #              fqn      , inter_gene_max_sapce
        model = Model("foo/T2SS", 11)

        core_genes = []
        model_genes = []
        for g_name in ('gspD', 'sctC', 'sctJ', 'sctN', 'abc',  'flgB'):
            core_gene = CoreGene(self.model_location, g_name, self.profile_factory)
            core_genes.append(core_gene)
            model_genes.append(ModelGene(core_gene, model))
        model_genes[4]._loner = True
        model_genes[5]._loner = True
        model_genes[5]._multi_system = True

        model.add_mandatory_gene(model_genes[0])
        model.add_mandatory_gene(model_genes[1])
        model.add_accessory_gene(model_genes[2])
        model.add_accessory_gene(model_genes[3])
        model.add_neutral_gene(model_genes[4])
        model.add_accessory_gene(model_genes[5])


        #     CoreHit(gene, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        #                                                     pos     score
        h11 = CoreHit(core_genes[0], "h11", 10, "replicon_1", 10, 1.0, 11.0, 1.0, 1.0, 10, 20)
        m_h11 = ModelHit(h11, gene_ref=model_genes[0], gene_status=GeneStatus.MANDATORY)

        h21 = CoreHit(core_genes[2], "h21", 10, "replicon_1", 20, 1.0, 21.0, 1.0, 1.0, 10, 20)
        m_h21 = ModelHit(h21, gene_ref=model_genes[1], gene_status=GeneStatus.MANDATORY)

        h31 = CoreHit(core_genes[1], "h31", 10, "replicon_1", 30, 1.0, 31.0, 1.0, 1.0, 10, 20)
        m_h31 = ModelHit(h31, gene_ref=model_genes[1], gene_status=GeneStatus.MANDATORY)

        h51 = CoreHit(core_genes[2], "h51", 10, "replicon_1", 50, 1.0, 51.0, 1.0, 1.0, 10, 20)
        m_h51 = ModelHit(h51, gene_ref=model_genes[2], gene_status=GeneStatus.ACCESSORY)

        h61 = CoreHit(core_genes[3], "h61", 10, "replicon_1", 60, 1.0, 61.0, 1.0, 1.0, 10, 20)
        m_h61 = ModelHit(h61, gene_ref=model_genes[3], gene_status=GeneStatus.ACCESSORY)

        # case replicon is linear with a single hit (not loner) between 2 clusters
        h70 = CoreHit(core_genes[3], "h70", 10, "replicon_1", 70, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h70 = ModelHit(h70, gene_ref=model_genes[3], gene_status=GeneStatus.ACCESSORY)
        h80 = CoreHit(core_genes[4], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[4], gene_status=GeneStatus.NEUTRAL)
        c0 = Cluster([m_h11, m_h21], model, self.hit_weights)
        c1 = Cluster([m_h70, m_h80], model, self.hit_weights)
        true_loners, true_clusters = _get_true_loners([c0, c1])
        self.assertEqual(len(true_clusters), 2)
        self.assertListEqual(true_clusters[0].hits, [m_h11, m_h21])
        self.assertListEqual(true_clusters[1].hits, [m_h70, m_h80])
        self.assertEqual(true_loners, {})

        # replicon is linear, 3 clusters, the last one contains only one hit (loner h80)
        h80 = CoreHit(core_genes[4], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[4], gene_status=GeneStatus.NEUTRAL)

        c0 = Cluster([m_h11, m_h21, m_h31], model, self.hit_weights)
        c1 = Cluster([m_h51, m_h61], model, self.hit_weights)
        c2 = Cluster([m_h80], model, self.hit_weights)
        true_loners, true_clusters = _get_true_loners([c0, c1, c2])
        self.assertEqual(len(true_clusters), 2)
        self.assertListEqual(true_clusters[0].hits, [m_h11, m_h21, m_h31])
        self.assertListEqual(true_clusters[1].hits, [m_h51, m_h61])
        self.assertEqual(len(true_loners), 1)
        self.assertListEqual(true_loners['abc'].hits, [m_h80])

        # replicon is linear, 2 clusters, the last one contains 3 hits including one loner
        # although the gene is mark as loner as the hit is in cluster
        # it is not considered as True Loner
        # so it's type should be a ModelHit not Loner
        h80 = CoreHit(core_genes[4], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[4], gene_status=GeneStatus.NEUTRAL)

        c0 = Cluster([m_h11, m_h21, m_h31], model, self.hit_weights)
        c1 = Cluster([m_h51, m_h61, m_h80], model, self.hit_weights)
        true_loners, true_clusters = _get_true_loners([c0, c1])
        self.assertEqual(len(true_clusters), 2)
        self.assertListEqual(true_clusters[0].hits, [m_h11, m_h21, m_h31])
        self.assertListEqual(true_clusters[1].hits, [m_h51, m_h61, m_h80])
        self.assertEqual(len(true_loners), 0)
        self.assertTrue(isinstance(true_clusters[1].hits[-1], ModelHit))
        self.assertFalse(isinstance(true_clusters[1].hits[-1], Loner))

        # case replicon is linear
        # one cluster with one hit loner
        h80 = CoreHit(core_genes[4], "h80", 10, "replicon_1", 80, 1.0, 80.0, 1.0, 1.0, 10, 20)
        m_h80 = ModelHit(h80, gene_ref=model_genes[4], gene_status=GeneStatus.NEUTRAL)
        c0 = Cluster([m_h80], model, self.hit_weights)
        true_loners, true_clusters = _get_true_loners([c0])
        self.assertEqual(len(true_clusters), 0)
        self.assertTrue(isinstance(true_loners['abc'], Cluster))
        self.assertTrue(isinstance(true_loners['abc'][0], Loner))
        self.assertListEqual(true_loners['abc'].hits, [m_h80])

        # replicon is linear,
        # 2 clusters,
        #  - one regular cluster
        #  - 3 loner multisystem (same gene)
        h90 = CoreHit(core_genes[5], "h90", 10, "replicon_1", 90, 1.0, 90.0, 1.0, 1.0, 10, 20)
        h100 = CoreHit(core_genes[5], "h90", 10, "replicon_1", 90, 1.0, 100.0, 1.0, 1.0, 10, 20)
        h110 = CoreHit(core_genes[5], "h90", 10, "replicon_1", 90, 1.0, 110.0, 1.0, 1.0, 10, 20)
        m_h90 = MultiSystem(h90, gene_ref=model_genes[5], gene_status=GeneStatus.ACCESSORY, counterpart=[h100, h110])
        m_h100 = MultiSystem(h100, gene_ref=model_genes[5], gene_status=GeneStatus.ACCESSORY, counterpart=[h90, h110])
        m_h110 = MultiSystem(h110, gene_ref=model_genes[5], gene_status=GeneStatus.ACCESSORY, counterpart=[h90, h100])
        c0 = Cluster([m_h11, m_h21, m_h31], model, self.hit_weights)
        c1 = Cluster([m_h90], model, self.hit_weights)
        c2 = Cluster([m_h100], model, self.hit_weights)
        c3 = Cluster([m_h110], model, self.hit_weights)

        true_loners, true_clusters = _get_true_loners([c0, c1,c2, c3])
        self.assertEqual(len(true_clusters), 1)
        self.assertEqual(len(true_loners), 1)
        # check that the loner is a LonerMultiSystem
        # that is the hit with best score mh110
        # and hold the right counterpart
        self.assertTrue(isinstance(true_loners['flgB'], Cluster))
        self.assertTrue(isinstance(true_loners['flgB'][0], LonerMultiSystem))
        self.assertListEqual(true_loners['flgB'].hits, [m_h110])
        self.assertListEqual(true_loners['flgB'][0].counterpart, [h90, h100])


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

        h10 = CoreHit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        m_h10 = ModelHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = CoreHit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        m_h20 = ModelHit(h20, gene_1, GeneStatus.MANDATORY)
        h30 = CoreHit(c_gene_3, "h30", 10, "replicon_2", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        m_h30 = ModelHit(h30, gene_1, GeneStatus.ACCESSORY)
        h50 = CoreHit(c_gene_3, "h50", 10, "replicon_2", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        m_h50 = ModelHit(h50, gene_1, GeneStatus.ACCESSORY)

        with self.assertRaises(MacsypyError) as ctx:
            with self.catch_log():
                Cluster([m_h10, m_h20, m_h30, m_h50], model_1, self.hit_weights)
        msg = "Cannot build a cluster from hits coming from different replicons"
        self.assertEqual(str(ctx.exception), msg)


    def test_replicon_name(self):
        model = Model("foo/T2SS", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)

        gene_1 = ModelGene(c_gene_1, model)
        gene_2 = ModelGene(c_gene_2, model)

        replicon_name = "replicon_1"
        #     CoreHit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = CoreHit(c_gene_1, "h10", 10, replicon_name, 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        m_h10 = ModelHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = CoreHit(c_gene_2, "h20", 10, replicon_name, 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        m_h20 = ModelHit(h20, gene_2, GeneStatus.MANDATORY)

        c1 = Cluster([m_h10, m_h20], model, self.hit_weights)
        self.assertEqual(c1.replicon_name, replicon_name)


    def test_len(self):
        model = Model("foo/T2SS", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)

        gene_1 = ModelGene(c_gene_1, model)
        gene_2 = ModelGene(c_gene_2, model)

        #     CoreHit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = CoreHit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        m_h10 = ModelHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = CoreHit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        m_h20 = ModelHit(h20, gene_2, GeneStatus.MANDATORY)

        c1 = Cluster([m_h10, m_h20], model, self.hit_weights)
        self.assertEqual(len(c1), 2)

    def test_loner(self):
        model = Model("foo/bar", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)

        gene_1 = ModelGene(c_gene_1, model, loner=True)
        gene_2 = ModelGene(c_gene_2, model)

        #     CoreHit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = CoreHit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        m_h10 = ModelHit(h10, gene_1, GeneStatus.MANDATORY)
        l_h10 = Loner(h10, gene_1, GeneStatus.MANDATORY)
        h20 = CoreHit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        m_h20 = ModelHit(h20, gene_2, GeneStatus.MANDATORY)

        c1 = Cluster([l_h10], model, self.hit_weights)
        c2 = Cluster([m_h20], model, self.hit_weights)
        c3 = Cluster([m_h10, m_h20], model, self.hit_weights)
        self.assertTrue(c1.loner)
        self.assertFalse(c2.loner)
        self.assertFalse(c3.loner)


    def test_multi_system(self):
        model = Model("foo/bar", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)

        gene_1 = ModelGene(c_gene_1, model, multi_system=True)
        gene_2 = ModelGene(c_gene_2, model)

        #     CoreHit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = CoreHit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        m_h10 = ModelHit(h10, gene_1, GeneStatus.MANDATORY)
        l_h10 = MultiSystem(h10, gene_1, GeneStatus.MANDATORY)
        h20 = CoreHit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        m_h20 = ModelHit(h20, gene_2, GeneStatus.MANDATORY)

        c1 = Cluster([l_h10], model, self.hit_weights)
        c2 = Cluster([m_h20], model, self.hit_weights)
        c3 = Cluster([m_h10, m_h20], model, self.hit_weights)
        self.assertTrue(c1.multi_system)
        self.assertFalse(c2.multi_system)
        self.assertFalse(c3.multi_system)


    def test_contains(self):
        model = Model("foo/T2SS", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)
        c_gene_3 = CoreGene(self.model_location, "sctJ", self.profile_factory)

        gene_1 = ModelGene(c_gene_1, model)
        gene_2 = ModelGene(c_gene_2, model)
        gene_3 = ModelGene(c_gene_3, model)

        #     CoreHit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = CoreHit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        m_h10 = ModelHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = CoreHit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        m_h20 = ModelHit(h20, gene_2, GeneStatus.MANDATORY)
        h30 = CoreHit(c_gene_3, "h30", 10, "replicon_1", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        m_h30 = ModelHit(h30, gene_3, GeneStatus.ACCESSORY)
        h50 = CoreHit(c_gene_3, "h50", 10, "replicon_1", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        m_h50 = ModelHit(h50, gene_3, GeneStatus.ACCESSORY)
        c1 = Cluster([m_h10, m_h20, m_h50], model, self.hit_weights)

        self.assertTrue(m_h10 in c1)
        self.assertFalse(m_h30 in c1)


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

        #     CoreHit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = CoreHit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        m_h10 = ModelHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = CoreHit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        m_h20 = ModelHit(h20, gene_2, GeneStatus.MANDATORY)

        c = Cluster([m_h10, m_h20], model, self.hit_weights)

        self.assertTrue(c.fulfilled_function(gene_1))
        self.assertFalse(c.fulfilled_function(gene_3))

        # test with several genes
        self.assertTrue(c.fulfilled_function(gene_3, gene_1))

        # The cluster contains exchangeable
        h50 = CoreHit(c_gene_4, "h50", 10, "replicon_1", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        m_h50 = ModelHit(h50, gene_4, GeneStatus.ACCESSORY)
        c = Cluster([m_h10, m_h50], model, self.hit_weights)
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

        #     CoreHit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h_gspd = CoreHit(c_gene_gspd, "h_gspd", 10, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        m_h_gspd = ModelHit(h_gspd, gene_gspd, GeneStatus.MANDATORY)
        h_tadz = CoreHit(c_gene_tadZ, "h_tadz", 20, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        m_h_tadz = ModelHit(h_tadz, gene_tadZ, GeneStatus.MANDATORY)

        h_sctj = CoreHit(c_gene_sctj, "h_sctj", 30, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        m_h_sctj = ModelHit(h_sctj, gene_sctj, GeneStatus.ACCESSORY)
        h_sctj_an = CoreHit(c_gene_sctJ_FLG, "h_sctj_an", 30, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        m_h_sctj_an = ModelHit(h_sctj_an, analog_sctJ_FLG, GeneStatus.ACCESSORY)

        h_sctn = CoreHit(c_gene_sctn, "sctn", 40, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        m_h_sctn = ModelHit(h_sctn, gene_sctn, GeneStatus.ACCESSORY)
        h_sctn_hom = CoreHit(c_gene_sctn_FLG, "h_scth_hom", 30, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        m_h_sctn_hom = ModelHit(h_sctn_hom, homolog_sctn_FLG, GeneStatus.ACCESSORY)

        h_toto = CoreHit(c_gene_sctn, "toto", 50, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        m_h_toto = ModelHit(h_toto, gene_toto, GeneStatus.NEUTRAL)

        h_flie = CoreHit(c_gene_flie, "h_flie", 100, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        l_flie = Loner(h_flie, gene_flie, GeneStatus.MANDATORY)
        ms_flie = MultiSystem(h_flie, gene_flie, GeneStatus.MANDATORY)
        lms_flie = LonerMultiSystem(h_flie, gene_flie, GeneStatus.MANDATORY)

        # 2 mandatory, 2 accessory no analog/homolog
        c1 = Cluster([m_h_gspd, m_h_tadz, m_h_sctj, m_h_sctn], model, self.hit_weights)
        self.assertEqual(c1.score, 3.0)

        # 2 mandatory, 2 accessory 1 neutral, no analog/homolog
        c1 = Cluster([m_h_gspd, m_h_tadz, m_h_sctj, m_h_sctn, m_h_toto], model, self.hit_weights)
        self.assertEqual(c1.score, 3.0)

        # 1 mandatory + 1 mandatory duplicated 1 time
        # 1 accessory + 1 accessory duplicated 1 times
        # no analog/homolog
        c1 = Cluster([m_h_gspd, m_h_tadz, m_h_sctj, m_h_sctn, m_h_gspd, m_h_sctn], model, self.hit_weights)
        self.assertEqual(c1.score, 3.0)

        # 2 mandatory
        # 1 accessory + 1 accessory homolog
        c1 = Cluster([m_h_gspd, m_h_tadz, m_h_sctj, m_h_sctn_hom], model, self.hit_weights)
        self.assertEqual(c1.score, 2.9)

        # 2 mandatory
        # 1 accessory homolog + 2 accessory (check that the score used is the accessory not homolg)
        c1 = Cluster([m_h_sctn_hom, m_h_gspd, m_h_tadz, m_h_sctj, m_h_sctn], model, self.hit_weights)
        self.assertEqual(c1.score, 3.0)

        # test true loners
        c1 = Cluster([l_flie], model, self.hit_weights)
        self.assertEqual(c1.score, 0.7)

        # test multi system out of cluster
        c1 = Cluster([ms_flie], model, self.hit_weights)
        self.assertEqual(c1.score, 0.7)

        # test multi system out of cluster
        c1 = Cluster([lms_flie], model, self.hit_weights)
        self.assertEqual(c1.score, 0.7)

        # test the cache score
        c1 = Cluster([ms_flie], model, self.hit_weights)
        self.assertEqual(c1.score, 0.7)

        non_valid_hit = ModelHit(h_sctn, gene_sctn, GeneStatus.FORBIDDEN)
        c1 = Cluster([m_h_gspd, non_valid_hit, m_h_tadz], model, self.hit_weights)
        with self.assertRaises(MacsypyError) as ctx:
            c1.score
        self.assertEqual(str(ctx.exception),
                         "a Cluster contains hit sctN 1 which is neither mandatory nor accessory: forbidden")


    def test_merge(self):
        model = Model("foo/T2SS", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)
        c_gene_3 = CoreGene(self.model_location, "sctJ", self.profile_factory)

        gene_1 = ModelGene(c_gene_1, model)
        gene_2 = ModelGene(c_gene_2, model)
        gene_3 = ModelGene(c_gene_3, model)

        #     CoreHit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = CoreHit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        m_h10 = ModelHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = CoreHit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        m_h20 = ModelHit(h20, gene_2, GeneStatus.MANDATORY)
        h30 = CoreHit(c_gene_3, "h30", 10, "replicon_1", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        m_h30 = ModelHit(h30, gene_3, GeneStatus.ACCESSORY)
        h50 = CoreHit(c_gene_3, "h50", 10, "replicon_1", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        m_h50 = ModelHit(h50, gene_3, GeneStatus.ACCESSORY)

        c1 = Cluster([m_h10, m_h20], model, self.hit_weights)
        c2 = Cluster([m_h30, m_h50], model, self.hit_weights)
        c1.merge(c2)
        self.assertListEqual(c1.hits, [m_h10, m_h20, m_h30, m_h50])

        c1 = Cluster([m_h10, m_h20], model, self.hit_weights)
        c2 = Cluster([m_h30, m_h50], model, self.hit_weights)
        c2.merge(c1)
        self.assertListEqual(c2.hits, [m_h30, m_h50, m_h10, m_h20])

        c1 = Cluster([m_h10, m_h20], model, self.hit_weights)
        c2 = Cluster([m_h30, m_h50], model, self.hit_weights)
        c1.merge(c2, before=True)
        self.assertListEqual(c1.hits, [m_h30, m_h50, m_h10, m_h20])

        model_2 = Model("foo/T3SS", 11)
        c_gene_3 = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_3 = ModelGene(c_gene_3, model)

        h30 = CoreHit(c_gene_3, "h30", 10, "replicon_2", 30, 1.0, 30.0, 1.0, 1.0, 10, 20)
        m_h30 = ModelHit(h30, gene_3, GeneStatus.ACCESSORY)
        h50 = CoreHit(c_gene_3, "h50", 10, "replicon_2", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)
        m_h50 = ModelHit(h50, gene_3, GeneStatus.ACCESSORY)
        c3 = Cluster([m_h30, m_h50], model_2, self.hit_weights)
        with self.assertRaises(MacsypyError) as ctx:
            c1.merge(c3)
        self.assertEqual(str(ctx.exception), "Try to merge Clusters from different model")


    def test_str(self):
        model = Model("foo/T2SS", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)

        gene_1 = ModelGene(c_gene_1, model)
        gene_2 = ModelGene(c_gene_2, model)

        #     CoreHit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = CoreHit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        m_h10 = ModelHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = CoreHit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        m_h20 = ModelHit(h20, gene_2, GeneStatus.MANDATORY)
        c1 = Cluster([m_h10, m_h20], model, self.hit_weights)
        s ="""Cluster:
- model = T2SS
- replicon = replicon_1
- hits = (h10, gspD, 10), (h20, sctC, 20)"""
        self.assertEqual(str(c1), s)
