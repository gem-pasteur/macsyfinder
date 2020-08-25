#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2020  Institut Pasteur (Paris) and CNRS.           #
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

from macsypy.hit import Hit, ValidHit
from macsypy.config import Config, MacsyDefaults
from macsypy.gene import CoreGene, ModelGene, Exchangeable, GeneStatus
from macsypy.profile import ProfileFactory
from macsypy.model import Model
from macsypy.registries import ModelLocation
from macsypy.cluster import Cluster
from macsypy.system import OrderedMatchMaker, UnorderedMatchMaker
from macsypy.system import System, RejectedClusters, LikelySystem, UnlikelySystem

from tests import MacsyTest


class MatchMakerTest(MacsyTest):

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        self.cfg = Config(MacsyDefaults(), args)

        self.model_name = 'foo'
        self.model_location = ModelLocation(path=os.path.join(args.models_dir, self.model_name))
        self.profile_factory = ProfileFactory(self.cfg)

        self.model = Model("foo/model_A", 10)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, self.model)

        c_gene_sctn_flg = CoreGene(self.model_location, "sctN_FLG", self.profile_factory)
        gene_sctn_flg = Exchangeable(c_gene_sctn_flg, gene_sctn)
        gene_sctn.add_exchangeable(gene_sctn_flg)

        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, self.model)
        c_gene_sctj_flg = CoreGene(self.model_location, "sctJ_FLG", self.profile_factory)
        gene_sctj_flg = Exchangeable(c_gene_sctj_flg, gene_sctj)
        gene_sctj.add_exchangeable(gene_sctj_flg)

        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, self.model)
        c_gene_flgb = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd_an = Exchangeable(c_gene_flgb, gene_gspd)
        gene_gspd.add_exchangeable(gene_gspd_an)

        c_gene_abc = CoreGene(self.model_location, "abc", self.profile_factory)
        gene_abc = ModelGene(c_gene_abc, self.model)
        c_gene_tadz = CoreGene(self.model_location, "tadZ", self.profile_factory)
        gene_abc_ho = Exchangeable(c_gene_tadz, gene_abc)
        gene_abc.add_exchangeable(gene_abc_ho)

        c_gene_toto = CoreGene(self.model_location, "toto", self.profile_factory)
        gene_toto = ModelGene(c_gene_toto, self.model)
        c_gene_totote = CoreGene(self.model_location, "totote", self.profile_factory)
        gene_toto_ho = Exchangeable(c_gene_totote, gene_toto)
        gene_toto.add_exchangeable(gene_toto_ho)

        self.model.add_mandatory_gene(gene_sctn)
        self.model.add_mandatory_gene(gene_sctj)
        self.model.add_accessory_gene(gene_gspd)
        self.model.add_neutral_gene(gene_toto)
        self.model.add_forbidden_gene(gene_abc)

        self.hits = {}
        h_sctj = Hit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.hits['h_sctj'] = ValidHit(h_sctj, gene_sctj, GeneStatus.MANDATORY)
        h_sctj_flg = Hit(c_gene_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.hits['h_sctj_flg'] = ValidHit(h_sctj_flg, gene_sctj_flg, GeneStatus.MANDATORY)
        h_sctn = Hit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.hits['h_sctn'] = ValidHit(h_sctn, gene_sctn, GeneStatus.MANDATORY)
        h_sctn_flg = Hit(c_gene_sctn_flg, "hit_sctn_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.hits['h_sctn_flg'] = ValidHit(h_sctn_flg, gene_sctn_flg, GeneStatus.MANDATORY)
        h_gspd = Hit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.hits['h_gspd'] = ValidHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
        h_gspd_an = Hit(c_gene_flgb, "hit_gspd_an", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.hits['h_gspd_an'] = ValidHit(h_gspd_an, gene_gspd_an, GeneStatus.ACCESSORY)
        h_abc = Hit(c_gene_abc, "hit_abc", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.hits['h_abc'] = ValidHit(h_abc, gene_abc, GeneStatus.FORBIDDEN)
        h_abc_ho = Hit(c_gene_tadz, "hit_abc_ho", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.hits['h_abc_ho'] = ValidHit(h_abc_ho, gene_abc_ho, GeneStatus.FORBIDDEN)
        h_toto = Hit(c_gene_toto, "hit_toto", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.hits['h_toto'] = ValidHit(h_toto, gene_toto, GeneStatus.NEUTRAL)
        h_toto_ho = Hit(c_gene_totote, "hit_toto_ho", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.hits['h_toto_ho'] = ValidHit(h_toto_ho, gene_toto_ho, GeneStatus.NEUTRAL)


    def bla_test_ordered_match(self):

        #####################
        # test single locus #
        #####################

        # it lack one mandatory gene
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 3
        c1 = Cluster([self.hits['h_sctj'], self.hits['h_gspd']], self.model)
        ordered_match_maker = OrderedMatchMaker(self.model)
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, RejectedClusters)
        self.assertEqual(res.reason,
                         "The quorum of mandatory genes required (2) is not reached: 1\n"
                         "The quorum of genes required (3) is not reached: 2")

        # all quorum are reached
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        c1 = Cluster([self.hits['h_sctj'], self.hits['h_sctn'], self.hits['h_gspd']], self.model)
        ordered_match_maker = OrderedMatchMaker(self.model)
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, System)

        # with one mandatory analog
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        c1 = Cluster([self.hits['h_sctj_flg'], self.hits['h_sctn'], self.hits['h_gspd']], self.model)
        ordered_match_maker = OrderedMatchMaker(self.model)
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, System)

        # with one accessory analog
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        c1 = Cluster([self.hits['h_sctj'], self.hits['h_sctn'], self.hits['h_gspd_an']], self.model)
        ordered_match_maker = OrderedMatchMaker(self.model)
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, System)

        # the min_gene_required quorum is not reached
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 4
        c1 = Cluster([self.hits['h_sctj'], self.hits['h_sctn_flg'], self.hits['h_gspd']], self.model)
        ordered_match_maker = OrderedMatchMaker(self.model)
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, RejectedClusters)
        self.assertEqual(res.reason,
                         "The quorum of genes required (4) is not reached: 3")

        # the min_gene_required quorum is not reached even there is a neutral
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 4
        c1 = Cluster([self.hits['h_sctj'], self.hits['h_sctn_flg'], self.hits['h_gspd'], self.hits['h_toto']], self.model)
        ordered_match_maker = OrderedMatchMaker(self.model)
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, RejectedClusters)
        self.assertEqual(res.reason,
                         "The quorum of genes required (4) is not reached: 3")

        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 4
        c1 = Cluster([self.hits['h_sctj'], self.hits['h_sctn_flg'], self.hits['h_gspd'], self.hits['h_toto_ho']], self.model)
        ordered_match_maker = OrderedMatchMaker(self.model)
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, RejectedClusters)
        self.assertEqual(res.reason,
                         "The quorum of genes required (4) is not reached: 3")

        # the cluster contain a forbidden gene
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        c1 = Cluster([self.hits['h_sctj'], self.hits['h_sctn'], self.hits['h_gspd'], self.hits['h_abc']], self.model)
        ordered_match_maker = OrderedMatchMaker(self.model)
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, RejectedClusters)
        self.assertEqual(res.reason, "There is 1 forbidden genes occurrence(s): abc")

        # the cluster contain a forbidden gene homolog
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        c1 = Cluster([self.hits['h_sctj'], self.hits['h_sctn'], self.hits['h_gspd'], self.hits['h_abc_ho']], self.model)
        ordered_match_maker = OrderedMatchMaker(self.model)
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, RejectedClusters)
        self.assertEqual(res.reason, "There is 1 forbidden genes occurrence(s): tadZ")

        #####################
        # test multi loci   #
        #####################
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        c1 = Cluster([self.hits['h_sctj'], self.hits['h_sctn']], self.model)
        c2 = Cluster([self.hits['h_gspd']], self.model)
        ordered_match_maker = OrderedMatchMaker(self.model)
        res = ordered_match_maker.match([c1, c2])
        self.assertIsInstance(res, System)

        # with one analog an one homolog
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        c1 = Cluster([self.hits['h_sctj_flg'], self.hits['h_sctn_flg']], self.model)
        c2 = Cluster([self.hits['h_gspd']], self.model)
        ordered_match_maker = OrderedMatchMaker(self.model)
        res = ordered_match_maker.match([c1, c2])
        self.assertIsInstance(res, System)

        # with one analog an one homolog and one forbidden in 3 clusters
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        c1 = Cluster([self.hits['h_sctj_flg'], self.hits['h_sctn_flg']], self.model)
        c2 = Cluster([self.hits['h_gspd']], self.model)
        c3 = Cluster([self.hits['h_abc']], self.model)
        ordered_match_maker = OrderedMatchMaker(self.model)
        res = ordered_match_maker.match([c1, c2, c3])
        self.assertEqual(res.reason, "There is 1 forbidden genes occurrence(s): abc")

    def test_unordered_match(self):

        # # it lack one mandatory gene
        # self.model._min_mandatory_genes_required = 2
        # self.model._min_genes_required = 3
        # hits = [self.hits['h_sctj'], self.hits['h_gspd']]
        # unordered_match_maker = UnorderedMatchMaker(self.model)
        # res = unordered_match_maker.match(hits)
        # self.assertIsInstance(res, UnlikelySystem)
        # self.assertEqual(res.reason,
        #                  "The quorum of mandatory genes required (2) is not reached: 1\n"
        #                  "The quorum of genes required (3) is not reached: 2")
        #
        # # all quorum are reached
        # self.model._min_mandatory_genes_required = 2
        # self.model._min_genes_required = 1
        # hits = [self.hits['h_sctj'], self.hits['h_sctn'], self.hits['h_gspd']]
        # unordered_match_maker = UnorderedMatchMaker(self.model)
        # res = unordered_match_maker.match(hits)
        # self.assertIsInstance(res, LikelySystem)
        #
        # # with one mandatory analog
        # self.model._min_mandatory_genes_required = 2
        # self.model._min_genes_required = 1
        # hits = [self.hits['h_sctj_flg'], self.hits['h_sctn'], self.hits['h_gspd']]
        # unordered_match_maker = UnorderedMatchMaker(self.model)
        # res = unordered_match_maker.match(hits)
        # self.assertIsInstance(res, LikelySystem)
        #
        # # with one accessory analog
        # self.model._min_mandatory_genes_required = 2
        # self.model._min_genes_required = 1
        # hits = [self.hits['h_sctj'], self.hits['h_sctn'], self.hits['h_gspd_an']]
        # unordered_match_maker = UnorderedMatchMaker(self.model)
        # res = unordered_match_maker.match(hits)
        # self.assertIsInstance(res, LikelySystem)
        #
        # # the min_gene_required quorum is not reached
        # self.model._min_mandatory_genes_required = 2
        # self.model._min_genes_required = 4
        # hits = [self.hits['h_sctj'], self.hits['h_sctn_flg'], self.hits['h_gspd']]
        # unordered_match_maker = UnorderedMatchMaker(self.model)
        # res = unordered_match_maker.match(hits)
        # self.assertIsInstance(res, UnlikelySystem)
        # self.assertEqual(res.reason,
        #                  "The quorum of genes required (4) is not reached: 3")
        #
        # # the min_gene_required quorum is not reached even there is a neutral
        # self.model._min_mandatory_genes_required = 2
        # self.model._min_genes_required = 4
        # hits = [self.hits['h_sctj'], self.hits['h_sctn_flg'], self.hits['h_gspd'], self.hits['h_toto']]
        # unordered_match_maker = UnorderedMatchMaker(self.model)
        # res = unordered_match_maker.match(hits)
        # self.assertIsInstance(res, UnlikelySystem)
        # self.assertEqual(res.reason,
        #                  "The quorum of genes required (4) is not reached: 3")
        #
        # self.model._min_mandatory_genes_required = 2
        # self.model._min_genes_required = 4
        # hits = [self.hits['h_sctj'], self.hits['h_sctn_flg'], self.hits['h_gspd'], self.hits['h_toto_ho']]
        # unordered_match_maker = UnorderedMatchMaker(self.model)
        # res = unordered_match_maker.match(hits)
        # self.assertIsInstance(res, UnlikelySystem)
        # self.assertEqual(res.reason,
        #                  "The quorum of genes required (4) is not reached: 3")

        # the hits contain a forbidden gene
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        hits = [self.hits['h_sctj'], self.hits['h_sctn'], self.hits['h_gspd'], self.hits['h_abc']]
        unordered_match_maker = UnorderedMatchMaker(self.model)
        res = unordered_match_maker.match(hits)
        self.assertIsInstance(res, UnlikelySystem)
        self.assertEqual(res.reason, "There is 1 forbidden genes occurrence(s): abc")

        # # the cluster contain a forbidden gene homolog
        # self.model._min_mandatory_genes_required = 2
        # self.model._min_genes_required = 1
        # c1 = Cluster([self.hits['h_sctj'], self.hits['h_sctn'], self.hits['h_gspd'], self.hits['h_abc_ho']], self.model)
        # ordered_match_maker = OrderedMatchMaker(self.model)
        # res = ordered_match_maker.match([c1])
        # self.assertIsInstance(res, RejectedClusters)
        # self.assertEqual(res.reason, "There is 1 forbidden genes occurrence(s): tadZ")
