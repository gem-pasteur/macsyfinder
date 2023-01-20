#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2023  Institut Pasteur (Paris) and CNRS.           #
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

from macsypy.hit import CoreHit, ModelHit
from macsypy.config import Config, MacsyDefaults
from macsypy.gene import CoreGene, ModelGene, Exchangeable, GeneStatus
from macsypy.profile import ProfileFactory
from macsypy.model import Model
from macsypy.registries import ModelLocation
from macsypy.cluster import Cluster
from macsypy.system import OrderedMatchMaker, UnorderedMatchMaker
from macsypy.system import System, RejectedCandidate, LikelySystem, UnlikelySystem
from macsypy.error import MacsypyError

from tests import MacsyTest


class MatchMakerTest(MacsyTest):

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
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

        c_gene_flgb = CoreGene(self.model_location, "flgB", self.profile_factory)
        gene_gspd_ex = Exchangeable(c_gene_flgb, gene_gspd)
        gene_gspd.add_exchangeable(gene_gspd_ex)

        c_gene_abc = CoreGene(self.model_location, "abc", self.profile_factory)
        gene_abc = ModelGene(c_gene_abc, self.model)
        c_gene_tadz = CoreGene(self.model_location, "tadZ", self.profile_factory)
        gene_abc_ex = Exchangeable(c_gene_tadz, gene_abc)
        gene_abc.add_exchangeable(gene_abc_ex)

        c_gene_toto = CoreGene(self.model_location, "toto", self.profile_factory)
        gene_toto = ModelGene(c_gene_toto, self.model)
        c_gene_totote = CoreGene(self.model_location, "totote", self.profile_factory)
        gene_toto_ex = Exchangeable(c_gene_totote, gene_toto)
        gene_toto.add_exchangeable(gene_toto_ex)

        self.model.add_mandatory_gene(gene_sctn)
        self.model.add_mandatory_gene(gene_sctj)
        self.model.add_accessory_gene(gene_gspd)
        self.model.add_neutral_gene(gene_toto)
        self.model.add_forbidden_gene(gene_abc)

        self.c_hits = {
                    'ch_sctj': CoreHit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20),
                    'ch_sctj_flg': CoreHit(c_gene_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20),
                    'ch_sctn': CoreHit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20),
                    'ch_sctn_flg': CoreHit(c_gene_sctn_flg, "hit_sctn_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20),
                    'ch_gspd': CoreHit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20),
                    'ch_gspd_ex': CoreHit(c_gene_flgb, "hit_gspd_an", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20),
                    'ch_abc': CoreHit(c_gene_abc, "hit_abc", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20),
                    'ch_abc_ex': CoreHit(c_gene_tadz, "hit_abc_ho", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20),
                    'ch_toto': CoreHit(c_gene_toto, "hit_toto", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20),
                    'ch_toto_ex': CoreHit(c_gene_totote, "hit_toto_ho", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20),
                    }
        self.m_hits = {
            'mh_sctj': ModelHit(self.c_hits['ch_sctj'], gene_sctj, GeneStatus.MANDATORY),
            'mh_sctj_flg': ModelHit(self.c_hits['ch_sctj_flg'], gene_sctj_flg, GeneStatus.MANDATORY),
            'mh_sctn': ModelHit(self.c_hits['ch_sctn'], gene_sctn, GeneStatus.MANDATORY),
            'mh_sctn_flg': ModelHit(self.c_hits['ch_sctn_flg'], gene_sctn_flg, GeneStatus.MANDATORY),
            'mh_gspd': ModelHit(self.c_hits['ch_gspd'], gene_gspd, GeneStatus.ACCESSORY),
            'mh_gspd_ex': ModelHit(self.c_hits['ch_gspd_ex'], gene_gspd_ex, GeneStatus.ACCESSORY),
            'mh_abc': ModelHit(self.c_hits['ch_abc'], gene_abc, GeneStatus.FORBIDDEN),
            'mh_abc_ex': ModelHit(self.c_hits['ch_abc_ex'], gene_abc_ex, GeneStatus.FORBIDDEN),
            'mh_toto': ModelHit(self.c_hits['ch_toto'], gene_toto, GeneStatus.NEUTRAL),
            'mh_toto_ex': ModelHit(self.c_hits['ch_toto_ex'], gene_toto_ex, GeneStatus.NEUTRAL)
        }


    def test_sort_hits_by_status(self):
        ordered_match_maker = OrderedMatchMaker(self.model, self.cfg.redundancy_penalty())
        mandatory_exp = [self.m_hits['mh_sctn'], self.m_hits['mh_sctj']]
        accessory_exp = [self.m_hits['mh_gspd']]
        neutral_exp = [self.m_hits['mh_toto']]
        forbidden_exp = [self.m_hits['mh_abc']]

        mandatory, accessory, neutral, forbidden = ordered_match_maker.sort_hits_by_status(mandatory_exp + accessory_exp + neutral_exp + forbidden_exp)
        self.assertListEqual([h.gene.name for h in mandatory_exp], [h.gene.name for h in mandatory])
        self.assertListEqual([h.gene.name for h in accessory_exp], [h.gene.name for h in accessory])
        self.assertListEqual([h.gene.name for h in neutral_exp], [h.gene.name for h in neutral])
        self.assertListEqual([h.gene.name for h in forbidden_exp], [h.gene.name for h in forbidden])

        # do the same but with exchangeable
        mandatory_exp_exch = [self.m_hits['mh_sctn_flg'], self.m_hits['mh_sctj_flg']]
        accessory_exp_exch = [self.m_hits['mh_gspd_ex']]
        neutral_exp_exch = [self.m_hits['mh_toto_ex']]
        forbidden_exp_exch = [self.m_hits['mh_abc_ex']]

        mandatory, accessory, neutral, forbidden = ordered_match_maker.sort_hits_by_status(mandatory_exp_exch +
                                                                                           accessory_exp_exch +
                                                                                           neutral_exp_exch +
                                                                                           forbidden_exp_exch)
        self.assertListEqual([h.gene.name for h in mandatory_exp_exch], [h.gene.name for h in mandatory])
        self.assertListEqual([h.gene.name for h in accessory_exp_exch], [h.gene.name for h in accessory])
        self.assertListEqual([h.gene.name for h in neutral_exp_exch], [h.gene.name for h in neutral])
        self.assertListEqual([h.gene.name for h in forbidden_exp_exch], [h.gene.name for h in forbidden])

        # test if gene_ref is the ModelGene
        # alternate_of return the ModelGene of the function
        self.assertListEqual([h.gene.name for h in mandatory_exp], [h.gene_ref.alternate_of().name for h in mandatory])
        self.assertListEqual([h.gene.name for h in accessory_exp], [h.gene_ref.alternate_of().name for h in accessory])
        self.assertListEqual([h.gene.name for h in neutral_exp], [h.gene_ref.alternate_of().name for h in neutral])
        self.assertListEqual([h.gene.name for h in forbidden_exp], [h.gene_ref.alternate_of().name for h in forbidden])

        # test if the hit does not refer to gene belonging to the model
        model2 = Model("foo/model_B", 10)
        cg_fliE = CoreGene(self.model_location, "fliE", self.profile_factory)
        ch_fliE = CoreHit(cg_fliE, "hit_fliE", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        mg_fliE = ModelGene(cg_fliE, model2)
        mh_fliE = ModelHit(ch_fliE, mg_fliE, GeneStatus.NEUTRAL)
        with self.assertRaises(MacsypyError) as ctx:
            with self.catch_log():
                ordered_match_maker.sort_hits_by_status([mh_fliE])
        self.assertEqual(str(ctx.exception),
                         "Gene 'fliE' not found in model 'foo/model_B'")

    def test_ordered_match(self):

        #####################
        # test single locus #
        #####################

        # it lack one mandatory gene
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 3
        c1 = Cluster([self.m_hits['mh_sctj'], self.m_hits['mh_gspd']], self.model, self.cfg.hit_weights())
        ordered_match_maker = OrderedMatchMaker(self.model, self.cfg.redundancy_penalty())
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, RejectedCandidate)
        self.assertEqual(res.reasons,
                         ["The quorum of mandatory genes required (2) is not reached: 1",
                          "The quorum of genes required (3) is not reached: 2"])

        # all quorum are reached
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        c1 = Cluster([self.m_hits['mh_sctj'], self.m_hits['mh_sctn'], self.m_hits['mh_gspd']], self.model, self.cfg.hit_weights())
        ordered_match_maker = OrderedMatchMaker(self.model, self.cfg.redundancy_penalty())
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, System)

        # with one mandatory analog
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        c1 = Cluster([self.m_hits['mh_sctj_flg'], self.m_hits['mh_sctn'], self.m_hits['mh_gspd']], self.model, self.cfg.hit_weights())
        ordered_match_maker = OrderedMatchMaker(self.model, self.cfg.redundancy_penalty())
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, System)

        # with one accessory analog
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        c1 = Cluster([self.m_hits['mh_sctj'], self.m_hits['mh_sctn'], self.m_hits['mh_gspd_ex']], self.model, self.cfg.hit_weights())
        ordered_match_maker = OrderedMatchMaker(self.model, self.cfg.redundancy_penalty())
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, System)

        # the min_gene_required quorum is not reached
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 4
        c1 = Cluster([self.m_hits['mh_sctj'], self.m_hits['mh_sctn_flg'], self.m_hits['mh_gspd']], self.model, self.cfg.hit_weights())
        ordered_match_maker = OrderedMatchMaker(self.model, self.cfg.redundancy_penalty())
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, RejectedCandidate)
        self.assertListEqual(res.reasons,
                             ["The quorum of genes required (4) is not reached: 3"])

        # the min_gene_required quorum is not reached even there is a neutral
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 4
        c1 = Cluster([self.m_hits['mh_sctj'], self.m_hits['mh_sctn_flg'], self.m_hits['mh_gspd'], self.m_hits['mh_toto']],
                     self.model, self.cfg.hit_weights())
        ordered_match_maker = OrderedMatchMaker(self.model, self.cfg.redundancy_penalty())
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, RejectedCandidate)
        self.assertEqual(res.reasons,
                         ["The quorum of genes required (4) is not reached: 3"])

        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 4
        c1 = Cluster([self.m_hits['mh_sctj'], self.m_hits['mh_sctn_flg'], self.m_hits['mh_gspd'], self.m_hits['mh_toto_ex']],
                     self.model, self.cfg.hit_weights())
        ordered_match_maker = OrderedMatchMaker(self.model, self.cfg.redundancy_penalty())
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, RejectedCandidate)
        self.assertEqual(res.reasons,
                         ["The quorum of genes required (4) is not reached: 3"])

        # the cluster contain a forbidden gene
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        c1 = Cluster([self.m_hits['mh_sctj'], self.m_hits['mh_sctn'], self.m_hits['mh_gspd'], self.m_hits['mh_abc']],
                     self.model, self.cfg.hit_weights())
        ordered_match_maker = OrderedMatchMaker(self.model, self.cfg.redundancy_penalty())
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, RejectedCandidate)
        self.assertEqual(res.reasons, ["There is 1 forbidden genes occurrence(s): abc"])

        # the cluster contain a forbidden gene homolog
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        c1 = Cluster([self.m_hits['mh_sctj'], self.m_hits['mh_sctn'], self.m_hits['mh_gspd'], self.m_hits['mh_abc_ex']],
                     self.model, self.cfg.hit_weights())
        ordered_match_maker = OrderedMatchMaker(self.model, self.cfg.redundancy_penalty())
        res = ordered_match_maker.match([c1])
        self.assertIsInstance(res, RejectedCandidate)
        self.assertEqual(res.reasons, ["There is 1 forbidden genes occurrence(s): tadZ"])

        #####################
        # test multi loci   #
        #####################
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        c1 = Cluster([self.m_hits['mh_sctj'], self.m_hits['mh_sctn']], self.model, self.cfg.hit_weights())
        c2 = Cluster([self.m_hits['mh_gspd']], self.model, self.cfg.hit_weights())
        ordered_match_maker = OrderedMatchMaker(self.model, self.cfg.redundancy_penalty())
        res = ordered_match_maker.match([c1, c2])
        self.assertIsInstance(res, System)

        # with one analog an one homolog
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        c1 = Cluster([self.m_hits['mh_sctj_flg'], self.m_hits['mh_sctn_flg']], self.model, self.cfg.hit_weights())
        c2 = Cluster([self.m_hits['mh_gspd_ex']], self.model, self.cfg.hit_weights())
        ordered_match_maker = OrderedMatchMaker(self.model, self.cfg.redundancy_penalty())
        res = ordered_match_maker.match([c1, c2])
        self.assertIsInstance(res, System)

        # with one analog an one homolog and one forbidden in 3 clusters
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        c1 = Cluster([self.m_hits['mh_sctj_flg'], self.m_hits['mh_sctn_flg']], self.model, self.cfg.hit_weights())
        c2 = Cluster([self.m_hits['mh_gspd']], self.model, self.cfg.hit_weights())
        c3 = Cluster([self.m_hits['mh_abc']], self.model, self.cfg.hit_weights())
        ordered_match_maker = OrderedMatchMaker(self.model, self.cfg.redundancy_penalty())
        res = ordered_match_maker.match([c1, c2, c3])
        self.assertEqual(res.reasons, ["There is 1 forbidden genes occurrence(s): abc"])

    def test_unordered_match(self):

        # it lack one mandatory gene
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 3
        hits = [self.m_hits['mh_sctj'], self.m_hits['mh_gspd']]
        unordered_match_maker = UnorderedMatchMaker(self.model)
        res = unordered_match_maker.match(hits)
        self.assertIsInstance(res, UnlikelySystem)
        self.assertEqual(res.reasons,
                         ["The quorum of mandatory genes required (2) is not reached: 1",
                          "The quorum of genes required (3) is not reached: 2"])

        # all quorum are reached
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        hits = [self.m_hits['mh_sctj'], self.m_hits['mh_sctn'], self.m_hits['mh_gspd']]
        unordered_match_maker = UnorderedMatchMaker(self.model)
        res = unordered_match_maker.match(hits)
        self.assertIsInstance(res, LikelySystem)

        # with one mandatory analog
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        hits = [self.m_hits['mh_sctj_flg'], self.m_hits['mh_sctn'], self.m_hits['mh_gspd']]
        unordered_match_maker = UnorderedMatchMaker(self.model)
        res = unordered_match_maker.match(hits)
        self.assertIsInstance(res, LikelySystem)

        # with one accessory analog
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        hits = [self.m_hits['mh_sctj'], self.m_hits['mh_sctn'], self.m_hits['mh_gspd_ex']]
        unordered_match_maker = UnorderedMatchMaker(self.model)
        res = unordered_match_maker.match(hits)
        self.assertIsInstance(res, LikelySystem)

        # the min_gene_required quorum is not reached
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 4
        hits = [self.m_hits['mh_sctj'], self.m_hits['mh_sctn_flg'], self.m_hits['mh_gspd']]
        unordered_match_maker = UnorderedMatchMaker(self.model)
        res = unordered_match_maker.match(hits)
        self.assertIsInstance(res, UnlikelySystem)
        self.assertEqual(res.reasons,
                         ["The quorum of genes required (4) is not reached: 3"])

        # the min_gene_required quorum is not reached even there is a neutral
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 4
        hits = [self.m_hits['mh_sctj'], self.m_hits['mh_sctn_flg'], self.m_hits['mh_gspd'], self.m_hits['mh_toto']]
        unordered_match_maker = UnorderedMatchMaker(self.model)
        res = unordered_match_maker.match(hits)
        self.assertIsInstance(res, UnlikelySystem)
        self.assertEqual(res.reasons,
                         ["The quorum of genes required (4) is not reached: 3"])

        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 4
        hits = [self.m_hits['mh_sctj'], self.m_hits['mh_sctn_flg'], self.m_hits['mh_gspd'], self.m_hits['mh_toto_ex']]
        unordered_match_maker = UnorderedMatchMaker(self.model)
        res = unordered_match_maker.match(hits)
        self.assertIsInstance(res, UnlikelySystem)
        self.assertEqual(res.reasons,
                         ["The quorum of genes required (4) is not reached: 3"])

        # the hits contain a forbidden gene
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        allowed_hits = [self.m_hits['mh_sctj'], self.m_hits['mh_sctn'], self.m_hits['mh_gspd']]
        forbidden_hits = [self.m_hits['mh_abc']]
        unordered_match_maker = UnorderedMatchMaker(self.model)
        res = unordered_match_maker.match(allowed_hits + forbidden_hits)
        self.assertIsInstance(res, LikelySystem)
        self.assertListEqual([(h.id, h.position) for h in res.hits],
                             [(h.id, h.position) for h in allowed_hits + forbidden_hits])
        self.assertListEqual(res._forbidden_hits, [self.m_hits['mh_abc']])

        # the hits contain a forbidden gene homolog
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        hits = [self.m_hits['mh_sctj'], self.m_hits['mh_sctn'], self.m_hits['mh_gspd'], self.m_hits['mh_abc_ex']]

        unordered_match_maker = UnorderedMatchMaker(self.model)
        res = unordered_match_maker.match(hits)
        self.assertIsInstance(res, LikelySystem)
        self.assertListEqual(res._forbidden_hits, [self.m_hits['mh_abc_ex']])

        # only one forbidden hit (no mandatory, accessory, neutral)
        self.model._min_mandatory_genes_required = 2
        self.model._min_genes_required = 1
        allowed_hits = []
        forbidden_hits = [self.m_hits['mh_abc']]
        unordered_match_maker = UnorderedMatchMaker(self.model)
        res = unordered_match_maker.match(allowed_hits + forbidden_hits)
        self.assertIsNone(res)

