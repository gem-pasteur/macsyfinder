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
import itertools

from macsypy.hit import CoreHit, ModelHit, Loner, HitWeight
from macsypy.config import Config, MacsyDefaults
from macsypy.gene import CoreGene, ModelGene, Exchangeable, GeneStatus
from macsypy.profile import ProfileFactory
from macsypy.model import Model
from macsypy.registries import ModelLocation
from macsypy.cluster import Cluster
from macsypy.system import System, HitSystemTracker, LikelySystem, UnlikelySystem, AbstractUnordered, RejectedCandidate
from macsypy.solution import Solution
from macsypy.serialization import TxtSystemSerializer, TsvSystemSerializer, TsvSolutionSerializer, \
    TxtLikelySystemSerializer, TxtUnikelySystemSerializer, TsvSpecialHitSerializer, TsvRejectedCandidatesSerializer

from tests import MacsyTest


class SerializationTest(MacsyTest):

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        self.cfg = Config(MacsyDefaults(), args)

        self.model_name = 'foo'
        self.model_location = ModelLocation(path=os.path.join(args.models_dir, self.model_name))
        self.profile_factory = ProfileFactory(self.cfg)
        self.hit_weights = HitWeight(**self.cfg.hit_weights())
        # reset the uniq id number for AbstractUnordered
        # to have predictable results for (Likely/Unlikely)Systems
        System._id = itertools.count(1)
        AbstractUnordered._id = itertools.count(1)

    def test_SystemSerializer_str(self):
        model_name = 'foo'
        model_location = ModelLocation(path=os.path.join(self.cfg.models_dir()[0], model_name))
        model_A = Model("foo/A", 10)
        model_B = Model("foo/B", 10)

        c_gene_sctn_flg = CoreGene(model_location, "sctN_FLG", self.profile_factory)
        gene_sctn_flg = ModelGene(c_gene_sctn_flg, model_B)
        c_gene_sctj_flg = CoreGene(model_location, "sctJ_FLG", self.profile_factory)
        gene_sctj_flg = ModelGene(c_gene_sctj_flg, model_B)
        c_gene_flgB = CoreGene(model_location, "flgB", self.profile_factory)
        c_gene_tadZ = CoreGene(model_location, "tadZ", self.profile_factory)
        gene_tadZ = ModelGene(c_gene_tadZ, model_B)

        c_gene_sctn = CoreGene(model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model_A)
        gene_sctn_hom = Exchangeable(c_gene_sctn_flg, gene_sctn)
        gene_sctn.add_exchangeable(gene_sctn_hom)

        c_gene_sctj = CoreGene(model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model_A)
        gene_sctj_an = Exchangeable(c_gene_sctj_flg, gene_sctj)
        gene_sctj.add_exchangeable(gene_sctj_an)

        c_gene_gspd = CoreGene(model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model_A)
        gene_gspd_an = Exchangeable(c_gene_flgB, gene_gspd)
        gene_gspd.add_exchangeable(gene_gspd_an)

        c_gene_abc = CoreGene(model_location, "abc", self.profile_factory)
        gene_abc = ModelGene(c_gene_abc, model_A)
        gene_abc_ho = Exchangeable(c_gene_tadZ, gene_abc)
        gene_abc.add_exchangeable(gene_abc_ho)

        model_A.add_mandatory_gene(gene_sctn)
        model_A.add_mandatory_gene(gene_sctj)
        model_A.add_accessory_gene(gene_gspd)
        model_A.add_forbidden_gene(gene_abc)

        model_B.add_mandatory_gene(gene_sctn_flg)
        model_B.add_mandatory_gene(gene_sctj_flg)
        model_B.add_accessory_gene(gene_gspd)
        model_B.add_accessory_gene(gene_tadZ)

        h_sctj = CoreHit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = CoreHit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd = CoreHit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_sctj_flg = CoreHit(c_gene_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_tadZ = CoreHit(c_gene_tadZ, "hit_tadZ", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        model_A._min_mandatory_genes_required = 2
        model_A._min_genes_required = 2
        c1 = Cluster([ModelHit(h_sctj, gene_sctj, GeneStatus.MANDATORY),
                      ModelHit(h_sctn, gene_sctn, GeneStatus.MANDATORY),
                      ModelHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
                      ],
                     model_A, self.hit_weights)

        c2 = Cluster([ModelHit(h_sctj, gene_sctj, GeneStatus.MANDATORY),
                      ModelHit(h_sctn, gene_sctn, GeneStatus.MANDATORY)],
                     model_A, self.hit_weights)

        model_B._min_mandatory_genes_required = 1
        model_B._min_genes_required = 2
        c3 = Cluster([ModelHit(h_sctj_flg, gene_sctj_flg, GeneStatus.MANDATORY),
                      ModelHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY),
                      ModelHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)],
                     model_B, self.hit_weights)

        sys_A = System(model_A, [c1, c2], self.cfg.redundancy_penalty())
        sys_A.id = "sys_id_A"
        sys_B = System(model_B, [c3], self.cfg.redundancy_penalty())
        sys_B.id = "sys_id_B"
        hit_multi_sys_tracker = HitSystemTracker([sys_A, sys_B])
        system_serializer = TxtSystemSerializer()

        sys_str = f"""system id = {sys_A.id}
model = foo/A
replicon = replicon_id
clusters = [('hit_sctj', 'sctJ', 1), ('hit_sctn', 'sctN', 1), ('hit_gspd', 'gspD', 1)], [('hit_sctj', 'sctJ', 1), ('hit_sctn', 'sctN', 1)]
occ = 2
wholeness = 1.000
loci nb = 2
score = 1.500

mandatory genes:
\t- sctN: 2 (sctN, sctN)
\t- sctJ: 2 (sctJ, sctJ)

accessory genes:
\t- gspD: 1 (gspD [sys_id_B])

neutral genes:
"""
        self.assertEqual(sys_str, system_serializer.serialize(sys_A, hit_multi_sys_tracker))


    def test_SystemSerializer_tsv(self):
        model = Model("foo/T2SS", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model, loner=True)
        c_gene_sctn_flg = CoreGene(self.model_location, "sctN_FLG", self.profile_factory)
        gene_sctn_flg = Exchangeable(c_gene_sctn_flg, gene_sctn)
        gene_sctn.add_exchangeable(gene_sctn_flg)
        model.add_accessory_gene(gene_sctn)

        #CoreHit(gene, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        #                                                           pos      score
        ch_gspd = CoreHit(c_gene_gspd, "h_gspd", 803, "replicon_id", 10, 1.0, 1.0, 1.0, 1.0, 10, 20)
        mh_gspd = ModelHit(ch_gspd, gene_ref=gene_gspd, gene_status=GeneStatus.MANDATORY)
        ch_sctj = CoreHit(c_gene_sctj, "h_sctj", 803, "replicon_id", 20, 1.0, 1.0, 1.0, 1.0, 20, 30)
        mh_sctj = ModelHit(ch_sctj, gene_ref=gene_sctj, gene_status=GeneStatus.ACCESSORY)

        ch_sctn_flg = CoreHit(c_gene_sctn_flg, "h_sctn_flg", 803, "replicon_id", 40, 1.0, 1.0, 1.0, 1.0, 30, 40)
        mh_sctn_flg = ModelHit(ch_sctn_flg, gene_ref=gene_sctn_flg, gene_status=GeneStatus.ACCESSORY)
        ch_sctn = CoreHit(c_gene_sctn, "h_sctn", 803, "replicon_id", 80, 1.0, 1.0, 1.0, 1.0, 30, 40)
        mh_sctn = Loner(ch_sctn, gene_ref=gene_sctn, gene_status=GeneStatus.ACCESSORY, counterpart=[mh_sctn_flg])

        c1 = Cluster([mh_gspd, mh_sctj], model, self.hit_weights)
        c2 = Cluster([mh_sctn], model, self.hit_weights)
        sys_multi_loci = System(model, [c1, c2], self.cfg.redundancy_penalty())
        # score                         1.5 .35 = 1.85
        hit_multi_sys_tracker = HitSystemTracker([sys_multi_loci])
        system_serializer = TsvSystemSerializer()

        sys_tsv = "\t".join(["replicon_id", "h_gspd", "gspD", "10", "foo/T2SS", sys_multi_loci.id, "1", "1",
                             "1.000", "1.850", "1", "gspD", "mandatory", "803",
                             "1.0", "1.000", "1.000", "1.000", "10", "20", "", ""])
        sys_tsv += "\n"
        sys_tsv += "\t".join(["replicon_id", "h_sctj", "sctJ", "20", "foo/T2SS", sys_multi_loci.id, "1", "1",
                              "1.000", "1.850", "1", "sctJ", "accessory", "803",
                              "1.0", "1.000", "1.000", "1.000", "20", "30", "", ""])
        sys_tsv += "\n"
        sys_tsv += "\t".join(["replicon_id", "h_sctn", "sctN", "80", "foo/T2SS", sys_multi_loci.id, "1", "-1",
                              "1.000", "1.850", "1", "sctN", "accessory", "803",
                              "1.0", "1.000", "1.000", "1.000", "30", "40", "h_sctn_flg", ""])
        sys_tsv += "\n"
        self.maxDiff = None
        self.assertEqual(sys_tsv, system_serializer.serialize(sys_multi_loci, hit_multi_sys_tracker))


    def test_SolutionSerializer_tsv(self):
        model_name = 'foo'
        model_location = ModelLocation(path=os.path.join(self.cfg.models_dir()[0], model_name))

        ###########
        # Model B #
        ###########
        model_B = Model("foo/B", 10)
        c_gene_sctn_flg = CoreGene(model_location, "sctN_FLG", self.profile_factory)
        gene_sctn_flg = ModelGene(c_gene_sctn_flg, model_B)
        c_gene_sctj_flg = CoreGene(model_location, "sctJ_FLG", self.profile_factory)
        gene_sctj_flg = ModelGene(c_gene_sctj_flg, model_B)
        c_gene_flgB = CoreGene(model_location, "flgB", self.profile_factory)
        gene_flgB = ModelGene(c_gene_flgB, model_B)
        c_gene_tadZ = CoreGene(model_location, "tadZ", self.profile_factory)
        gene_tadZ = ModelGene(c_gene_tadZ, model_B)

        model_B.add_mandatory_gene(gene_sctn_flg)
        model_B.add_mandatory_gene(gene_sctj_flg)
        model_B.add_accessory_gene(gene_flgB)
        model_B.add_accessory_gene(gene_tadZ)

        ###########
        # Model A #
        ###########
        model_A = Model("foo/A", 10)
        c_gene_sctn = CoreGene(model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model_A)
        gene_sctn_hom = Exchangeable(c_gene_sctn_flg, gene_sctn)
        gene_sctn.add_exchangeable(gene_sctn_hom)

        c_gene_sctj = CoreGene(model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model_A)
        gene_sctj_an = Exchangeable(c_gene_sctj_flg, gene_sctj)
        gene_sctj.add_exchangeable(gene_sctj_an)

        c_gene_gspd = CoreGene(model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model_A)
        gene_gspd_an = Exchangeable(c_gene_flgB, gene_gspd)
        gene_gspd.add_exchangeable(gene_gspd_an)

        c_gene_abc = CoreGene(model_location, "abc", self.profile_factory)
        gene_abc = ModelGene(c_gene_abc, model_A, loner=True)
        gene_abc_ho = Exchangeable(c_gene_tadZ, gene_abc)
        gene_abc.add_exchangeable(gene_abc_ho)

        model_A.add_mandatory_gene(gene_sctn)
        model_A.add_mandatory_gene(gene_sctj)
        model_A.add_accessory_gene(gene_gspd)
        model_A.add_accessory_gene(gene_abc)



        #       CoreHit(gene, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        #                                                           pos      score
        h_sctj = CoreHit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        mh_sctj = ModelHit(h_sctj, gene_sctj, GeneStatus.MANDATORY)
        h_sctn = CoreHit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 2, 1.0, 1.0, 1.0, 1.0, 10, 20)
        mh_sctn = ModelHit(h_sctn, gene_sctn, GeneStatus.MANDATORY)
        h_gspd = CoreHit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 3, 1.0, 1.0, 1.0, 1.0, 10, 20)
        mh_gspd = ModelHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)

        h_sctj_flg = CoreHit(c_gene_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 10, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_flgB = CoreHit(c_gene_flgB, "hit_flgB", 803, "replicon_id", 11, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_abc = CoreHit(c_gene_abc, "hit_abc", 803, "replicon_id", 20, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_abc2 = CoreHit(c_gene_abc, "hit_abc2", 803, "replicon_id", 50, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_tadZ = CoreHit(c_gene_tadZ, "hit_tadZ", 803, "replicon_id", 40, 1.0, 1.0, 1.0, 1.0, 10, 20)
        mh_sctj_flg = ModelHit(h_sctj_flg, gene_sctj_flg, GeneStatus.MANDATORY)
        mh_flgB = ModelHit(h_flgB, gene_flgB, GeneStatus.ACCESSORY)
        mh_abc = ModelHit(h_abc, gene_abc, GeneStatus.ACCESSORY)
        mh_abc2 = ModelHit(h_abc2, gene_abc, GeneStatus.ACCESSORY)
        mh_tadZ = ModelHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY)

        model_A._min_mandatory_genes_required = 2
        model_A._min_genes_required = 2
        c1 = Cluster([mh_sctj, mh_sctn, mh_gspd], model_A, self.hit_weights)
        c2 = Cluster([mh_sctj, mh_sctn], model_A, self.hit_weights)
        c3 = Cluster([Loner(h_abc, gene_ref=gene_abc, gene_status=GeneStatus.ACCESSORY, counterpart=[mh_abc2])],
                     model_A, self.hit_weights)

        model_B._min_mandatory_genes_required = 1
        model_B._min_genes_required = 2
        c5 = Cluster([mh_sctj_flg, mh_tadZ, mh_flgB], model_B, self.hit_weights)

        sys_A = System(model_A, [c1, c2, c3], self.cfg.redundancy_penalty())
        # score =               2.5, 2 , 0.35 = 4.85 - (2 * 1.5) = 1.85

        sys_A.id = "sys_id_A"
        sys_B = System(model_B, [c5], self.cfg.redundancy_penalty())
        # score =                2.0
        sys_B.id = "sys_id_B"

        sol = Solution([sys_A, sys_B])
        sol_id = '12'

        hit_multi_sys_tracker = HitSystemTracker([sys_A, sys_B])
        sol_serializer = TsvSolutionSerializer()

        sol_tsv = '\t'.join([sol_id, 'replicon_id', 'hit_sctj', 'sctJ', '1', 'foo/A', 'sys_id_A',
                            '2', '1', '1.000', '1.850', '2', 'sctJ', 'mandatory',
                            '803', '1.0', '1.000', '1.000', '1.000', '10', '20', '', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id, 'replicon_id', 'hit_sctn', 'sctN', '2', 'foo/A', 'sys_id_A',
                             '2', '1', '1.000', '1.850', '2', 'sctN', 'mandatory',
                             '803', '1.0', '1.000', '1.000', '1.000', '10', '20', '', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id, 'replicon_id', 'hit_gspd', 'gspD', '3', 'foo/A', 'sys_id_A',
                             '2', '1', '1.000', '1.850', '2', 'gspD', 'accessory',
                             '803', '1.0', '1.000', '1.000', '1.000', '10', '20', '', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id, 'replicon_id', 'hit_sctj', 'sctJ', '1', 'foo/A', 'sys_id_A',
                             '2', '2', '1.000', '1.850', '2', 'sctJ', 'mandatory',
                             '803', '1.0', '1.000', '1.000', '1.000', '10', '20', '', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id, 'replicon_id', 'hit_sctn', 'sctN', '2', 'foo/A', 'sys_id_A',
                             '2', '2', '1.000', '1.850', '2', 'sctN', 'mandatory',
                             '803', '1.0', '1.000', '1.000', '1.000', '10', '20', '', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id, 'replicon_id', 'hit_abc', 'abc', '20', 'foo/A', 'sys_id_A',
                             '2', '-1', '1.000', '1.850', '2', 'abc', 'accessory',
                             '803', '1.0', '1.000', '1.000', '1.000', '10', '20', 'hit_abc2', ''])
        sol_tsv += "\n"
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id, 'replicon_id', 'hit_sctj_flg', 'sctJ_FLG', '10', 'foo/B', 'sys_id_B',
                             '1', '1', '0.750', '2.000', '1', 'sctJ_FLG', 'mandatory',
                             '803', '1.0', '1.000', '1.000', '1.000', '10', '20', '', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id, 'replicon_id', 'hit_flgB', 'flgB', '11', 'foo/B', 'sys_id_B',
                              '1', '1', '0.750', '2.000', '1', 'flgB', 'accessory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', '', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id, 'replicon_id', 'hit_tadZ', 'tadZ', '40', 'foo/B', 'sys_id_B',
                              '1', '1', '0.750', '2.000', '1', 'tadZ', 'accessory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', '', ''])
        sol_tsv += "\n"
        sol_tsv += "\n"
        ser = sol_serializer.serialize(sol, sol_id, hit_multi_sys_tracker)
        self.maxDiff = None
        self.assertEqual(ser, sol_tsv)


    def test_LikelySystemSerializer_txt(self):
        model = Model("foo/FOO", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model)
        model.add_accessory_gene(gene_sctn)
        c_gene_abc = CoreGene(self.model_location, "abc", self.profile_factory)
        gene_abc = ModelGene(c_gene_abc, model)
        model.add_forbidden_gene(gene_abc)

        hit_1 = CoreHit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ModelHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = CoreHit(c_gene_sctj, "hit_2", 803, "replicon_id", 2, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ModelHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = CoreHit(c_gene_sctn, "hit_3", 803, "replicon_id", 3, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ModelHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        hit_4 = CoreHit(c_gene_abc, "hit_4", 803, "replicon_id", 4, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_4 = ModelHit(hit_4, gene_abc, GeneStatus.FORBIDDEN)

        ls_1 = LikelySystem(model, [v_hit_1], [v_hit_2, v_hit_3], [], [v_hit_4])
        hit_multi_sys_tracker = HitSystemTracker([ls_1])
        ser = TxtLikelySystemSerializer()

        txt = ser.serialize(ls_1, hit_multi_sys_tracker)
        expected_txt = """This replicon contains genetic materials needed for system foo/FOO
WARNING there quorum is reached but there is also some forbidden genes.

system id = replicon_id_FOO_1
model = foo/FOO
replicon = replicon_id
hits = [('hit_1', 'gspD', 1), ('hit_2', 'sctJ', 2), ('hit_3', 'sctN', 3), ('hit_4', 'abc', 4)]
wholeness = 1.000

mandatory genes:
\t- gspD: 1 (gspD)

accessory genes:
\t- sctJ: 1 (sctJ)
\t- sctN: 1 (sctN)

neutral genes:

forbidden genes:
\t- abc: 1 (abc)

Use ordered replicon to have better prediction.
"""
        self.assertEqual(txt, expected_txt)


    def test_UnlikelySystemSerializer_txt(self):
        model = Model("foo/FOO", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model)
        model.add_accessory_gene(gene_sctn)
        c_gene_abc = CoreGene(self.model_location, "abc", self.profile_factory)
        gene_abc = ModelGene(c_gene_abc, model)
        model.add_forbidden_gene(gene_abc)

        hit_1 = CoreHit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ModelHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = CoreHit(c_gene_sctj, "hit_2", 803, "replicon_id", 2, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ModelHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = CoreHit(c_gene_sctn, "hit_3", 803, "replicon_id", 3, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ModelHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        hit_4 = CoreHit(c_gene_abc, "hit_4", 803, "replicon_id", 4, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_4 = ModelHit(hit_4, gene_abc, GeneStatus.FORBIDDEN)
        ser = TxtUnikelySystemSerializer()

        ls_1 = UnlikelySystem(model, [v_hit_1], [v_hit_2, v_hit_3], [], [v_hit_4], ["the reason why"])
        txt = ser.serialize(ls_1)
        expected_txt = """This replicon probably not contains a system foo/FOO:
the reason why

system id = replicon_id_FOO_1
model = foo/FOO
replicon = replicon_id
hits = [('hit_1', 'gspD', 1), ('hit_2', 'sctJ', 2), ('hit_3', 'sctN', 3), ('hit_4', 'abc', 4)]
wholeness = 1.000

mandatory genes:
\t- gspD: 1 (gspD)

accessory genes:
\t- sctJ: 1 (sctJ)
\t- sctN: 1 (sctN)

neutral genes:

forbidden genes:
\t- abc: 1 (abc)

Use ordered replicon to have better prediction.
"""
        self.assertEqual(txt, expected_txt)


    def test_SpecialHitSerializer_tsv(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        cfg = Config(MacsyDefaults(), args)

        model_name = 'foo'
        models_location = ModelLocation(path=os.path.join(args.models_dir, model_name))

        # we need to reset the ProfileFactory
        # because it's a like a singleton
        # so other tests are influenced by ProfileFactory and it's configuration
        # for instance search_genes get profile without hmmer_exe
        profile_factory = ProfileFactory(cfg)
        model = Model("foo/T2SS", 10)

        gene_name = "gspD"
        cg_gspd = CoreGene(models_location, gene_name, profile_factory)
        mg_gspd = ModelGene(cg_gspd, model, loner=True)

        gene_name = "sctJ"
        cg_sctj = CoreGene(models_location, gene_name, profile_factory)
        mg_sctj = ModelGene(cg_sctj, model)

        gene_name = "abc"
        cg_abc = CoreGene(models_location, gene_name, profile_factory)
        mg_abc = ModelGene(cg_abc, model)

        model.add_mandatory_gene(mg_gspd)
        model.add_accessory_gene(mg_sctj)
        model.add_accessory_gene(mg_abc)

        chit_abc = CoreHit(cg_abc, "hit_abc", 803, "replicon_id", 3, 1.0, 1.0, 1.0, 1.0, 10, 20)
        chit_sctj = CoreHit(cg_sctj, "hit_sctj", 803, "replicon_id", 4, 1.0, 1.0, 1.0, 1.0, 10, 20)
        chit_gspd1 = CoreHit(cg_gspd, "hit_gspd1", 803, "replicon_id", 20, 1.0, 2.0, 1.0, 1.0, 10, 20)
        chit_gspd2 = CoreHit(cg_gspd, "hit_gspd2", 803, "replicon_id", 30, 1.0, 3.0, 1.0, 1.0, 10, 20)
        mhit_abc = ModelHit(chit_abc, mg_abc, GeneStatus.ACCESSORY)
        mhit_sctj = ModelHit(chit_sctj, mg_sctj, GeneStatus.ACCESSORY)
        mhit_gspd1 = ModelHit(chit_gspd1, mg_gspd, GeneStatus.MANDATORY)
        mhit_gspd2 = ModelHit(chit_gspd2, mg_gspd, GeneStatus.MANDATORY)
        l_gspd1 = Loner(mhit_gspd1, counterpart=[mhit_gspd2])
        l_gspd2 = Loner(mhit_gspd2, counterpart=[mhit_gspd1])
        ser = TsvSpecialHitSerializer()
        txt = ser.serialize([l_gspd1, l_gspd2])

        expected_txt = "\t".join(['replicon', 'model_fqn', 'function', 'gene_name', 'hit_id', 'hit_pos', 'hit_status',
                                'hit_seq_len', 'hit_i_eval', 'hit_score', 'hit_profile_cov',
                                'hit_seq_cov', 'hit_begin_match', 'hit_end_match'])
        expected_txt += "\n"
        expected_txt += "\t".join(['replicon_id', 'foo/T2SS', 'gspD', 'gspD', 'hit_gspd1', '20', 'mandatory', '803',
                                '1.000e+00', '2.000', '1.000', '1.000', '10', '20'])
        expected_txt += "\n"
        expected_txt += "\t".join(['replicon_id', 'foo/T2SS', 'gspD', 'gspD', 'hit_gspd2', '30', 'mandatory', '803',
                                '1.000e+00', '3.000', '1.000', '1.000', '10', '20'])
        expected_txt += "\n"
        self.maxDiff = None
        self.assertEqual(txt, expected_txt)


    def test_RejectedClusters_tsv(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.res_search_dir = "blabla"

        cfg = Config(MacsyDefaults(), args)
        model_name = 'foo'
        models_location = ModelLocation(path=os.path.join(args.models_dir, model_name))
        profile_factory = ProfileFactory(cfg)

        model = Model("foo/T2SS", 11)

        gene_name = "gspD"
        c_gene_gspd = CoreGene(models_location, gene_name, profile_factory)
        gene_1 = ModelGene(c_gene_gspd, model)
        gene_name = "sctC"
        c_gene_sctc = CoreGene(models_location, gene_name, profile_factory)
        gene_2 = ModelGene(c_gene_sctc, model)
        gene_name = 'tadZ'
        c_gene_tadz = CoreGene(models_location, gene_name, profile_factory)
        gene_3 = ModelGene(c_gene_tadz, model)
        gene_name = 'abc'
        c_gene_abc = CoreGene(models_location, gene_name, profile_factory)
        gene_4 = Exchangeable(c_gene_abc, gene_3)
        gene_3.add_exchangeable(gene_4)

        model.add_mandatory_gene(gene_1)
        model.add_accessory_gene(gene_2)
        model.add_accessory_gene(gene_3)

        #     CoreHit(gene, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = CoreHit(c_gene_gspd, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        mh_10 = ModelHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = CoreHit(c_gene_sctc, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        mh_20 = ModelHit(h20, gene_2, GeneStatus.ACCESSORY)
        h40 = CoreHit(c_gene_gspd, "h10", 10, "replicon_1", 40, 1.0, 10.0, 1.0, 1.0, 10, 20)
        mh_40 = ModelHit(h40, gene_1, GeneStatus.MANDATORY)
        h50 = CoreHit(c_gene_sctc, "h20", 10, "replicon_1", 50, 1.0, 20.0, 1.0, 1.0, 10, 20)
        mh_50 = ModelHit(h50, gene_2, GeneStatus.ACCESSORY)
        h60 = CoreHit(c_gene_gspd, "h60", 10, "replicon_1", 60, 1.0, 10.0, 1.0, 1.0, 10, 20)
        mh_60 = ModelHit(h60, gene_1, GeneStatus.MANDATORY)
        h70 = CoreHit(c_gene_abc, "h70", 10, "replicon_1", 70, 1.0, 10.0, 1.0, 1.0, 10, 20)
        mh_70 = ModelHit(h70, gene_4, GeneStatus.ACCESSORY)

        hit_weights = HitWeight(**cfg.hit_weights())
        c1 = Cluster([mh_10, mh_20], model, hit_weights)
        c2 = Cluster([mh_40, mh_50], model, hit_weights)
        c3 = Cluster([mh_60, mh_70], model, hit_weights)

        r_c1 = RejectedCandidate(model, [c1, c2], ["The reasons to reject this candidate"])

        r_c2 = RejectedCandidate(model, [c2, c3], ["reason One", "reason Two"])

        model_fam_name = 'foo'
        model_vers = '0.0b2'

        ser = TsvRejectedCandidatesSerializer()
        tsv = ser.serialize([r_c1, r_c2])

        expected_tsv = '\t'.join(['candidate_id', 'replicon', 'model_fqn', 'cluster_id', 'hit_id', 'hit_pos', 'gene_name', 'function', 'reasons'])
        expected_tsv += '\n'
        expected_tsv += '\t'.join(['replicon_1_T2SS_1', 'replicon_1', 'foo/T2SS', c1.id, 'h10', '10', 'gspD', 'gspD', 'The reasons to reject this candidate'])
        expected_tsv += '\n'
        expected_tsv += '\t'.join(['replicon_1_T2SS_1', 'replicon_1', 'foo/T2SS', c1.id, 'h20', '20', 'sctC', 'sctC', 'The reasons to reject this candidate'])
        expected_tsv += '\n'
        expected_tsv += '\t'.join(['replicon_1_T2SS_1', 'replicon_1', 'foo/T2SS', c2.id, 'h10', '40', 'gspD', 'gspD', 'The reasons to reject this candidate'])
        expected_tsv += '\n'
        expected_tsv += '\t'.join(['replicon_1_T2SS_1', 'replicon_1', 'foo/T2SS', c2.id, 'h20', '50', 'sctC', 'sctC', 'The reasons to reject this candidate'])
        expected_tsv += '\n'
        expected_tsv += '\n'
        expected_tsv += '\t'.join(['replicon_1_T2SS_2', 'replicon_1', 'foo/T2SS', c2.id, 'h10', '40', 'gspD', 'gspD', 'reason One/reason Two'])
        expected_tsv += '\n'
        expected_tsv += '\t'.join(['replicon_1_T2SS_2', 'replicon_1', 'foo/T2SS', c2.id, 'h20', '50', 'sctC', 'sctC', 'reason One/reason Two'])
        expected_tsv += '\n'
        expected_tsv += '\t'.join(['replicon_1_T2SS_2', 'replicon_1', 'foo/T2SS', c3.id, 'h60', '60', 'gspD', 'gspD', 'reason One/reason Two'])
        expected_tsv += '\n'
        expected_tsv += '\t'.join(['replicon_1_T2SS_2', 'replicon_1', 'foo/T2SS', c3.id, 'h70', '70', 'abc', 'tadZ', 'reason One/reason Two'])
        expected_tsv += '\n'
        expected_tsv += '\n'

        self.maxDiff = None
        self.assertEqual(expected_tsv, tsv)
