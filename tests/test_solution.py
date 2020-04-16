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
from macsypy.system import System
from macsypy.solution import Solution, compute_max_bound, find_best_solution
from tests import MacsyTest


def _build_systems(cfg, profile_factory):
    model_name = 'foo'
    model_location = ModelLocation(path=os.path.join(cfg.models_dir(), model_name))
    model_1 = Model("foo/A", 10)
    model_2 = Model("foo/B", 10)
    model_3 = Model("foo/C", 10)
    model_4 = Model("foo/D", 10)
    model_5 = Model("foo/E", 10)

    c_gene_sctn_flg = CoreGene(model_location, "sctN_FLG", profile_factory)
    gene_sctn_flg = ModelGene(c_gene_sctn_flg, model_2)
    c_gene_sctj_flg = CoreGene(model_location, "sctJ_FLG", profile_factory)
    gene_sctj_flg = ModelGene(c_gene_sctj_flg, model_2)
    c_gene_flgB = CoreGene(model_location, "flgB", profile_factory)
    gene_flgB = ModelGene(c_gene_flgB, model_2)
    c_gene_tadZ = CoreGene(model_location, "tadZ", profile_factory)
    gene_tadZ = ModelGene(c_gene_tadZ, model_2)

    c_gene_sctn = CoreGene(model_location, "sctN", profile_factory)
    gene_sctn = ModelGene(c_gene_sctn, model_1)
    gene_sctn_hom = Exchangeable(c_gene_sctn_flg, gene_sctn)
    gene_sctn.add_exchangeable(gene_sctn_hom)

    c_gene_sctj = CoreGene(model_location, "sctJ", profile_factory)
    gene_sctj = ModelGene(c_gene_sctj, model_1)
    gene_sctj_an = Exchangeable(c_gene_sctj_flg, gene_sctj)
    gene_sctj.add_exchangeable(gene_sctj_an)

    c_gene_gspd = CoreGene(model_location, "gspD", profile_factory)
    gene_gspd = ModelGene(c_gene_gspd, model_1)
    gene_gspd_an = Exchangeable(c_gene_flgB, gene_gspd)
    gene_gspd.add_exchangeable(gene_gspd_an)

    c_gene_abc = CoreGene(model_location, "abc", profile_factory)
    gene_abc = ModelGene(c_gene_abc, model_1)
    gene_abc_ho = Exchangeable(c_gene_tadZ, gene_abc)
    gene_abc.add_exchangeable(gene_abc_ho)

    model_1.add_mandatory_gene(gene_sctn)
    model_1.add_mandatory_gene(gene_sctj)
    model_1.add_accessory_gene(gene_gspd)
    model_1.add_forbidden_gene(gene_abc)

    model_2.add_mandatory_gene(gene_sctn_flg)
    model_2.add_mandatory_gene(gene_sctj_flg)
    model_2.add_accessory_gene(gene_flgB)
    model_2.add_accessory_gene(gene_tadZ)

    model_3.add_mandatory_gene(gene_sctn_flg)
    model_3.add_mandatory_gene(gene_sctj_flg)
    model_3.add_mandatory_gene(gene_flgB)
    model_3.add_accessory_gene(gene_tadZ)
    model_3.add_accessory_gene(gene_gspd)

    model_4.add_mandatory_gene(gene_abc)
    model_4.add_accessory_gene(gene_sctn)

    model_5.add_accessory_gene(gene_gspd)

    h_sctj = Hit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
    h_sctn = Hit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
    h_gspd = Hit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

    h_sctj_flg = Hit(c_gene_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
    h_flgB = Hit(c_gene_flgB, "hit_flgB", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
    h_tadZ = Hit(c_gene_tadZ, "hit_tadZ", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

    h_abc = Hit(c_gene_abc, "hit_abc", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

    model_1._min_mandatory_genes_required = 2
    model_1._min_genes_required = 2
    c1 = Cluster([ValidHit(h_sctj, gene_sctj, GeneStatus.MANDATORY),
                  ValidHit(h_sctn, gene_sctn, GeneStatus.MANDATORY),
                  ValidHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
                  ],
                 model_1)

    c2 = Cluster([ValidHit(h_sctj, gene_sctj, GeneStatus.MANDATORY),
                  ValidHit(h_sctn, gene_sctn, GeneStatus.MANDATORY)],
                 model_1)

    model_2._min_mandatory_genes_required = 1
    model_2._min_genes_required = 2
    c3 = Cluster([ValidHit(h_sctj_flg, gene_sctj_flg, GeneStatus.MANDATORY),
                  ValidHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY),
                  ValidHit(h_flgB, gene_flgB, GeneStatus.ACCESSORY)],
                 model_2)
    model_3._min_mandatory_genes_required = 1
    model_3._min_genes_required = 2
    c4 = Cluster([ValidHit(h_sctj_flg, gene_sctj_flg, GeneStatus.MANDATORY),
                  ValidHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY),
                  ValidHit(h_flgB, gene_flgB, GeneStatus.MANDATORY),
                  ValidHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)],
                 model_3)
    model_4._min_mandatory_genes_required = 1
    model_4._min_genes_required = 1
    c5 = Cluster([ValidHit(h_abc, gene_abc, GeneStatus.MANDATORY),
                  ValidHit(h_sctn, gene_sctn, GeneStatus.ACCESSORY)],
                 model_4)

    model_5._min_mandatory_genes_required = 0
    model_5._min_genes_required = 1
    c6 = Cluster([ValidHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)],
                 model_5)
    s0 = System(model_1, [c1, c2])
    # we need to tweek the replicon_id to have stable ressults
    # whatever the number of tests ran
    # or the tests order
    s0.id = "replicon_id_A_0"
    s1 = System(model_2, [c3])
    s1.id = "replicon_id_B_1"
    s2 = System(model_3, [c4])
    s2.id = "replicon_id_C_2"
    s3 = System(model_4, [c5])
    s3.id = "replicon_id_D_3"
    s4 = System(model_5, [c6])
    s4.id = "replicon_id_E_4"

    return s0, s1, s2, s3, s4


class SolutionTest(MacsyTest):

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        self.cfg = Config(MacsyDefaults(), args)

        # we need to reset the ProfileFactory
        # because it's a like a singleton
        # so other tests are influenced by ProfileFactory and it's configuration
        # for instance search_genes get profile without hmmer_exe
        self.profile_factory = ProfileFactory(self.cfg)
        self.systems = _build_systems(self.cfg, self.profile_factory)

    def test_hits(self):
        systems = self.systems[0:2]
        sol = Solution(systems)
        hits = set([h.hit for h in systems[0].hits] + [h.hit for h in systems[1].hits])
        self.assertSetEqual(hits, sol.hits)

    def test_iadd(self):
        sol = Solution([self.systems[0]])
        sol1 = sol
        sol += self.systems[1]
        self.assertTrue(isinstance(sol, Solution))
        sol2 = Solution(self.systems[0:2])
        self.assertSetEqual(sol.hits, sol2.hits)
        self.assertEqual(sol.score, sol2.score)

    def test_eq(self):
        sol_0 = Solution([self.systems[0]])
        sol_0_bis = Solution([self.systems[0]])
        sol_1 = Solution([self.systems[1]])
        self.assertEqual(sol_0, sol_0_bis)
        self.assertNotEqual(sol_0, sol_1)

    def test_str(self):
        systems = self.systems[0:2]
        sol = Solution(systems)
        expec_sol_str = """Score of the solution = 3.5
Sys_ID=replicon_id_A_0 Score=1.5
Sys_ID=replicon_id_B_1 Score=2.0"""
        self.assertEqual(str(sol),  expec_sol_str)

    def test_score(self):
        systems = self.systems[0:2]
        sol = Solution(systems)
        self.assertEqual(sum(s.score for s in systems), sol.score)

    def test_is_compatible(self):
        sol = Solution([self.systems[0]])
        self.assertTrue(sol.is_compatible(self.systems[1]))
        self.assertFalse(sol.is_compatible(self.systems[2]))


class SolutionFuncTest(MacsyTest):

    @classmethod
    def setUpClass(cls) -> None:
        # to turn on debugging
        # uncomment the 3 following lines
        # import macsypy
        # macsypy.init_logger()
        # macsypy.logger_set_level('DEBUG')

        pass

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        self.cfg = Config(MacsyDefaults(), args)
        # we need to reset the ProfileFactory
        # because it's a like a singleton
        # so other tests are influenced by ProfileFactory and it's configuration
        # for instance search_genes get profile without hmmer_exe
        self.profile_factory = ProfileFactory(self.cfg)
        self.systems = _build_systems(self.cfg, self.profile_factory)


    def test_compute_max_bound(self):
        self.assertEqual(compute_max_bound([self.systems[0]]), self.systems[0].score)
        self.assertEqual(compute_max_bound(self.systems[:-1]), 8.0)


    def test_find_best_solution(self):
        sorted_systems = sorted(self.systems, key=lambda s: (- s.score, s.id))

        sorted_syst = sorted_systems[:-1]
        # sorted_syst = [('replicon_id_C_2', 3.0), ('replicon_id_B_1', 2.0), ('replicon_id_A_0', 1.5), ('replicon_id_D_3', 1.5)]
        # replicon_id_C_2 ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_B_1 ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB']
        # replicon_id_A_0 ['hit_sctj', 'hit_sctn', 'hit_gspd', 'hit_sctj', 'hit_sctn']
        # replicon_id_D_3 ['hit_abc', 'hit_sctn']
        # C and D are compatible
        # B and A are compatible
        # So the best Solution expected is C D
        best_sol = find_best_solution(sorted_syst, Solution([]), Solution([]))
        self.assertEqual(best_sol, Solution([sorted_syst[0], sorted_syst[-1]]))

        sorted_syst = sorted_systems[:-2]
        # sorted_syst = [('replicon_id_C_2', 3.0), ('replicon_id_B_1', 2.0), ('replicon_id_A_0', 1.5)]
        # replicon_id_C_2 ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_B_1 ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB']
        # replicon_id_A_0 ['hit_sctj', 'hit_sctn', 'hit_gspd', 'hit_sctj', 'hit_sctn']
        # C is alone
        # B and A are compatible
        # So the best Solution expected is B and A
        best_sol = find_best_solution(sorted_syst, Solution([]), Solution([]))
        self.assertEqual(best_sol, Solution(sorted_syst[1:]))

        sorted_syst = [sorted_systems[0],
                       sorted_systems[1],
                       sorted_systems[-1]]
        # sorted_syst = [('replicon_id_C_2', 3.0), ('replicon_id_B_1', 2.0), ('replicon_id_E_4', 0.5)]
        # replicon_id_C_2 ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_B_1 ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB']
        # replicon_id_E_4 ['hit_gspd']
        # C is alone
        # B and E are compatible
        # So the best Solution expected is C
        # but the branch B E is stoped as the max bound is lower than best_solution C
        best_sol = find_best_solution(sorted_syst, Solution([]), Solution([]))
        self.assertEqual(best_sol, Solution([sorted_syst[0]]))
