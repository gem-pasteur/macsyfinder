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

from macsypy.hit import Hit, ValidHit, HitWeight
from macsypy.config import Config, MacsyDefaults
from macsypy.gene import CoreGene, ModelGene, Exchangeable, GeneStatus
from macsypy.profile import ProfileFactory
from macsypy.model import Model
from macsypy.registries import ModelLocation
from macsypy.cluster import Cluster
from macsypy.system import System
from macsypy.solution import find_best_solutions
from tests import MacsyTest


def _build_systems(cfg, profile_factory):
    model_name = 'foo'
    model_location = ModelLocation(path=os.path.join(cfg.models_dir(), model_name))
    model_A = Model("foo/A", 10)
    model_B = Model("foo/B", 10)
    model_C = Model("foo/C", 10)
    model_D = Model("foo/D", 10)
    model_E = Model("foo/E", 10)
    model_F = Model("foo/F", 10)
    model_G = Model("foo/G", 10)
    model_H = Model("foo/H", 10)
    model_I = Model("foo/I", 10)
    model_J = Model("foo/J", 10)
    model_K = Model("foo/K", 10)

    c_gene_sctn_flg = CoreGene(model_location, "sctN_FLG", profile_factory)
    gene_sctn_flg = ModelGene(c_gene_sctn_flg, model_B)
    c_gene_sctj_flg = CoreGene(model_location, "sctJ_FLG", profile_factory)
    gene_sctj_flg = ModelGene(c_gene_sctj_flg, model_B)
    c_gene_flgB = CoreGene(model_location, "flgB", profile_factory)
    gene_flgB = ModelGene(c_gene_flgB, model_B)
    c_gene_tadZ = CoreGene(model_location, "tadZ", profile_factory)
    gene_tadZ = ModelGene(c_gene_tadZ, model_B)

    c_gene_sctn = CoreGene(model_location, "sctN", profile_factory)
    gene_sctn = ModelGene(c_gene_sctn, model_A)
    gene_sctn_hom = Exchangeable(c_gene_sctn_flg, gene_sctn)
    gene_sctn.add_exchangeable(gene_sctn_hom)

    c_gene_sctj = CoreGene(model_location, "sctJ", profile_factory)
    gene_sctj = ModelGene(c_gene_sctj, model_A)
    gene_sctj_an = Exchangeable(c_gene_sctj_flg, gene_sctj)
    gene_sctj.add_exchangeable(gene_sctj_an)

    c_gene_gspd = CoreGene(model_location, "gspD", profile_factory)
    gene_gspd = ModelGene(c_gene_gspd, model_A)
    gene_gspd_an = Exchangeable(c_gene_flgB, gene_gspd)
    gene_gspd.add_exchangeable(gene_gspd_an)

    c_gene_abc = CoreGene(model_location, "abc", profile_factory)
    gene_abc = ModelGene(c_gene_abc, model_A)
    gene_abc_ho = Exchangeable(c_gene_tadZ, gene_abc)
    gene_abc.add_exchangeable(gene_abc_ho)

    c_gene_sctc = CoreGene(model_location, "sctC", profile_factory)
    gene_sctc = ModelGene(c_gene_sctc, model_J)

    model_A.add_mandatory_gene(gene_sctn)
    model_A.add_mandatory_gene(gene_sctj)
    model_A.add_accessory_gene(gene_gspd)
    model_A.add_forbidden_gene(gene_abc)

    model_B.add_mandatory_gene(gene_sctn_flg)
    model_B.add_mandatory_gene(gene_sctj_flg)
    model_B.add_accessory_gene(gene_flgB)
    model_B.add_accessory_gene(gene_tadZ)

    model_C.add_mandatory_gene(gene_sctn_flg)
    model_C.add_mandatory_gene(gene_sctj_flg)
    model_C.add_mandatory_gene(gene_flgB)
    model_C.add_accessory_gene(gene_tadZ)
    model_C.add_accessory_gene(gene_gspd)

    model_D.add_mandatory_gene(gene_abc)
    model_D.add_accessory_gene(gene_sctn)

    model_E.add_accessory_gene(gene_gspd)

    model_F.add_mandatory_gene(gene_abc)

    # model G idem as C
    model_G.add_mandatory_gene(gene_sctn_flg)
    model_G.add_mandatory_gene(gene_sctj_flg)
    model_G.add_mandatory_gene(gene_flgB)
    model_G.add_accessory_gene(gene_tadZ)
    model_G.add_accessory_gene(gene_gspd)

    # mode lH idem as D
    model_H.add_mandatory_gene(gene_abc)
    model_H.add_accessory_gene(gene_sctn)

    # model I
    model_I.add_mandatory_gene(gene_flgB)
    model_I.add_accessory_gene(gene_tadZ)

    # model J
    model_J.add_mandatory_gene(gene_abc)
    model_J.add_mandatory_gene(gene_gspd)
    model_J.add_accessory_gene(gene_tadZ)
    model_J.add_accessory_gene(gene_sctc)

    # model K
    model_K.add_mandatory_gene(gene_flgB)
    model_K.add_mandatory_gene(gene_sctn_flg)
    model_K.add_accessory_gene(gene_sctj_flg)
    model_K.add_accessory_gene(gene_sctn)


    h_sctj = Hit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
    h_sctn = Hit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 2, 1.0, 1.0, 1.0, 1.0, 10, 20)
    h_gspd = Hit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 3, 1.0, 1.0, 1.0, 1.0, 10, 20)

    h_sctj_flg = Hit(c_gene_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 4, 1.0, 1.0, 1.0, 1.0, 10, 20)
    h_flgB = Hit(c_gene_flgB, "hit_flgB", 803, "replicon_id", 5, 1.0, 1.0, 1.0, 1.0, 10, 20)
    h_tadZ = Hit(c_gene_tadZ, "hit_tadZ", 803, "replicon_id", 6, 1.0, 1.0, 1.0, 1.0, 10, 20)

    h_abc = Hit(c_gene_abc, "hit_abc", 803, "replicon_id", 7, 1.0, 1.0, 1.0, 1.0, 10, 20)

    model_A._min_mandatory_genes_required = 2
    model_A._min_genes_required = 2
    hit_weights = HitWeight(**cfg.hit_weights())
    c1 = Cluster([ValidHit(h_sctj, gene_sctj, GeneStatus.MANDATORY),
                  ValidHit(h_sctn, gene_sctn, GeneStatus.MANDATORY),
                  ValidHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
                  ],
                 model_A, hit_weights)

    c2 = Cluster([ValidHit(h_sctj, gene_sctj, GeneStatus.MANDATORY),
                  ValidHit(h_sctn, gene_sctn, GeneStatus.MANDATORY)],
                 model_A, hit_weights)

    model_B._min_mandatory_genes_required = 1
    model_B._min_genes_required = 2
    c3 = Cluster([ValidHit(h_sctj_flg, gene_sctj_flg, GeneStatus.MANDATORY),
                  ValidHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY),
                  ValidHit(h_flgB, gene_flgB, GeneStatus.ACCESSORY)],
                 model_B, hit_weights)
    model_C._min_mandatory_genes_required = 1
    model_C._min_genes_required = 2
    c4 = Cluster([ValidHit(h_sctj_flg, gene_sctj_flg, GeneStatus.MANDATORY),
                  ValidHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY),
                  ValidHit(h_flgB, gene_flgB, GeneStatus.MANDATORY),
                  ValidHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)],
                 model_C, hit_weights)
    model_D._min_mandatory_genes_required = 1
    model_D._min_genes_required = 1
    c5 = Cluster([ValidHit(h_abc, gene_abc, GeneStatus.MANDATORY),
                  ValidHit(h_sctn, gene_sctn, GeneStatus.ACCESSORY)],
                 model_D, hit_weights)

    model_E._min_mandatory_genes_required = 0
    model_E._min_genes_required = 1
    c6 = Cluster([ValidHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)],
                 model_E, hit_weights)

    model_F._min_mandatory_genes_required = 1
    model_F._min_genes_required = 1
    c7 = Cluster([ValidHit(h_abc, gene_abc, GeneStatus.MANDATORY)],
                 model_F, hit_weights)

    model_H._min_mandatory_genes_required = 1
    model_H._min_genes_required = 1

    model_I._min_mandatory_genes_required = 1
    model_I._min_genes_required = 1
    c8 = Cluster([ValidHit(h_flgB, gene_flgB, GeneStatus.MANDATORY),
                  ValidHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY)],
                 model_I, hit_weights)

    model_J._min_mandatory_genes_required = 1
    model_J._min_genes_required = 1
    c9 = Cluster([ValidHit(h_abc, gene_abc, GeneStatus.MANDATORY),
                  ValidHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY)],
                 model_I, hit_weights)

    model_K._min_mandatory_genes_required = 1
    model_K._min_genes_required = 1
    c10 = Cluster([ValidHit(h_flgB, gene_flgB, GeneStatus.MANDATORY),
                  ValidHit(h_sctn, gene_sctn, GeneStatus.ACCESSORY)],
                  model_K, hit_weights)

    systems = {}

    systems['A'] = System(model_A, [c1, c2], cfg.redundancy_penalty())  # 5 hits
    # we need to tweek the replicon_id to have stable ressults
    # whatever the number of tests ran
    # or the tests order
    systems['A'].id = "replicon_id_A"
    systems['B'] = System(model_B, [c3], cfg.redundancy_penalty())  # 3 hits
    systems['B'].id = "replicon_id_B"
    systems['C'] = System(model_C, [c4], cfg.redundancy_penalty())  # 4 hits
    systems['C'].id = "replicon_id_C"
    systems['D'] = System(model_D, [c5], cfg.redundancy_penalty())  # 2 hits
    systems['D'].id = "replicon_id_D"
    systems['E'] = System(model_E, [c6], cfg.redundancy_penalty())  # 1 hit
    systems['E'].id = "replicon_id_E"
    systems['F'] = System(model_F, [c7], cfg.redundancy_penalty())  # 1 hit
    systems['F'].id = "replicon_id_F"
    systems['G'] = System(model_G, [c4], cfg.redundancy_penalty())  # 4 hits
    systems['G'].id = "replicon_id_G"
    systems['H'] = System(model_H, [c5], cfg.redundancy_penalty())  # 2 hits
    systems['H'].id = "replicon_id_H"
    systems['I'] = System(model_I, [c8], cfg.redundancy_penalty())  # 2 hits
    systems['I'].id = "replicon_id_I"
    systems['J'] = System(model_J, [c9], cfg.redundancy_penalty())  # 2 hits
    systems['J'].id = "replicon_id_J"
    systems['K'] = System(model_K, [c10], cfg.redundancy_penalty())  # 2 hits
    systems['K'].id = "replicon_id_K"

    return systems


class SolutionExplorerTest(MacsyTest):

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
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        self.cfg = Config(MacsyDefaults(), args)
        # we need to reset the ProfileFactory
        # because it's a like a singleton
        # so other tests are influenced by ProfileFactory and it's configuration
        # for instance search_genes get profile without hmmer_exe
        self.profile_factory = ProfileFactory(self.cfg)
        self.systems = _build_systems(self.cfg, self.profile_factory)


    def test_find_best_solution(self):
        systems = [self.systems[k] for k in 'ABCD']
        sorted_syst = sorted(systems, key=lambda s: (- s.score, s.id))
        # sorted_syst = [('replicon_id_C', 3.0), ('replicon_id_B', 2.0), ('replicon_id_A', 1.5), ('replicon_id_D', 1.5)]
        # replicon_id_C ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_B ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB']
        # replicon_id_A ['hit_sctj', 'hit_sctn', 'hit_gspd', 'hit_sctj', 'hit_sctn']
        # replicon_id_D ['hit_abc', 'hit_sctn']
        # C and D are compatible 4.5
        # B and A are compatible 3.5
        # B and D are compatible 3.5
        # So the best Solution expected is C D 4.5
        best_sol, score = find_best_solutions(sorted_syst)
        expected_sol = [[self.systems[k] for k in 'CD']]
        # The order of solutions are not relevant
        # The order of systems in each solutions are not relevant
        # transform list in set to compare them
        best_sol = {frozenset(sol) for sol in best_sol}
        expected_sol = {frozenset(sol) for sol in expected_sol}
        self.assertEqual(score, 4.5)
        self.assertSetEqual(best_sol, expected_sol)

        systems = [self.systems[k] for k in 'ABC']
        sorted_syst = sorted(systems, key=lambda s: (- s.score, s.id))
        # sorted_syst = [('replicon_id_C', 3.0), ('replicon_id_B', 2.0), ('replicon_id_A', 1.5)]
        # replicon_id_C ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_B ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB']
        # replicon_id_A ['hit_sctj', 'hit_sctn', 'hit_gspd', 'hit_sctj', 'hit_sctn']
        # C is alone 3.0
        # B and A are compatible 3.5
        # So the best Solution expected is B and A
        best_sol, score = find_best_solutions(sorted_syst)
        expected_sol = [[self.systems[k] for k in 'BA']]
        best_sol = {frozenset(sol) for sol in best_sol}
        expected_sol = {frozenset(sol) for sol in expected_sol}
        self.assertEqual(score, 3.5)
        self.assertSetEqual(best_sol, expected_sol)

        systems = [self.systems[k] for k in 'BCE']
        sorted_syst = sorted(systems, key=lambda s: (- s.score, s.id))
        # sorted_syst = [('replicon_id_C', 3.0), ('replicon_id_B', 2.0), ('replicon_id_E', 0.5)]
        # replicon_id_C ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_B ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB']
        # replicon_id_E ['hit_gspd']
        # C is alone 3.0
        # B and E are compatible 2.5
        # So the best Solution expected is C
        best_sol, score = find_best_solutions(sorted_syst)
        expected_sol = [[self.systems[k] for k in 'C']]
        best_sol = {frozenset(sol) for sol in best_sol}
        expected_sol = {frozenset(sol) for sol in expected_sol}
        self.assertEqual(score, 3.0)
        self.assertSetEqual(best_sol, expected_sol)

        # systems = [('replicon_id_C', 3.0), ('replicon_id_B', 2.0), ('replicon_id_A', 1.5),
        #            ('replicon_id_D', 1.5), ('replicon_id_E', 0.5)]
        # replicon_id_C ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_B ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB']
        # replicon_id_A ['hit_sctj', 'hit_sctn', 'hit_gspd', 'hit_sctj', 'hit_sctn']
        # replicon_id_D ['hit_abc', 'hit_sctn']
        # replicon_id_E ['hit_gspd']
        # C and D are compatible 4.5
        # B and A are compatible 3.5
        # B and E are compatible 2.5
        # D and E are compatible 2.0
        systems = [self.systems[k] for k in 'ABCDE']
        sorted_syst = sorted(systems, key=lambda s: (- s.score, s.id))
        best_sol, score = find_best_solutions(sorted_syst)
        expected_sol = [[self.systems[k] for k in 'CD']]
        best_sol = {frozenset(sol) for sol in best_sol}
        expected_sol = {frozenset(sol) for sol in expected_sol}
        self.assertEqual(score, 4.5)
        self.assertSetEqual(best_sol, expected_sol)

        # systems = [('replicon_id_C', 3.0), ('replicon_id_B', 2.0), ('replicon_id_A', 1.5),
        #            ('replicon_id_D', 1.5), ('replicon_id_E', 0.5), ('replicon_id_F', 1.0)]
        # replicon_id_C ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_B ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB']
        # replicon_id_A ['hit_sctj', 'hit_sctn', 'hit_gspd', 'hit_sctj', 'hit_sctn']
        # replicon_id_D ['hit_abc', 'hit_sctn']
        # replicon_id_E ['hit_gspd']
        # replicon_id_F ['hit_abc']
        # C and D are compatible 4.5
        # C and F are compatible 4.0
        # B and A and F are compatible 4.5
        # B and D and E are compatible 4.0
        # B and E and F are compatible 3.5
        # So the best Solution expected are C D / B A F
        systems = [self.systems[k] for k in 'ABCDEF']
        sorted_syst = sorted(systems, key=lambda s: (- s.score, s.id))
        best_sol, score = find_best_solutions(sorted_syst)
        expected_sol = [[self.systems[k] for k in 'BAF'],  # 3 + 5 + 1 = 9 hits
                        [self.systems[k] for k in 'CD']]   # 4 + 2 = 7 hits
        best_sol_unordered = {frozenset(sol) for sol in best_sol}
        expected_sol_unordered = {frozenset(sol) for sol in expected_sol}
        self.assertEqual(score, 4.5)
        # test if the composition is right
        self.assertSetEqual(best_sol_unordered, expected_sol_unordered)
        # test if solution order is right
        best_sol_ordered = [frozenset(sol) for sol in best_sol]
        expected_sol_ordered = [frozenset(sol) for sol in expected_sol]
        self.assertListEqual(best_sol_ordered, expected_sol_ordered)

        systems = [self.systems[k] for k in 'ABCDGH']
        sorted_syst = sorted(systems, key=lambda s: (- s.score, s.id))
        # sorted_syst = [('replicon_id_C', 3.0), ('replicon_id_B', 2.0), ('replicon_id_A', 1.5), ('replicon_id_D', 1.5)
        #                ('replicon_id_G', 3.0), ('replicon_id_H', 1.5)]
        # replicon_id_A ['hit_sctj', 'hit_sctn', 'hit_gspd', 'hit_sctj', 'hit_sctn']
        # replicon_id_B ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB']
        # replicon_id_C ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_D ['hit_abc', 'hit_sctn']
        # replicon_id_G ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_H ['hit_abc', 'hit_sctn']
        # C and D are compatible 4.5   wholeness = 0.8  + 1.0 = 1.8
        # C and H are compatible 4.5               0.8  + 1.0 = 1.8
        # G and D are compatible 4.5               0.8  + 1.0 = 1.8
        # G and H are compatible 4.5               0.8  + 1.0 = 1.8
        # B and A are compatible 3.5               0.75 + 1.0 = 1.75
        # So the best Solution expected are C D / C H / G D / G H with score 4.5

        best_sol, score = find_best_solutions(sorted_syst)
        expected_sol = [[self.systems[k] for k in 'CD'],  # 4 + 2 hits
                        [self.systems[k] for k in 'CH'],  # 4 + 2
                        [self.systems[k] for k in 'GD'],  # 4 + 2
                        [self.systems[k] for k in 'GH']]  # 4 + 2
        best_sol = {frozenset(sol) for sol in best_sol}
        expected_sol = {frozenset(sol) for sol in expected_sol}
        self.assertEqual(score, 4.5)
        self.assertSetEqual(best_sol, expected_sol)

        systems = [self.systems[k] for k in 'HJKI']
        sorted_syst = sorted(systems, key=lambda s: (- s.score, s.id))
        best_sol, score = find_best_solutions(sorted_syst)

        # check if solution is ordered by woleness average (3rd criterion)
        # first criterion nb of hits
        # second citerion nb of systems
        # replicon_id_H ['hit_abc', 'hit_sctn']
        # replicon_id_I ['hit_flgB', 'hit_tadZ']
        # replicon_id_J ['hit_abc', 'hit_tadZ']
        # replicon_id_K ['hit_flgB', 'hit_sctn']
        #                                                             score  Nb hits  nb sys wholeness
        expected_sol = [[self.systems[k] for k in 'HI'],  # 1.5 + 1.5 = 3.0    4        2       1.0
                        [self.systems[k] for k in 'JK']]  # 1.5 + 1.5 = 3.0    4        2       0.5
        best_sol = {frozenset(sol) for sol in best_sol}
        expected_sol = {frozenset(sol) for sol in expected_sol}
        self.assertEqual(score, 3.0)
        self.assertSetEqual(best_sol, expected_sol)
