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
from macsypy.system import System, HitSystemTracker
from macsypy.score import ComposedScore, BestSystemSelector
from macsypy.error import MacsypyError
from tests import MacsyTest


class ComposedScoreTest(MacsyTest):

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
        self.systems = self.build_systems()
        self.hit_system_tracker = HitSystemTracker(self.systems)


    def build_systems(self):
        model_name = 'foo'
        model_location = ModelLocation(path=os.path.join(self.cfg.models_dir(), model_name))
        model_1 = Model("foo/T2SS", 10)
        model_2 = Model("foo/T3SS", 10)
        model_3 = Model("foo/T4SS", 10)

        c_gene_sctn_flg = CoreGene(model_location, "sctN_FLG", self.profile_factory)
        gene_sctn_flg = ModelGene(c_gene_sctn_flg, model_2)
        c_gene_sctj_flg = CoreGene(model_location, "sctJ_FLG", self.profile_factory)
        gene_sctj_flg = ModelGene(c_gene_sctj_flg, model_2)
        c_gene_flgB = CoreGene(model_location, "flgB", self.profile_factory)
        gene_flgB = ModelGene(c_gene_flgB, model_2)
        c_gene_tadZ = CoreGene(model_location, "tadZ", self.profile_factory)
        gene_tadZ = ModelGene(c_gene_tadZ, model_2)

        c_gene_sctn = CoreGene(model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model_1)
        gene_sctn_hom = Exchangeable(c_gene_sctn_flg, gene_sctn)
        gene_sctn.add_exchangeable(gene_sctn_hom)

        c_gene_sctj = CoreGene(model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model_1)
        gene_sctj_an = Exchangeable(c_gene_sctj_flg, gene_sctj)
        gene_sctj.add_exchangeable(gene_sctj_an)

        c_gene_gspd = CoreGene(model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model_1)
        gene_gspd_an = Exchangeable(c_gene_flgB, gene_gspd)
        gene_gspd.add_exchangeable(gene_gspd_an)

        c_gene_abc = CoreGene(model_location, "abc", self.profile_factory)
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
        model_3.add_accessory_gene(gene_flgB)
        model_3.add_accessory_gene(gene_tadZ)
        model_3.add_accessory_gene(gene_gspd)
        
        h_sctj = Hit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = Hit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd = Hit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_sctj_flg = Hit(c_gene_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_flgB = Hit(c_gene_flgB, "hit_flgB", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_tadZ = Hit(c_gene_tadZ, "hit_tadZ", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        model_1._min_mandatory_genes_required = 2
        model_1._min_genes_required = 2
        c1 = Cluster([ValidHit(h_sctj, gene_sctj, GeneStatus.MANDATORY),
                      ValidHit(h_sctn, gene_sctn, GeneStatus.MANDATORY),
                      ValidHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
                      ],
                     model_1)

        c2 = Cluster([ValidHit(h_sctj, gene_sctj, GeneStatus.MANDATORY),
                      ValidHit(h_sctn, gene_sctn, GeneStatus.MANDATORY),
                      ValidHit(h_flgB, gene_gspd, GeneStatus.ACCESSORY)],
                     model_1)

        model_2._min_mandatory_genes_required = 1
        model_2._min_genes_required = 2
        c3 = Cluster([ValidHit(h_sctj_flg, gene_sctj_flg, GeneStatus.MANDATORY),
                      ValidHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY),
                      ValidHit(h_flgB, gene_flgB, GeneStatus.ACCESSORY)],
                     model_2)
        c4 = Cluster([ValidHit(h_sctj_flg, gene_sctj_flg, GeneStatus.MANDATORY),
                      ValidHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY),
                      ValidHit(h_flgB, gene_flgB, GeneStatus.ACCESSORY),
                      ValidHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)],
                     model_3)
        s0 = System(model_1, [c1])
        s1 = System(model_1, [c1, c2])
        s2 = System(model_2, [c3])
        s3 = System(model_3, [c4])
        return s0, s1, s2, s3


    def test_system(self):
        cs = ComposedScore(self.systems[0], self.hit_system_tracker)
        self.assertEqual(cs.system, self.systems[0])

    def test_sys_score(self):
        cs = ComposedScore(self.systems[0], self.hit_system_tracker)
        self.assertEqual(cs.sys_score, self.systems[0].score)


    def test_overlapping_genes(self):
        # h_gspd is also in s3
        cs = ComposedScore(self.systems[0], self.hit_system_tracker)
        self.assertEqual(cs.overlapping_genes, 1)

        # sctJ_FLG, tadZ, flgB are also in s3
        cs = ComposedScore(self.systems[2], self.hit_system_tracker)
        self.assertEqual(cs.overlapping_genes, 3)

        # sctJ_FLG, tadZ, flgB are also in s2 and gspd in s1 and s0
        cs = ComposedScore(self.systems[3], self.hit_system_tracker)
        self.assertEqual(cs.overlapping_genes, 4)


    def test_overlapping_length(self):
        # h_gspd is also in s3
        cs = ComposedScore(self.systems[0], self.hit_system_tracker)
        self.assertEqual(cs.overlapping_length, 1)

        # sctJ_FLG, tadZ, flgB are also in s3 and flgB is also in S1
        cs = ComposedScore(self.systems[2], self.hit_system_tracker)
        self.assertEqual(cs.overlapping_length, 4)

        # sctJ_FLG, tadZ, are also in s2
        # flgB is in s2 and s1
        # and gspd in s1 and s0
        cs = ComposedScore(self.systems[3], self.hit_system_tracker)
        self.assertEqual(cs.overlapping_length, 6)


class BestSystemSelectorTest(MacsyTest):

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
        self.systems = self.build_systems()
        self.hit_system_tracker = HitSystemTracker(self.systems)


    def build_systems(self):
        """
        :return: systems

         --- sctn_flg - sctj_flg - flgB ------ tadZ -- sctN --    s0 m1 T2SS
         --- sctn_flg - sctj_flg - flgB ------ tadZ ----------    s1 m1 T2SS
         --- sctn_flg - sctj_flg - flgB -------------- sctN --    s2 m1 T2SS
         ------------------------------- gspD--tadZ-----------    s3 m2 T3SS
         ------------------------------------abc------ sctN --    s4 m3 T4SS
        """
        model_name = 'foo'
        model_location = ModelLocation(path=os.path.join(self.cfg.models_dir(), model_name))
        model_1 = Model("foo/T2SS", 10)
        model_2 = Model("foo/T3SS", 10)
        model_3 = Model("foo/T4SS", 10)

        c_gene_sctn_flg = CoreGene(model_location, "sctN_FLG", self.profile_factory)
        gene_sctn_flg = ModelGene(c_gene_sctn_flg, model_2)

        c_gene_sctj_flg = CoreGene(model_location, "sctJ_FLG", self.profile_factory)
        gene_sctj_flg = ModelGene(c_gene_sctj_flg, model_2)

        c_gene_flgB = CoreGene(model_location, "flgB", self.profile_factory)
        gene_flgB = ModelGene(c_gene_flgB, model_2)

        c_gene_tadZ = CoreGene(model_location, "tadZ", self.profile_factory)
        gene_tadZ = ModelGene(c_gene_tadZ, model_2)

        c_gene_sctn = CoreGene(model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model_1)

        c_gene_gspd = CoreGene(model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model_1)

        c_gene_abc = CoreGene(model_location, "abc", self.profile_factory)
        gene_abc = ModelGene(c_gene_abc, model_1)

        model_1.add_mandatory_gene(gene_sctn_flg)
        model_1.add_mandatory_gene(gene_sctj_flg)
        model_1.add_accessory_gene(gene_flgB)
        model_1.add_accessory_gene(gene_tadZ)
        model_1.add_accessory_gene(gene_sctn)

        model_2.add_mandatory_gene(gene_gspd)
        model_2.add_accessory_gene(gene_tadZ)
        model_2.add_accessory_gene(gene_sctn)

        model_3.add_mandatory_gene(gene_abc)
        model_3.add_accessory_gene(gene_tadZ)
        model_3.add_accessory_gene(gene_sctn)

        h_sctn_flg = Hit(c_gene_sctn_flg, "hit_sctn_flg", 803, "replicon_id", 10, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctj_flg = Hit(c_gene_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 11, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_flgB = Hit(c_gene_flgB, "hit_flgB", 803, "replicon_id", 12, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_tadZ = Hit(c_gene_tadZ, "hit_tadZ", 803, "replicon_id", 20, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = Hit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 21, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_gspd = Hit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 15, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_abc = Hit(c_gene_abc, "hit_abc", 803, "replicon_id", 17, 1.0, 1.0, 1.0, 1.0, 10, 20)

        model_1._min_mandatory_genes_required = 2
        model_1._min_genes_required = 2
        m1_c1 = Cluster([ValidHit(h_sctn_flg, gene_sctn_flg, GeneStatus.MANDATORY),
                         ValidHit(h_sctj_flg, gene_sctj_flg, GeneStatus.MANDATORY),
                         ValidHit(h_flgB, gene_flgB, GeneStatus.ACCESSORY)
                         ],
                        model_1)

        m1_c2 = Cluster([ValidHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY)],
                        model_1)

        m1_c3 = Cluster([ValidHit(h_sctn, gene_sctn, GeneStatus.ACCESSORY)],
                        model_1)
        m1_c4 = Cluster([ValidHit(h_sctn_flg, gene_sctn_flg, GeneStatus.MANDATORY),
                         ValidHit(h_sctj_flg, gene_sctj_flg, GeneStatus.MANDATORY),
                         ],
                        model_1)

        model_2._min_mandatory_genes_required = 1
        model_2._min_genes_required = 2
        m2_c1 = Cluster([ValidHit(h_gspd, gene_gspd, GeneStatus.MANDATORY)],
                         model_2)
        m2_c2 = Cluster([ValidHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY),
                         ValidHit(h_sctn, gene_sctn, GeneStatus.ACCESSORY)],
                         model_2)

        model_3._min_mandatory_genes_required = 1
        model_3._min_genes_required = 2
        m3_c1 = Cluster([ValidHit(h_abc, gene_abc, GeneStatus.MANDATORY)],
                         model_3)
        m3_c2 = Cluster([
                         ValidHit(h_sctn, gene_sctn, GeneStatus.ACCESSORY)],
                         model_3)

        s0 = System(model_1, [m1_c1, m1_c2, m1_c3])
        s1 = System(model_1, [m1_c1, m1_c2])
        s2 = System(model_1, [m1_c1, m1_c3])
        s3 = System(model_1, [m1_c4, m1_c2, m1_c3])
        s4 = System(model_2, [m2_c1, m2_c2])
        s5 = System(model_3, [m3_c1, m3_c2])
        return s0, s1, s2, s3, s4, s5

    def test_BestSystemSelector(self):
        bs = BestSystemSelector(self.systems[:3], self.hit_system_tracker)
        self.assertTrue(isinstance(bs, BestSystemSelector))
        with self.assertRaises(MacsypyError) as ctx:
            BestSystemSelector(self.systems, self.hit_system_tracker)
        self.assertTrue(str(ctx.exception).startswith(
            "Cannot build Score with system from different models:"))

    def test_best_system(self):
        bs = BestSystemSelector(self.systems[:4], self.hit_system_tracker)
        best_systems = bs.best_system()
        self.assertEqual(len(best_systems), 1)
        self.assertEqual(best_systems[0], self.systems[0])

        bs = BestSystemSelector(self.systems[1:4], self.hit_system_tracker)
        best_systems = bs.best_system()
        self.assertEqual(len(best_systems), 1)
        self.assertEqual(best_systems[0], self.systems[1])

        bs = BestSystemSelector([self.systems[1], self.systems[3]], self.hit_system_tracker)
        best_systems = bs.best_system()
        self.assertEqual(len(best_systems), 1)
        self.assertEqual(best_systems[0], self.systems[1])