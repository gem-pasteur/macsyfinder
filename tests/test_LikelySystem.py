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

from macsypy.hit import Hit, ValidHit
from macsypy.config import Config, MacsyDefaults
from macsypy.gene import CoreGene, ModelGene, Exchangeable, GeneStatus
from macsypy.profile import ProfileFactory
from macsypy.model import Model
from macsypy.registries import ModelLocation
from macsypy.system import LikelySystem, UnlikelySystem, HitSystemTracker, ClusterSystemTracker

from tests import MacsyTest


class LikelySystemTest(MacsyTest):

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        self.cfg = Config(MacsyDefaults(), args)

        self.model_name = 'foo'
        self.model_location = ModelLocation(path=os.path.join(args.models_dir, self.model_name))
        self.profile_factory = ProfileFactory(self.cfg)


    def test_init(self):
        model = Model("foo/model_A", 10)
        # test if id is well incremented
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        ls_1 = LikelySystem(model, [v_hit_1], [v_hit_2], [], [])
        self.assertTrue(ls_1.id.startswith('replicon_id_model_A_'))

        ls_2 = LikelySystem(model, [v_hit_1, v_hit_2], [], [], [])
        # check if the id of the second likelysystem is well increased
        self.assertEqual(int(ls_2.id.split('_')[-1]), int(ls_1.id.split('_')[-1]) + 1)

    def test_hits(self):
        model = Model("foo/T2SS", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model)
        model.add_accessory_gene(gene_sctn)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        ls_1 = LikelySystem(model, [v_hit_1], [v_hit_2, v_hit_3], [], [])

        self.assertListEqual(ls_1.hits, [v_hit_1, v_hit_2, v_hit_3])

    def test_str(self):
        model = Model("foo/T2SS", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model)
        model.add_accessory_gene(gene_sctn)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        ls_1 = LikelySystem(model, [v_hit_1], [v_hit_2, v_hit_3], [], [])
        expected_str = ', '.join([f"({h.id}, {h.gene.name}, {h.position})" for h in (v_hit_1, v_hit_2, v_hit_3)])
        self.assertEqual(str(ls_1), expected_str)

    def test_forbidden_hits(self):
        model = Model("foo/T2SS", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_forbidden_gene(gene_sctj)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model)
        model.add_forbidden_gene(gene_sctn)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.FORBIDDEN)
        hit_3 = Hit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.FORBIDDEN)
        ls_1 = LikelySystem(model, [v_hit_1], [], [], [v_hit_2, v_hit_3])
        self.assertListEqual(ls_1.forbidden_hits, [v_hit_2, v_hit_3])

    def test_allowed_hits(self):
        model = Model("foo/T2SS", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model)
        model.add_accessory_gene(gene_sctn)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        ls_1 = LikelySystem(model, [v_hit_1], [v_hit_2], [], [v_hit_3])
        self.assertListEqual(ls_1.allowed_hits, [v_hit_1, v_hit_2])


class UnlikelySystemTest(MacsyTest):

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        self.cfg = Config(MacsyDefaults(), args)

        self.model_name = 'foo'
        self.model_location = ModelLocation(path=os.path.join(args.models_dir, self.model_name))
        self.profile_factory = ProfileFactory(self.cfg)

    def test_init(self):
        model = Model("foo/model_A", 10)
        # test if id is well incremented
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        uls_1 = UnlikelySystem(model, [v_hit_1], [v_hit_2], [], [], ["reason"])
        self.assertTrue(uls_1.id.startswith('replicon_id_model_A_'))

        uls_2 = UnlikelySystem(model, [v_hit_1], [v_hit_2], [], [], ["reason"])
        self.assertEqual(int(uls_2.id.split('_')[-1]), int(uls_1.id.split('_')[-1]) + 1)

    def test_reason(self):
        model = Model("foo/model_A", 10)
        # test if id is well incremented
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_forbidden_gene(gene_sctj)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.FORBIDDEN)
        reason_2 = ["forbidden gene"]
        uls_2 = UnlikelySystem(model, [v_hit_1], [], [], [v_hit_2], reason_2)
        self.assertEqual(uls_2.reasons, reason_2)


    def test_hits(self):
        model = Model("foo/T2SS", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model)
        model.add_accessory_gene(gene_sctn)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        uls_1 = UnlikelySystem(model, [v_hit_1], [v_hit_2, v_hit_3], [], [], ["reason"])

        self.assertListEqual(uls_1.hits, [v_hit_1, v_hit_2, v_hit_3])


    def test_str(self):
        model = Model("foo/T2SS", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model)
        model.add_accessory_gene(gene_sctn)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 2, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(c_gene_sctn, "hit_3", 803, "replicon_id", 3, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        uls_1 = UnlikelySystem(model, [v_hit_1], [v_hit_2, v_hit_3], [], [], ["reason"])

        expected_str = """(hit_1, gspD, 1), (hit_2, sctJ, 2), (hit_3, sctN, 3): These hits does not probably constitute a system because:
reason"""
        self.assertEqual(str(uls_1), expected_str)

