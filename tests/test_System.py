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
from macsypy.system import System, HitSystemTracker, ClusterSystemTracker
from macsypy.error import MacsypyError
from tests import MacsyTest


class SystemTest(MacsyTest):

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.verbosity = 3
        self.cfg = Config(MacsyDefaults(), args)

        self.model_name = 'foo'
        self.model_location = ModelLocation(path=os.path.join(args.models_dir, self.model_name))
        self.profile_factory = ProfileFactory(self.cfg)
        self.hit_weights = HitWeight(**self.cfg.hit_weights())


    def test_init(self):
        model = Model("foo/T2SS", 10)
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
        system_1 = System(model,
                          [Cluster([v_hit_1, v_hit_2], model, self.hit_weights)],
                          self.cfg.redundancy_penalty())
        self.assertTrue(system_1.id.startswith('replicon_id_T2SS_'))

        system_2 = System(model,
                          [Cluster([v_hit_1, v_hit_2], model, self.hit_weights)],
                          self.cfg.redundancy_penalty())
        self.assertEqual(int(system_2.id.split('_')[-1]), int(system_1.id.split('_')[-1]) + 1)

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
        system_1 = System(model, [Cluster([v_hit_1, v_hit_2], model, self.hit_weights),
                                  Cluster([v_hit_3], model, self.hit_weights)],
                          self.cfg.redundancy_penalty())

        self.assertEqual(system_1.hits, [v_hit_1, v_hit_2, v_hit_3])


    def test_position(self):
        model = Model("foo/mod1", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model, loner=True)
        model.add_accessory_gene(gene_sctn)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 10, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 20, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        system_1 = System(model, [Cluster([v_hit_1, v_hit_2], model, self.hit_weights),
                                  Cluster([v_hit_3], model, self.hit_weights)],
                          self.cfg.redundancy_penalty())
        # loner are not to take in account to compute position if system contains none loner hit
        self.assertEqual(system_1.position, (10, 20))

        model = Model("foo/mod2", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model, loner=True)
        model.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model, loner=True)
        model.add_accessory_gene(gene_sctj)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model, loner=True)
        model.add_accessory_gene(gene_sctn)
        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 10, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 20, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        system_1 = System(model, [Cluster([v_hit_1, v_hit_2], model, self.hit_weights),
                                  Cluster([v_hit_3], model, self.hit_weights)],
                          self.cfg.redundancy_penalty())
        # loner are not to take in account to compute position if system contains none loner hit
        self.assertEqual(system_1.position, (1, 20))


    def test_multi_loci(self):
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
        c1 = Cluster([v_hit_1, v_hit_2], model, self.hit_weights)
        c2 = Cluster([v_hit_1, v_hit_3], model, self.hit_weights)
        sys_single_locus = System(model, [c1])
        self.assertFalse(sys_single_locus.multi_loci)
        sys_multi_loci = System(model, [c1, c2])
        self.assertTrue(sys_multi_loci.multi_loci)
        c1 = Cluster([v_hit_1, v_hit_2], model, self.hit_weights)
        c3 = Cluster([v_hit_3], model, self.hit_weights)
        sys_single_locus_plus_loner = System(model, [c1, c3], self.cfg.redundancy_penalty())
        self.assertFalse(sys_single_locus_plus_loner.multi_loci)

    def test_loci(self):
        model = Model("foo/T2SS", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model, loner=True)
        model.add_accessory_gene(gene_sctn)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        c1 = Cluster([v_hit_1, v_hit_2], model, self.hit_weights)
        c2 = Cluster([v_hit_1, v_hit_3], model, self.hit_weights)
        sys_single_locus = System(model, [c1], self.cfg.redundancy_penalty())
        self.assertEqual(sys_single_locus.loci, 1)
        sys_multi_loci = System(model, [c1, c2], self.cfg.redundancy_penalty())
        self.assertEqual(sys_multi_loci.loci, 2)
        c1 = Cluster([v_hit_1, v_hit_2], model, self.hit_weights)
        c3 = Cluster([v_hit_3], model, self.hit_weights)
        sys_single_locus_plus_loner = System(model, [c1, c3], self.cfg.redundancy_penalty())
        self.assertEqual(sys_single_locus_plus_loner.loci, 1)

    def test_wholeness(self):
        model_1 = Model("foo/T2SS", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model_1)
        model_1.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model_1)
        model_1.add_accessory_gene(gene_sctj)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model_1, loner=True)
        model_1.add_accessory_gene(gene_sctn)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        m1_v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        m1_v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        m1_v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        c1 = Cluster([m1_v_hit_1, m1_v_hit_2], model_1, self.hit_weights)
        c2 = Cluster([m1_v_hit_1, m1_v_hit_3], model_1, self.hit_weights)
        s = System(model_1, [c1])
        self.assertEqual(s.wholeness, 2 / 3)
        s = System(model_1, [c1, c2])
        self.assertEqual(s.wholeness, 3 / 3)

        model_2 = Model("foo/T2SS", 10, max_nb_genes=2)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        m2_gene_gspd = ModelGene(c_gene_gspd, model_2)
        model_2.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        m2_gene_sctj = ModelGene(c_gene_sctj, model_2)
        model_2.add_accessory_gene(gene_sctj)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        m2_gene_sctn = ModelGene(c_gene_sctn, model_2, loner=True)
        model_2.add_accessory_gene(gene_sctn)

        m2_v_hit_1 = ValidHit(hit_1, m2_gene_gspd, GeneStatus.MANDATORY)
        m2_v_hit_2 = ValidHit(hit_2, m2_gene_sctj, GeneStatus.ACCESSORY)
        c3 = Cluster([m2_v_hit_1, m2_v_hit_2], model_2, self.hit_weights)
        s = System(model_2, [c3])
        self.assertEqual(s.wholeness, 1)

    def test_occurrence(self):
        model = Model("foo/T2SS", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model, loner=True)
        model.add_accessory_gene(gene_sctn)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        c = Cluster([v_hit_1, v_hit_2, v_hit_3], model, self.hit_weights)
        s = System(model, [c])
        self.assertEqual(s.occurrence(), 1)
        c1 = Cluster([v_hit_1, v_hit_2, v_hit_3], model, self.hit_weights)
        c2 = Cluster([v_hit_2, v_hit_3], model, self.hit_weights)
        s = System(model, [c1, c2], self.cfg.redundancy_penalty())
        # The estimation of occurrence number is based on mandatory only
        self.assertEqual(s.occurrence(), 1)
        c1 = Cluster([v_hit_1, v_hit_2, v_hit_3], model, self.hit_weights)
        c2 = Cluster([v_hit_1, v_hit_3], model, self.hit_weights)
        s = System(model, [c1, c2], self.cfg.redundancy_penalty())
        self.assertEqual(s.occurrence(), 2)

        ##########################
        # with multi system gene #
        ##########################
        # the multi_system genes should not count in the occurence
        model = Model("foo/T2SS", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_tadz = CoreGene(self.model_location, "tadZ", self.profile_factory)
        gene_tadz = ModelGene(c_gene_tadz, model)
        model.add_mandatory_gene(gene_tadz)
        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_mandatory_gene(gene_sctj)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model, multi_system=True)
        model.add_mandatory_gene(gene_sctn)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_tadz, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_tadz, GeneStatus.MANDATORY)
        hit_3_1 = Hit(c_gene_sctj, "hit_3_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3_1 = ValidHit(hit_3_1, gene_sctj, GeneStatus.MANDATORY)
        hit_3_2 = Hit(c_gene_sctj, "hit_3_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3_2 = ValidHit(hit_3_2, gene_sctj, GeneStatus.MANDATORY)
        hit_4_1 = Hit(c_gene_sctn, "hit_4_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_4_1 = ValidHit(hit_4_1, gene_sctn, GeneStatus.MANDATORY)
        hit_4_2 = Hit(c_gene_sctn, "hit_4_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_4_2 = ValidHit(hit_4_2, gene_sctn, GeneStatus.MANDATORY)
        hit_4_3 = Hit(c_gene_sctn, "hit_4_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_4_3 = ValidHit(hit_4_3, gene_sctn, GeneStatus.MANDATORY)

        c = Cluster([v_hit_1, v_hit_2, v_hit_3_1, v_hit_3_2, v_hit_4_1, v_hit_4_2, v_hit_4_3], model, self.hit_weights)
        s = System(model, [c], self.cfg.redundancy_penalty())
        self.assertEqual(s.occurrence(), 1)


    def test_score(self):
        model = Model("foo/T2SS", 10)
        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        c_gene_tadZ = CoreGene(self.model_location, "tadZ", self.profile_factory)
        gene_tadZ = ModelGene(c_gene_tadZ, model)
        model.add_mandatory_gene(gene_tadZ)

        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        c_gene_sctJ_FLG = CoreGene(self.model_location, "sctJ_FLG", self.profile_factory)
        sctJ_FLG = Exchangeable(c_gene_sctJ_FLG, gene_sctj)
        gene_sctj.add_exchangeable(sctJ_FLG)
        model.add_accessory_gene(gene_sctj)

        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model, loner=True)
        c_gene_sctn_FLG = CoreGene(self.model_location, "sctN_FLG", self.profile_factory)
        sctn_FLG = Exchangeable(c_gene_sctn_FLG, gene_sctj)
        gene_sctn.add_exchangeable(sctn_FLG)
        model.add_accessory_gene(gene_sctn)

        c_gene_abc = CoreGene(self.model_location, "abc", self.profile_factory)
        gene_abc = ModelGene(c_gene_abc, model)
        model.add_neutral_gene(gene_abc)

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
        v_h_sctj_an = ValidHit(h_sctj_an, sctJ_FLG, GeneStatus.ACCESSORY)

        h_sctn = Hit(c_gene_sctn, "sctn", 40, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_h_sctn = ValidHit(h_sctn, gene_sctn, GeneStatus.ACCESSORY)
        h_sctn_hom = Hit(c_gene_sctn_FLG, "h_scth_hom", 30, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_h_sctn_hom = ValidHit(h_sctn_hom, sctn_FLG, GeneStatus.ACCESSORY)

        h_abc = Hit(c_gene_abc, "h_abc", 50, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_h_abc = ValidHit(h_abc, gene_abc, GeneStatus.NEUTRAL)

        h_flie_1 = Hit(c_gene_flie, "h_flie_1", 150, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_h_flie_1 = ValidHit(h_flie_1, gene_flie, GeneStatus.MANDATORY)

        h_flie_2 = Hit(c_gene_flie, "h_flie_2", 300, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_h_flie_2 = ValidHit(h_flie_2, gene_flie, GeneStatus.MANDATORY)

        h_flie_3 = Hit(c_gene_flie, "h_flie_3", 25, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_h_flie_3 = ValidHit(h_flie_3, gene_flie, GeneStatus.MANDATORY)

        # system with
        # 1 cluster
        # 2 mandatory, 2 accessory no analog/homolog
        s = System(model,
                   [Cluster([v_h_gspd, v_h_tadz, v_h_sctj, v_h_sctn], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 3)

        # system with
        # 2 clusters
        # 2 mandatory, 2 accessory no analog/homolog no duplicates
        s = System(model, [Cluster([v_h_gspd, v_h_tadz], model, self.hit_weights),
                           Cluster([v_h_sctj, v_h_sctn], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 3)

        # system with
        # 1 cluster
        # 1 mandatory + 1 mandatory duplicated 1 time
        # 1 accessory + 1 accessory duplicated 1 times
        # no analog/homolog
        s = System(model,
                   [Cluster([v_h_gspd, v_h_tadz, v_h_sctj, v_h_sctn, v_h_gspd, v_h_sctn], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 3)

        # system with 2 clusters
        # 1rst cluster 1 mandatory + 1 accessory
        #                  1       +   0.5
        # 2nd cluster 1 mandatory + 1 accessory already in first cluster
        #                  1      +    0.5
        # 3 - (2 - 1) * 1.5 = 1.5
        s = System(model,
                   [Cluster([v_h_gspd, v_h_sctj], model, self.hit_weights),
                    Cluster([v_h_tadz, v_h_sctj], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 1.5)

        # system with 1 cluster
        # 2 mandatory
        # 1 accessory + 1 accessory exchangeable
        s = System(model,
                   [Cluster([v_h_gspd, v_h_tadz, v_h_sctj_an, v_h_sctn], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 2.9)

        # system with 2 cluster
        # 1 mandatory + 1 accessory
        #    1        +      0.5
        # 1 mandatory + 1 accessory exchangeable same role as cluster_1 accessory
        #    1        +      0.4
        # system penalty due to 2 genes with same role in 2 clusters: -1.5
        #    2.9 - 1.5 = 1.8
        s = System(model,
                   [Cluster([v_h_gspd, v_h_sctj], model, self.hit_weights),
                    Cluster([v_h_tadz, v_h_sctj_an], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 1.4)

        # system with 2 cluster
        # 1 mandatory + 1 accessory + 1 neutral + same gene as accessory
        #    1        +      0.5        0             0
        # 2 clusters of loner multi systems
        #    + 0.7
        # system penalty: 0
        #    2.2
        s = System(model,
                   [Cluster([v_h_gspd, v_h_sctj, v_h_abc, v_h_sctj], model, self.hit_weights),
                    Cluster([v_h_flie_1], model, self.hit_weights),
                    Cluster([v_h_flie_2], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 2.2)

        # system with 2 cluster
        # 1 mandatory + 1 accessory + 1 neutral + same gene as accessory + mandatory loner
        #    1        +      0.5        0             0                  +   1
        # 2 clusters of loner multi systems one is already in first cluster
        #    + 0
        # system penalty: 0
        #    2.5
        s = System(model,
                   [Cluster([v_h_gspd, v_h_sctj, v_h_abc, v_h_sctj, v_h_flie_1], model, self.hit_weights),
                    Cluster([v_h_flie_2], model, self.hit_weights),
                    Cluster([v_h_flie_3], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 2.5)


    def test_is_compatible(self):
        model_A = Model("foo/A", 10)
        model_B = Model("foo/B", 10)
        model_C = Model("foo/C", 10)
        model_D = Model("foo/D", 10)

        c_gene_sctn_flg = CoreGene(self.model_location, "sctN_FLG", self.profile_factory)
        gene_sctn_flg = ModelGene(c_gene_sctn_flg, model_B)
        c_gene_sctj_flg = CoreGene(self.model_location, "sctJ_FLG", self.profile_factory)
        gene_sctj_flg = ModelGene(c_gene_sctj_flg, model_B)
        c_gene_flgB = CoreGene(self.model_location, "flgB", self.profile_factory)
        gene_flgB = ModelGene(c_gene_flgB, model_B)
        c_gene_tadZ = CoreGene(self.model_location, "tadZ", self.profile_factory)
        gene_tadZ = ModelGene(c_gene_tadZ, model_B)

        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model_A, multi_system=True)
        gene_sctn_hom = Exchangeable(c_gene_sctn_flg, gene_sctn)
        gene_sctn.add_exchangeable(gene_sctn_hom)

        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model_A)
        gene_sctj_an = Exchangeable(c_gene_sctj_flg, gene_sctj)
        gene_sctj.add_exchangeable(gene_sctj_an)

        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model_A)
        gene_gspd_an = Exchangeable(c_gene_flgB, gene_gspd)
        gene_gspd.add_exchangeable(gene_gspd_an)

        c_gene_abc = CoreGene(self.model_location, "abc", self.profile_factory)
        gene_abc = ModelGene(c_gene_abc, model_A)
        gene_abc_ho = Exchangeable(c_gene_tadZ, gene_abc)
        gene_abc.add_exchangeable(gene_abc_ho)

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

        h_sctj = Hit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = Hit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd = Hit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_sctj_flg = Hit(c_gene_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_flgB = Hit(c_gene_flgB, "hit_flgB", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_tadZ = Hit(c_gene_tadZ, "hit_tadZ", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_abc = Hit(c_gene_abc, "hit_abc", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        model_A._min_mandatory_genes_required = 2
        model_A._min_genes_required = 2
        c1 = Cluster([ValidHit(h_sctj, gene_sctj, GeneStatus.MANDATORY),
                      ValidHit(h_sctn, gene_sctn, GeneStatus.MANDATORY),
                      ValidHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
                      ],
                     model_A, self.hit_weights)

        c2 = Cluster([ValidHit(h_sctj, gene_sctj, GeneStatus.MANDATORY),
                      ValidHit(h_sctn, gene_sctn, GeneStatus.MANDATORY)],
                     model_A, self.hit_weights)

        model_B._min_mandatory_genes_required = 1
        model_B._min_genes_required = 2
        c3 = Cluster([ValidHit(h_sctj_flg, gene_sctj_flg, GeneStatus.MANDATORY),
                      ValidHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY),
                      ValidHit(h_flgB, gene_flgB, GeneStatus.ACCESSORY)],
                     model_B, self.hit_weights)
        model_C._min_mandatory_genes_required = 1
        model_C._min_genes_required = 2
        c4 = Cluster([ValidHit(h_sctj_flg, gene_sctj_flg, GeneStatus.MANDATORY),
                      ValidHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY),
                      ValidHit(h_flgB, gene_flgB, GeneStatus.MANDATORY),
                      ValidHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)],
                     model_C, self.hit_weights)
        model_D._min_mandatory_genes_required = 1
        model_D._min_genes_required = 1
        c5 = Cluster([ValidHit(h_abc, gene_abc, GeneStatus.MANDATORY),
                      ValidHit(h_sctn, gene_sctn, GeneStatus.ACCESSORY)],
                     model_D, self.hit_weights)

        sys_A = System(model_A, [c1, c2], self.cfg.redundancy_penalty())
        # we need to tweek the replicon_id to have stable results
        # whatever the number of tests ran
        # or the tests order
        sys_A.id = "replicon_id_A"
        sys_B = System(model_B, [c3], self.cfg.redundancy_penalty())
        sys_B.id = "replicon_id_B"
        sys_C = System(model_C, [c4], self.cfg.redundancy_penalty())
        sys_C.id = "replicon_id_C"
        sys_D = System(model_D, [c5], self.cfg.redundancy_penalty())
        sys_D.id = "replicon_id_D"

        self.assertTrue(sys_A.is_compatible(sys_B))
        self.assertFalse(sys_A.is_compatible(sys_C))  # share h_gspd
        self.assertTrue(sys_A.is_compatible(sys_D))   # share h_sctn but sctn is defined as multi_system im model A


    def test_HitSystemTracker(self):
        model_1 = Model("foo/T2SS", 10)
        model_2 = Model("foo/T3SS", 10)

        c_gene_sctn_flg = CoreGene(self.model_location, "sctN_FLG", self.profile_factory)
        gene_sctn_flg = ModelGene(c_gene_sctn_flg, model_2)
        c_gene_sctj_flg = CoreGene(self.model_location, "sctJ_FLG", self.profile_factory)
        gene_sctj_flg = ModelGene(c_gene_sctj_flg, model_2)
        c_gene_flgB = CoreGene(self.model_location, "flgB", self.profile_factory)
        gene_flgB = ModelGene(c_gene_flgB, model_2)
        c_gene_tadZ = CoreGene(self.model_location, "tadZ", self.profile_factory)
        gene_tadZ = ModelGene(c_gene_tadZ, model_2)

        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model_1)
        gene_sctn_hom = Exchangeable(c_gene_sctn_flg, gene_sctn)
        gene_sctn.add_exchangeable(gene_sctn_hom)

        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model_1)
        gene_sctj_an = Exchangeable(c_gene_sctj_flg, gene_sctj)
        gene_sctj.add_exchangeable(gene_sctj_an)

        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model_1)
        gene_gspd_an = Exchangeable(c_gene_flgB, gene_gspd)
        gene_gspd.add_exchangeable(gene_gspd_an)

        c_gene_abc = CoreGene(self.model_location, "abc", self.profile_factory)
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
                     model_1,
                     self.hit_weights)

        c2 = Cluster([ValidHit(h_sctj, gene_sctj, GeneStatus.MANDATORY),
                      ValidHit(h_sctn, gene_sctn, GeneStatus.MANDATORY),
                      ValidHit(h_flgB, gene_gspd, GeneStatus.ACCESSORY)],
                     model_1,
                     self.hit_weights)

        model_2._min_mandatory_genes_required = 1
        model_2._min_genes_required = 2
        c3 = Cluster([ValidHit(h_sctj_flg, gene_sctj_flg, GeneStatus.MANDATORY),
                      ValidHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY),
                      ValidHit(h_flgB, gene_flgB, GeneStatus.ACCESSORY)],
                     model_2,
                     self.hit_weights)
        s1 = System(model_1, [c1], self.cfg.redundancy_penalty())
        s2 = System(model_1, [c1, c2], self.cfg.redundancy_penalty())
        s3 = System(model_2, [c3], self.cfg.redundancy_penalty())

        track_multi_systems_hit = HitSystemTracker([s1, s2, s3])
        self.assertSetEqual({s1, s2}, track_multi_systems_hit[h_sctj])
        self.assertSetEqual({s1, s2}, track_multi_systems_hit[h_gspd])
        self.assertSetEqual({s3}, track_multi_systems_hit[h_tadZ])
        self.assertSetEqual({s2, s3}, track_multi_systems_hit[h_flgB])


    def test_ClusteSystemTracker(self):
        model_1 = Model("foo/mod1", 10)
        model_2 = Model("foo/mod2", 10)

        c_gene_sctn_flg = CoreGene(self.model_location, "sctN_FLG", self.profile_factory)
        gene_sctn_flg = ModelGene(c_gene_sctn_flg, model_2)
        c_gene_sctj_flg = CoreGene(self.model_location, "sctJ_FLG", self.profile_factory)
        gene_sctj_flg = ModelGene(c_gene_sctj_flg, model_2)
        c_gene_flgB = CoreGene(self.model_location, "flgB", self.profile_factory)
        gene_flgB = ModelGene(c_gene_flgB, model_2)
        c_gene_tadZ = CoreGene(self.model_location, "tadZ", self.profile_factory)
        gene_tadZ = ModelGene(c_gene_tadZ, model_2)

        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model_1)
        gene_sctn_hom = Exchangeable(c_gene_sctn, gene_sctn_flg)
        gene_sctn_flg.add_exchangeable(gene_sctn_hom)

        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model_1)
        gene_sctj_an = Exchangeable(c_gene_sctj, gene_sctj_flg)
        gene_sctj_flg.add_exchangeable(gene_sctj_an)

        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model_1)
        gene_gspd_an = Exchangeable(c_gene_gspd, gene_tadZ)
        gene_tadZ.add_exchangeable(gene_gspd_an)

        c_gene_abc = CoreGene(self.model_location, "abc", self.profile_factory)
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
                     model_1, self.hit_weights)

        c2 = Cluster([ValidHit(h_sctj, gene_sctj, GeneStatus.MANDATORY),
                      ValidHit(h_sctn, gene_sctn, GeneStatus.MANDATORY),
                      ValidHit(h_flgB, gene_gspd, GeneStatus.ACCESSORY)],
                     model_1, self.hit_weights)

        model_2._min_mandatory_genes_required = 1
        model_2._min_genes_required = 2
        c3 = Cluster([ValidHit(h_sctj_flg, gene_sctj_flg, GeneStatus.MANDATORY),
                      ValidHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY),
                      ValidHit(h_flgB, gene_flgB, GeneStatus.ACCESSORY)],
                     model_2, self.hit_weights)
        s1 = System(model_1, [c1], self.cfg.redundancy_penalty())
        s2 = System(model_1, [c1, c2], self.cfg.redundancy_penalty())
        s3 = System(model_2, [c3], self.cfg.redundancy_penalty())

        track_multi_systems_cluster = ClusterSystemTracker([s1, s2, s3])
        self.assertSetEqual({s1, s2}, track_multi_systems_cluster[c1])
        self.assertSetEqual({s2}, track_multi_systems_cluster[c2])
        self.assertSetEqual({s3}, track_multi_systems_cluster[c3])


    def test_count(self):
        model = Model("foo/T2SS", 10)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model)
        c_gene_sctn_flg = CoreGene(self.model_location, "sctN_FLG", self.profile_factory)
        gene_sctn_flg = Exchangeable(c_gene_sctn_flg, gene_sctn)
        gene_sctn.add_exchangeable(gene_sctn_flg)

        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        c_gene_sctj_flg = CoreGene(self.model_location, "sctJ_FLG", self.profile_factory)
        gene_sctj_flg = Exchangeable(c_gene_sctj_flg, gene_sctj)
        gene_sctj.add_exchangeable(gene_sctj_flg)

        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        c_gene_flgB = CoreGene(self.model_location, "flgB", self.profile_factory)
        gene_gspd_an = Exchangeable(c_gene_flgB, gene_gspd)
        gene_gspd.add_exchangeable(gene_gspd_an)

        c_gene_abc = CoreGene(self.model_location, "abc", self.profile_factory)
        gene_abc = ModelGene(c_gene_abc, model)
        c_gene_tadZ = CoreGene(self.model_location, "tadZ", self.profile_factory)
        gene_abc_ho = Exchangeable(c_gene_tadZ, gene_abc)
        gene_abc.add_exchangeable(gene_abc_ho)

        c_gene_toto = CoreGene(self.model_location, "toto", self.profile_factory)
        gene_toto = ModelGene(c_gene_toto, model)
        c_gene_totote = CoreGene(self.model_location, "totote", self.profile_factory)
        gene_toto_ho = Exchangeable(c_gene_totote, gene_toto)
        gene_toto.add_exchangeable(gene_toto_ho)

        c_gene_not_in_model = CoreGene(self.model_location, "toto", self.profile_factory)
        gene_not_in_model = ModelGene(c_gene_not_in_model, model)

        model.add_mandatory_gene(gene_sctn)
        model.add_mandatory_gene(gene_sctj)
        model.add_accessory_gene(gene_gspd)
        model.add_neutral_gene(gene_toto)
        model.add_forbidden_gene(gene_abc)

        h_sctj = Hit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctj_flg = Hit(c_gene_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = Hit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn_flg = Hit(c_gene_sctn_flg, "hit_sctn_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd = Hit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd_an = Hit(c_gene_flgB, "hit_gspd_an", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_abc = Hit(c_gene_abc, "hit_abc", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_abc_ho = Hit(c_gene_tadZ, "hit_abc_ho", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_toto = Hit(c_gene_toto, "hit_toto", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_toto_ho = Hit(c_gene_totote, "hit_toto_ho", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_not_in_model = Hit(c_gene_not_in_model, "hit_not_in_model", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        model._min_mandatory_genes_required = 2
        model._min_genes_required = 1

        v_h_sctj = ValidHit(h_sctj, gene_sctj, GeneStatus.MANDATORY)
        v_h_sctn = ValidHit(h_sctn, gene_sctn, GeneStatus.MANDATORY)
        v_h_gspd = ValidHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
        v_h_abc = ValidHit(h_abc, gene_abc, GeneStatus.FORBIDDEN)
        v_h_toto = ValidHit(h_toto, gene_toto, GeneStatus.NEUTRAL)
        c1 = Cluster([v_h_sctj,
                      v_h_sctn,
                      v_h_gspd,
                      v_h_abc,
                      v_h_toto],
                     model, self.hit_weights)
        s1 = System(model, [c1], self.cfg.redundancy_penalty())

        self.assertDictEqual(s1.mandatory_occ, {'sctJ': [v_h_sctj], 'sctN': [v_h_sctn]})
        self.assertDictEqual(s1.accessory_occ, {'gspD': [v_h_gspd]})
        self.assertDictEqual(s1.neutral_occ, {'toto': [v_h_toto]})

        # test with homolog and analog
        v_h_sctj = ValidHit(h_sctj, gene_sctj, GeneStatus.MANDATORY)
        v_h_sctn = ValidHit(h_sctn, gene_sctn, GeneStatus.MANDATORY)
        v_h_gspd = ValidHit(h_gspd_an, gene_gspd, GeneStatus.ACCESSORY)
        v_h_abc = ValidHit(h_abc, gene_abc, GeneStatus.FORBIDDEN)
        v_h_toto = ValidHit(h_toto_ho, gene_toto, GeneStatus.NEUTRAL)
        c1 = Cluster([v_h_sctj,
                      v_h_sctn,
                      v_h_gspd,
                      v_h_abc,
                      v_h_toto],
                     model, self.hit_weights)
        s1 = System(model, [c1], self.cfg.redundancy_penalty())

        self.assertDictEqual(s1.mandatory_occ, {'sctJ': [v_h_sctj], 'sctN': [v_h_sctn]})
        self.assertDictEqual(s1.accessory_occ, {'gspD': [v_h_gspd]})
        self.assertDictEqual(s1.neutral_occ, {'toto': [v_h_toto]})

        v_h_sctj = ValidHit(h_sctj, gene_sctj, GeneStatus.MANDATORY)
        v_h_not_in_model = ValidHit(h_not_in_model, gene_not_in_model, GeneStatus.MANDATORY)
        c1 = Cluster([v_h_sctj, v_h_not_in_model], model, self.hit_weights)
        with self.assertRaises(MacsypyError) as ctx:
            s1 = System(model, [c1], self.cfg.redundancy_penalty())
        self.assertEqual(str(ctx.exception),
                         "gene 'toto' does not belong to 'mandatory' genes in model 'T2SS'")

    def test_multi_loci(self):
        model = Model("foo/T2SS", 10)
        c_gene_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model)
        c_gene_sctn_flg = CoreGene(self.model_location, "sctN_FLG", self.profile_factory)
        gene_sctn_flg = Exchangeable(c_gene_sctn_flg, gene_sctn)
        gene_sctn.add_exchangeable(gene_sctn_flg)

        c_gene_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        c_gene_sctj_flg = CoreGene(self.model_location, "sctJ_FLG", self.profile_factory)
        gene_sctj_flg = Exchangeable(c_gene_sctj_flg, gene_sctj)
        gene_sctj.add_exchangeable(gene_sctj_flg)

        c_gene_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        c_gene_flgB = CoreGene(self.model_location, "flgB", self.profile_factory)
        gene_gspd_an = Exchangeable(c_gene_flgB, gene_gspd)
        gene_gspd.add_exchangeable(gene_gspd_an)

        c_gene_abc = CoreGene(self.model_location, "abc", self.profile_factory)
        gene_abc = ModelGene(c_gene_abc, model)
        c_gene_tadZ = CoreGene(self.model_location, "tadZ", self.profile_factory)
        gene_abc_ho = Exchangeable(c_gene_tadZ, gene_abc)
        gene_abc.add_exchangeable(gene_abc_ho)

        c_gene_toto = CoreGene(self.model_location, "toto", self.profile_factory)
        gene_toto = ModelGene(c_gene_toto, model)
        c_gene_totote = CoreGene(self.model_location, "totote", self.profile_factory)
        gene_toto_ho = Exchangeable(c_gene_totote, gene_toto)
        gene_toto.add_exchangeable(gene_toto_ho)

        model.add_mandatory_gene(gene_sctn)
        model.add_mandatory_gene(gene_sctj)
        model.add_accessory_gene(gene_gspd)
        model.add_neutral_gene(gene_toto)
        model.add_forbidden_gene(gene_abc)

        h_sctj = Hit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctj_flg = Hit(c_gene_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = Hit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn_flg = Hit(c_gene_sctn_flg, "hit_sctn_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd = Hit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd_an = Hit(c_gene_flgB, "hit_gspd_an", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_abc = Hit(c_gene_abc, "hit_abc", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_abc_ho = Hit(c_gene_tadZ, "hit_abc_ho", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_toto = Hit(c_gene_toto, "hit_toto", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_toto_ho = Hit(c_gene_totote, "hit_toto_ho", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        model._min_mandatory_genes_required = 2
        model._min_genes_required = 1

        v_h_sctj = ValidHit(h_sctj, gene_sctj, GeneStatus.MANDATORY)
        v_h_sctn = ValidHit(h_sctn, gene_sctn, GeneStatus.MANDATORY)
        v_h_gspd = ValidHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
        v_h_toto = ValidHit(h_toto, gene_toto, GeneStatus.NEUTRAL)
        c1 = Cluster([v_h_sctj,
                      v_h_sctn,
                      v_h_gspd,
                      v_h_toto],
                     model, self.hit_weights)
        s1 = System(model, [c1], self.cfg.redundancy_penalty())

        self.assertFalse(s1.multi_loci)
        c2 = Cluster([v_h_sctj,
                      v_h_sctn],
                     model, self.hit_weights)
        c3 = Cluster([v_h_gspd,
                      v_h_toto],
                     model, self.hit_weights)
        s2 = System(model, [c2, c3], self.cfg.redundancy_penalty())
        self.assertTrue(s2.multi_loci)
