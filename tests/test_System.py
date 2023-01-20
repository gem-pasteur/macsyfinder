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

from macsypy.hit import CoreHit, ModelHit, Loner, MultiSystem, LonerMultiSystem, HitWeight
from macsypy.config import Config, MacsyDefaults
from macsypy.gene import CoreGene, ModelGene, Exchangeable, GeneStatus
from macsypy.profile import ProfileFactory
from macsypy.model import Model
from macsypy.registries import ModelLocation
from macsypy.cluster import Cluster
from macsypy.system import System, HitSystemTracker
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

        hit_1 = CoreHit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ModelHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = CoreHit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ModelHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
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

        hit_1 = CoreHit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ModelHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = CoreHit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ModelHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = CoreHit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ModelHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
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

        hit_1 = CoreHit(c_gene_gspd, "hit_1", 803, "replicon_id", 10, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ModelHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = CoreHit(c_gene_sctj, "hit_2", 803, "replicon_id", 20, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ModelHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = CoreHit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ModelHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
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
        hit_1 = CoreHit(c_gene_gspd, "hit_1", 803, "replicon_id", 10, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ModelHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = CoreHit(c_gene_sctj, "hit_2", 803, "replicon_id", 20, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ModelHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = CoreHit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ModelHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
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

        hit_1 = CoreHit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ModelHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = CoreHit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ModelHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = CoreHit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ModelHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
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

    def test_loci_nb(self):
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

        hit_1 = CoreHit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        m_hit_1 = ModelHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = CoreHit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        m_hit_2 = ModelHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = CoreHit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        m_hit_3 = ModelHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        c1 = Cluster([m_hit_1, m_hit_2], model, self.hit_weights)
        c2 = Cluster([m_hit_1, m_hit_3], model, self.hit_weights)
        sys_single_locus = System(model, [c1], self.cfg.redundancy_penalty())
        self.assertEqual(sys_single_locus.loci_nb, 1)
        sys_multi_loci = System(model, [c1, c2], self.cfg.redundancy_penalty())
        self.assertEqual(sys_multi_loci.loci_nb, 2)
        c1 = Cluster([m_hit_1, m_hit_2], model, self.hit_weights)
        h_loner = Loner(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        c3 = Cluster([h_loner], model, self.hit_weights)
        sys_single_locus_plus_loner = System(model, [c1, c3], self.cfg.redundancy_penalty())
        self.assertEqual(sys_single_locus_plus_loner.loci_nb, 1)
        # if max_gene_requird == 1 we authorize cluster with one gene
        c4 = Cluster([m_hit_1], model, self.hit_weights)
        sys_single_locus_of_one_hit_not_loner = System(model, [c4], self.cfg.redundancy_penalty())
        self.assertEqual(sys_single_locus_of_one_hit_not_loner.loci_nb, 1)

    def test_loci_num(self):
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

        hit_1 = CoreHit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ModelHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = CoreHit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ModelHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = CoreHit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ModelHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        c1 = Cluster([v_hit_1, v_hit_2], model, self.hit_weights)
        c2 = Cluster([v_hit_1, v_hit_3], model, self.hit_weights)
        sys_single_locus = System(model, [c1], self.cfg.redundancy_penalty())
        self.assertListEqual(sys_single_locus.loci_num, [1])
        sys_multi_loci = System(model, [c1, c2], self.cfg.redundancy_penalty())
        self.assertListEqual(sys_multi_loci.loci_num, [1, 2])
        c1 = Cluster([v_hit_1, v_hit_2], model, self.hit_weights)
        h_loner = Loner(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        c3 = Cluster([h_loner], model, self.hit_weights)
        sys_single_locus_plus_loner = System(model, [c1, c3], self.cfg.redundancy_penalty())
        self.assertListEqual(sys_single_locus_plus_loner.loci_num, [1, -1])

        c4 = Cluster([v_hit_1], model, self.hit_weights)
        sys_single_locus_of_one_hit_not_loner = System(model, [c4], self.cfg.redundancy_penalty())
        self.assertEqual(sys_single_locus_of_one_hit_not_loner.loci_num, [1])

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

        hit_1 = CoreHit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        m1_v_hit_1 = ModelHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = CoreHit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        m1_v_hit_2 = ModelHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = CoreHit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        m1_v_hit_3 = ModelHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
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

        m2_v_hit_1 = ModelHit(hit_1, m2_gene_gspd, GeneStatus.MANDATORY)
        m2_v_hit_2 = ModelHit(hit_2, m2_gene_sctj, GeneStatus.ACCESSORY)
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

        hit_1 = CoreHit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ModelHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = CoreHit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ModelHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = CoreHit(c_gene_sctn, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ModelHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
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

        hit_1 = CoreHit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ModelHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = CoreHit(c_gene_tadz, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ModelHit(hit_2, gene_tadz, GeneStatus.MANDATORY)
        hit_3_1 = CoreHit(c_gene_sctj, "hit_3_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3_1 = ModelHit(hit_3_1, gene_sctj, GeneStatus.MANDATORY)
        hit_3_2 = CoreHit(c_gene_sctj, "hit_3_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3_2 = ModelHit(hit_3_2, gene_sctj, GeneStatus.MANDATORY)
        hit_4_1 = CoreHit(c_gene_sctn, "hit_4_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_4_1 = ModelHit(hit_4_1, gene_sctn, GeneStatus.MANDATORY)
        hit_4_2 = CoreHit(c_gene_sctn, "hit_4_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_4_2 = ModelHit(hit_4_2, gene_sctn, GeneStatus.MANDATORY)
        hit_4_3 = CoreHit(c_gene_sctn, "hit_4_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_4_3 = ModelHit(hit_4_3, gene_sctn, GeneStatus.MANDATORY)

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
        gene_flie = ModelGene(c_gene_flie, model, loner=True)
        model.add_mandatory_gene(gene_flie)

        c_gene_flgb = CoreGene(self.model_location, "flgB", self.profile_factory)
        gene_flgb = ModelGene(c_gene_flgb, model, multi_system=True)
        model.add_accessory_gene(gene_flgb)

        c_gene_flgc = CoreGene(self.model_location, "flgC", self.profile_factory)
        gene_flgc = ModelGene(c_gene_flgb, model, loner=True, multi_system=True)
        model.add_accessory_gene(gene_flgc)

        h_gspd = CoreHit(c_gene_gspd, "h_gspd", 10, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        mh_gspd = ModelHit(h_gspd, gene_gspd, GeneStatus.MANDATORY)

        h_tadz = CoreHit(c_gene_tadZ, "h_tadz", 20, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        mh_tadz = ModelHit(h_tadz, gene_tadZ, GeneStatus.MANDATORY)

        h_sctj = CoreHit(c_gene_sctj, "h_sctj", 30, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        mh_sctj = ModelHit(h_sctj, gene_sctj, GeneStatus.ACCESSORY)
        h_sctj_an = CoreHit(c_gene_sctJ_FLG, "h_sctj_an", 30, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        mh_sctj_an = ModelHit(h_sctj_an, sctJ_FLG, GeneStatus.ACCESSORY)

        h_sctn = CoreHit(c_gene_sctn, "sctn", 40, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        mh_sctn = ModelHit(h_sctn, gene_sctn, GeneStatus.ACCESSORY)
        l_sctn = Loner(h_sctn, gene_sctn, GeneStatus.ACCESSORY)
        h_sctn_hom = CoreHit(c_gene_sctn_FLG, "h_scth_hom", 30, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        mh_sctn_hom = ModelHit(h_sctn_hom, sctn_FLG, GeneStatus.ACCESSORY)
        l_sctn_hom = ModelHit(h_sctn_hom, sctn_FLG, GeneStatus.ACCESSORY)

        h_abc = CoreHit(c_gene_abc, "h_abc", 50, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        mh_abc = ModelHit(h_abc, gene_abc, GeneStatus.NEUTRAL)

        h_flie_1 = CoreHit(c_gene_flie, "h_flie_1", 150, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_flie_2 = CoreHit(c_gene_flie, "h_flie_2", 300, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        # mh_flie_1 is a false loner (gene tag as loner hit in cluster)
        mh_flie_1 = ModelHit(h_flie_1, gene_flie, GeneStatus.MANDATORY)
        mh_flie_2 = ModelHit(h_flie_2, gene_flie, GeneStatus.MANDATORY)

        # l_h_flie_1 is a True Loner gene tag as loner, hit  outside cluster
        # mh_flie_2 is a counterpart of l_h_flie_1
        l_h_flie_1 = Loner(h_flie_1, gene_flie, GeneStatus.MANDATORY, counterpart=[mh_flie_2])
        l_h_flie_2 = Loner(h_flie_2, gene_flie, GeneStatus.MANDATORY, counterpart=[mh_flie_1])

        # mh_flie_3 is used as False loner (loner in cluster)
        h_flie_3 = CoreHit(c_gene_flie, "h_flie_3", 25, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        mh_flie_3 = ModelHit(h_flie_3, gene_flie, GeneStatus.MANDATORY)

        # mh_flgb is used as multi_system in cluster
        h_flgb = CoreHit(c_gene_flgb, "h_flgb", 60, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        mh_flgb = MultiSystem(h_flgb, gene_flgb, GeneStatus.ACCESSORY)

        # mh_flgc is used as multi_system in cluster
        h_flgc = CoreHit(c_gene_flgc, "h_flgc", 80, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        mh_flgc = LonerMultiSystem(h_flgc, gene_flgc, GeneStatus.ACCESSORY)

        # system with
        # 1 cluster
        # 2 mandatory, 2 accessory no analog/homolog sctn is a loner included in cluster
        #    2 * 1   +     2 * 0.5 = 3
        s = System(model,
                   [Cluster([mh_gspd, mh_tadz, mh_sctj, mh_sctn], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 3)

        # system with
        # 2 clusters
        # 2 mandatory, 2 accessory no analog/homolog no duplicates
        #     [1 + 0.5] + [1 + 0.5] = 3
        s = System(model, [Cluster([mh_gspd, mh_tadz], model, self.hit_weights),
                           Cluster([mh_sctj, mh_sctn], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 3)

        # system with
        # 1 cluster
        # 2 mandatory gspd + tadZ + 1 mandatory duplicated 1 time gspd
        # 1 accessory sctj + sctn + 1 accessory duplicated 1 times sctn
        # no analog/homolog
        #  1 + 1 + 0 + 0.5 + 0.5 + 0 = 3
        s = System(model,
                   [Cluster([mh_gspd, mh_tadz, mh_sctj, mh_sctn, mh_gspd, mh_sctn], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 3)

        # system with 2 clusters
        # 1rst cluster 1 mandatory gspd + 1 accessory sctj
        #                  1       +   0.5
        # 2nd cluster 1 mandatory tadZ  + 1 accessory sctj already in first cluster
        #                  1      +    0.5
        # [1 + .5] + [1 + .5]
        # 3 - (2 - 1) * 1.5 = 1.5
        s = System(model,
                   [Cluster([mh_gspd, mh_sctj], model, self.hit_weights),
                    Cluster([mh_tadz, mh_sctj], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 1.5)

        # system with 1 cluster
        # 2 mandatory gspd + tadZ
        # 1 accessory sctn + 1 accessory exchangeable sctj_an
        #  1 + 1 + .5 + (.5 * .8) = 2.9
        s = System(model,
                   [Cluster([mh_gspd, mh_tadz, mh_sctj_an, mh_sctn], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 2.9)

        # system with 2 clusters
        # 1 mandatory + 1 accessory
        #    1        +      0.5
        # 1 mandatory + 1 accessory exchangeable same role as cluster_1 accessory
        #    1        +      0.4
        # system penalty due to 2 genes with same role in 2 clusters: -1.5
        #    2.9 - 1.5 = 1.8
        s = System(model,
                   [Cluster([mh_gspd, mh_sctj], model, self.hit_weights),
                    Cluster([mh_tadz, mh_sctj_an], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 1.4)

        # system with 2 cluster
        # 1 mandatory + 1 accessory + 1 neutral + same gene as accessory
        #    gspd            sctj        abc           sctj
        #    1        +      0.5          0             0
        # 1 mandatory true loner (loner alone)
        #       flie_1
        #    + (1 * 0.7)
        # system penalty: 0
        #    2.2
        s = System(model,
                   [Cluster([mh_gspd, mh_sctj, mh_abc, mh_sctj], model, self.hit_weights),
                    Cluster([l_h_flie_1], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 2.2)

        # system with 2 cluster
        # 1 mandatory + 1 accessory + 1 neutral + same gene as accessory + mandatory "faux loner"
        #    1        +      0.5        0             0                  +   1
        # 1 cluster of true loner already in first cluster
        #    + 0
        # system penalty: 0
        #    2.5
        s = System(model,
                   [Cluster([mh_gspd, mh_sctj, mh_abc, mh_sctj, mh_flie_1], model, self.hit_weights),
                    Cluster([l_h_flie_2], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 2.5)

        # system with 2 clusters
        # 1 mandatory gspd + 1 accessory sctj
        #        1                 0.5
        # 1 mandatory tadz + 1 accessory multisystems (in this cluster)
        #        1                 0.5
        # system penalty: 0
        #      3.0
        s = System(model,
                   [Cluster([mh_gspd, mh_sctj], model, self.hit_weights),
                    Cluster([mh_tadz, mh_flgb], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 3.0)

        # system with 1 cluster + 1 multi system out of cluster
        # 1 mandatory gspd + 1 accessory sctj
        #        1                 0.5
        # 1 accessory multisystems flgb (out of cluster)
        #      0.5   *   0.7  = 0.35
        # system penalty: 0
        #      1.85
        s = System(model,
                   [Cluster([mh_gspd, mh_sctj], model, self.hit_weights),
                    Cluster([mh_flgb], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 1.85)

        # system with 1 cluster + 1 Loner multi system (out of cluster)
        # 1 mandatory gspd + 1 accessory sctj
        #        1                 0.5
        # 1 accessory multisystems flgc (out of cluster)
        #      0.5   *   0.7  = 0.35
        # system penalty: 0
        #      1.85
        s = System(model,
                   [Cluster([mh_gspd, mh_sctj], model, self.hit_weights),
                    Cluster([mh_flgc], model, self.hit_weights)],
                   self.cfg.redundancy_penalty())
        self.assertEqual(s.score, 1.85)

    def test_is_compatible(self):
        model_A = Model("foo/A", 10)
        model_B = Model("foo/B", 10)
        model_C = Model("foo/C", 10)
        model_D = Model("foo/D", 10)
        model_E = Model("foo/E", 10)
        model_F = Model("foo/F", 10)

        cg_sctn_flg = CoreGene(self.model_location, "sctN_FLG", self.profile_factory)
        cg_sctj_flg = CoreGene(self.model_location, "sctJ_FLG", self.profile_factory)
        cg_flgB = CoreGene(self.model_location, "flgB", self.profile_factory)
        cg_tadZ = CoreGene(self.model_location, "tadZ", self.profile_factory)
        cg_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)
        cg_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        cg_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        cg_abc = CoreGene(self.model_location, "abc", self.profile_factory)

        mg_A_sctn = ModelGene(cg_sctn, model_A, multi_system=True)
        mg_A_sctn_hom = Exchangeable(cg_sctn_flg, mg_A_sctn, multi_system=True)
        mg_A_sctn.add_exchangeable(mg_A_sctn_hom)
        mg_A_sctj = ModelGene(cg_sctj, model_A)
        mg_A_sctj_an = Exchangeable(cg_sctj_flg, mg_A_sctj)
        mg_A_sctj.add_exchangeable(mg_A_sctj_an)
        mg_A_gspd = ModelGene(cg_gspd, model_A)
        mg_A_gspd_an = Exchangeable(cg_flgB, mg_A_gspd)
        mg_A_gspd.add_exchangeable(mg_A_gspd_an)
        mg_A_abc = ModelGene(cg_abc, model_A)
        mg_A_abc_ho = Exchangeable(cg_tadZ, mg_A_abc)
        mg_A_abc.add_exchangeable(mg_A_abc_ho)

        mg_B_sctn_flg = ModelGene(cg_sctn_flg, model_B)
        mg_B_sctj_flg = ModelGene(cg_sctj_flg, model_B)
        mg_B_flgB = ModelGene(cg_flgB, model_B)
        mg_B_tadZ = ModelGene(cg_tadZ, model_B, loner=True)

        mg_C_sctn_flg = ModelGene(cg_sctn_flg, model_C)
        mg_C_sctj_flg = ModelGene(cg_sctj_flg, model_C)
        mg_C_flgB = ModelGene(cg_flgB, model_C)
        mg_C_tadZ = ModelGene(cg_tadZ, model_C)
        mg_C_gspd = ModelGene(cg_gspd, model_C)

        mg_D_abc = ModelGene(cg_abc, model_D)
        mg_D_sctn = ModelGene(cg_sctn, model_D, multi_system=True)

        mg_E_abc = ModelGene(cg_abc, model_E)
        mg_E_sctn = ModelGene(cg_sctn, model_E, multi_model=True)

        mg_F_abc = ModelGene(cg_abc, model_F)
        mg_F_sctn = ModelGene(cg_sctn, model_F, multi_model=True)

        model_A.add_mandatory_gene(mg_A_sctn)  # multi system
        model_A.add_mandatory_gene(mg_A_sctj)
        model_A.add_accessory_gene(mg_A_gspd)
        model_A.add_forbidden_gene(mg_A_abc)

        model_B.add_mandatory_gene(mg_B_sctn_flg)
        model_B.add_mandatory_gene(mg_B_sctj_flg)
        model_B.add_accessory_gene(mg_B_flgB)
        model_B.add_accessory_gene(mg_B_tadZ)  # loner

        model_C.add_mandatory_gene(mg_C_sctn_flg)
        model_C.add_mandatory_gene(mg_C_sctj_flg)
        model_C.add_mandatory_gene(mg_C_flgB)
        model_C.add_accessory_gene(mg_C_tadZ)
        model_C.add_accessory_gene(mg_C_gspd)

        model_D.add_mandatory_gene(mg_D_abc)
        model_D.add_accessory_gene(mg_D_sctn)

        model_E.add_mandatory_gene(mg_E_abc)
        model_E.add_accessory_gene(mg_E_sctn)

        model_F.add_mandatory_gene(mg_F_abc)
        model_F.add_accessory_gene(mg_F_sctn)

        ch_sctj = CoreHit(cg_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        ch_sctn = CoreHit(cg_sctn, "hit_sctn", 803, "replicon_id", 2, 1.0, 1.0, 1.0, 1.0, 10, 20)
        ch_gspd = CoreHit(cg_gspd, "hit_gspd", 803, "replicon_id", 3, 1.0, 1.0, 1.0, 1.0, 10, 20)
        ch_sctj_flg = CoreHit(cg_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 4, 1.0, 1.0, 1.0, 1.0, 10, 20)
        ch_flgB = CoreHit(cg_flgB, "hit_flgB", 803, "replicon_id", 5, 1.0, 1.0, 1.0, 1.0, 10, 20)
        ch_tadZ = CoreHit(cg_tadZ, "hit_tadZ", 803, "replicon_id", 6, 1.0, 1.0, 1.0, 1.0, 10, 20)
        ch_abc = CoreHit(cg_abc, "hit_abc", 803, "replicon_id", 7, 1.0, 1.0, 1.0, 1.0, 10, 20)

        ch_sctj_2 = CoreHit(cg_sctj, "hit_sctj_2", 803, "replicon_id", 20, 1.0, 1.0, 1.0, 1.0, 10, 20)
        ch_sctn_2 = CoreHit(cg_sctn, "hit_sctn_2", 803, "replicon_id", 21, 1.0, 1.0, 1.0, 1.0, 10, 20)
        ch_gspd_2 = CoreHit(cg_gspd, "hit_gspd_2", 803, "replicon_id", 22, 1.0, 1.0, 1.0, 1.0, 10, 20)
        ch_abc_2 = CoreHit(cg_abc, "hit_abc_2", 803, "replicon_id", 23, 1.0, 1.0, 1.0, 1.0, 10, 20)

        model_A._min_mandatory_genes_required = 2
        model_A._min_genes_required = 2
        c1_A = Cluster([ModelHit(ch_sctj, mg_A_sctj, GeneStatus.MANDATORY),
                        MultiSystem(ch_sctn, mg_A_sctn, GeneStatus.MANDATORY),
                        ModelHit(ch_gspd, mg_A_gspd, GeneStatus.ACCESSORY)
                        ],
                       model_A, self.hit_weights)

        c2_A = Cluster([ModelHit(ch_sctj, mg_A_sctj, GeneStatus.MANDATORY),
                        MultiSystem(ch_sctn, mg_A_sctn, GeneStatus.MANDATORY)],
                       model_A, self.hit_weights)

        c3_A = Cluster([ModelHit(ch_sctj_2, mg_A_sctj, GeneStatus.MANDATORY),
                        MultiSystem(ch_sctn_2, mg_A_sctn, GeneStatus.MANDATORY),
                        ModelHit(ch_gspd_2, mg_A_gspd, GeneStatus.ACCESSORY)
                        ],
                       model_A, self.hit_weights)

        c4_A = Cluster([ModelHit(ch_sctj, mg_A_sctj, GeneStatus.MANDATORY),
                        MultiSystem(ch_sctn_2, mg_A_sctn, GeneStatus.MANDATORY),
                        ModelHit(ch_gspd_2, mg_A_gspd, GeneStatus.ACCESSORY)
                        ],
                       model_A, self.hit_weights)
        c5_A = Cluster([ModelHit(ch_sctj_2, mg_A_sctj, GeneStatus.MANDATORY),
                        MultiSystem(ch_sctn_2, mg_A_sctn, GeneStatus.MANDATORY),
                        ModelHit(ch_gspd, mg_A_gspd, GeneStatus.ACCESSORY)
                        ],
                       model_A, self.hit_weights)

        model_B._min_mandatory_genes_required = 1
        model_B._min_genes_required = 2
        c6_B = Cluster([ModelHit(ch_sctj_flg, mg_B_sctj_flg, GeneStatus.MANDATORY),
                        ModelHit(ch_tadZ, mg_B_tadZ, GeneStatus.ACCESSORY),
                        ModelHit(ch_flgB, mg_B_flgB, GeneStatus.ACCESSORY)],
                       model_B, self.hit_weights)

        model_C._min_mandatory_genes_required = 1
        model_C._min_genes_required = 2
        c7_C = Cluster([ModelHit(ch_sctj_flg, mg_C_sctj_flg, GeneStatus.MANDATORY),
                        ModelHit(ch_tadZ, mg_C_tadZ, GeneStatus.ACCESSORY),
                        ModelHit(ch_flgB, mg_C_flgB, GeneStatus.MANDATORY),
                        ModelHit(ch_gspd, mg_C_gspd, GeneStatus.ACCESSORY)],
                       model_C, self.hit_weights)

        model_D._min_mandatory_genes_required = 1
        model_D._min_genes_required = 1
        c8_D = Cluster([ModelHit(ch_abc, mg_D_abc, GeneStatus.MANDATORY),
                        MultiSystem(ch_sctn, mg_D_sctn, GeneStatus.ACCESSORY)],
                       model_D, self.hit_weights)

        model_E._min_mandatory_genes_required = 1
        model_E._min_genes_required = 1
        c9_E = Cluster([ModelHit(ch_abc, mg_E_abc, GeneStatus.MANDATORY),
                        ModelHit(ch_sctn, mg_E_sctn, GeneStatus.ACCESSORY)],
                       model_E, self.hit_weights)

        model_F._min_mandatory_genes_required = 1
        model_F._min_genes_required = 1
        c10_F = Cluster([ModelHit(ch_abc_2, mg_F_abc, GeneStatus.MANDATORY),
                        ModelHit(ch_sctn, mg_F_sctn, GeneStatus.ACCESSORY)],
                       model_F, self.hit_weights)

        sys_A1 = System(model_A, [c1_A, c2_A], self.cfg.redundancy_penalty())
        # we need to tweek the replicon_id to have stable results
        # whatever the number of tests ran
        # or the tests order
        sys_A1.id = "replicon_id_A1"
        sys_A2 = System(model_A, [c3_A], self.cfg.redundancy_penalty())
        sys_A2.id = "replicon_id_A2"
        sys_A3 = System(model_A, [c4_A], self.cfg.redundancy_penalty())
        sys_A3.id = "replicon_id_A3"
        sys_A4 = System(model_A, [c5_A], self.cfg.redundancy_penalty())
        sys_A4.id = "replicon_id_A4"
        sys_B = System(model_B, [c6_B], self.cfg.redundancy_penalty())
        sys_B.id = "replicon_id_B"
        sys_C = System(model_C, [c7_C], self.cfg.redundancy_penalty())
        sys_C.id = "replicon_id_C"
        sys_D = System(model_D, [c8_D], self.cfg.redundancy_penalty())
        sys_D.id = "replicon_id_D"
        sys_E = System(model_E, [c9_E], self.cfg.redundancy_penalty())
        sys_E.id = "replicon_id_E"
        sys_F = System(model_F, [c10_F], self.cfg.redundancy_penalty())
        sys_F.id = "replicon_id_F"
        self.assertFalse(sys_A1.is_compatible(sys_A3))  # False: same model: share sctj (not multi_system)
        self.assertTrue(sys_A3.is_compatible(sys_A4))  # True: same model: share sctn (multi_system)
        self.assertTrue(sys_A1.is_compatible(sys_A2))   # True: same model: do not share any hit
        self.assertTrue(sys_A1.is_compatible(sys_B))    # True: different models: do not share any hit
        self.assertFalse(sys_A1.is_compatible(sys_C))   # False: diffrent models: share h_gspd
        self.assertFalse(sys_A1.is_compatible(sys_D))   # False: different models: share h_sctn (which is multi_system im model A but not multi model)
        self.assertFalse(sys_D.is_compatible(sys_E))    # False: different modesls: share sctn but sctn is multi_model only in E
        self.assertTrue(sys_E.is_compatible(sys_F))     # True: different Model: share sctn which is multi_model in E and F


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

        h_sctj = CoreHit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = CoreHit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd = CoreHit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_sctj_flg = CoreHit(c_gene_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_flgB = CoreHit(c_gene_flgB, "hit_flgB", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_tadZ = CoreHit(c_gene_tadZ, "hit_tadZ", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        model_1._min_mandatory_genes_required = 2
        model_1._min_genes_required = 2
        c1 = Cluster([ModelHit(h_sctj, gene_sctj, GeneStatus.MANDATORY),
                      ModelHit(h_sctn, gene_sctn, GeneStatus.MANDATORY),
                      ModelHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
                      ],
                     model_1,
                     self.hit_weights)

        c2 = Cluster([ModelHit(h_sctj, gene_sctj, GeneStatus.MANDATORY),
                      ModelHit(h_sctn, gene_sctn, GeneStatus.MANDATORY),
                      ModelHit(h_flgB, gene_gspd, GeneStatus.ACCESSORY)],
                     model_1,
                     self.hit_weights)

        model_2._min_mandatory_genes_required = 1
        model_2._min_genes_required = 2
        c3 = Cluster([ModelHit(h_sctj_flg, gene_sctj_flg, GeneStatus.MANDATORY),
                      ModelHit(h_tadZ, gene_tadZ, GeneStatus.ACCESSORY),
                      ModelHit(h_flgB, gene_flgB, GeneStatus.ACCESSORY)],
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

        h_sctj = CoreHit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = CoreHit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd = CoreHit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd_an = CoreHit(c_gene_flgB, "hit_gspd_an", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_abc = CoreHit(c_gene_abc, "hit_abc", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_toto = CoreHit(c_gene_toto, "hit_toto", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_toto_ho = CoreHit(c_gene_totote, "hit_toto_ho", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_not_in_model = CoreHit(c_gene_not_in_model, "hit_not_in_model", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        model._min_mandatory_genes_required = 2
        model._min_genes_required = 1

        mh_sctj = ModelHit(h_sctj, gene_sctj, GeneStatus.MANDATORY)
        mh_sctn = ModelHit(h_sctn, gene_sctn, GeneStatus.MANDATORY)
        mh_gspd = ModelHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
        mh_abc = ModelHit(h_abc, gene_abc, GeneStatus.FORBIDDEN)
        mh_toto = ModelHit(h_toto, gene_toto, GeneStatus.NEUTRAL)
        c1 = Cluster([mh_sctj,
                      mh_sctn,
                      mh_gspd,
                      mh_abc,
                      mh_toto],
                     model, self.hit_weights)
        s1 = System(model, [c1], self.cfg.redundancy_penalty())

        self.assertDictEqual(s1.mandatory_occ, {'sctJ': [mh_sctj], 'sctN': [mh_sctn]})
        self.assertDictEqual(s1.accessory_occ, {'gspD': [mh_gspd]})
        self.assertDictEqual(s1.neutral_occ, {'toto': [mh_toto]})

        # test with homolog and analog
        mh_sctj = ModelHit(h_sctj, gene_sctj, GeneStatus.MANDATORY)
        mh_sctn = ModelHit(h_sctn, gene_sctn, GeneStatus.MANDATORY)
        mh_gspd = ModelHit(h_gspd_an, gene_gspd, GeneStatus.ACCESSORY)
        mh_abc = ModelHit(h_abc, gene_abc, GeneStatus.FORBIDDEN)
        mh_toto = ModelHit(h_toto_ho, gene_toto, GeneStatus.NEUTRAL)
        c1 = Cluster([mh_sctj,
                      mh_sctn,
                      mh_gspd,
                      mh_abc,
                      mh_toto],
                     model, self.hit_weights)
        s1 = System(model, [c1], self.cfg.redundancy_penalty())

        self.assertDictEqual(s1.mandatory_occ, {'sctJ': [mh_sctj], 'sctN': [mh_sctn]})
        self.assertDictEqual(s1.accessory_occ, {'gspD': [mh_gspd]})
        self.assertDictEqual(s1.neutral_occ, {'toto': [mh_toto]})

        mh_sctj = ModelHit(h_sctj, gene_sctj, GeneStatus.MANDATORY)
        mh_not_in_model = ModelHit(h_not_in_model, gene_not_in_model, GeneStatus.MANDATORY)
        c1 = Cluster([mh_sctj, mh_not_in_model], model, self.hit_weights)
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
        gene_gspd = ModelGene(c_gene_gspd, model, loner=True)
        c_gene_flgB = CoreGene(self.model_location, "flgB", self.profile_factory)
        gene_gspd_an = Exchangeable(c_gene_flgB, gene_gspd, loner=True)
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

        h_sctj = CoreHit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = CoreHit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd = CoreHit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_toto = CoreHit(c_gene_toto, "hit_toto", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        model._min_mandatory_genes_required = 2
        model._min_genes_required = 1

        mh_sctj = ModelHit(h_sctj, gene_sctj, GeneStatus.MANDATORY)
        mh_sctn = ModelHit(h_sctn, gene_sctn, GeneStatus.MANDATORY)
        mh_gspd = ModelHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
        l_h_gspd = Loner(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
        mh_toto = ModelHit(h_toto, gene_toto, GeneStatus.NEUTRAL)
        c1 = Cluster([mh_sctj,
                      mh_sctn,
                      mh_gspd,
                      mh_toto],
                     model, self.hit_weights)
        # 1 cluster
        s1 = System(model, [c1], self.cfg.redundancy_penalty())
        self.assertFalse(s1.multi_loci)

        c2 = Cluster([mh_sctj,
                      mh_sctn],
                     model, self.hit_weights)

        c3 = Cluster([mh_gspd,
                      mh_toto],
                     model, self.hit_weights)
        # 2 clusters of 2 hits each
        s2 = System(model, [c2, c3], self.cfg.redundancy_penalty())
        self.assertTrue(s2.multi_loci)

        # c2 contains 2 hits
        # c4 one hit  one loner
        c4 = Cluster([l_h_gspd], model, self.hit_weights)
        s3 = System(model, [c2, c4], self.cfg.redundancy_penalty())
        self.assertFalse(s3.multi_loci)

        # c2 contains (2 hits)
        # c5 contains one hit not loner
        model._min_mandatory_genes_required = 1
        model._min_genes_required = 1
        c5 = Cluster([mh_toto], model, self.hit_weights)
        s4 = System(model, [c2, c5], self.cfg.redundancy_penalty())
        self.assertTrue(s4.multi_loci)


    def test_get_loners(self):
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
        gene_gspd = ModelGene(c_gene_gspd, model, loner=True)
        c_gene_flgB = CoreGene(self.model_location, "flgB", self.profile_factory)
        gene_gspd_an = Exchangeable(c_gene_flgB, gene_gspd)
        gene_gspd.add_exchangeable(gene_gspd_an)

        c_gene_abc = CoreGene(self.model_location, "abc", self.profile_factory)
        gene_abc = ModelGene(c_gene_abc, model, loner=True, multi_system=True)
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
        model.add_accessory_gene(gene_abc)
        model.add_neutral_gene(gene_toto)

        h_sctj = CoreHit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctj_flg = CoreHit(c_gene_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = CoreHit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn_flg = CoreHit(c_gene_sctn_flg, "hit_sctn_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd = CoreHit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd_an = CoreHit(c_gene_flgB, "hit_gspd_an", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_abc = CoreHit(c_gene_abc, "hit_abc", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_abc_ho = CoreHit(c_gene_tadZ, "hit_abc_ho", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_toto = CoreHit(c_gene_toto, "hit_toto", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_toto_ho = CoreHit(c_gene_totote, "hit_toto_ho", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        model._min_mandatory_genes_required = 2
        model._min_genes_required = 1

        mh_sctj = ModelHit(h_sctj, gene_sctj, GeneStatus.MANDATORY)
        mh_sctn = ModelHit(h_sctn, gene_sctn, GeneStatus.MANDATORY)
        mh_gspd = ModelHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
        l_h_gspd = Loner(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
        l_ms_abc = LonerMultiSystem(h_abc, gene_abc, GeneStatus.ACCESSORY)
        mh_toto = ModelHit(h_toto, gene_toto, GeneStatus.NEUTRAL)

        c1 = Cluster([mh_sctj,
                      mh_sctn,
                      mh_gspd,
                      mh_toto],
                     model, self.hit_weights)
        # 1 cluster
        s1 = System(model, [c1], self.cfg.redundancy_penalty())
        self.assertSetEqual(set(),
                            s1.get_loners())

        c2 = Cluster([mh_sctj,
                      mh_sctn],
                     model, self.hit_weights)
        c3 = Cluster([mh_gspd,
                      mh_toto],
                     model, self.hit_weights)
        # 2 clusters of 2 hits each
        s2 = System(model, [c2, c3], self.cfg.redundancy_penalty())
        self.assertTrue(s2.multi_loci)

        # c2 contains 2 hits
        # c4 one hit one loner
        # c5 one hit Loner Multi Systems
        c4 = Cluster([l_h_gspd], model, self.hit_weights)
        c5 = Cluster([l_ms_abc], model, self.hit_weights)
        s3 = System(model, [c2, c4, c5], self.cfg.redundancy_penalty())
        self.assertSetEqual({l_h_gspd, l_ms_abc},
                            s3.get_loners())


    def test_get_hits_encoding_multisystem(self):
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
        gene_gspd = ModelGene(c_gene_gspd, model, loner=True)
        c_gene_flgB = CoreGene(self.model_location, "flgB", self.profile_factory)
        gene_gspd_an = Exchangeable(c_gene_flgB, gene_gspd, loner=True)
        gene_gspd.add_exchangeable(gene_gspd_an)

        c_gene_abc = CoreGene(self.model_location, "abc", self.profile_factory)
        gene_abc = ModelGene(c_gene_abc, model, loner=True, multi_system=True)
        c_gene_tadZ = CoreGene(self.model_location, "tadZ", self.profile_factory)
        gene_abc_ho = Exchangeable(c_gene_tadZ, gene_abc, loner=True, multi_system=True)
        gene_abc.add_exchangeable(gene_abc_ho)

        c_gene_toto = CoreGene(self.model_location, "toto", self.profile_factory)
        gene_toto = ModelGene(c_gene_toto, model)
        c_gene_totote = CoreGene(self.model_location, "totote", self.profile_factory)
        gene_toto_ho = Exchangeable(c_gene_totote, gene_toto)
        gene_toto.add_exchangeable(gene_toto_ho)

        model.add_mandatory_gene(gene_sctn)
        model.add_mandatory_gene(gene_sctj)
        model.add_accessory_gene(gene_gspd)
        model.add_accessory_gene(gene_abc)
        model.add_neutral_gene(gene_toto)

        h_sctj = CoreHit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = CoreHit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd = CoreHit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_abc = CoreHit(c_gene_abc, "hit_abc", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_abc_ho = CoreHit(c_gene_tadZ, "hit_abc_ho", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_toto = CoreHit(c_gene_toto, "hit_toto", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        model._min_mandatory_genes_required = 2
        model._min_genes_required = 1

        mh_sctj = ModelHit(h_sctj, gene_sctj, GeneStatus.MANDATORY)
        mh_sctn = ModelHit(h_sctn, gene_sctn, GeneStatus.MANDATORY)
        mh_gspd = ModelHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
        l_h_gspd = Loner(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
        l_ms_abc = Loner(h_abc, gene_abc, GeneStatus.ACCESSORY)
        ms_tadZ = ModelHit(h_abc_ho, gene_abc_ho, GeneStatus.ACCESSORY)
        mh_toto = ModelHit(h_toto, gene_toto, GeneStatus.NEUTRAL)

        c1 = Cluster([mh_sctj,
                      mh_sctn,
                      mh_gspd,
                      mh_toto],
                     model, self.hit_weights)
        # 1 cluster
        s1 = System(model, [c1], self.cfg.redundancy_penalty())
        self.assertSetEqual(set(),
                            s1.get_hits_encoding_multisystem())

        c2 = Cluster([mh_sctj,
                      mh_sctn],
                     model, self.hit_weights)

        c3 = Cluster([mh_sctj,
                      mh_sctn,
                      ms_tadZ,
                      mh_toto],
                     model, self.hit_weights)
        c4 = Cluster([l_h_gspd], model, self.hit_weights)
        c5 = Cluster([l_ms_abc], model, self.hit_weights)

        # c2 2 regular hits
        # c3 4 hits with 1 multi_systems (homolog)
        s2 = System(model, [c2, c3], self.cfg.redundancy_penalty())
        self.assertSetEqual({ms_tadZ},
                            s2.get_hits_encoding_multisystem())

        # c2 contains 2 hits
        # c4 one hit one loner
        # c5 one hit Loner MultiSystems
        s3 = System(model, [c2, c4, c5], self.cfg.redundancy_penalty())
        self.assertSetEqual({l_ms_abc},
                            s3.get_hits_encoding_multisystem())


    def test_get_multisystems(self):
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
        gene_gspd = ModelGene(c_gene_gspd, model, loner=True)
        c_gene_flgB = CoreGene(self.model_location, "flgB", self.profile_factory)
        gene_gspd_an = Exchangeable(c_gene_flgB, gene_gspd, loner=True)
        gene_gspd.add_exchangeable(gene_gspd_an)

        c_gene_abc = CoreGene(self.model_location, "abc", self.profile_factory)
        gene_abc = ModelGene(c_gene_abc, model, loner=True, multi_system=True)
        c_gene_tadZ = CoreGene(self.model_location, "tadZ", self.profile_factory)
        gene_abc_ho = Exchangeable(c_gene_tadZ, gene_abc, loner=True, multi_system=True)
        gene_abc.add_exchangeable(gene_abc_ho)

        c_gene_toto = CoreGene(self.model_location, "toto", self.profile_factory)
        gene_toto = ModelGene(c_gene_toto, model)
        c_gene_totote = CoreGene(self.model_location, "totote", self.profile_factory)
        gene_toto_ho = Exchangeable(c_gene_totote, gene_toto)
        gene_toto.add_exchangeable(gene_toto_ho)

        model.add_mandatory_gene(gene_sctn)
        model.add_mandatory_gene(gene_sctj)
        model.add_accessory_gene(gene_gspd)
        model.add_accessory_gene(gene_abc)
        model.add_neutral_gene(gene_toto)

        h_sctj = CoreHit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = CoreHit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd = CoreHit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_abc = CoreHit(c_gene_abc, "hit_abc", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_abc_ho = CoreHit(c_gene_tadZ, "hit_abc_ho", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_toto = CoreHit(c_gene_toto, "hit_toto", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        model._min_mandatory_genes_required = 2
        model._min_genes_required = 1

        mh_sctj = ModelHit(h_sctj, gene_sctj, GeneStatus.MANDATORY)
        mh_sctn = ModelHit(h_sctn, gene_sctn, GeneStatus.MANDATORY)
        mh_gspd = ModelHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
        l_h_gspd = Loner(h_gspd, gene_gspd, GeneStatus.ACCESSORY)
        l_ms_abc = LonerMultiSystem(h_abc, gene_abc, GeneStatus.ACCESSORY)
        ms_tadZ = MultiSystem(h_abc_ho, gene_abc_ho, GeneStatus.ACCESSORY)
        mh_toto = ModelHit(h_toto, gene_toto, GeneStatus.NEUTRAL)

        c1 = Cluster([mh_sctj,
                      mh_sctn,
                      mh_gspd,
                      mh_toto],
                     model, self.hit_weights)
        # 1 cluster
        s1 = System(model, [c1], self.cfg.redundancy_penalty())
        self.assertSetEqual(set(),
                            s1.get_multisystems())

        c2 = Cluster([mh_sctj,
                      mh_sctn],
                     model, self.hit_weights)

        c3 = Cluster([mh_sctj,
                      mh_sctn,
                      ms_tadZ,
                      mh_toto],
                     model, self.hit_weights)
        c4 = Cluster([l_h_gspd], model, self.hit_weights)
        c5 = Cluster([l_ms_abc], model, self.hit_weights)

        # c2 2 regular hits
        # c3 4 hits with 1 multi_systems (homolog)
        s2 = System(model, [c2, c3], self.cfg.redundancy_penalty())
        self.assertSetEqual({ms_tadZ},
                            s2.get_multisystems())

        # c2 contains 2 hits
        # c4 one hit one loner
        # c5 one hit Loner MultiSystems
        s3 = System(model, [c2, c4, c5], self.cfg.redundancy_penalty())
        self.assertSetEqual({l_ms_abc},
                            s3.get_multisystems())

    def test_fulfilled_function(self):
        model = Model("foo/T2SS", 11)

        cg_gspd = CoreGene(self.model_location, "gspD", self.profile_factory)
        cg_sctc = CoreGene(self.model_location, "sctC", self.profile_factory)
        cg_sctj = CoreGene(self.model_location, "sctJ", self.profile_factory)
        cg_sctj_flg = CoreGene(self.model_location, "sctJ_FLG", self.profile_factory)
        cg_tadz = CoreGene(self.model_location, "tadZ", self.profile_factory)
        cg_sctn = CoreGene(self.model_location, "sctN", self.profile_factory)

        mg_gspd = ModelGene(cg_gspd, model)
        mg_sctc = ModelGene(cg_sctc, model)
        mg_sctj = ModelGene(cg_sctj, model)
        mg_sctj_flg = Exchangeable(cg_sctj_flg, mg_sctj)
        mg_sctj.add_exchangeable(mg_sctj_flg)
        mg_tadz = ModelGene(cg_tadz, model)
        mg_sctn = ModelGene(cg_sctn, model)

        model.add_mandatory_gene(mg_gspd)
        model.add_mandatory_gene(mg_sctc)
        model.add_accessory_gene(mg_sctj)
        model.add_accessory_gene(mg_tadz)
        model.add_accessory_gene(mg_sctn)

        #     CoreHit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        ch_gspd = CoreHit(cg_gspd, "gspD", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        ch_sctc = CoreHit(cg_sctc, "sctC", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        ch_sctj = CoreHit(cg_sctj, "sctJ", 10, "replicon_1", 30, 1.0, 20.0, 1.0, 1.0, 10, 20)
        ch_sctj_flg = CoreHit(cg_sctj_flg, "sctJ_FLG", 10, "replicon_1", 40, 1.0, 50.0, 1.0, 1.0, 10, 20)
        ch_tadz = CoreHit(cg_tadz, "tadZ", 10, "replicon_1", 50, 1.0, 50.0, 1.0, 1.0, 10, 20)

        mh_gspd = ModelHit(ch_gspd, mg_gspd, GeneStatus.MANDATORY)
        mh_sctc = ModelHit(ch_sctc, mg_sctc, GeneStatus.MANDATORY)
        mh_sctj = ModelHit(ch_sctj, mg_sctj, GeneStatus.ACCESSORY)
        mh_sctj_flg = ModelHit(ch_sctj_flg, mg_sctj_flg, GeneStatus.ACCESSORY)
        mh_tadz = ModelHit(ch_tadz, mg_tadz, GeneStatus.ACCESSORY)

        c1 = Cluster([mh_gspd, mh_sctc], model, self.hit_weights)
        c2 = Cluster([mh_sctj, mh_tadz], model, self.hit_weights)
        c3 = Cluster([mh_tadz, mh_sctj_flg], model, self.hit_weights)

        # One system with 2 regular clusters, No exchangeable
        s1 = System(model, [c1, c2])
        self.assertSetEqual(s1.fulfilled_function(mg_gspd),
                            {mg_gspd.name})
        self.assertFalse(s1.fulfilled_function(mg_sctn))

        # test with several genes
        self.assertSetEqual(s1.fulfilled_function(mg_gspd, mg_sctn, mg_tadz),
                            {mg_gspd.name, mg_tadz.name})

        # One cluster contains exchangeable
        s2 = System(model, [c1, c3])
        self.assertSetEqual(s2.fulfilled_function(mg_sctj),
                            {mg_sctj.name})

        # Test with function as string
        self.assertSetEqual(s1.fulfilled_function(mg_gspd.name),
                            {mg_gspd.name})
