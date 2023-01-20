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

from macsypy.hit import CoreHit, ModelHit, HitWeight
from macsypy.config import Config, MacsyDefaults
from macsypy.gene import CoreGene, ModelGene, GeneStatus
from macsypy.profile import ProfileFactory
from macsypy.model import Model
from macsypy.registries import ModelLocation
from macsypy.cluster import Cluster
from macsypy.system import RejectedCandidate

from tests import MacsyTest


class TestRejectedCandidate(MacsyTest):

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
        self.hit_weights = HitWeight(self.cfg.hit_weights())


    def test_init(self):
        model = Model("foo/T2SS", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_1 = ModelGene(c_gene_1, model)
        model.add_mandatory_gene(gene_1)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)
        gene_2 = ModelGene(c_gene_2, model)
        model.add_accessory_gene(gene_2)

        #     CoreHit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = CoreHit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        v_h10 = ModelHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = CoreHit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        v_h20 = ModelHit(h20, gene_2, GeneStatus.ACCESSORY)
        c1 = Cluster([v_h10, v_h20], model, self.hit_weights)
        r_c = RejectedCandidate(model, c1, ["bla"])
        self.assertListEqual(r_c.clusters, [c1])
        self.assertEqual(r_c.reasons, ['bla'])


    def test_str(self):
        model = Model("foo/T2SS", 11)

        c_gene_1 = CoreGene(self.model_location, "gspD", self.profile_factory)
        gene_1 = ModelGene(c_gene_1, model)
        model.add_mandatory_gene(gene_1)
        c_gene_2 = CoreGene(self.model_location, "sctC", self.profile_factory)
        gene_2 = ModelGene(c_gene_2, model)
        model.add_accessory_gene(gene_2)

        #     CoreHit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = CoreHit(c_gene_1, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        v_h10 = ModelHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = CoreHit(c_gene_2, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        v_h20 = ModelHit(h20, gene_2, GeneStatus.ACCESSORY)
        h40 = CoreHit(c_gene_1, "h40", 10, "replicon_1", 40, 1.0, 10.0, 1.0, 1.0, 10, 20)
        v_h40 = ModelHit(h40, gene_1, GeneStatus.MANDATORY)
        h50 = CoreHit(c_gene_2, "h50", 10, "replicon_1", 50, 1.0, 20.0, 1.0, 1.0, 10, 20)
        v_h50 = ModelHit(h50, gene_2, GeneStatus.ACCESSORY)
        c1 = Cluster([v_h10, v_h20], model, self.hit_weights)
        c2 = Cluster([v_h40, v_h50], model, self.hit_weights)
        r_c = RejectedCandidate(model, [c1, c2], ["bla"])

        expected_str = """Cluster:
- model = T2SS
- replicon = replicon_1
- hits = (h10, gspD, 10), (h20, sctC, 20)
Cluster:
- model = T2SS
- replicon = replicon_1
- hits = (h40, gspD, 40), (h50, sctC, 50)
This candidate has been rejected because:
\t- bla
"""
        self.assertEqual(expected_str, str(r_c))


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
        rc = RejectedCandidate(model, [Cluster([v_hit_1, v_hit_2], model, self.hit_weights),
                                      Cluster([v_hit_3], model, self.hit_weights)],
                                      ["bla bla"])

        self.assertEqual(rc.hits, [v_hit_1, v_hit_2, v_hit_3])
        self.assertEqual(rc.reasons, ["bla bla"])
