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
import itertools

from macsypy.hit import Hit, ValidHit, HitWeight
from macsypy.config import Config, MacsyDefaults
from macsypy.gene import CoreGene, ModelGene, Exchangeable, GeneStatus
from macsypy.profile import ProfileFactory
from macsypy.model import Model
from macsypy.registries import ModelLocation
from macsypy.cluster import Cluster
from macsypy.system import System, HitSystemTracker, LikelySystem, UnlikelySystem, AbstractSetOfHits
from macsypy.serialization import TxtSystemSerializer, TsvSystemSerializer, TsvSolutionSerializer, \
    TxtLikelySystemSerializer, TxtUnikelySystemSerializer

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
        # reset the uniq id number for AbstractSetOfHits
        # to have predictable results
        AbstractSetOfHits._id = itertools.count(1)

    def test_SystemSerializer_str(self):
        model_name = 'foo'
        model_location = ModelLocation(path=os.path.join(self.cfg.models_dir(), model_name))
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

        h_sctj = Hit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = Hit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd = Hit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_sctj_flg = Hit(c_gene_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_tadZ = Hit(c_gene_tadZ, "hit_tadZ", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

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
                      ValidHit(h_gspd, gene_gspd, GeneStatus.ACCESSORY)],
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
        gene_sctn = ModelGene(c_gene_sctn, model)
        c_gene_sctn_flg = CoreGene(self.model_location, "sctN_FLG", self.profile_factory)
        gene_sctn_flg = Exchangeable(c_gene_sctn_flg, gene_sctn)
        gene_sctn.add_exchangeable(gene_sctn_flg)
        model.add_accessory_gene(gene_sctn)

        h_gspd = Hit(c_gene_gspd, "h_gspd", 803, "replicon_id", 10, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_h_gspd = ValidHit(h_gspd, gene_gspd, GeneStatus.MANDATORY)
        h_sctj = Hit(c_gene_sctj, "h_sctj", 803, "replicon_id", 20, 1.0, 1.0, 1.0, 1.0, 20, 30)
        v_h_sctj = ValidHit(h_sctj, gene_sctj, GeneStatus.ACCESSORY)
        h_sctn_flg = Hit(c_gene_sctn_flg, "h_sctn_flg", 803, "replicon_id", 30, 1.0, 1.0, 1.0, 1.0, 30, 40)
        v_h_sctn_flg = ValidHit(h_sctn_flg, gene_sctn_flg, GeneStatus.ACCESSORY)
        c1 = Cluster([v_h_gspd, v_h_sctj], model, self.hit_weights)
        c2 = Cluster([v_h_sctn_flg], model, self.hit_weights)
        sys_multi_loci = System(model, [c1, c2], self.cfg.redundancy_penalty())
        hit_multi_sys_tracker = HitSystemTracker([sys_multi_loci])
        system_serializer = TsvSystemSerializer()

        sys_tsv = "\t".join(["replicon_id", "h_gspd", "gspD", "10", "foo/T2SS", sys_multi_loci.id, "1",
                             "1.000", "1.900", "1", "gspD", "mandatory", "803",
                             "1.0", "1.000", "1.000", "1.000", "10", "20", ""])
        sys_tsv += "\n"
        sys_tsv += "\t".join(["replicon_id", "h_sctj", "sctJ", "20", "foo/T2SS", sys_multi_loci.id, "1",
                              "1.000", "1.900", "1", "sctJ", "accessory", "803",
                              "1.0", "1.000", "1.000", "1.000", "20", "30", ""])
        sys_tsv += "\n"
        sys_tsv += "\t".join(["replicon_id", "h_sctn_flg", "sctN_FLG", "30", "foo/T2SS", sys_multi_loci.id, "1",
                              "1.000", "1.900", "1", "sctN", "accessory", "803",
                              "1.0", "1.000", "1.000", "1.000", "30", "40", ""])
        sys_tsv += "\n"
        self.assertEqual(sys_tsv, system_serializer.serialize(sys_multi_loci, hit_multi_sys_tracker))


    def test_SolutionSerializer_tsv(self):
        model_name = 'foo'
        model_location = ModelLocation(path=os.path.join(self.cfg.models_dir(), model_name))
        model_A = Model("foo/A", 10)
        model_B = Model("foo/B", 10)

        c_gene_sctn_flg = CoreGene(model_location, "sctN_FLG", self.profile_factory)
        gene_sctn_flg = ModelGene(c_gene_sctn_flg, model_B)
        c_gene_sctj_flg = CoreGene(model_location, "sctJ_FLG", self.profile_factory)
        gene_sctj_flg = ModelGene(c_gene_sctj_flg, model_B)
        c_gene_flgB = CoreGene(model_location, "flgB", self.profile_factory)
        gene_flgB = ModelGene(c_gene_flgB, model_B)
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
        model_B.add_accessory_gene(gene_flgB)
        model_B.add_accessory_gene(gene_tadZ)

        h_sctj = Hit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = Hit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd = Hit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_sctj_flg = Hit(c_gene_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_flgB = Hit(c_gene_flgB, "hit_flgB", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_tadZ = Hit(c_gene_tadZ, "hit_tadZ", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

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

        sys_A = System(model_A, [c1, c2], self.cfg.redundancy_penalty())
        sys_A.id = "sys_id_A"
        sys_B = System(model_B, [c3], self.cfg.redundancy_penalty())
        sys_B.id = "sys_id_B"

        sol = [sys_A, sys_B]
        sol_id = '12'

        hit_multi_sys_tracker = HitSystemTracker([sys_A, sys_B])
        system_serializer = TsvSolutionSerializer()

        sol_tsv = '\t'.join([sol_id, 'replicon_id', 'hit_sctj', 'sctJ', '1', 'foo/A', 'sys_id_A',
                            '2', '1.000', '1.500', '2', 'sctJ', 'mandatory',
                            '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id, 'replicon_id', 'hit_sctn', 'sctN', '1', 'foo/A', 'sys_id_A',
                             '2', '1.000', '1.500', '2', 'sctN', 'mandatory',
                             '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id, 'replicon_id', 'hit_gspd', 'gspD', '1', 'foo/A', 'sys_id_A',
                             '2', '1.000', '1.500', '2', 'gspD', 'accessory',
                             '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id, 'replicon_id', 'hit_sctj', 'sctJ', '1', 'foo/A', 'sys_id_A',
                             '2', '1.000', '1.500', '2', 'sctJ', 'mandatory',
                             '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id, 'replicon_id', 'hit_sctn', 'sctN', '1', 'foo/A', 'sys_id_A',
                             '2', '1.000', '1.500', '2', 'sctN', 'mandatory',
                             '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id, 'replicon_id', 'hit_sctj_flg', 'sctJ_FLG', '1', 'foo/B', 'sys_id_B',
                             '1', '0.750', '2.000', '1', 'sctJ_FLG', 'mandatory',
                             '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id, 'replicon_id', 'hit_tadZ', 'tadZ', '1', 'foo/B', 'sys_id_B',
                             '1', '0.750', '2.000', '1', 'tadZ', 'accessory',
                             '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id, 'replicon_id', 'hit_flgB', 'flgB', '1', 'foo/B', 'sys_id_B',
                             '1', '0.750', '2.000', '1', 'flgB', 'accessory',
                             '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += "\n"
        ser = system_serializer.serialize(sol, sol_id, hit_multi_sys_tracker)
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

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 2, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(c_gene_sctn, "hit_3", 803, "replicon_id", 3, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        hit_4 = Hit(c_gene_abc, "hit_4", 803, "replicon_id", 4, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_4 = ValidHit(hit_4, gene_abc, GeneStatus.FORBIDDEN)

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

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 2, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(c_gene_sctn, "hit_3", 803, "replicon_id", 3, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        hit_4 = Hit(c_gene_abc, "hit_4", 803, "replicon_id", 4, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_4 = ValidHit(hit_4, gene_abc, GeneStatus.FORBIDDEN)
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