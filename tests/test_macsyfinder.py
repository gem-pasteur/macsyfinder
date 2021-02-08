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
import sys
import shutil
import tempfile
import argparse
from io import StringIO
import logging
import unittest
import itertools

from macsypy.config import Config, MacsyDefaults
from macsypy.gene import CoreGene, ModelGene, Exchangeable, GeneStatus, GeneBank
from macsypy.profile import ProfileFactory
from macsypy.registries import ModelLocation
from macsypy.hit import Hit, ValidHit, HitWeight
from macsypy.model import Model, ModelBank
from macsypy.system import System, HitSystemTracker, RejectedClusters, AbstractSetOfHits, LikelySystem, UnlikelySystem
from macsypy.cluster import Cluster

from macsypy.scripts.macsyfinder import systems_to_txt, systems_to_tsv, rejected_clst_to_txt, solutions_to_tsv, \
    likely_systems_to_txt, likely_systems_to_tsv, unlikely_systems_to_txt
from macsypy.scripts.macsyfinder import list_models, parse_args, search_systems

import macsypy
from tests import MacsyTest, which


class TestMacsyfinder(MacsyTest):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        AbstractSetOfHits._id = itertools.count(1)

    def tearDown(self):
        try:
            shutil.rmtree(self.tmp_dir)
        except:
            pass


    def test_list_models(self):
        cmd_args = argparse.Namespace()
        cmd_args.models_dir = os.path.join(self._data_dir, 'fake_model_dir')
        cmd_args.list_models = True
        rcv_list_models = list_models(cmd_args)
        exp_list_models = """set_1
      /def_1_1
      /def_1_2
      /def_1_3
      /def_1_4
set_2
      /level_1
              /def_1_1
              /def_1_2
              /level_2
                      /def_2_3
                      /def_2_4
"""
        self.assertEqual(exp_list_models, rcv_list_models)


    def test_systems_to_txt(self):
        system_str = f"""# macsyfinder {macsypy.__version__}
# {' '.join(sys.argv)}
# No Systems found
"""
        f_out = StringIO()
        track_multi_systems_hit = HitSystemTracker([])
        systems_to_txt([], track_multi_systems_hit, f_out)
        self.assertMultiLineEqual(system_str, f_out.getvalue())

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
        # test if id is well incremented
        gene_name = "gspD"
        c_gene_gspd = CoreGene(models_location, gene_name, profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        gene_name = "sctJ"
        c_gene_sctj = CoreGene(models_location, gene_name, profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        system_1 = System(model,
                          [Cluster([v_hit_1, v_hit_2], model, HitWeight(**cfg.hit_weights()))],
                          cfg.redundancy_penalty())

        system_str = f"""# macsyfinder {macsypy.__version__}
# {' '.join(sys.argv)}
# Systems found:

system id = replicon_id_T2SS_{next(System._id) - 1}
model = foo/T2SS
replicon = replicon_id
clusters = [('hit_1', 'gspD', 1), ('hit_2', 'sctJ', 1)]
occ = 1
wholeness = 1.000
loci nb = 1
score = 1.500

mandatory genes:
\t- gspD: 1 (gspD)

accessory genes:
\t- sctJ: 1 (sctJ)

neutral genes:

============================================================
"""

        f_out = StringIO()
        track_multi_systems_hit = HitSystemTracker([system_1])
        systems_to_txt([system_1], track_multi_systems_hit, f_out)
        self.assertMultiLineEqual(system_str, f_out.getvalue())


    def test_systems_to_tsv(self):
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
            c_gene_gspd = CoreGene(models_location, gene_name, profile_factory)
            gene_gspd = ModelGene(c_gene_gspd, model)
            model.add_mandatory_gene(gene_gspd)
            gene_name = "sctJ"
            c_gene_sctj = CoreGene(models_location, gene_name, profile_factory)
            gene_sctj = ModelGene(c_gene_sctj, model)
            model.add_accessory_gene(gene_sctj)

            hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
            v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
            hit_2 = Hit(c_gene_sctj, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
            v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
            system_1 = System(model,
                              [Cluster([v_hit_1, v_hit_2], model, HitWeight(**cfg.hit_weights()))],
                              cfg.redundancy_penalty())

            system_tsv = f"""# macsyfinder {macsypy.__version__}
# {' '.join(sys.argv)}
# Systems found:
"""
            system_tsv += "\t".join(["replicon", "hit_id", "gene_name", "hit_pos", "model_fqn", "sys_id", "sys_loci",
                                     "sys_wholeness", "sys_score", "sys_occ", "hit_gene_ref", "hit_status",
                                     "hit_seq_len", "hit_i_eval", "hit_score", "hit_profile_cov", "hit_seq_cov",
                                     "hit_begin_match", "hit_end_match", "used_in"])
            system_tsv += "\n"
            system_tsv += "\t".join([ "replicon_id", "hit_1", "gspD", "1", "foo/T2SS", system_1.id,
                                     "1", "1.000", "1.500", "1", "gspD", "mandatory", "803", "1.0", "1.000", "1.000",
                                     "1.000", "10", "20", ""])
            system_tsv += "\n"
            system_tsv += "\t".join(["replicon_id", "hit_2", "sctJ", "1", "foo/T2SS", system_1.id, "1", "1.000",
                                     "1.500", "1", "sctJ", "accessory", "803", "1.0", "1.000", "1.000", "1.000", "10",
                                     "20", ""])
            system_tsv += "\n\n"

            f_out = StringIO()
            track_multi_systems_hit = HitSystemTracker([system_1])
            systems_to_tsv([system_1], track_multi_systems_hit, f_out)
            self.assertMultiLineEqual(system_tsv, f_out.getvalue())

            # test No system found
            system_str = f"""# macsyfinder {macsypy.__version__}
# {' '.join(sys.argv)}
# No Systems found
"""
            f_out = StringIO()
            track_multi_systems_hit = HitSystemTracker([])
            systems_to_tsv([], track_multi_systems_hit, f_out)
            self.assertMultiLineEqual(system_str, f_out.getvalue())


    def test_solutions_to_tsv(self):
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

        model_A = Model("foo/A", 10)
        model_B = Model("foo/B", 10)
        model_C = Model("foo/C", 10)

        c_gene_sctn_flg = CoreGene(models_location, "sctN_FLG", profile_factory)
        gene_sctn_flg = ModelGene(c_gene_sctn_flg, model_B)
        c_gene_sctj_flg = CoreGene(models_location, "sctJ_FLG", profile_factory)
        gene_sctj_flg = ModelGene(c_gene_sctj_flg, model_B)
        c_gene_flgB = CoreGene(models_location, "flgB", profile_factory)
        gene_flgB = ModelGene(c_gene_flgB, model_B)
        c_gene_tadZ = CoreGene(models_location, "tadZ", profile_factory)
        gene_tadZ = ModelGene(c_gene_tadZ, model_B)

        c_gene_sctn = CoreGene(models_location, "sctN", profile_factory)
        gene_sctn = ModelGene(c_gene_sctn, model_A)
        gene_sctn_hom = Exchangeable(c_gene_sctn_flg, gene_sctn)
        gene_sctn.add_exchangeable(gene_sctn_hom)

        c_gene_sctj = CoreGene(models_location, "sctJ", profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model_A)
        gene_sctj_an = Exchangeable(c_gene_sctj_flg, gene_sctj)
        gene_sctj.add_exchangeable(gene_sctj_an)

        c_gene_gspd = CoreGene(models_location, "gspD", profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model_A)
        gene_gspd_an = Exchangeable(c_gene_flgB, gene_gspd)
        gene_gspd.add_exchangeable(gene_gspd_an)

        c_gene_abc = CoreGene(models_location, "abc", profile_factory)
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

        h_sctj = Hit(c_gene_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = Hit(c_gene_sctn, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd = Hit(c_gene_gspd, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        h_sctj_flg = Hit(c_gene_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_flgB = Hit(c_gene_flgB, "hit_flgB", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_tadZ = Hit(c_gene_tadZ, "hit_tadZ", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

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

        sys_A = System(model_A, [c1, c2], cfg.redundancy_penalty())
        sys_A.id = "sys_id_A"
        sys_B = System(model_B, [c3], cfg.redundancy_penalty())
        sys_B.id = "sys_id_B"
        sys_C = System(model_C, [c4], cfg.redundancy_penalty())
        sys_C.id = "sys_id_C"

        sol_1 = [sys_A, sys_B]
        sol_2 = [sys_A, sys_C]
        sol_id_1 = '1'
        sol_id_2 = '2'

        sol_tsv = f"""# macsyfinder {macsypy.__version__}
# {' '.join(sys.argv)}
# Systems found:
"""
        sol_tsv += "\t".join(["sol_id", "replicon", "hit_id", "gene_name", "hit_pos", "model_fqn", "sys_id", "sys_loci",
                              "sys_wholeness", "sys_score", "sys_occ", "hit_gene_ref", "hit_status",
                              "hit_seq_len", "hit_i_eval", "hit_score", "hit_profile_cov", "hit_seq_cov",
                              "hit_begin_match", "hit_end_match", "used_in"])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_1, 'replicon_id', 'hit_sctj', 'sctJ', '1', 'foo/A', 'sys_id_A',
                              '2', '1.000', '1.500', '2', 'sctJ', 'mandatory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_1, 'replicon_id', 'hit_sctn', 'sctN', '1', 'foo/A', 'sys_id_A',
                              '2', '1.000', '1.500', '2', 'sctN', 'mandatory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_1, 'replicon_id', 'hit_gspd', 'gspD', '1', 'foo/A', 'sys_id_A',
                              '2', '1.000', '1.500', '2', 'gspD', 'accessory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_1, 'replicon_id', 'hit_sctj', 'sctJ', '1', 'foo/A', 'sys_id_A',
                              '2', '1.000', '1.500', '2', 'sctJ', 'mandatory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_1, 'replicon_id', 'hit_sctn', 'sctN', '1', 'foo/A', 'sys_id_A',
                              '2', '1.000', '1.500', '2', 'sctN', 'mandatory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_1, 'replicon_id', 'hit_sctj_flg', 'sctJ_FLG', '1', 'foo/B', 'sys_id_B',
                              '1', '0.750', '2.000', '1', 'sctJ_FLG', 'mandatory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_1, 'replicon_id', 'hit_tadZ', 'tadZ', '1', 'foo/B', 'sys_id_B',
                              '1', '0.750', '2.000', '1', 'tadZ', 'accessory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_1, 'replicon_id', 'hit_flgB', 'flgB', '1', 'foo/B', 'sys_id_B',
                              '1', '0.750', '2.000', '1', 'flgB', 'accessory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_2, 'replicon_id', 'hit_sctj', 'sctJ', '1', 'foo/A', 'sys_id_A',
                              '2', '1.000', '1.500', '2', 'sctJ', 'mandatory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_2, 'replicon_id', 'hit_sctn', 'sctN', '1', 'foo/A', 'sys_id_A',
                              '2', '1.000', '1.500', '2', 'sctN', 'mandatory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_2, 'replicon_id', 'hit_gspd', 'gspD', '1', 'foo/A', 'sys_id_A',
                              '2', '1.000', '1.500', '2', 'gspD', 'accessory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_2, 'replicon_id', 'hit_sctj', 'sctJ', '1', 'foo/A', 'sys_id_A',
                              '2', '1.000', '1.500', '2', 'sctJ', 'mandatory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_2, 'replicon_id', 'hit_sctn', 'sctN', '1', 'foo/A', 'sys_id_A',
                              '2', '1.000', '1.500', '2', 'sctN', 'mandatory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', ''])
        sol_tsv += "\n"
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_2, 'replicon_id', 'hit_sctj_flg', 'sctJ_FLG', '1', 'foo/C', 'sys_id_C',
                              '1', '0.800', '3.000', '1', 'sctJ_FLG', 'mandatory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', 'sys_id_B'])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_2, 'replicon_id', 'hit_tadZ', 'tadZ', '1', 'foo/C', 'sys_id_C',
                              '1', '0.800', '3.000', '1', 'tadZ', 'accessory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', 'sys_id_B'])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_2, 'replicon_id', 'hit_flgB', 'flgB', '1', 'foo/C', 'sys_id_C',
                              '1', '0.800', '3.000', '1', 'flgB', 'mandatory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', 'sys_id_B'])
        sol_tsv += "\n"
        sol_tsv += '\t'.join([sol_id_2, 'replicon_id', 'hit_gspd', 'gspD', '1', 'foo/C', 'sys_id_C',
                              '1', '0.800', '3.000', '1', 'gspD', 'accessory',
                              '803', '1.0', '1.000', '1.000', '1.000', '10', '20', 'sys_id_A'])
        sol_tsv += "\n"
        sol_tsv += "\n"

        f_out = StringIO()
        hit_multi_sys_tracker = HitSystemTracker([sys_A, sys_B])
        solutions_to_tsv([sol_1, sol_2], hit_multi_sys_tracker, f_out)
        self.assertMultiLineEqual(sol_tsv, f_out.getvalue())


    def test_rejected_clst_to_txt(self):
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
        model.add_mandatory_gene(gene_1)
        model.add_accessory_gene(gene_2)

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(c_gene_gspd, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        v_h10 = ValidHit(h10, gene_1, GeneStatus.MANDATORY)
        h20 = Hit(c_gene_sctc, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        v_h20 = ValidHit(h20, gene_2, GeneStatus.ACCESSORY)
        h40 = Hit(c_gene_gspd, "h10", 10, "replicon_1", 40, 1.0, 10.0, 1.0, 1.0, 10, 20)
        v_h40 = ValidHit(h40, gene_1, GeneStatus.MANDATORY)
        h50 = Hit(c_gene_sctc, "h20", 10, "replicon_1", 50, 1.0, 20.0, 1.0, 1.0, 10, 20)
        v_h50 = ValidHit(h50, gene_2, GeneStatus.ACCESSORY)
        hit_weights = HitWeight(**cfg.hit_weights())
        c1 = Cluster([v_h10, v_h20], model, hit_weights)
        c2 = Cluster([v_h40, v_h50], model, hit_weights)
        r_c = RejectedClusters(model, [c1, c2], ["The reasons to reject this clusters"])

        rej_clst_str = f"""# macsyfinder {macsypy.__version__}
# {' '.join(sys.argv)}
# Rejected clusters:

Cluster:
- model = T2SS
- replicon = replicon_1
- hits = (h10, gspD, 10), (h20, sctC, 20)
Cluster:
- model = T2SS
- replicon = replicon_1
- hits = (h10, gspD, 40), (h20, sctC, 50)
These clusters have been rejected because:
\t- The reasons to reject this clusters
============================================================
"""

        f_out = StringIO()
        rejected_clst_to_txt([r_c], f_out)
        self.maxDiff = None
        self.assertMultiLineEqual(rej_clst_str, f_out.getvalue())

        rej_clst_str = f"""# macsyfinder {macsypy.__version__}
# {' '.join(sys.argv)}
# No Rejected clusters
"""
        f_out = StringIO()
        rejected_clst_to_txt([], f_out)
        self.assertMultiLineEqual(rej_clst_str, f_out.getvalue())


    def test_likely_systems_to_txt(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'unordered'
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
        # test if id is well incremented
        gene_name = "gspD"
        c_gene_gspd = CoreGene(models_location, gene_name, profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        gene_name = "sctJ"
        c_gene_sctj = CoreGene(models_location, gene_name, profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)
        gene_name = "sctC"
        c_gene_sctc = CoreGene(models_location, gene_name, profile_factory)
        gene_sctc = ModelGene(c_gene_sctc, model)
        model.add_neutral_gene(gene_sctc)
        gene_name = "tadZ"
        c_gene_tadz = CoreGene(models_location, gene_name, profile_factory)
        gene_tadz = ModelGene(c_gene_tadz, model)
        model.add_forbidden_gene(gene_tadz)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 804, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(c_gene_sctc, "hit_3", 805, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctc, GeneStatus.NEUTRAL)
        hit_4 = Hit(c_gene_tadz, "hit_4", 806, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_4 = ValidHit(hit_4, gene_tadz, GeneStatus.FORBIDDEN)

        system_1 = LikelySystem(model, [v_hit_1], [v_hit_2], [v_hit_3], [v_hit_4])

        system_str = f"""# macsyfinder {macsypy.__version__}
# {' '.join(sys.argv)}
# Systems found:

This replicon contains genetic materials needed for system foo/T2SS
WARNING there quorum is reached but there is also some forbidden genes.

system id = replicon_id_T2SS_1
model = foo/T2SS
replicon = replicon_id
hits = [('hit_1', 'gspD', 1), ('hit_2', 'sctJ', 1), ('hit_3', 'sctC', 1), ('hit_4', 'tadZ', 1)]
wholeness = 1.000

mandatory genes:
\t- gspD: 1 (gspD)

accessory genes:
\t- sctJ: 1 (sctJ)

neutral genes:
\t- sctC: 1 (sctC)

forbidden genes:
\t- tadZ: 1 (tadZ)

Use ordered replicon to have better prediction.

"""

        f_out = StringIO()
        track_multi_systems_hit = HitSystemTracker([system_1])
        likely_systems_to_txt([system_1], track_multi_systems_hit, f_out)
        self.assertMultiLineEqual(system_str, f_out.getvalue())

        f_out = StringIO()
        likely_systems_to_txt([], track_multi_systems_hit, f_out)
        expected_out = f"""# macsyfinder {macsypy.__version__}
# {' '.join(sys.argv)}
# No Likely Systems found
"""
        self.assertEqual(expected_out, f_out.getvalue())


    def test_likely_systems_to_tsv(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'unordered'
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
        # test if id is well incremented
        gene_name = "gspD"
        c_gene_gspd = CoreGene(models_location, gene_name, profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        gene_name = "sctJ"
        c_gene_sctj = CoreGene(models_location, gene_name, profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)
        gene_name = "sctC"
        c_gene_sctc = CoreGene(models_location, gene_name, profile_factory)
        gene_sctc = ModelGene(c_gene_sctc, model)
        model.add_neutral_gene(gene_sctc)
        gene_name = "tadZ"
        c_gene_tadz = CoreGene(models_location, gene_name, profile_factory)
        gene_tadz = ModelGene(c_gene_tadz, model)
        model.add_forbidden_gene(gene_tadz)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 804, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(c_gene_sctc, "hit_3", 805, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctc, GeneStatus.NEUTRAL)
        hit_4 = Hit(c_gene_tadz, "hit_4", 806, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_4 = ValidHit(hit_4, gene_tadz, GeneStatus.FORBIDDEN)

        system_1 = LikelySystem(model, [v_hit_1], [v_hit_2], [v_hit_3], [v_hit_4])

        sol_tsv = f"""# macsyfinder {macsypy.__version__}
# {' '.join(sys.argv)}
# Likely Systems found:"""
        sol_tsv += "\n\n"
        sol_tsv += "\t".join(["replicon", "hit_id", "gene_name", "hit_pos", "model_fqn", "sys_id", "sys_wholeness",
                              "hit_gene_ref", "hit_status", "hit_seq_len", "hit_i_eval", "hit_score",
                              "hit_profile_cov", "hit_seq_cov", "hit_begin_match", "hit_end_match", "used_in"])
        sol_tsv += "\n"
        sol_tsv += '\t'.join(["replicon_id", "hit_1", "gspD", "1", "foo/T2SS", "replicon_id_T2SS_1", "1.000",
                              "gspD", "mandatory", "803", "1.0", "1.000", "1.000", "1.000", "10", "20", ""])
        sol_tsv += "\n"
        sol_tsv += '\t'.join(["replicon_id", "hit_2", "sctJ", "1", "foo/T2SS", "replicon_id_T2SS_1", "1.000",
                              "sctJ", "accessory", "804", "1.0", "1.000", "1.000", "1.000", "10", "20", ""])
        sol_tsv += "\n"
        sol_tsv += '\t'.join(["replicon_id", "hit_4", "tadZ", "1", "foo/T2SS", "replicon_id_T2SS_1", "1.000",
                              "tadZ", "forbidden", "806", "1.0", "1.000", "1.000", "1.000", "10", "20", ""])
        sol_tsv += "\n"
        sol_tsv += '\t'.join(["replicon_id", "hit_3", "sctC", "1", "foo/T2SS", "replicon_id_T2SS_1", "1.000",
                              "sctC", "neutral", "805", "1.0", "1.000", "1.000", "1.000", "10", "20", ""])
        sol_tsv += "\n"
        sol_tsv += "\n"

        f_out = StringIO()
        track_multi_systems_hit = HitSystemTracker([system_1])
        likely_systems_to_tsv([system_1], track_multi_systems_hit, f_out)
        self.assertMultiLineEqual(sol_tsv, f_out.getvalue())

        f_out = StringIO()
        likely_systems_to_tsv([], track_multi_systems_hit, f_out)
        expected_out = f"""# macsyfinder {macsypy.__version__}
# {' '.join(sys.argv)}
# No Likely Systems found
"""
        self.assertEqual(expected_out, f_out.getvalue())


    def test_unnlikely_systems_to_txt(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'unordered'
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
        # test if id is well incremented
        gene_name = "gspD"
        c_gene_gspd = CoreGene(models_location, gene_name, profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)
        model.add_mandatory_gene(gene_gspd)
        gene_name = "sctJ"
        c_gene_sctj = CoreGene(models_location, gene_name, profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model)
        model.add_accessory_gene(gene_sctj)
        gene_name = "sctC"
        c_gene_sctc = CoreGene(models_location, gene_name, profile_factory)
        gene_sctc = ModelGene(c_gene_sctc, model)
        model.add_neutral_gene(gene_sctc)
        gene_name = "tadZ"
        c_gene_tadz = CoreGene(models_location, gene_name, profile_factory)
        gene_tadz = ModelGene(c_gene_tadz, model)
        model.add_forbidden_gene(gene_tadz)

        hit_1 = Hit(c_gene_gspd, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(c_gene_sctj, "hit_2", 804, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(c_gene_sctc, "hit_3", 805, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctc, GeneStatus.NEUTRAL)
        hit_4 = Hit(c_gene_tadz, "hit_4", 806, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_4 = ValidHit(hit_4, gene_tadz, GeneStatus.FORBIDDEN)
        reason = "why it not a system"
        system_1 = UnlikelySystem(model, [v_hit_1], [v_hit_2], [v_hit_3], [v_hit_4], reason)

        exp_txt = f"""# macsyfinder {macsypy.__version__}
# {' '.join(sys.argv)}
# Unlikely Systems found:

This replicon probably not contains a system foo/T2SS:
{reason}

system id = replicon_id_T2SS_1
model = foo/T2SS
replicon = replicon_id
hits = [('hit_1', 'gspD', 1), ('hit_2', 'sctJ', 1), ('hit_3', 'sctC', 1), ('hit_4', 'tadZ', 1)]
wholeness = 1.000

mandatory genes:
\t- gspD: 1 (gspD)

accessory genes:
\t- sctJ: 1 (sctJ)

neutral genes:
\t- sctC: 1 (sctC)

forbidden genes:
\t- tadZ: 1 (tadZ)

Use ordered replicon to have better prediction.

============================================================
"""

        f_out = StringIO()
        unlikely_systems_to_txt([system_1], f_out)
        self.assertMultiLineEqual(exp_txt, f_out.getvalue())

        f_out = StringIO()
        unlikely_systems_to_txt([], f_out)
        expected_out = f"""# macsyfinder {macsypy.__version__}
# {' '.join(sys.argv)}
# No Unlikely Systems found
"""
        self.assertEqual(expected_out, f_out.getvalue())


    def test_parse_args(self):
        command_line = "macsyfinder --sequence-db test_1.fasta --db-type=gembase --models-dir data/models/ " \
                       "--models functional all -w 4 --out test_1-all"
        parser, args = parse_args(command_line.split()[1:])
        self.assertIsNone(args.cfg_file)
        self.assertIsNone(args.coverage_profile)
        self.assertIsNone(args.hmmer)
        self.assertIsNone(args.i_evalue_sel)
        self.assertIsNone(args.inter_gene_max_space)
        self.assertIsNone(args.max_nb_genes)
        self.assertIsNone(args.min_genes_required)
        self.assertIsNone(args.min_mandatory_genes_required)
        self.assertIsNone(args.multi_loci)
        self.assertIsNone(args.previous_run)
        self.assertIsNone(args.profile_suffix)
        self.assertIsNone(args.replicon_topology)
        self.assertIsNone(args.res_extract_suffix)
        self.assertIsNone(args.res_search_suffix)
        self.assertIsNone(args.topology_file)
        self.assertFalse(args.idx)
        self.assertFalse(args.list_models)
        self.assertFalse(args.mute)
        self.assertFalse(args.relative_path)
        self.assertEqual(args.db_type, 'gembase')
        self.assertEqual(args.models_dir, 'data/models/')
        self.assertEqual(args.out_dir, 'test_1-all')
        self.assertEqual(args.sequence_db, 'test_1.fasta')
        self.assertEqual(args.verbosity, 0)
        self.assertEqual(args.worker, 4)

        self.assertListEqual(args.models, ['functional', 'all'])

        command_line = "macsyfinder --sequence-db test_!.fasta " \
                       "--db-type=ordered_replicon --models-dir data/models/ " \
                       "--models functional all -w 4 --out test_1-all " \
                       "--mute --multi-loci TXSscan/T2SS,TXSScan/T3SS --relative-path"
        parser, args = parse_args(command_line.split()[1:])
        self.assertEqual(args.db_type, 'ordered_replicon')
        self.assertEqual(args.multi_loci, "TXSscan/T2SS,TXSScan/T3SS")
        self.assertTrue(args.relative_path)
        self.assertTrue(args.mute)

        command_line = "macsyfinder --sequence-db test_1.dasta " \
                       "--db-type=ordered_replicon --models-dir data/models/ " \
                       "--i-evalue-sel=0.5 " \
                       "--min-genes-required TXSScan/T2SS 15 --min-genes-required TXSScan/Flagellum 10"
        parser, args = parse_args(command_line.split()[1:])
        self.assertEqual(args.i_evalue_sel, 0.5)
        self.assertListEqual(args.min_genes_required, [['TXSScan/T2SS', '15'], ['TXSScan/Flagellum', '10']])


    @unittest.skipIf(not which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_search_systems(self):
        logger = logging.getLogger('macsypy.macsyfinder')
        macsypy.logger_set_level(level='ERROR')
        defaults = MacsyDefaults()

        out_dir = os.path.join(self.tmp_dir, 'macsyfinder_test_search_systems')
        os.mkdir(out_dir)

        # test gembase replicon
        seq_db = self.find_data('base', 'VICH001.B.00001.C001.prt')
        model_dir = self.find_data('data_set', 'models')
        args = f"--sequence-db {seq_db} --db-type=gembase --models-dir {model_dir} --models set_1 all -w 4 -o {out_dir}"

        _, parsed_args = parse_args(args.split())
        config = Config(defaults, parsed_args)
        model_bank = ModelBank()
        gene_bank = GeneBank()
        profile_factory = ProfileFactory(config)

        systems, rejected_clst = search_systems(config, model_bank, gene_bank, profile_factory, logger)
        expected_sys_id = ['VICH001.B.00001.C001_MSH_5', 'VICH001.B.00001.C001_MSH_7',
                           'VICH001.B.00001.C001_T4P_25', 'VICH001.B.00001.C001_T4P_23',
                           'VICH001.B.00001.C001_T4P_21', 'VICH001.B.00001.C001_T4P_22',
                           'VICH001.B.00001.C001_T4P_17', 'VICH001.B.00001.C001_T4P_16',
                           'VICH001.B.00001.C001_T4bP_26', 'VICH001.B.00001.C001_T4P_24',
                           'VICH001.B.00001.C001_T4P_18', 'VICH001.B.00001.C001_T4P_19',
                           'VICH001.B.00001.C001_T4P_20',
                           'VICH001.B.00001.C001_T2SS_10', 'VICH001.B.00001.C001_T2SS_9'
                           ]
        self.assertListEqual([s.id for s in systems], expected_sys_id)

        expected_scores = [10.5, 10.0, 12.0, 9.5, 9.0, 8.5, 6.0, 5.0, 5.5, 10.5, 7.5, 7.0, 8.0, 8.3, 7.5]
        self.assertListEqual([s.score for s in systems], expected_scores)
        self.assertEqual(len(rejected_clst), 11)

        # test hits but No Systems
        args = f"--sequence-db {seq_db} --db-type=gembase --models-dir {model_dir} --models set_1 Tad -w 4 -o {out_dir}"
        _, parsed_args = parse_args(args.split())
        config = Config(defaults, parsed_args)
        model_bank = ModelBank()
        gene_bank = GeneBank()
        profile_factory = ProfileFactory(config)
        systems, rejected_clst = search_systems(config, model_bank, gene_bank, profile_factory, logger)
        self.assertEqual(systems, [])

        # test No hits
        seq_db = self.find_data('base', 'test_1.fasta')
        args = f"--sequence-db {seq_db} --db-type=gembase --models-dir {model_dir} --models set_1 T4bP -w 4 -o {out_dir}"
        _, parsed_args = parse_args(args.split())
        config = Config(defaults, parsed_args)
        model_bank = ModelBank()
        gene_bank = GeneBank()
        profile_factory = ProfileFactory(config)
        systems, rejected_clst = search_systems(config, model_bank, gene_bank, profile_factory, logger)
        self.assertEqual(systems, [])
        self.assertEqual(rejected_clst, [])


    def test_search_systems_unordered(self):
        logger = logging.getLogger('macsypy.macsyfinder')
        macsypy.logger_set_level(level='ERROR')
        defaults = MacsyDefaults()

        out_dir = os.path.join(self.tmp_dir, 'macsyfinder_test_search_systems')
        os.mkdir(out_dir)
        seq_db = self.find_data('base', 'VICH001.B.00001.C001.prt')
        model_dir = self.find_data('data_set', 'models')
        # test unordered replicon
        args = f"--sequence-db {seq_db} --db-type=unordered --models-dir {model_dir} --models set_1 all -w 4 -o {out_dir}"

        _, parsed_args = parse_args(args.split())
        config = Config(defaults, parsed_args)
        model_bank = ModelBank()
        gene_bank = GeneBank()
        profile_factory = ProfileFactory(config)

        systems, uncomplete_sys = search_systems(config, model_bank, gene_bank, profile_factory, logger)
        expected_sys_id = ['Unordered_T2SS_4', 'Unordered_MSH_3', 'Unordered_T4P_5', 'Unordered_T4bP_6']
        self.assertListEqual([s.id for s in systems], expected_sys_id)

        expected_uncomplete_sys_id = ['Unordered_Archaeal-T4P_1', 'Unordered_ComM_2', 'Unordered_Tad_7']
        self.assertListEqual([s.id for s in uncomplete_sys], expected_uncomplete_sys_id)


    def test_search_systems_model_unknown(self):
        logger = logging.getLogger('macsypy.macsyfinder')
        macsypy.logger_set_level(level='ERROR')
        defaults = MacsyDefaults()

        out_dir = os.path.join(self.tmp_dir, 'macsyfinder_test_search_systems')
        os.mkdir(out_dir)
        seq_db = self.find_data('base', 'test_1.fasta')
        model_dir = self.find_data('data_set', 'models')
        args = f"--sequence-db {seq_db} --db-type=gembase --models-dir {model_dir} --models nimporaoik -w 4 -o {out_dir}"

        _, parsed_args = parse_args(args.split())
        config = Config(defaults, parsed_args)
        model_bank = ModelBank()
        gene_bank = GeneBank()
        profile_factory = ProfileFactory(config)

        exit_ori = sys.exit
        sys.exit = self.fake_exit
        try:
            with self.assertRaises(TypeError) as ctx:
                _ = search_systems(config, model_bank, gene_bank, profile_factory, logger)
            self.assertEqual(str(ctx.exception),
                             "macsyfinder: \"No such model definition: 'nimporaoik'\"")
        finally:
            sys.exit = exit_ori
