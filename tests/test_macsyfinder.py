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
import sys
import shutil
import tempfile
import argparse
from io import StringIO
import logging
import unittest

from macsypy.config import Config, MacsyDefaults
from macsypy.gene import CoreGene, ModelGene, GeneStatus, GeneBank
from macsypy.profile import ProfileFactory
from macsypy.registries import ModelRegistry, scan_models_dir, ModelLocation
from macsypy.hit import Hit, ValidHit
from macsypy.model import Model, ModelBank
from macsypy.system import System, HitSystemTracker
from macsypy.cluster import Cluster, RejectedClusters
from macsypy.error import MacsypyError

from macsypy.scripts.macsyfinder import systems_to_txt, systems_to_tsv, rejected_clst_to_txt, parse_args, search_systems
from macsypy.scripts.macsyfinder import list_models

import macsypy
from tests import MacsyTest, which


class TestMacsyfinder(MacsyTest):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()


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
        args.sequence_db = self.find_data("base", "test_base.fa")
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
        system_1 = System(model, [Cluster([v_hit_1, v_hit_2], model)])

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
            args.sequence_db = self.find_data("base", "test_base.fa")
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
            system_1 = System(model, [Cluster([v_hit_1, v_hit_2], model)])

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


    def test_rejected_clst_to_txt(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
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
        c1 = Cluster([v_h10, v_h20], model)
        c2 = Cluster([v_h40, v_h50], model)
        r_c = RejectedClusters(model, [c1, c2], "The reasons to reject this clusters")

        rej_clst_str = f"""# macsyfinder {macsypy.__version__}
# {' '.join(sys.argv)}
# Rejected clusters:

Cluster:
    - model: T2SS
    - hits: (h10, gspD, 10), (h20, sctC, 20)
Cluster:
    - model: T2SS
    - hits: (h10, gspD, 40), (h20, sctC, 50)
These clusters has been rejected because:
The reasons to reject this clusters
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


    def test_parse_args(self):
        command_line = "macsyfinder --sequence-db VICH001.B.00001.C001.prt --db-type=gembase --models-dir data/models/ " \
                       "--models TFF-SF_final all -w 4 --out VICH001-all"
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
        self.assertEqual(args.out_dir, 'VICH001-all')
        self.assertEqual(args.sequence_db, 'VICH001.B.00001.C001.prt')
        self.assertEqual(args.verbosity, 0)
        self.assertEqual(args.worker, 4)

        self.assertListEqual(args.models, [['TFF-SF_final', 'all']])

        command_line = "macsyfinder --sequence-db VICH001.B.00001.C001.prt " \
                       "--db-type=ordered_replicon --models-dir data/models/ " \
                       "--models TFF-SF_final all -w 4 --out VICH001-all " \
                       "--mute --multi-loci TXSscan/T2SS,TXSScan/T3SS --relative-path"
        parser, args = parse_args(command_line.split()[1:])
        self.assertEqual(args.db_type, 'ordered_replicon')
        self.assertEqual(args.multi_loci, "TXSscan/T2SS,TXSScan/T3SS")
        self.assertTrue(args.relative_path)
        self.assertTrue(args.mute)

        command_line = "macsyfinder --sequence-db VICH001.B.00001.C001.prt " \
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
        seq_db = self.find_data('base', 'VICH001.B.00001.C001.prt')
        model_dir = self.find_data('data_set', 'models')
        args = f"--sequence-db {seq_db} --db-type=gembase --models-dir {model_dir} --models set_1 all -w 4 -o {out_dir}"

        _, parsed_args = parse_args(args.split())
        config = Config(defaults, parsed_args)
        model_bank = ModelBank()
        gene_bank = GeneBank()
        profile_factory = ProfileFactory(config)
        systems, rejected_clst = search_systems(config, model_bank, gene_bank, profile_factory, logger)
        expected_sys_id = ['VICH001.B.00001.C001_MSH_1', 'VICH001.B.00001.C001_MSH_2',
                           'VICH001.B.00001.C001_T4P_14', 'VICH001.B.00001.C001_T4P_12',
                           'VICH001.B.00001.C001_T4P_10', 'VICH001.B.00001.C001_T4P_11',
                           'VICH001.B.00001.C001_T4P_6', 'VICH001.B.00001.C001_T4P_5',
                           'VICH001.B.00001.C001_T4bP_15', 'VICH001.B.00001.C001_T4P_13',
                           'VICH001.B.00001.C001_T4P_7', 'VICH001.B.00001.C001_T4P_8',
                           'VICH001.B.00001.C001_T4P_9',
                           'VICH001.B.00001.C001_T2SS_4', 'VICH001.B.00001.C001_T2SS_3'
                           ]

        self.assertListEqual([s.id for s in systems], expected_sys_id)

        expected_scores = [10.5, 10.0, 12.0, 9.5, 9.0, 8.5, 6.0, 5.0, 5.5, 10.5, 7.5, 7.0, 8.0, 8.25, 7.5]
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
        seq_db = self.find_data('base', 'test_aesu.fa')
        args = f"--sequence-db {seq_db} --db-type=gembase --models-dir {model_dir} --models set_1 T4bP -w 4 -o {out_dir}"
        _, parsed_args = parse_args(args.split())
        config = Config(defaults, parsed_args)
        model_bank = ModelBank()
        gene_bank = GeneBank()
        profile_factory = ProfileFactory(config)
        systems, rejected_clst = search_systems(config, model_bank, gene_bank, profile_factory, logger)
        self.assertEqual(systems, [])
        self.assertEqual(rejected_clst, [])

    def test_search_systems_model_unknown(self):
        logger = logging.getLogger('macsypy.macsyfinder')
        macsypy.logger_set_level(level='ERROR')
        defaults = MacsyDefaults()

        out_dir = os.path.join(self.tmp_dir, 'macsyfinder_test_search_systems')
        os.mkdir(out_dir)
        seq_db = self.find_data('base', 'VICH001.B.00001.C001.prt')
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
