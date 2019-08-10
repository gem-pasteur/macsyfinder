# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Neron                                         #
# Copyright (c) 2014-2019  Institut Pasteur (Paris) and CNRS.                  #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################

import os
import sys
import shutil
import tempfile
import argparse
from io import StringIO

from macsypy.config import Config, MacsyDefaults
from macsypy.gene import ProfileFactory, Gene, GeneStatus
from macsypy.registries import ModelRegistry, scan_models_dir, ModelLocation
from macsypy.hit import Hit, ValidHit
from macsypy.model import Model
from macsypy.system import System, HitSystemTracker
from macsypy.cluster import Cluster, RejectedClusters
from macsypy.scripts.macsyfinder import systems_to_file, rejected_clst_to_file, parse_args
import macsypy
from tests import MacsyTest


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
        cmd_args.models_dir = os.path.join(self._data_dir, 'data_set_1', 'models')
        cmd_args.list_models = True
        registry = ModelRegistry()
        models_location = scan_models_dir(cmd_args.models_dir)
        for ml in models_location:
            registry.add(ml)
        list_models = """set_1
      /CONJ
      /Flagellum
      /T2SS
      /T3SS
      /T4P
      /T4SS_typeF
      /T4SS_typeI
      /T9SS
      /Tad
set_2
      /CONJ
      /Flagellum
      /T2SS
      /T3SS
      /T4P
      /T4SS_typeF
      /T4SS_typeI
      /T9SS
      /Tad
"""
        self.assertEqual(str(registry), list_models)

    def test_systems_to_file(self):
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
        gene_gspd = Gene(profile_factory, "gspD", model, models_location)
        model.add_mandatory_gene(gene_gspd)
        gene_sctj = Gene(profile_factory, "sctJ", model, models_location)
        model.add_accessory_gene(gene_sctj)

        hit_1 = Hit(gene_gspd, model, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(gene_sctj, model, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        system_1 = System(model, [Cluster([v_hit_1, v_hit_2], model)])

        system_str = """# macsyfinder {}
# {}
# Systems found:

system id = replicon_id_T2SS_{}
model = foo/T2SS
replicon = replicon_id
clusters = [('gspD', 1), ('sctJ', 1)]
occ = 1
wholeness = 1.000
loci nb = 1
score = 1.500

mandatory genes:
\t- gspD: 1 (gspD)

accessory genes:
\t- sctJ: 1 (sctJ)

============================================================
""".format(macsypy.__version__, ' '.join(sys.argv), next(System._id) - 1)

        f_out = StringIO()
        track_multi_systems_hit = HitSystemTracker([system_1])
        systems_to_file([system_1], track_multi_systems_hit, f_out)
        self.assertMultiLineEqual(system_str, f_out.getvalue())

    def test_rejected_clst_to_file(self):
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

        gene_1 = Gene(profile_factory, "gspD", model, models_location)
        gene_2 = Gene(profile_factory, "sctC", model, models_location)

        #     Hit(gene, model, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #         profile_coverage, sequence_coverage, begin_match, end_match
        h10 = Hit(gene_1, model, "h10", 10, "replicon_1", 10, 1.0, 10.0, 1.0, 1.0, 10, 20)
        h20 = Hit(gene_2, model, "h20", 10, "replicon_1", 20, 1.0, 20.0, 1.0, 1.0, 10, 20)
        h40 = Hit(gene_1, model, "h10", 10, "replicon_1", 40, 1.0, 10.0, 1.0, 1.0, 10, 20)
        h50 = Hit(gene_2, model, "h20", 10, "replicon_1", 50, 1.0, 20.0, 1.0, 1.0, 10, 20)
        c1 = Cluster([h10, h20], model)
        c2 = Cluster([h40, h50], model)
        r_c = RejectedClusters(model, [c1, c2], "The reasons to reject this clusters")

        rej_clst_str = """# macsyfinder {}
# {}
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
""".format(macsypy.__version__, ' '.join(sys.argv))

        f_out = StringIO()
        rejected_clst_to_file([r_c], f_out)
        self.maxDiff = None
        self.assertMultiLineEqual(rej_clst_str, f_out.getvalue())


    def test_parse_args(self):
        command_line = "macsyfinder --sequence-db VICH001.B.00001.C001.prt --db-type=gembase --models-dir data/models/ " \
                       "--models TFF-SF_final all -w 4 --out VICH001-all"
        args = parse_args(command_line.split()[1:])
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
        args = parse_args(command_line.split()[1:])
        self.assertEqual(args.db_type, 'ordered_replicon')
        self.assertEqual(args.multi_loci, "TXSscan/T2SS,TXSScan/T3SS")
        self.assertTrue(args.relative_path)
        self.assertTrue(args.mute)

        command_line = "macsyfinder --sequence-db VICH001.B.00001.C001.prt " \
                       "--db-type=ordered_replicon --models-dir data/models/ " \
                       "--i-evalue-sel=0.5 " \
                       "--min-genes-required TXSScan/T2SS 15 --min-genes-required TXSScan/Flagellum 10"
        args = parse_args(command_line.split()[1:])
        self.assertEqual(args.i_evalue_sel, 0.5)
        self.assertListEqual(args.min_genes_required, [['TXSScan/T2SS', '15'], ['TXSScan/Flagellum', '10']])
