#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2025  Institut Pasteur (Paris) and CNRS.           #
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
import logging
import unittest
import platform
import itertools

from pandas import __version__ as pd_vers
import macsylib
from macsylib.config import Config, MacsyDefaults
from macsylib.registries import ModelRegistry, scan_models_dir
from macsylib.system import System, RejectedCandidate, AbstractUnordered
from macsylib.utils import get_def_to_detect

from macsyfinder.scripts.msf import list_models, parse_args, search_systems, get_version_message
from macsyfinder import __version__ as msf_vers, __commit__ as msf_commit, __citation__ as msf_citation
from tests import MacsyTest


class TestMacsyfinder(MacsyTest):

    def setUp(self):
        self._tmp_dir = tempfile.TemporaryDirectory(prefix='test_msf_macsyfinder_')
        self.tmp_dir = self._tmp_dir.name
        self._reset_id()
        AbstractUnordered._id = itertools.count(1)


    def tearDown(self):
        self._tmp_dir.cleanup()


    def _fill_model_registry(self, config):
        model_registry = ModelRegistry()

        for model_dir in config.models_dir():
            models_loc_available = scan_models_dir(model_dir,
                                                   profile_suffix=config.profile_suffix(),
                                                   relative_path=config.relative_path())
            for model_loc in models_loc_available:
                model_registry.add(model_loc)
        return model_registry


    def _reset_id(self):
        """
        reset System._id and RejectedCluster._id to get predictable ids
        """
        System._id = itertools.count(1)
        RejectedCandidate._id = itertools.count(1)


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


    @unittest.skipIf(platform.system() == 'Windows' or os.getuid() == 0, 'Skip test on Windows or if run as root')
    def test_list_models_no_permissions(self):
        # on gitlab it is not allowed to change the permission of a directory
        # located in tests/data
        # So I need to copy it in /tmp
        tmp_dir = tempfile.TemporaryDirectory(prefix='test_msl_Config_')
        model_dir_name = 'fake_model_dir'
        src_model_dir = self.find_data(model_dir_name)
        dst_model_dir = os.path.join(tmp_dir.name, 'fake_model_dir')
        shutil.copytree(src_model_dir, dst_model_dir)

        log = logging.getLogger('macsylib')
        log.setLevel(logging.WARNING)
        cmd_args = argparse.Namespace()
        cmd_args.models_dir = dst_model_dir
        cmd_args.list_models = True
        models_dir_perm = os.stat(cmd_args.models_dir).st_mode
        try:
            os.chmod(cmd_args.models_dir, 0o110)
            with self.catch_log(log_name='macsylib') as log:
                rcv_list_models = list_models(cmd_args)
                log_msg = log.get_value().strip()
            self.assertEqual(rcv_list_models, '')
            self.assertEqual(log_msg, f"{cmd_args.models_dir} is not readable: [Errno 13] Permission denied: '{cmd_args.models_dir}' : skip it.")
        finally:
            os.chmod(cmd_args.models_dir, models_dir_perm)
            tmp_dir.cleanup()


    def test_get_version_message(self):
        py_vers = sys.version.replace('\n', ' ')
        exp_vers_msg = f"""MacSyFinder {msf_vers} {msf_commit}
using:
- Python {py_vers}
- MacSyLib {macsylib.__version__} {macsylib.__commit__}
- NetworkX {macsylib.solution.nx.__version__}
- Pandas {pd_vers}

MacSyFinder is distributed under the terms of the GNU General Public License (GPLv3).
See the COPYING file for details.

If you use this software please cite:
{msf_citation}

and don't forget to cite the used models:
msf_data cite <model>
"""
        self.maxDiff = None
        self.assertEqual(exp_vers_msg, get_version_message())


    def test_parse_args(self):
        command_line = "macsyfinder --sequence-db test_1.fasta --db-type=gembase --models-dir data/models/ " \
                       "--models functional all -w 4 --out-dir test_1-all"
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
        self.assertIsNone(args.index_dir)
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

        command_line = "macsyfinder --sequence-db test_1.fasta " \
                       "--db-type=ordered_replicon --models-dir data/models/ " \
                       "--models functional all -w 4 --out-dir test_1-all " \
                       "--mute --multi-loci TXSscan/T2SS,TXSScan/T3SS --relative-path --index-dir the_idx_dir"
        parser, args = parse_args(command_line.split()[1:])
        self.assertEqual(args.db_type, 'ordered_replicon')
        self.assertEqual(args.index_dir, 'the_idx_dir')
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


    @unittest.skipIf(not shutil.which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_search_systems(self):
        logger = logging.getLogger('macsylib.macsyfinder')
        macsylib.logger_set_level(level='ERROR')
        defaults = MacsyDefaults()

        out_dir = os.path.join(self.tmp_dir, 'macsyfinder_test_search_systems')
        os.mkdir(out_dir)

        # test gembase replicon
        seq_db = self.find_data('base', 'VICH001.B.00001.C001.prt')
        model_dir = self.find_data('data_set', 'models')
        args = f"--sequence-db {seq_db} --db-type=gembase --models-dir {model_dir} --models set_1 all -w 4" \
               f" -o {out_dir} --index-dir {out_dir}"

        _, parsed_args = parse_args(args.split())
        config = Config(defaults, parsed_args)
        model_registry = self._fill_model_registry(config)
        def_to_detect, models_fam_name, models_version = get_def_to_detect(config.models(), model_registry)

        self._reset_id()
        systems, rejected_clst = search_systems(config, model_registry, def_to_detect, logger)
        expected_sys_id = ['VICH001.B.00001.C001_MSH_1',
                           'VICH001.B.00001.C001_T4P_13', 'VICH001.B.00001.C001_T4P_11', 'VICH001.B.00001.C001_T4P_9',
                           'VICH001.B.00001.C001_T4P_10', 'VICH001.B.00001.C001_T4P_5', 'VICH001.B.00001.C001_T4P_4',
                           'VICH001.B.00001.C001_T4bP_14', 'VICH001.B.00001.C001_T4P_12', 'VICH001.B.00001.C001_T4P_6',
                           'VICH001.B.00001.C001_T4P_7', 'VICH001.B.00001.C001_T4P_8',
                           'VICH001.B.00001.C001_T2SS_3', 'VICH001.B.00001.C001_T2SS_2']

        self.assertListEqual([s.id for s in systems], expected_sys_id)

        expected_scores = [10.5, 12.0, 9.5, 9.0, 8.5, 6.0, 5.0, 5.5, 10.5, 7.5, 7.0, 8.0, 8.06, 7.5]
        self.assertListEqual([s.score for s in systems], expected_scores)
        self.assertEqual(len(rejected_clst), 10)

        # test hits but No Systems
        args = f"--sequence-db {seq_db} --db-type=gembase --models-dir {model_dir} --models set_1 Tad -w 4" \
               f" -o {out_dir} --index-dir {out_dir}"
        _, parsed_args = parse_args(args.split())
        config = Config(defaults, parsed_args)
        model_registry = self._fill_model_registry(config)
        def_to_detect, models_fam_name, models_version = get_def_to_detect(config.models(), model_registry)
        self._reset_id()
        systems, rejected_clst = search_systems(config, model_registry, def_to_detect, logger)
        self.assertEqual(systems, [])

        # test No hits
        seq_db = self.find_data('base', 'test_1.fasta')
        args = f"--sequence-db {seq_db} --db-type=gembase --models-dir {model_dir} --models set_1 T4bP -w 4" \
               f" -o {out_dir} --index-dir {out_dir}"
        _, parsed_args = parse_args(args.split())
        config = Config(defaults, parsed_args)
        model_registry = self._fill_model_registry(config)
        def_to_detect, models_fam_name, models_version = get_def_to_detect(config.models(), model_registry)
        self._reset_id()
        systems, rejected_clst = search_systems(config, model_registry, def_to_detect, logger)
        self.assertEqual(systems, [])
        self.assertEqual(rejected_clst, [])

        # test multisystems
        # multisytem hit are not in System (to small cluster)
        # no system
        seq_db = self.find_data('base', 'test_12.fasta')
        model_dir = self.find_data('models')
        args = f"--sequence-db {seq_db} --db-type=gembase --models-dir {model_dir} " \
               f"--models functional T12SS-multisystem -w 4 -o {out_dir} --index-dir {out_dir}"
        _, parsed_args = parse_args(args.split())
        config = Config(defaults, parsed_args)
        model_registry = self._fill_model_registry(config)
        def_to_detect, models_fam_name, models_version = get_def_to_detect(config.models(), model_registry)
        self._reset_id()
        systems, rejected_clst = search_systems(config, model_registry, def_to_detect, logger)

        self.assertEqual(systems, [])
        self.assertEqual([r.id for r in rejected_clst],
                         ['VICH001.B.00001.C001_T12SS-multisystem_1', 'VICH001.B.00001.C001_T12SS-multisystem_2'])

        # multisystem is in System, so it can play role for other cluster
        # 2 systems found
        seq_db = self.find_data('base', 'test_13.fasta')
        model_dir = self.find_data('models')
        args = f"--sequence-db {seq_db} --db-type=gembase --models-dir {model_dir} " \
               f"--models functional T12SS-multisystem -w 4 -o {out_dir} --index-dir {out_dir}"
        _, parsed_args = parse_args(args.split())
        config = Config(defaults, parsed_args)
        model_registry = self._fill_model_registry(config)
        def_to_detect, models_fam_name, models_version = get_def_to_detect(config.models(), model_registry)
        self._reset_id()
        systems, rejected_clst = search_systems(config, model_registry, def_to_detect, logger)
        self.assertEqual({s.id for s in systems},
                         {'VICH001.B.00001.C001_T12SS-multisystem_3',
                          'VICH001.B.00001.C001_T12SS-multisystem_2',
                          'VICH001.B.00001.C001_T12SS-multisystem_1'})
        self.assertEqual([r.id for r in rejected_clst],
                         ['VICH001.B.00001.C001_T12SS-multisystem_1'])

    @unittest.skipIf(not shutil.which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_db_type_set_to_gembase(self):
        logger = logging.getLogger('macsylib.macsyfinder')
        macsylib.logger_set_level(level='WARNING')
        defaults = MacsyDefaults()

        out_dir = os.path.join(self.tmp_dir, 'macsyfinder_test_search_systems')
        os.mkdir(out_dir)

        # test gembase replicon
        seq_db = self.find_data('base', 'ordered_replicon_base.fasta')
        model_dir = self.find_data('data_set', 'models')
        args = f"--sequence-db {seq_db} --db-type=gembase --models-dir {model_dir} --models set_1 T4P -w 4" \
               f" -o {out_dir} --index-dir {out_dir}"

        _, parsed_args = parse_args(args.split())
        config = Config(defaults, parsed_args)
        model_registry = self._fill_model_registry(config)
        def_to_detect, models_fam_name, models_version = get_def_to_detect(config.models(), model_registry)

        self._reset_id()
        with self.catch_log() as log:
            systems, rejected_clst = search_systems(config, model_registry, def_to_detect, logger)
            log_msg = log.get_value().split('\n')[-2] # the message finish with empty line
        self.assertEqual(log_msg,
                         f"Most of replicons contains only ONE sequence are you sure that '{seq_db}' is a 'gembase'.")


    @unittest.skipIf(not shutil.which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_search_systems_unordered(self):
        logger = logging.getLogger('macsylib.macsyfinder')
        macsylib.logger_set_level(level='ERROR')
        defaults = MacsyDefaults()

        out_dir = os.path.join(self.tmp_dir, 'macsyfinder_test_search_systems')
        os.mkdir(out_dir)
        seq_db = self.find_data('base', 'VICH001.B.00001.C001.prt')
        model_dir = self.find_data('data_set', 'models')
        # test unordered replicon
        args = f"--sequence-db {seq_db} --db-type=unordered --models-dir {model_dir} --models set_1 all -w 4" \
               f" -o {out_dir} --index-dir {out_dir}"

        _, parsed_args = parse_args(args.split())
        config = Config(defaults, parsed_args)

        model_registry = self._fill_model_registry(config)
        def_to_detect, models_fam_name, models_version = get_def_to_detect(config.models(), model_registry)
        systems, uncomplete_sys = search_systems(config, model_registry, def_to_detect, logger)
        expected_sys_id = ['VICH001.B.00001.C001_T2SS_4', 'VICH001.B.00001.C001_MSH_3',
                           'VICH001.B.00001.C001_T4P_5', 'VICH001.B.00001.C001_T4bP_6']
        self.assertListEqual([s.id for s in systems], expected_sys_id)

        expected_uncomplete_sys_id = ['VICH001.B.00001.C001_Archaeal-T4P_1', 'VICH001.B.00001.C001_ComM_2',
                                      'VICH001.B.00001.C001_Tad_7']
        self.assertListEqual([s.id for s in uncomplete_sys], expected_uncomplete_sys_id)
