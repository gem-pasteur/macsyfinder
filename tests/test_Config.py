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
import shutil
import unittest
from argparse import Namespace
from configparser import ConfigParser, ParsingError
import tempfile
from time import strftime

from macsypy.config import MacsyDefaults, Config, NoneConfig

from tests import MacsyTest


class TestConfig(MacsyTest):

    def setUp(self):
        self._current_dir = os.getcwd()
        self.tmp_dir = os.path.join(tempfile.gettempdir(),
                                    'test_macsyfinder_Config')
        if os.path.exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
        os.mkdir(self.tmp_dir)
        self.defaults = MacsyDefaults()
        self.parsed_args = Namespace()


    def tearDown(self):
        os.chdir(self._current_dir)
        try:
            shutil.rmtree(self.tmp_dir)
            #pass
        except Exception as err:
            pass

    def test_str_2_tuple(self):
        s = 'set_1/Flagellum 12 set_1/t4ss 13'
        expected = [('set_1/Flagellum', '12'), ('set_1/t4ss', '13')]
        cfg = Config(self.defaults, self.parsed_args)
        self.assertListEqual(cfg._str_2_tuple(s), expected)

        with self.assertRaises(ValueError) as ctx:
            s = 'set_1/Flagellum 12 set_1/t4ss'
            cfg._str_2_tuple(s)
        self.assertEqual(str(ctx.exception),
                         f"You must provide a list of model name and value separated by spaces: {s}")

    def test_set_option(self):
        cfg = Config(self.defaults, self.parsed_args)
        opts = {'GOOD': "GOOD", "BAD": None}
        self.assertFalse('GOOD' in cfg._options)
        self.assertFalse('BAD' in cfg._options)
        cfg._set_options(opts)
        self.assertTrue('GOOD' in cfg._options)
        self.assertEqual(cfg._options['GOOD'], 'GOOD')
        self.assertFalse('BAD' in cfg._options)


    def test_set_log_level(self):
        cfg = Config(self.defaults, self.parsed_args)
        cfg._set_log_level(20)
        self.assertEqual(cfg.log_level(), 20)
        cfg._set_log_level('WARNING')
        self.assertEqual(cfg.log_level(), 30)
        with self.assertRaises(ValueError) as ctx:
            cfg._set_log_level('FOO')
        self.assertEqual(str(ctx.exception),
                         "Invalid value for log_level: FOO.")



    def test_config_file_2_dict(self):
        cfg = Config(self.defaults, self.parsed_args)
        res = cfg._config_file_2_dict('nimportnaoik')
        self.assertDictEqual({}, res)

        cfg_file = self.find_data(os.path.join('conf_files', 'macsy_test_conf.conf'))
        res = cfg._config_file_2_dict(cfg_file)
        expected = {'db_type': 'gembase',
                    'inter_gene_max_space': 'set_1/T2SS 2 set_1/Flagellum 4',
                    'min_mandatory_genes_required': 'set_1/T2SS 5 set_1/Flagellum 9',
                    'replicon_topology': 'circular',
                    'sequence_db': '/path/to/sequence/bank/fasta_file',
                    'topology_file': '/the/path/to/the/topology/to/use'}
        self.assertDictEqual(expected, res)

        bad_cfg_file = self.find_data(os.path.join('conf_files', 'macsy_test_bad_conf.conf'))
        with self.assertRaises(ParsingError):
            cfg._config_file_2_dict(bad_cfg_file)


    def test_Config(self):
        cfg = Config(self.defaults, self.parsed_args)
        methods_needing_args = {'inter_gene_max_space': None,
                                'max_nb_genes': None,
                                'min_genes_required': None,
                                'min_mandatory_genes_required': None,
                                'multi_loci': None
                                }

        for opt, val in self.defaults.items():
            if opt == 'out_dir':
                self.assertEqual(cfg.out_dir(),
                                 os.path.join(cfg.res_search_dir(), f'macsyfinder-{strftime("%Y%m%d_%H-%M-%S")}')
                                 )
            elif opt == 'multi_loci':
                self.assertFalse(cfg.multi_loci('whatever'))
            elif opt in methods_needing_args:
                self.assertEqual(getattr(cfg, opt)('whatever'), val,
                                 msg=f"test of '{opt}' failed : expected {getattr(cfg, opt)('whatever')} !=  got {val}")
            elif opt == 'models_dir':
                val = self.defaults['system_models_dir']
                self.assertEqual(getattr(cfg, opt)(), val,
                                 msg=f"test of '{opt}' failed : expected {getattr(cfg, opt)()} !=  got {val}")
            else:
                self.assertEqual(getattr(cfg, opt)(), val,
                                 msg=f"test of '{opt}' failed : expected {getattr(cfg, opt)()} !=  got {val}")


    def test_cmd_config_file(self):
        methods_needing_args = {'inter_gene_max_space': [('set_1/Flagellum', 4), ('set_1/T2SS', 2)],
                                'max_nb_genes':  [('set_1/Flagellum', 6), ('set_1/T3SS', 3)],
                                'min_genes_required': [('set_1/Flagellum', 8), ('set_1/T4SS', 4)],
                                'min_mandatory_genes_required': [('set_1/Flagellum', 12), ('set_1/T6SS', 6)],
                                'multi_loci': {'set_1/Flagellum', 'T4SS'}
                                }

        self.parsed_args.cfg_file = self.find_data(os.path.join('conf_files', 'macsy_models.conf'))
        cfg = Config(self.defaults, self.parsed_args)

        expected_values = {k: v for k, v in self.defaults.items()}
        expected_values['cfg_file'] = self.parsed_args.cfg_file
        expected_values.update(methods_needing_args)

        for opt, val in expected_values.items():
            if opt == 'out_dir':
                self.assertEqual(cfg.out_dir(),
                                 os.path.join(cfg.res_search_dir(),
                                              f'macsyfinder-{strftime("%Y%m%d_%H-%M-%S")}')
                                 )
            elif opt == 'multi_loci':
                self.assertTrue(cfg.multi_loci('set_1/Flagellum'))
                self.assertTrue(cfg.multi_loci('set_1/T4SS'))
                self.assertFalse(cfg.multi_loci('set_1/T6SS'))
            elif opt in methods_needing_args:
                for model, genes in expected_values[opt]:
                    self.assertEqual(getattr(cfg, opt)(model), genes)
            elif opt == 'models_dir':
                self.assertEqual(getattr(cfg, opt)(), self.defaults['system_models_dir'])
            else:
                self.assertEqual(getattr(cfg, opt)(), val, msg=f"error with {opt}")

        self.parsed_args.cfg_file = 'niportnaoik'
        with self.assertRaises(ValueError) as ctx:
            Config(self.defaults, self.parsed_args)
        self.assertEqual(str(ctx.exception),
                         "Config file niportnaoik not found.")


    def test_project_config_file(self):
        os.chdir(self.tmp_dir)
        methods_needing_args = {'inter_gene_max_space': [('set_1/Flagellum', 4), ('set_1/T2SS', 2)],
                                'max_nb_genes':  [('set_1/Flagellum', 6), ('set_1/T3SS', 3)],
                                'min_genes_required': [('set_1/Flagellum', 8), ('set_1/T4SS', 4)],
                                'min_mandatory_genes_required': [('set_1/Flagellum', 12), ('set_1/T6SS', 6)],
                                'multi_loci': {'set_1/Flagellum', 'T4SS'}
                                }
        hmmer_opts_in_file = {'coverage_profile': 0.8,
                              'e_value_search': 0.12}

        try:
            shutil.copyfile(self.find_data(os.path.join('conf_files', 'project.conf')),
                            os.path.join(self.tmp_dir, 'macsyfinder.conf')
                            )
            cfg = Config(self.defaults, self.parsed_args)

            expected_values = {k: v for k, v in self.defaults.items()}
            expected_values.update(methods_needing_args)
            expected_values.update(hmmer_opts_in_file)

            for opt, val in expected_values.items():
                if opt == 'out_dir':
                    self.assertEqual(cfg.out_dir(),
                                     os.path.join(cfg.res_search_dir(),
                                                  f'macsyfinder-{strftime("%Y%m%d_%H-%M-%S")}')
                                     )
                elif opt == 'multi_loci':
                    self.assertTrue(cfg.multi_loci('set_1/Flagellum'))
                    self.assertTrue(cfg.multi_loci('set_1/T4SS'))
                    self.assertFalse(cfg.multi_loci('set_1/T6SS'))
                elif opt in methods_needing_args:
                    for model, genes in expected_values[opt]:
                        self.assertEqual(getattr(cfg, opt)(model), genes)
                elif opt == 'models_dir':
                    self.assertEqual(getattr(cfg, opt)(), self.defaults['system_models_dir'])
                else:
                    self.assertEqual(getattr(cfg, opt)(), val)
        except:
            os.chdir(self._current_dir)


    def test_Config_file_bad_values(self):
        ori_conf_file = self.find_data(os.path.join('conf_files', 'macsy_models.conf'))
        config_parser = ConfigParser()
        config_parser.read(ori_conf_file)
        config_parser.add_section('general')
        config_parser.set('general', 'worker', 'foo')
        with tempfile.TemporaryDirectory() as tmpdirname:
            dest_conf_file = os.path.join(tmpdirname, 'macsyfinder.conf')
            with open(dest_conf_file, 'w') as cfg_file:
                config_parser.write(cfg_file)
            self.parsed_args.cfg_file = dest_conf_file
            with self.assertRaises(ValueError) as ctx:
                Config(self.defaults, self.parsed_args)
            self.assertEqual(str(ctx.exception),
                             "Invalid value in config_file for option 'worker': "
                             "invalid literal for int() with base 10: 'foo'")


    def test_Config_default_conf_file(self):
        methods_needing_args = {'inter_gene_max_space': [('set_1/Flagellum', 4), ('set_1/T2SS', 2)],
                                'max_nb_genes':  [('set_1/Flagellum', 6), ('set_1/T3SS', 3)],
                                'min_genes_required': [('set_1/Flagellum', 8), ('set_1/T4SS', 4)],
                                'min_mandatory_genes_required': [('set_1/Flagellum', 12), ('set_1/T6SS', 6)],
                                'multi_loci': {'set_1/Flagellum', 'set_1/T4SS'}
                                }
        with tempfile.TemporaryDirectory() as tmpdirname:
            ori_conf_file = self.find_data(os.path.join('conf_files', 'macsy_models.conf'))
            dest_conf_file = os.path.join(tmpdirname, 'macsyfinder.conf')
            shutil.copy(ori_conf_file, dest_conf_file)
            os.environ['MACSY_CONF'] = dest_conf_file
            virtual_env = os.environ.get("VIRTUAL_ENV")
            del os.environ["VIRTUAL_ENV"]
            try:
                cfg = Config(self.defaults, self.parsed_args)

                expected_values = {k: v for k, v in self.defaults.items()}
                expected_values.update(methods_needing_args)
                for opt, val in expected_values.items():
                    if opt == 'out_dir':
                        self.assertEqual(cfg.out_dir(),
                                         os.path.join(cfg.res_search_dir(),
                                                      f'macsyfinder-{strftime("%Y%m%d_%H-%M-%S")}')
                                         )
                    elif opt == 'multi_loci':
                        self.assertTrue(cfg.multi_loci('set_1/Flagellum'))
                        self.assertTrue(cfg.multi_loci('set_1/T4SS'))
                        self.assertFalse(cfg.multi_loci('set_1/T6SS'))
                    elif opt in methods_needing_args:
                        for model, genes in expected_values[opt]:
                            self.assertEqual(getattr(cfg, opt)(model), genes)
                    elif opt == 'models_dir':
                        self.assertEqual(getattr(cfg, opt)(), self.defaults['system_models_dir'])
                    else:
                        self.assertEqual(getattr(cfg, opt)(), val)
            finally:
                os.environ["VIRTUAL_ENV"] = virtual_env

    def test_Config_conf_file_virtualenv(self):
        methods_needing_args = {'inter_gene_max_space': [('set_1/Flagellum', 4), ('set_1/T2SS', 2)],
                                'min_mandatory_genes_required': [('set_1/Flagellum', 9), ('set_1/T2SS', 5)]
                               }
        modified_args = {'db_type': 'gembase',
                         'coverage_profile': 0.1,
                         'replicon_topology': 'circular'
                        }

        with tempfile.TemporaryDirectory() as tmpdirname:
            ori_conf_file = self.find_data(os.path.join('conf_files', 'macsy_virtualenv_test.conf'))
            conf_dir = os.path.join(tmpdirname, 'etc', 'macsyfinder')
            os.makedirs(conf_dir)
            dest_conf_file = os.path.join(conf_dir, 'macsyfinder.conf')
            shutil.copy(ori_conf_file, dest_conf_file)
            virtual_env = os.environ.get("VIRTUAL_ENV")

            os.environ['VIRTUAL_ENV'] = tmpdirname
            try:
                cfg = Config(self.defaults, self.parsed_args)

                expected_values = {k: v for k, v in self.defaults.items()}
                expected_values.update(methods_needing_args)
                expected_values.update(modified_args)
                for opt, val in expected_values.items():
                    if opt == 'out_dir':
                        self.assertEqual(cfg.out_dir(),
                                         os.path.join(cfg.res_search_dir(),
                                                      f'macsyfinder-{strftime("%Y%m%d_%H-%M-%S")}')
                                         )
                    elif opt in ('max_nb_genes', 'min_genes_required', 'multi_loci'):  # not set in cfg file
                        pass
                    elif opt in methods_needing_args:
                        for model, genes in expected_values[opt]:
                            self.assertEqual(getattr(cfg, opt)(model), genes)
                    elif opt == 'models_dir':
                        self.assertEqual(getattr(cfg, opt)(), self.defaults['system_models_dir'])
                    else:
                        self.assertEqual(getattr(cfg, opt)(), val)
            finally:
                if virtual_env:
                    os.environ["VIRTUAL_ENV"] = virtual_env


    def test_Config_args(self):
        methods_needing_args = {'inter_gene_max_space': [('set_1/Flagellum', '14'), ('set_1/T2SS', '12')],
                                'max_nb_genes': [('set_1/Flagellum', '16'), ('set_1/T3SS', '13')],
                                'min_genes_required': [('set_1/Flagellum', '18'), ('set_1/T4SS', '14')],
                                'min_mandatory_genes_required': [('set_1/Flagellum', '22'), ('set_1/T6SS', '16')],
                                'multi_loci': 'set_1/Flagellum, set_1/T4SS',
                                }
        for opt, value in methods_needing_args.items():
            setattr(self.parsed_args, opt, value)

        simple_opt = {'hmmer': 'foo',
                      'i_evalue_sel': 20,
                      'replicon_topology': 'linear',
                      'db_type': 'gembase',
                      'sequence_db': self.find_data(os.path.join('base', 'test_1.fasta')),
                      'topology_file': __file__  # test only the existence of a file
                      }

        for opt, val in simple_opt.items():
            setattr(self.parsed_args, opt, val)

        cfg = Config(self.defaults, self.parsed_args)

        expected_values = {k: v for k, v in self.defaults.items()}
        expected_values.update(methods_needing_args)
        expected_values.update(simple_opt)
        for opt, val in expected_values.items():
            if opt == 'out_dir':
                self.assertEqual(cfg.out_dir(),
                                 os.path.join(cfg.res_search_dir(), f'macsyfinder-{strftime("%Y%m%d_%H-%M-%S")}')
                                 )
            elif opt == 'multi_loci':
                self.assertTrue(cfg.multi_loci('set_1/Flagellum'))
                self.assertTrue(cfg.multi_loci('set_1/T4SS'))
                self.assertFalse(cfg.multi_loci('set_1/T6SS'))
            elif opt in methods_needing_args:
                for model, genes in expected_values[opt]:
                    self.assertEqual(getattr(cfg, opt)(model), int(genes))
            elif opt == 'models_dir':
                self.assertEqual(getattr(cfg, opt)(), self.defaults['system_models_dir'])
            else:
                self.assertEqual(getattr(cfg, opt)(), val,
                                 msg=f"{opt} failed: expected: val '{val}' != got '{getattr(cfg, opt)()}'")

    def test_Config_file_n_args(self):
        cfg_needing_args = {'inter_gene_max_space': [('set_1/Flagellum', '4'), ('set_1/T2SS', '2')],
                                'max_nb_genes': [('set_1/Flagellum', '6'), ('set_1/T3SS', '3')],
                                'min_genes_required': [('set_1/Flagellum', '8'), ('set_1/T4SS', '4')],
                                'min_mandatory_genes_required': [('set_1/Flagellum', '12'), ('set_1/T6SS', '6')],
                                'multi_loci': 'set_1/Flagellum, set_1/T4SS',
                                }

        self.parsed_args.cfg_file = self.find_data(os.path.join('conf_files', 'macsy_models.conf'))
        expected_values = {k: v for k, v in self.defaults.items()}
        expected_values['cfg_file'] = self.parsed_args.cfg_file
        expected_values.update(cfg_needing_args)

        cmd_needing_args = {'min_genes_required': [('set_1/Flagellum', 18), ('T4SS', 14)],
                            'min_mandatory_genes_required': [('set_1/Flagellum', 22), ('set_1/T6SS', 16)],
                            }
        for opt, value in cmd_needing_args.items():
            setattr(self.parsed_args, opt, ' '.join([f"{m} {v}" for m, v in value]))

        simple_opt = {'hmmer': 'foo',
                      'i_evalue_sel': 20,
                      'db_type': 'gembase'}
        for opt, val in simple_opt.items():
            setattr(self.parsed_args, opt, val)

        cfg = Config(self.defaults, self.parsed_args)

        expected_values.update(cmd_needing_args)
        expected_values.update(simple_opt)

        for opt, exp_val in expected_values.items():
            if opt == 'out_dir':
                self.assertEqual(cfg.out_dir(),
                                 os.path.join(cfg.res_search_dir(), f"macsyfinder-{strftime('%Y%m%d_%H-%M-%S')}")
                                 )
            elif opt == 'multi_loci':
                self.assertTrue(cfg.multi_loci('set_1/Flagellum'))
                self.assertTrue(cfg.multi_loci('set_1/T4SS'))
                self.assertFalse(cfg.multi_loci('set_1/T6SS'))
            elif opt in cfg_needing_args:
                for model, val in expected_values[opt]:
                    self.assertEqual(getattr(cfg, opt)(model), int(val))
            elif opt == 'models_dir':
                self.assertEqual(getattr(cfg, opt)(), self.defaults['system_models_dir'])
            else:
                self.assertEqual(getattr(cfg, opt)(), exp_val)


    def test_model_conf(self):
        self.parsed_args.models_dir = self.find_data('models')
        self.parsed_args.models = "Model_w_conf all"
        cfg = Config(self.defaults, self.parsed_args)
        expected_weights = {'mandatory': 13.0,
                            'accessory': 14.0,
                            'neutral': 0.0,
                            'itself': 11.0,
                            'exchangeable': 12.0,
                            'out_of_cluster': 10.0}
        self.assertDictEqual(cfg.hit_weights(), expected_weights)
        self.assertEqual(cfg.i_evalue_sel(), 0.012)
        self.assertEqual(cfg.e_value_search(), 0.12)
        self.assertEqual(cfg.coverage_profile(), 0.55)
        self.assertFalse(cfg.cut_ga())


    def test_bad_values(self):
        invalid_syntax = {'inter_gene_max_space': 'set_1/Flagellum 4 2',
                          'max_nb_genes': 'set_1/Flagellum set_1/T3SS 3',
                          'min_genes_required': 'set_1/Flagellum set_1/T4SS 4',
                          'min_mandatory_genes_required': '12 set_1/T6SS 6',
                          }

        for opt, val in invalid_syntax.items():
            args = Namespace()
            setattr(args, opt, val)
            with self.assertRaises(ValueError) as ctx:
                Config(self.defaults, args)
            self.assertEqual(str(ctx.exception), f"Invalid syntax for '{opt}': You must provide a list of model name "
                                                 f"and value separated by spaces: {val}.")

        int_error = {'inter_gene_max_space': 'set_1/Flagellum 4.2 set_1/T2SS 2',
                     'max_nb_genes': 'set_1/Flagellum 4 set_1/T3SS FOO',
                     'min_genes_required': 'set_1/Flagellum FOO set_1/T4SS 4',
                     'min_mandatory_genes_required': 'set_1/Flagellum 12 set_1/T6SS 6.4',
                     }

        for opt, val in int_error.items():
            args = Namespace()
            setattr(args, opt, val)
            with self.assertRaises(ValueError):
                Config(self.defaults, args)


    def test_bad_db_type(self):
        self.parsed_args.db_type = "FOO"
        with self.assertRaises(ValueError) as ctx:
            Config(self.defaults, self.parsed_args)
        self.assertEqual(str(ctx.exception),
                         "db_type as unauthorized value : 'FOO'.")

    def test_bad_topology(self):
        self.parsed_args.replicon_topology = "FOO"
        with self.assertRaises(ValueError) as ctx:
            Config(self.defaults, self.parsed_args)
        self.assertEqual(str(ctx.exception),
                         "replicon_topology as unauthorized value : 'FOO'.")

    def test_bad_topology_file(self):
        self.parsed_args.topology_file = "FOO"
        with self.assertRaises(ValueError) as ctx:
            Config(self.defaults, self.parsed_args)
        self.assertEqual(str(ctx.exception),
                         "topology_file 'FOO' does not exists or is not a file.")


    def test_bad_sequence_db(self):
        self.parsed_args.sequence_db = "FOO"
        with self.assertRaises(ValueError) as ctx:
            Config(self.defaults, self.parsed_args)
        self.assertEqual(str(ctx.exception),
                         "sequence_db 'FOO' does not exists or is not a file.")


    def test_bad_models_dir(self):
        self.parsed_args.models_dir = "FOO"
        with self.assertRaises(ValueError) as ctx:
            Config(self.defaults, self.parsed_args)
        self.assertEqual(str(ctx.exception),
                         "models_dir 'FOO' does not exists or is not a directory.")

    def test_save(self):
        self.parsed_args.max_nb_genes = [['Set_1/T2SS', 5], ['set_1/Flagelum', 12]]
        self.parsed_args.multi_loci = 'Set_1/T2SS,set_1/Flagelum'
        self.parsed_args.models = ['Set_1', 'T9SS', 'T3SS', 'T4SS_typeI']
        cfg = Config(self.defaults, self.parsed_args)
        expected = {k: v for k, v in cfg._options.items() if v is not None}
        expected['max_nb_genes'] = 'Set_1/T2SS 5 set_1/Flagelum 12'
        expected['models'] = 'Set_1 T9SS T3SS T4SS_typeI'
        # save in file 'macsyfinder.conf'
        with tempfile.TemporaryDirectory() as tmpdirname:
            cfg_path = os.path.join(tmpdirname, 'macsyfinder.conf')
            cfg.save(path_or_buf=cfg_path)
            cfg.save(path_or_buf='/tmp/macsyfinder.conf')
            new_args = Namespace()
            new_args.cfg_file = cfg_path
            restored_cfg = Config(self.defaults, new_args)
            self.maxDiff = None
            # the option cfg-file differ from the 2 configs
            # None in cfg
            # cfg_path in restored_cfg
            self.assertEqual(restored_cfg._options['cfg_file'], cfg_path)
            del(cfg._options['cfg_file'])
            del(restored_cfg._options['cfg_file'])
            self.assertDictEqual(cfg._options, restored_cfg._options)


    def test_out_dir(self):
        cfg = Config(self.defaults, self.parsed_args)
        self.assertEqual(cfg.out_dir(),
                         os.path.join(cfg.res_search_dir(), f'macsyfinder-{strftime("%Y%m%d_%H-%M-%S")}')
                         )
        self.parsed_args.out_dir = 'foo'
        cfg = Config(self.defaults, self.parsed_args)
        self.assertEqual(cfg.out_dir(), 'foo')

    def test_working_dir(self):
        cfg = Config(self.defaults, self.parsed_args)
        self.assertEqual(cfg.out_dir(), cfg.working_dir())

    def test_models_dir(self):
        self.parsed_args.models_dir = self.tmp_dir
        cfg = Config(self.defaults, self.parsed_args)
        self.assertEqual(cfg.models_dir(), [self.tmp_dir])

    def test_previous_n_sequence_db(self):
        self.parsed_args.previous_run = self.find_data(os.path.join('data_set', 'results_1'))
        sequence_db = self.find_data(os.path.join('base', 'test_1.fasta'))
        self.parsed_args.sequence_db = sequence_db
        with self.catch_log() as log:
            cfg = Config(self.defaults, self.parsed_args)
            # The config set the parsed_args.sequence_db to None
            catch_msg = log.get_value().strip()
        self.assertEqual(cfg.sequence_db(), 'tests/data/base/VICH001.B.00001.C001.prt')
        self.maxDiff = None
        self.assertEqual(f"ignore sequence_db '{sequence_db}' use sequence_db from previous_run "
                         f"'{os.path.abspath(cfg.previous_run())}'.",
                         catch_msg)

    def test_previous_wo_cfg(self):
        self.parsed_args.previous_run = self.find_data(os.path.join('data_set'))
        with self.assertRaises(ValueError) as ctx:
            Config(self.defaults, self.parsed_args)
        self.assertEqual(str(ctx.exception),
                         f"No config file found in dir {self.parsed_args.previous_run}")

    def test_cut_ga(self):
        # default value
        cfg = Config(self.defaults, self.parsed_args)
        self.assertTrue(cfg.cut_ga())
        # user specify --no-cu-ga
        self.parsed_args.no_cut_ga = True
        cfg = Config(self.defaults, self.parsed_args)
        self.assertFalse(cfg.cut_ga())

        # the user does not specify anything
        # the modeler set cut-ga false in model_conf
        self.parsed_args.cut_ga = None
        self.parsed_args.no_cut_ga = None
        self.parsed_args.models_dir = self.find_data('models')
        self.parsed_args.models = ['TFF-SF', 'all']
        cfg = Config(self.defaults, self.parsed_args)
        self.assertFalse(cfg.cut_ga())

        # the user specify --cut-ga
        # the modeler set cut-ga false in model_conf
        self.parsed_args.cut_ga = True
        self.parsed_args.no_cut_ga = None
        self.parsed_args.models_dir = self.find_data('models')
        self.parsed_args.models = ['TFF-SF', 'all']
        cfg = Config(self.defaults, self.parsed_args)
        self.assertTrue(cfg.cut_ga())

    @unittest.skipIf(not shutil.which('hmmsearch'), "'hmmsearch' not found in PATH")
    def test_hmmsearch(self):
        cfg = Config(self.defaults, self.parsed_args)
        self.assertEqual(cfg.hmmer(), shutil.which('hmmsearch'))

    def test_hmmsearch_set_by_user(self):
        self.parsed_args.hmmer = 'nimportnaoik'
        cfg = Config(self.defaults, self.parsed_args)
        self.assertEqual(cfg.hmmer(), self.parsed_args.hmmer)

    def test_no_hmmsearch(self):
        # if hmmsearch is not found in path the default is set to None (shutil.which('hmmsearch'))
        self.defaults.hmmer = None
        with self.assertRaises(ValueError) as ctx:
            with self.catch_log():
                Config(self.defaults, self.parsed_args)
        self.assertEqual(str(ctx.exception),
                         "'hmmsearch' NOT found in your PATH, Please specify hmmsearch path with --hmmer opt or"
                         " install 'hmmer' package.")

    def test_e_value_search(self):
        cfg = Config(self.defaults, self.parsed_args)
        self.assertEqual(self.defaults.e_value_search, cfg.e_value_search())

        self.parsed_args.e_value_search = 1.0
        cfg = Config(self.defaults, self.parsed_args)
        self.assertEqual(cfg.e_value_search(), 1.0)

    def test_hit_weights(self):
        cfg = Config(self.defaults, self.parsed_args)
        default = {k: self.defaults[f"{k}_weight"] for k in ('mandatory', 'accessory', 'neutral', 'itself',
                                                             'exchangeable', 'out_of_cluster')}
        self.assertDictEqual(default, cfg.hit_weights())

        self.parsed_args.mandatory_weight = 2
        self.parsed_args.accessory_weight = 1
        self.parsed_args.exchangeable_weight = 0.5
        self.parsed_args.out_of_cluster_weight = 0.2

        cfg = Config(self.defaults, self.parsed_args)
        expected = {'mandatory': 2,
                    'accessory': 1,
                    'neutral': self.defaults.neutral_weight,
                    'itself': self.defaults.itself_weight,
                    'exchangeable': .5,
                    'out_of_cluster': .2
                    }
        self.assertDictEqual(cfg.hit_weights(), expected)


class TestNoneConfig(MacsyTest):

    def test_NoneConfig(self):
        n_cfg = NoneConfig()

        self.assertIsNone(n_cfg.max_nb_genes("model_1"))
        self.assertIsNone(n_cfg.working_dir())
