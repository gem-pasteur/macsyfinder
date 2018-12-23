# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Néron                                         #
# Copyright © 2014  Institut Pasteur, Paris.                                   #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################

import os
import logging
from argparse import Namespace
from configparser import ParsingError

from macsypy.config_new import MacsyDefaults, Config

from tests import MacsyTest


class TestConfig(MacsyTest):

    def setUp(self):
        self.defaults = MacsyDefaults()
        self.parsed_args = Namespace()


    def test_str_2_tuple(self):
        s = 'Flagellum 12 t4ss 13'
        expected = [('Flagellum', '12'), ('t4ss', '13')]
        cfg = Config(self.defaults, self.parsed_args)
        self.assertListEqual(cfg._str_2_tuple(s), expected)

        with self.assertRaises(ValueError) as ctx:
            s = 'Flagellum 12 t4ss'
            cfg._str_2_tuple(s)
        self.assertEqual(str(ctx.exception),
                         "You must provide a list of model name and value separated by spaces")


    def test_config_file_2_dict(self):
        cfg = Config(self.defaults, self.parsed_args)
        res = cfg._config_file_2_dict(self.defaults, ['nimportnaoik'])
        self.assertDictEqual({}, res)

        cfg_file = self.find_data(os.path.join('conf_files', 'macsy_test_conf.conf'))
        res = cfg._config_file_2_dict(self.defaults, [cfg_file])
        expected = {'db_type': 'gembase',
                    'inter_gene_max_space': 'T2SS 2 Flagellum 4',
                    'min_mandatory_genes_required': 'T2SS 5 Flagellum 9',
                    'replicon_topology': 'circular',
                    'sequence_db': '/path/to/sequence/bank/fasta_file',
                    'topology_file': '/the/path/to/the/topology/to/use'}
        self.assertDictEqual(expected, res)

        bad_cfg_file = self.find_data(os.path.join('conf_files', 'macsy_test_bad_conf.conf'))
        with self.assertRaises(ParsingError):
            res = cfg._config_file_2_dict(self.defaults, [bad_cfg_file])


    def test_Config(self):
        cfg = Config(self.defaults, self.parsed_args)
        methods_needing_args = {'inter_gene_max_space': None,
                                'max_nb_genes': None,
                                'min_genes_required': None,
                                'min_mandatory_genes_required':None,
                                }
        for opt, val in self.defaults.items():
            if opt in methods_needing_args:
                self.assertEqual(getattr(cfg, opt)('whatever'), val)
            else:
                self.assertEqual(getattr(cfg, opt)(), val)


    def test_Config_file(self):
        methods_needing_args = {'inter_gene_max_space': [('Flagellum', 4), ('T2SS', 2)],
                                'max_nb_genes':  [('Flagellum', 6), ('T3SS', 3)],
                                'min_genes_required': [('Flagellum', 8), ('T4SS', 4)],
                                'min_mandatory_genes_required': [('Flagellum', 12), ('T6SS', 6)],
                                }

        self.parsed_args.cfg_file = self.find_data(os.path.join('conf_files', 'macsy_models.conf'))
        cfg = Config(self.defaults, self.parsed_args)

        expected_values = {k: v for k, v in self.defaults.items()}
        expected_values['cfg_file'] = self.parsed_args.cfg_file
        expected_values.update(methods_needing_args)

        for opt, val in expected_values.items():
            if opt in methods_needing_args:
                for model, genes in methods_needing_args[opt]:
                    self.assertEqual(getattr(cfg, opt)(model), genes)
            else:
                self.assertEqual(getattr(cfg, opt)(), val)


    def test_Config_args(self):
        methods_needing_args = {'inter_gene_max_space': [('Flagellum', 14), ('T2SS', 12)],
                                'max_nb_genes': [('Flagellum', 16), ('T3SS', 13)],
                                'min_genes_required': [('Flagellum', 18), ('T4SS', 14)],
                                'min_mandatory_genes_required': [('Flagellum', 22), ('T6SS', 16)],
                                }

        self.parsed_args.cfg_file = self.find_data(os.path.join('conf_files', 'macsy_models.conf'))
        cfg = Config(self.defaults, self.parsed_args)

        expected_values = {k: v for k, v in self.defaults.items()}
        expected_values['cfg_file'] = self.parsed_args.cfg_file
        expected_values.update(methods_needing_args)
        print()
        for opt, val in expected_values.items():
            print("======================")
            print(opt)
            if opt in methods_needing_args:
                print("methods with args")
                print(cfg._options[opt])
                for model, genes in methods_needing_args[opt]:
                    self.assertEqual(getattr(cfg, opt)(model), genes)
            else:
                self.assertEqual(getattr(cfg, opt)(), val)