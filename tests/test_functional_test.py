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


import shutil
import tempfile
import os
import inspect
import unittest
import json

from tests import MacsyTest, which
from macsypy.scripts import macsyfinder
import macsypy
from macsypy import model


@unittest.skip("skipping until macsyfinder api is not stable")
class Test(MacsyTest):

    def setUp(self):
        self.tmp_dir = tempfile.gettempdir()
        self.gene_bank = macsypy.gene.GeneBank()
        self.model_bank = model.ModelBank()
        self.profile_factory = macsypy.gene.ProfileFactory()
        search_systems.system_name_generator = search_systems.SystemNameGenerator()

    def tearDown(self):
        try:
            shutil.rmtree(self.out_dir)
        except:
            pass


    @unittest.skipIf(not which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_basic_run(self):
        """
        test if returncode of macsyfinder is 0 and
        test each element of the json
        macsyfinder is launched to search T9SS T3SS T4SS_typeI systems
        with test_aesu.fa sequence db in gembase format
        """
        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_test_basic_run')
        os.makedirs(self.out_dir)
        args = "--out-dir={out_dir} --sequence-db={seq_db} --db-type=gembase --models-dir={models_dir}" \
                  " --models {models}".format(out_dir=self.out_dir,
                                              models_dir=os.path.join(self._data_dir, 'data_set_1', 'models'),
                                              seq_db=os.path.join(self._data_dir, 'base', 'test_aesu.fa'),
                                              models="set_1 T9SS T3SS T4SS_typeI",
                                              )

        macsyfinder.main(args=args.split(), loglevel='ERROR')

        expected_result_path = self.find_data('data_set_1', 'basic_run_results', 'results.macsyfinder.json')
        with open(expected_result_path) as expected_result_file:
            expected_result_json = json.load(expected_result_file)

        test_result_path = os.path.join(self.out_dir, 'results.macsyfinder.json')
        with open(test_result_path) as test_result_file:
            test_result_json = json.load(test_result_file)

        # it should have only one occurrence of T9SS
        self.assertEqual(len(test_result_json), 1,
                         "different type of systems expected: 1  retrieved: {}".format(len(test_result_json)))
        expected_result_json = expected_result_json[0]
        test_result_json = test_result_json[0]
        self.assertEqual(expected_result_json['name'],
                         test_result_json['name'],
                         "type of system name expected: {} retrieved: {}".format(expected_result_json['name'],
                                                                                 test_result_json['name']))
        self.assertEqual(expected_result_json['occurrence_number'],
                         test_result_json['occurrence_number'],
                         "occurrence number expected {} retrieved: {}".format(expected_result_json['occurrence_number'],
                                                                              test_result_json['occurrence_number']))
        self.assertDictEqual(expected_result_json['replicon'],
                             test_result_json['replicon'],
                             "replicon expected {} retrieved: {}".format(expected_result_json['occurrence_number'],
                                                                         test_result_json['occurrence_number']))
        self.assertEqual(expected_result_json['id'],
                         test_result_json['id'],
                         "system occurrence id expected {} retrieved: {}".format(expected_result_json['id'],
                                                                                 test_result_json['id']))
        self.assertDictEqual(expected_result_json['summary']['mandatory'],
                             test_result_json['summary']['mandatory'],
                             "\nmandatory genes expected: {}"
                             "\nmandatory genes retrieved: {}".format(expected_result_json['summary']['mandatory'],
                                                                      test_result_json['summary']['mandatory']))
        self.assertDictEqual(expected_result_json['summary']['accessory'],
                             test_result_json['summary']['accessory'],
                             "\naccessory genes expected: {}"
                             "\naccessory genes retrieved: {}".format(expected_result_json['summary']['accessory'],
                                                                      test_result_json['summary']['accessory']))
        self.assertDictEqual(expected_result_json['summary']['forbidden'],
                             test_result_json['summary']['forbidden'],
                             "\nforbidden genes expected: {}"
                             "\nforbidden genes retrieved: {}".format(expected_result_json['summary']['forbidden'],
                                                                      test_result_json['summary']['forbidden']))
        self.assertListEqual(expected_result_json['genes'], test_result_json['genes'],
                             "\ngenes expected: {}"
                             "\ngenes retrieved: {}".format(expected_result_json['genes'], test_result_json['genes']))


    def test_T2SS_ordered_circular(self):
        args = "--sequence-db {seq_db} --db-type ordered_replicon --replicon-topology {topology}  " \
               "--models-dir {models_dir} " \
               "-m {models} -o {{out_dir}}".format(
                                                   models_dir=self.find_data('functional_tests', 'models'),
                                                   seq_db=self.find_data('functional_tests',
                                                                       'acav001_T2SS-Ori_T4P-multi.fasta'),
                                                   models="TXSS_ori T2SS",
                                                   topology="circular"
                                                   )
        self._test_macsyfinder_run(args)


    def test_T2SS_ordered_circular_single_locus(self):
        args = "--sequence-db {seq_db} --db-type ordered_replicon --replicon-topology {topology}  " \
               "--models-dir {models_dir} " \
               "-m {models} -o {{out_dir}}".format(models_dir=self.find_data('functional_tests', 'models'),
                                                   seq_db=self.find_data('functional_tests',
                                                                         'acav001_T2SS-Ori_T4P-multi.fasta'),
                                                   models="TXSS_modified T2SS-single-locus",
                                                   topology="circular"
                                                   )
        self._test_macsyfinder_run(args)


    def test_T2SS_ordered_linear(self):
        args = "--sequence-db {seq_db} --db-type ordered_replicon --replicon-topology {topology}  " \
               "--models-dir {models_dir} " \
               "-m {models} -o {{out_dir}}".format(models_dir=self.find_data('functional_tests', 'models'),
                                                   seq_db=self.find_data('functional_tests',
                                                                         'acav001_T2SS-Ori_T4P-multi.fasta'),
                                                   models="TXSS_ori T2SS",
                                                   topology="linear"
                                                   )
        self._test_macsyfinder_run(args)


    def test_T2SS_ordered_linear_multi_loci(self):
        args = "--sequence-db {seq_db} --db-type ordered_replicon --replicon-topology {topology}  " \
               "--models-dir {models_dir} " \
               "-m {models} -o {{out_dir}} " \
               "--multi-loci TXSS_modified/T2SS-single-locus".format(
                                                 models_dir=self.find_data('functional_tests', 'models'),
                                                 seq_db=self.find_data('functional_tests',
                                                                       'acav001_T2SS-Ori_T4P-multi.fasta'),
                                                 models="TXSS_modified T2SS-single-locus",
                                                 topology="linear"
                                                 )
        self._test_macsyfinder_run(args)


    def test_T2SS_ordered_linear_single_locus(self):
        args = "--sequence-db {seq_db} --db-type ordered_replicon --replicon-topology {topology}  " \
               "--models-dir {models_dir} " \
               "-m {models} -o {{out_dir}}".format(models_dir=self.find_data('functional_tests', 'models'),
                                                   seq_db=self.find_data('functional_tests',
                                                                         'acav001_T2SS-Ori_T4P-multi.fasta'),
                                                   models="TXSS_modified T2SS-single-locus",
                                                   topology="linear"
                                                   )
        self._test_macsyfinder_run(args)


    def test_T4P_ordered_linear(self):
        args = "--sequence-db {seq_db} --db-type ordered_replicon --replicon-topology {topology}  " \
               "--models-dir {models_dir} " \
               "-m {models} -o {{out_dir}}".format(models_dir=self.find_data('functional_tests', 'models'),
                                                   seq_db=self.find_data('functional_tests',
                                                                         'acav001_T2SS-Ori_T4P-multi.fasta'),
                                                   models="TXSS_modified T4P-single-locus",
                                                   topology="linear"
                                                   )
        self._test_macsyfinder_run(args)


    def test_T4P_ordered_linear_multi_loci(self):
        args = "--sequence-db {seq_db} --db-type ordered_replicon --replicon-topology {topology}  " \
               "--models-dir {models_dir} " \
               "-m {models} -o {{out_dir}} " \
               "--multi-loci TXSS_modified/T4P-single-locus".format(
                                                 models_dir=self.find_data('functional_tests', 'models'),
                                                 seq_db=self.find_data('functional_tests',
                                                                       'acav001_T2SS-Ori_T4P-multi.fasta'),
                                                 models="TXSS_modified T4P-single-locus",
                                                 topology="linear"
                                                 )
        self._test_macsyfinder_run(args)

    def _test_macsyfinder_run(self, args_tpl):
        # get the name of the calling function
        test_name = inspect.stack()[1].function
        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_{}'.format(test_name))
        os.makedirs(self.out_dir)
        args = args_tpl.format(out_dir=self.out_dir)
        macsyfinder.main(args=args.split(),
                         loglevel='ERROR',
                         models=self.model_bank,
                         genes=self.gene_bank,
                         profiles=self.profile_factory)

        expected_result_path = self.find_data('functional_tests', test_name[5:], 'results.macsyfinder.json')
        get_results = os.path.join(self.out_dir, 'results.macsyfinder.json')
        self.assertJsonEqual(expected_result_path, get_results)