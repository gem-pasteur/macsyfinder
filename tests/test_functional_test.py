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


import shutil
import tempfile
import os
import inspect
import unittest
import json

from tests import MacsyTest, which
from macsypy.scripts import macsyfinder
from macsypy.error import OptionError


class Test(MacsyTest):

    def setUp(self):
        self.tmp_dir = tempfile.gettempdir()

    def tearDown(self):
        try:
            shutil.rmtree(self.out_dir)
        except:
            pass

    @unittest.skip("skipping until macsyfinder api is not stable")
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

    @unittest.skip("skipping until macsyfinder api is not stable")
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

    @unittest.skip("skipping until macsyfinder api is not stable")
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

    @unittest.skip("skipping until macsyfinder api is not stable")
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

    @unittest.skip("skipping until macsyfinder api is not stable")
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

    @unittest.skip("skipping until macsyfinder api is not stable")
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

    @unittest.skip("skipping until macsyfinder api is not stable")
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

    @unittest.skip("skipping until macsyfinder api is not stable")
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


    def test_working_dir_exists(self):
        args = "--sequence-db {seq_db} --db-type ordered_replicon " \
               "--models-dir {models_dir} " \
               "-m {models} -o {{out_dir}} ".format(
                                                 models_dir=self.find_data('functional_tests', 'models'),
                                                 seq_db=self.find_data('functional_tests',
                                                                       'acav001_T2SS-Ori_T4P-multi.fasta'),
                                                 models="TXSS_modified T4P-single-locus",
                                                 topology="linear"
                                                 )
        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_working_dir_exists')
        os.makedirs(self.out_dir)
        open(os.path.join(self.out_dir, 'toto.empty'), 'w').close()

        args = args.format(out_dir=self.out_dir)
        with self.assertRaises(ValueError) as ctx:
            macsyfinder.main(args=args.split(), loglevel='ERROR')
        self.assertEqual(str(ctx.exception),
                         f"'{self.out_dir}' already exists and is not a empty")

    def test_working_dir_exists_and_not_dir(self):
        args = "--sequence-db {seq_db} --db-type ordered_replicon " \
               "--models-dir {models_dir} " \
               "-m {models} -o {{out_dir}} ".format(
                                                 models_dir=self.find_data('functional_tests', 'models'),
                                                 seq_db=self.find_data('functional_tests',
                                                                       'acav001_T2SS-Ori_T4P-multi.fasta'),
                                                 models="TXSS_modified T4P-single-locus",
                                                 topology="linear"
                                                 )
        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_working_dir_exists_and_not_dir')
        try:
            open(self.out_dir, 'w').close()

            args = args.format(out_dir=self.out_dir)
            with self.assertRaises(ValueError) as ctx:
                macsyfinder.main(args=args.split(), loglevel='ERROR')
            self.assertEqual(str(ctx.exception),
                             f"'{self.out_dir}' already exists and is not a directory")
        finally:
            os.unlink(self.out_dir)


    def test_no_models(self):
        args = "--sequence-db {seq_db} --db-type ordered_replicon " \
               "--models-dir {models_dir} " \
               "-o {{out_dir}} ".format(
                                         models_dir=self.find_data('functional_tests', 'models'),
                                         seq_db=self.find_data('functional_tests',
                                                               'acav001_T2SS-Ori_T4P-multi.fasta'),
                                         topology="linear"
                                        )
        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_no_models')
        args = args.format(out_dir=self.out_dir)
        with self.catch_io(out=True):
            with self.assertRaises(OptionError) as ctx:
                macsyfinder.main(args=args.split(), loglevel='ERROR')
        self.assertEqual(str(ctx.exception),
                         "argument --models or --previous-run is required.")

    def test_no_seq_db(self):
        args = "--db-type ordered_replicon " \
               "--models-dir {models_dir} " \
               "-o {{out_dir}} " \
               "-m {models} -o {{out_dir}} ".format(
                                                     models_dir=self.find_data('functional_tests', 'models'),
                                                     models="TXSS_modified T4P-single-locus",
                                                     topology="linear"
                                                     )
        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_no_seq_db')

        args = args.format(out_dir=self.out_dir)
        with self.catch_io(out=True):
            with self.assertRaises(OptionError) as ctx:
                macsyfinder.main(args=args.split(), loglevel='ERROR')
        self.assertEqual(str(ctx.exception),
                         "argument --sequence-db or --previous-run is required.")

    def test_no_db_type(self):
        args = "--sequence-db {seq_db} " \
               "--models-dir {models_dir} " \
               "-o {{out_dir}} "\
               "-m {models} -o {{out_dir}} ".format(seq_db=self.find_data('functional_tests',
                                                                          'acav001_T2SS-Ori_T4P-multi.fasta'),
                                                    models_dir=self.find_data('functional_tests', 'models'),
                                                    models="TXSS_modified T4P-single-locus",
                                                    topology="linear"
                                                    )

        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_no_db_type')
        args = args.format(out_dir=self.out_dir)
        with self.catch_io(out=True):
            with self.assertRaises(OptionError) as ctx:
                macsyfinder.main(args=args.split(), loglevel='ERROR')
        self.assertEqual(str(ctx.exception),
                         "argument --db-type or --previous-run is required.")

    def test_model_unkown(self):
        args = "--sequence-db {seq_db} " \
               "--db-type ordered_replicon " \
               "--models-dir {models_dir} " \
               "-o {{out_dir}} "\
               "-m {models} -o {{out_dir}} ".format(seq_db=self.find_data('functional_tests',
                                                                          'acav001_T2SS-Ori_T4P-multi.fasta'),
                                                    models_dir=self.find_data('functional_tests', 'models'),
                                                    models="TXSS_modified Unkown_model",
                                                    topology="linear"
                                                    )

        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_model_unkwon')
        os.makedirs(self.out_dir)

        args = args.format(out_dir=self.out_dir)
        with self.assertRaises(ValueError) as ctx:
            macsyfinder.main(args=args.split(), loglevel='ERROR')
        self.assertEqual(str(ctx.exception),
                         "Unkown_model does not match with any definitions")


    def _test_macsyfinder_run(self, args_tpl):
        # get the name of the calling function
        test_name = inspect.stack()[1].function
        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_{}'.format(test_name))
        os.makedirs(self.out_dir)
        args = args_tpl.format(out_dir=self.out_dir)

        macsyfinder.main(args=args.split(),
                         loglevel='ERROR'
                         )

        expected_result_path = self.find_data('functional_tests', test_name[5:], 'results.macsyfinder.json')
        get_results = os.path.join(self.out_dir, 'results.macsyfinder.json')
        self.assertJsonEqual(expected_result_path, get_results)