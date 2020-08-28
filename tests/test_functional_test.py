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
import itertools

from tests import MacsyTest, which
from macsypy.scripts import macsyfinder
from macsypy.error import OptionError
from macsypy.system import AbstractSetOfHits


class Test(MacsyTest):

    def setUp(self):
        self.tmp_dir = tempfile.gettempdir()
        # reset AbstractSetOfHits internal id to have predictable results (Systems, ...) id
        # it's works only if there is only one replicon
        # for gembase the order is not guarantee

        AbstractSetOfHits._id = itertools.count(1)

    def tearDown(self):
        try:
            shutil.rmtree(self.out_dir)
        except:
            pass

    @unittest.skipIf(not which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_gembase(self):
        """

        """
        args = f"--sequence-db={self.find_data('base', 'gembase.fasta')} " \
               "--db-type=gembase " \
               f"--models-dir={self.find_data('models')} " \
               "--models TFF-SF Archaeal-T4P ComM MSH T2SS T4bP T4P Tad " \
               "--out-dir={out_dir}"

        self._macsyfinder_run(args)
        expected_result_dir = self.find_data("functional_tests_gembase")
        for file_name in ('all_systems.tsv',
                          'all_best_systems.tsv',
                          'best_systems.tsv'):
            print(f"\n##################################### {file_name} ##########################################\n")
            expected_result = self.find_data(expected_result_dir, file_name)
            get_results = os.path.join(self.out_dir, file_name)
            self.assertTsvEqual(expected_result, get_results, comment="#")
        file_name = 'rejected_clusters.txt'
        expected_result = self.find_data(expected_result_dir, file_name)
        get_results = os.path.join(self.out_dir, file_name)
        self.assertFileEqual(expected_result, get_results, comment="#")


    # @unittest.skip("skipping until macsyfinder api is not stable")
    # def test_T2SS_ordered_circular(self):
    #     args = "--sequence-db {seq_db} --db-type ordered_replicon --replicon-topology {topology}  " \
    #            "--models-dir {models_dir} " \
    #            "-m {models} -o {{out_dir}}".format(
    #                                                models_dir=self.find_data('functional_tests', 'models'),
    #                                                seq_db=self.find_data('functional_tests',
    #                                                                    'acav001_T2SS-Ori_T4P-multi.fasta'),
    #                                                models="TXSS_ori T2SS",
    #                                                topology="circular"
    #                                                )
    #     self._test_macsyfinder_run(args)
    #
    # @unittest.skip("skipping until macsyfinder api is not stable")
    # def test_T2SS_ordered_circular_single_locus(self):
    #     args = "--sequence-db {seq_db} --db-type ordered_replicon --replicon-topology {topology}  " \
    #            "--models-dir {models_dir} " \
    #            "-m {models} -o {{out_dir}}".format(models_dir=self.find_data('functional_tests', 'models'),
    #                                                seq_db=self.find_data('functional_tests',
    #                                                                      'acav001_T2SS-Ori_T4P-multi.fasta'),
    #                                                models="TXSS_modified T2SS-single-locus",
    #                                                topology="circular"
    #                                                )
    #     self._test_macsyfinder_run(args)
    #
    # @unittest.skip("skipping until macsyfinder api is not stable")
    # def test_T2SS_ordered_linear(self):
    #     args = "--sequence-db {seq_db} --db-type ordered_replicon --replicon-topology {topology}  " \
    #            "--models-dir {models_dir} " \
    #            "-m {models} -o {{out_dir}}".format(models_dir=self.find_data('functional_tests', 'models'),
    #                                                seq_db=self.find_data('functional_tests',
    #                                                                      'acav001_T2SS-Ori_T4P-multi.fasta'),
    #                                                models="TXSS_ori T2SS",
    #                                                topology="linear"
    #                                                )
    #     self._test_macsyfinder_run(args)
    #
    # @unittest.skip("skipping until macsyfinder api is not stable")
    # def test_T2SS_ordered_linear_multi_loci(self):
    #     args = "--sequence-db {seq_db} --db-type ordered_replicon --replicon-topology {topology}  " \
    #            "--models-dir {models_dir} " \
    #            "-m {models} -o {{out_dir}} " \
    #            "--multi-loci TXSS_modified/T2SS-single-locus".format(
    #                                              models_dir=self.find_data('functional_tests', 'models'),
    #                                              seq_db=self.find_data('functional_tests',
    #                                                                    'acav001_T2SS-Ori_T4P-multi.fasta'),
    #                                              models="TXSS_modified T2SS-single-locus",
    #                                              topology="linear"
    #                                              )
    #     self._test_macsyfinder_run(args)
    #
    # @unittest.skip("skipping until macsyfinder api is not stable")
    # def test_T2SS_ordered_linear_single_locus(self):
    #     args = "--sequence-db {seq_db} --db-type ordered_replicon --replicon-topology {topology}  " \
    #            "--models-dir {models_dir} " \
    #            "-m {models} -o {{out_dir}}".format(models_dir=self.find_data('functional_tests', 'models'),
    #                                                seq_db=self.find_data('functional_tests',
    #                                                                      'acav001_T2SS-Ori_T4P-multi.fasta'),
    #                                                models="TXSS_modified T2SS-single-locus",
    #                                                topology="linear"
    #                                                )
    #     self._test_macsyfinder_run(args)


    @unittest.skipIf(not which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_only_loners(self):
        args = f"--sequence-db {self.find_data('base', 'MOBP1_twice.fasta')} " \
               "--db-type ordered_replicon " \
               "--replicon-topology linear  " \
               f"--models-dir {self.find_data('models')} " \
               "-m test_loners MOB_cf_T5SS " \
               "-o {out_dir}"
        self._macsyfinder_run(args)

        expected_result_dir = self.find_data("functional_tests_only_loners")
        for file_name in ('all_systems.tsv',
                          'all_best_systems.tsv',
                          'best_systems.tsv',
                          'rejected_clusters.txt'):
            expected_result = self.find_data(expected_result_dir, file_name)
            get_results = os.path.join(self.out_dir, file_name)
            self.assertFileEqual(expected_result, get_results, comment="#")


    @unittest.skipIf(not which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_T4P_ordered_linear(self):
        # TODO how to specify multi_loci = false when multi_loci =True is set in xml
        args = f"--sequence-db {self.find_data('base', 'one_replicon.fasta')} " \
               "--db-type ordered_replicon " \
               "--replicon-topology linear  " \
               f"--models-dir {self.find_data('models')} " \
               "-m TFF-SF T4P_single_locus " \
               "-o {out_dir}"
        self._macsyfinder_run(args)

        expected_result_dir = self.find_data("functional_tests_T4P_ordered_linear")
        for file_name in ('all_systems.tsv',
                          'all_best_systems.tsv',
                          'best_systems.tsv',
                          'rejected_clusters.txt'):
            expected_result = self.find_data(expected_result_dir, file_name)
            get_results = os.path.join(self.out_dir, file_name)
            self.assertFileEqual(expected_result, get_results, comment="#")


    @unittest.skipIf(not which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_T4P_ordered_linear_multi_loci(self):
        args = f"--sequence-db {self.find_data('base', 'one_replicon.fasta')} " \
               "--db-type ordered_replicon " \
               "--replicon-topology linear  " \
               f"--models-dir {self.find_data('models')} " \
               "-m TFF-SF T4P " \
               "-o {out_dir} " \
               "--multi-loci TFF-SF/T4P"

        self._macsyfinder_run(args)

        expected_result_dir = self.find_data("functional_tests_T4P_ordered_linear_multi_loci")
        for file_name in ('all_systems.tsv',
                          'all_best_systems.tsv',
                          'best_systems.tsv',
                          'rejected_clusters.txt'):
            expected_result = self.find_data(expected_result_dir, file_name)
            get_results = os.path.join(self.out_dir, file_name)
            self.assertFileEqual(expected_result, get_results, comment="#")


    def test_working_dir_exists(self):
        args = f"--sequence-db {self.find_data('base', 'one_replicon.fasta')} " \
               "--db-type ordered_replicon " \
               f"--models-dir {self.find_data('models')} " \
               "-m TFF-SF T4P " \
               "-o {out_dir}"

        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_working_dir_exists')
        os.makedirs(self.out_dir)
        open(os.path.join(self.out_dir, 'toto.empty'), 'w').close()

        args = args.format(out_dir=self.out_dir)
        with self.assertRaises(ValueError) as ctx:
            macsyfinder.main(args=args.split(), loglevel='ERROR')
        self.assertEqual(str(ctx.exception),
                         f"'{self.out_dir}' already exists and is not a empty")

    def test_working_dir_exists_and_not_dir(self):
        args = f"--sequence-db {self.find_data('base', 'one_replicon.fasta')} " \
               "--db-type ordered_replicon " \
               f"--models-dir {self.find_data('models')} " \
               "-m TFF-SF T4P " \
               "-o {out_dir} "

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
        args = f"--sequence-db {self.find_data('base', 'one_replicon.fasta')} " \
               "--db-type ordered_replicon " \
               f"--models-dir {self.find_data('models')} " \
               "-o {out_dir} "

        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_no_models')
        args = args.format(out_dir=self.out_dir)
        with self.catch_io(out=True):
            with self.assertRaises(OptionError) as ctx:
                macsyfinder.main(args=args.split(), loglevel='ERROR')
        self.assertEqual(str(ctx.exception),
                         "argument --models or --previous-run is required.")

    def test_no_seq_db(self):
        args = "--db-type ordered_replicon " \
               f"--models-dir {self.find_data('models')} " \
               "-m TFF-SF T4P " \
               "-o {out_dir} "

        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_no_seq_db')

        args = args.format(out_dir=self.out_dir)
        with self.catch_io(out=True):
            with self.assertRaises(OptionError) as ctx:
                macsyfinder.main(args=args.split(), loglevel='ERROR')
        self.assertEqual(str(ctx.exception),
                         "argument --sequence-db or --previous-run is required.")

    def test_no_db_type(self):
        args = f"--sequence-db {self.find_data('base', 'one_replicon.fasta')} " \
               f"--models-dir {self.find_data('models')} " \
               "-m TFF-SF T4P " \
               "-o {out_dir} "

        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_no_db_type')
        args = args.format(out_dir=self.out_dir)
        with self.catch_io(out=True):
            with self.assertRaises(OptionError) as ctx:
                macsyfinder.main(args=args.split(), loglevel='ERROR')
        self.assertEqual(str(ctx.exception),
                         "argument --db-type or --previous-run is required.")

    def test_model_unknown(self):
        args = f"--sequence-db {self.find_data('base', 'one_replicon.fasta')} " \
               "--db-type ordered_replicon " \
               f"--models-dir {self.find_data('models')} " \
               "-m TFF-SF Unknown_model " \
               "-o {out_dir}"

        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_model_unkwon')
        os.makedirs(self.out_dir)

        args = args.format(out_dir=self.out_dir)
        with self.assertRaises(ValueError) as ctx:
            macsyfinder.main(args=args.split(), loglevel='ERROR')
        self.assertEqual(str(ctx.exception),
                         "Unknown_model does not match with any definitions")


    def _macsyfinder_run(self, args_tpl):
        # get the name of the calling function
        test_name = inspect.stack()[1].function
        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_{}'.format(test_name))
        os.makedirs(self.out_dir)
        args = args_tpl.format(out_dir=self.out_dir)
        print("\n############################################3")
        print(args)
        print("############################################3")
        macsyfinder.main(args=args.split(),
                         loglevel='ERROR'
                         )
