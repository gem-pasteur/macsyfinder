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
        self.all_systems_tsv = "all_systems.tsv"
        self.all_systems_txt = "all_systems.txt"
        self.all_best_solutions = "all_best_solutions.tsv"
        self.best_solution = "best_solution.tsv"
        self.rejected_clusters = "rejected_clusters.txt"
        self.uncomplete_systems = "uncomplete_systems.txt"


    def tearDown(self):
        try:
            pass
            # self.out_dir is set in self._macsyfinder_run
            shutil.rmtree(self.out_dir)
        except:
            pass

    @unittest.skipIf(not which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_gembase(self):
        """

        """
        expected_result_dir = self.find_data("functional_test_gembase")
        args = "--db-type=gembase " \
               f"--models-dir={self.find_data('models')} " \
               "--models TFF-SF Archaeal-T4P ComM MSH T2SS T4bP T4P Tad " \
               "--out-dir={out_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"

        self._macsyfinder_run(args)
        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#")
        expected_result = self.find_data(expected_result_dir, self.rejected_clusters)
        get_results = os.path.join(self.out_dir, self.rejected_clusters)
        self.assertFileEqual(expected_result, get_results, comment="#")


    def test_only_loners(self):
        expected_result_dir = self.find_data("functional_tests_only_loners")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear  " \
               f"--models-dir {self.find_data('models')} " \
               "-m test_loners MOB_cf_T5SS " \
               "-o {out_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.rejected_clusters):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertFileEqual(expected_result, get_results, comment="#")


    def test_ordered_circular(self):
        # genetic organization of test_3.fasta
        # gene       abc    mfp    omf    omf    abc    gspd
        # gene id   01397  01398  01548  01562  01399  01400
        # pos        8      9      19     27     37     38
        # clst                 ]               [
        # syst (abc,37),  (gspd, 38), (abc,2), (mfp,3)

        expected_result_dir = self.find_data("functional_test_ordered_circular")
        # TODO how to specify multi_loci = false when multi_loci =True is set in xml
        args = "--db-type ordered_replicon " \
               "--replicon-topology circular  " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-simple-exch " \
               "-o {out_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.rejected_clusters):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertFileEqual(expected_result, get_results, comment="#")


    def test_ordered_linear(self):
        # genetic organization of test_3.fasta
        # gene       abc    mfp    omf    omf    abc    gspd
        # gene id   01397  01398  01548  01562  01399  01400
        # pos        8      9      19     27     37     38
        # clst    [            ]               [           ]
        # syst  no system

        expected_result_dir = self.find_data("functional_test_ordered_linear")
        # TODO how to specify multi_loci = false when multi_loci =True is set in xml
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear  " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-simple-exch " \
               "-o {out_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.rejected_clusters):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertFileEqual(expected_result, get_results, comment="#")


    def test_ordered_multi_system(self):
        # genetic organization of test_1.fasta
        #
        # gene       omf    mfp    abc    mfp    abc    gspd   omf    omf    omf
        # gene id   01360  01361  01397  01398  01399  01400  01506  01548  01562
        # pos         2      3     11     12     13      14    23     32     46
        # clst      [         ]   [                        ]  [  ]   [  ]   [  ]
        # syst                    [abc    mfp    abc    gspd   omf    omf    omf]

        expected_result_dir = self.find_data("functional_test_ordered_multi_system")
        # TODO how to specify multi_loci = false when multi_loci =True is set in xml
        args = "--db-type ordered_replicon " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-multi-syst-exch " \
               "-o {out_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.rejected_clusters):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertFileEqual(expected_result, get_results, comment="#")


    def test_ordered_multi_system_loner_in_clust(self):
        # genetic organization of test_2.fasta
        #
        # gene       abc    mfp    abc    gspd   omf    omf    omf
        # gene id   01397  01398  01399  01400  01506  01548  01562
        # pos        8      9      10      11    13     29     43
        # clst     [                               ]   [  ]   [  ]
        # syst     [abc    mfp    abc    gspd   omf]

        expected_result_dir = self.find_data("functional_test_ordered_multi_system_loner_in_clust")
        # TODO how to specify multi_loci = false when multi_loci =True is set in xml
        args = "--db-type ordered_replicon " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-multi-syst-exch " \
               "-o {out_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.rejected_clusters):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertFileEqual(expected_result, get_results, comment="#")


    def test_ordered_multi_loci(self):
        # genetic organization of test_4.fasta
        #
        # gene       abc    mfp    abc    gspd   omf    omf
        # gene id   01397  01398  01399  01400  01548  01562
        # pos        6      7      14      15    26     40
        # clst     [         ]   [           ]
        # syst    abc, mfp, abc, gspd, omf, omf

        expected_result_dir = self.find_data("functional_test_ordered_multi_loci")
        args = "--db-type ordered_replicon " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-simple-exch " \
               "-o {out_dir} " \
               "--multi-loci functional/T12SS-simple-exch " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"

        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.rejected_clusters):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertFileEqual(expected_result, get_results, comment="#")


    def test_ordered_single_loci(self):
        # genetic organization of test_4.fasta
        #
        # gene       abc    mfp    abc    gspd   omf    omf
        # gene id   01397  01398  01399  01400  01548  01562
        # pos        6      7      14      15    26     40
        # clst     [         ]   [           ]
        # syst    no system

        expected_result_dir = self.find_data("functional_test_ordered_single_loci")
        args = "--db-type ordered_replicon " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-simple-exch " \
               "-o {out_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"

        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.rejected_clusters):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertFileEqual(expected_result, get_results, comment="#")

    def test_degenerated_systems(self):
        # genetic organization of test_4.fasta
        #
        # gene      mfp   gspd
        # gene id  01398  01400
        # pos        7      15
        # syst    abc    mfp
        # inter_gene_max_space="8"
        expected_result_dir = self.find_data("functional_test_degenerated_systems")
        args = "--db-type ordered_replicon " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional  degenerated_systems " \
               "-o {out_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.rejected_clusters):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertFileEqual(expected_result, get_results, comment="#")


    def test_uncomplete_degenerated_systems(self):
        # genetic organization of test_4.fasta
        #
        # gene      mfp   gspd
        # gene id  01398  01400
        # pos        7      15
        # syst    abc    mfp
        # inter_gene_max_space="5"
        expected_result_dir = self.find_data("functional_test_uncomplete_degenerated_systems")
        args = "--db-type ordered_replicon " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional  uncomplete_degenerated_systems " \
               "-o {out_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.rejected_clusters):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertFileEqual(expected_result, get_results, comment="#")


    def test_unordered(self):
        # genetic organization of test_4.fasta
        #
        # gene       abc    mfp    abc    gspd   omf    omf
        # gene id   01397  01398  01399  01400  01548  01562
        # pos        6      7      14      15    26     40
        # syst    abc    mfp    abc    gspd   omf    omf
        expected_result_dir = self.find_data("functional_test_unordered")
        args = "--db-type unordered " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-simple-exch " \
               "-o {out_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_systems_txt,
                          self.uncomplete_systems):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertFileEqual(expected_result, get_results, comment="#")


    def test_working_dir_exists(self):
        args = f"--sequence-db {self.find_data('base', 'one_replicon.fasta')} " \
               "--db-type ordered_replicon " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS " \
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
               "-m functional T12SS " \
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
               "-m functional T12SS " \
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
               "-m functional T12SS " \
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
               "-m functional Unknown_model " \
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
        # print("\n############################################")
        # print(args)
        # print("##############################################")
        macsyfinder.main(args=args.split(),
                         loglevel='ERROR'
                         )
