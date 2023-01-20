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


import shutil
import tempfile
import os
import sys
import inspect
import unittest
import itertools

import colorlog

from tests import MacsyTest, which
from macsypy.scripts import macsyfinder
from macsypy.error import OptionError
from macsypy.system import System, AbstractUnordered, RejectedCandidate


class Test(MacsyTest):

    @classmethod
    def setUpClass(cls) -> None:
        cls._index_dir = os.path.join(tempfile.gettempdir(), 'test_macsyfinder_index')
        if not os.path.exists(cls._index_dir):
            os.makedirs(cls._index_dir)

    @classmethod
    def tearDownClass(cls) -> None:
        if os.path.exists(cls._index_dir):
            shutil.rmtree(cls._index_dir)

        logger = colorlog.getLogger('macsypy')
        for h in logger.handlers[:]:
            logger.removeHandler(h)

    def setUp(self):
        self.tmp_dir = tempfile.gettempdir()
        # reset System, AbstractUnordered internal id to have predictable results (Systems, ...) id
        # it's works only if there is only one replicon
        # for gembase the order is not guarantee

        System._id = itertools.count(1)
        AbstractUnordered._id = itertools.count(1)

        self.all_systems_tsv = "all_systems.tsv"
        self.all_systems_txt = "all_systems.txt"
        self.all_best_solutions = "all_best_solutions.tsv"
        self.best_solution = "best_solution.tsv"
        self.summary = "best_solution_summary.tsv"
        self.rejected_candidates_txt = "rejected_candidates.txt"
        self.rejected_candidates_tsv = "rejected_candidates.tsv"
        self.uncomplete_systems = "uncomplete_systems.txt"
        self.loners = "best_solution_loners.tsv"
        self.multisystems = "best_solution_multisystems.tsv"


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
               "--models TFF-SF all " \
               "--out-dir={out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"

        self._macsyfinder_run(args)
        #  all_systems_tsv, all_best_solutions, best_solution
        # provides system with non predictable id
        # and the order between equivalent solutions is not predictable
        # so the output is not predictable
        # even the biological meaning is good
        # so I disabled them from test until I found a fix
        # we steel check the overall system found (summary)
        # the loners and multisystems hits
        for file_name in (# self.all_systems_tsv,
                          # self.all_best_solutions,
                          # self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        expected_result = self.find_data(expected_result_dir, self.rejected_candidates_txt)
        get_results = os.path.join(self.out_dir, self.rejected_candidates_txt)
        self.assertFileEqual(expected_result, get_results, comment="#")

    @unittest.skipIf(not which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_timeout(self):
        """

        """
        expected_result_dir = self.find_data("functional_test_timeout")
        args = "--db-type=gembase " \
               f"--models-dir={self.find_data('models')} " \
               "--models TFF-SF all " \
               "--out-dir={out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path " \
               "--timeout 2 "

        with self.catch_io(out=True, err=True):
            with self.catch_log() as log:
                self._macsyfinder_run(args)
                log_msg = log.get_value().strip()
        self.assertEqual(log_msg,
                         """Timeout is over. Aborting
The THHY002.0321.00001.C001 cannot be solved in time skip it!
The replicon THHY002.0321.00001.C001 cannot be solved before timeout. SKIP IT.""")

        # test only if THHY002.0321.00001.C001 has_been skipped
        for file_name in (self.summary,):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

                with open(get_results) as results:
                    lines = results.readlines()
                self.assertEqual(lines[4],
                                 "# WARNING: The replicon 'THHY002.0321.00001.C001' has been SKIPPED. Cannot be solved before timeout.\n")

    def test_only_loners(self):
        # genetic organization of MOBP1_twice.fast
        # gene        MOBP1          MOBP1
        # gene_id     0832           0885
        # pos          8               19
        #             [ ]             [  ]
        expected_result_dir = self.find_data("functional_test_only_loners")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m test_loners MOB_cf_T5SS " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")


    def test_ordered_circular(self):
        # genetic organization of test_3.fasta
        # gene       abc    mfp    omf    omf    abc    gspd
        # gene id   01397  01398  01548  01562  01399  01400
        # pos        2      3      19     27     37     38
        # clst                 ]               [
        # syst (abc,2), (mfp,3), (abc,37), (gspd, 38)
        # in T12SS-simple-exch omf is not a loner

        expected_result_dir = self.find_data("functional_test_ordered_circular")
        args = "--db-type ordered_replicon " \
               "--replicon-topology circular " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-simple-exch " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")


    def test_ordered_linear(self):
        # genetic organization of test_3.fasta
        # gene       abc    mfp    omf    omf    abc    gspd
        # gene id   01397  01398  01548  01562  01399  01400
        # pos        2      3      19     27     37     38
        # clst    [  M      A  ]    M     M     [ M      A   ]
        # syst  no system
        # in T12SS-simple-exch omf is not a loner
        expected_result_dir = self.find_data("functional_test_ordered_linear")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-simple-exch " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")


    def test_ordered_1_cluster_3_loners(self):
        # genetic organization of test_1.fasta
        #
        # gene       omf    mfp    abc    mfp    abc    gspd   omf    omf    omf
        # gene id   01360  01361  01397  01398  01399  01400  01506  01548  01562
        # pos         2      3     11     12     13      14    23     32     46
        # clst      [ML      M]   [ M      M      M      ME]  [ML]   [ML]   [ML]
        # syst                    [abc    mfp    abc    gspd                omf] with equivalent for omf46 [omf23 omf32]
        # score                     1      1      0      0                   .7  = 2.7
        # loners      X                                        omf    omf    omf
        # omf2 colocate with mfp3  => not considerde as Loner
        expected_result_dir = self.find_data("functional_test_ordered_1_cluster_3_loners")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-loner " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")

    def test_ordered_1_cluster_and_clusters_of_loners(self):
        # genetic organization of test_15.fasta
        #
        # gene       abc    mfp    abc   gspd    omf   omf    omf    omf
        # gene id   01397  01398  01399  01400  01506 01360  01548  01562
        # pos        9      10     11     12      21   22     44     45
        # clst     [ M      M      M      ME ]   [ML   ML]   [ML     ML]
        # syst      [abc    mfp    abc    gspd]    2 clusters of 2 loners but considered as 4 loners
        # score       1      1      0      0                   .7  = 2.7
        # loners                                  omf   omf    omf    omf
        # omf 23 24 and 30 31 are clusters of loners with same gene  => consider as Loner
        expected_result_dir = self.find_data("functional_test_ordered_1_cluster_and_clusters_of_loners")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-loner " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")

    def test_ordered_2_clusters_3_loners(self):
        # genetic organization of test_5.fasta
        #
        # gene       omf    mfp    abc    mfp    abc    gspd   omf    omf    omf     abc    mfp   gspd
        # gene id   01360  01361  01397  01398  01399  01400  01506  01548  01562   01150  01361  0409
        # pos         2      3     11     12     13      14    23     32     46       55     56    57
        # clst      [ ML     ML]  [M      M      M       ME]  [ML]   [ML]   [ML]    [ M       M    ME]
        # syst                    [abc    mfp    abc    gspd] [omf]  [omf]  [omf]   [abc     mfp  gspd ]
        # 2 systems [abc    mfp    abc    gspd  omf46] with equivalent for omf46 [omf23 omf32]
        # score       1      1      0      0      .7  = 2.7
        #           [abc     mfp    gspd   omf46] with equivalent for omf46 [omf23 omf32]
        #             1       1      0     .7        =  2.7
        expected_result_dir = self.find_data("functional_test_ordered_2_clusters_3_loners")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-loner " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")


    def test_ordered_2_clusters_1_loner(self):
        # genetic organization of test_6.fasta
        #
        # gene       mfp    abc    mfp    abc    gspd   omf    abc    mfp   gspd
        # gene id   01361  01397  01398  01399  01400  01506  01150  01361  00409
        # pos         2      10     11     12     13     22     51    52     53
        # clst        M    [ M       M      M     ME]   [ML]   [M      M     ME]
        # syst              [abc    mfp    abc   gspd]  [omf]  [abc   mfp   gspd ]
        # 2 systems [abc    mfp    abc    gspd  omf22] with warning 1 loner 2 systems
        # score       1      1      0      0     .7   = 2.7
        #           [abc     mfp    gspd   omf22] with with warning 1 loner 2 systems
        #             1       1      0      .7   =  2.7
        expected_result_dir = self.find_data("functional_test_ordered_2_clusters_1_loner")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-loner-exch " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")


    def test_ordered_1_loner_in_clust(self):
        # genetic organization of test_2.fasta
        #
        # gene       abc    mfp    abc    gspd   omf    omf    omf
        # gene id   01397  01398  01399  01400  01506  01548  01562
        # pos        8      9      10      11    13     29     43
        # clst     [ M      M       M      ME     M]   [ML]   [ML]
        # syst     [abc    mfp    abc    gspd   omf]
        # score      1      1      0      0      1   = 3.0

        expected_result_dir = self.find_data("functional_test_ordered_1_loner_in_clust")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-loner-exch " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")


    def test_ordered_1_loner_exch_in_clust(self):
        # genetic organization of test_8.fasta
        #
        # gene       abc    mfp    abc    gspd   gspf    omf    omf
        # gene id   01397  01398  01399  01400  02599  01548  01562
        # pos        8      9      10      11    13     29     43
        # clst     [ M      M      M       ME    MLE]   [ML]  [ML]
        # syst     [abc    mfp    abc    gspd   gspf]
        # score      1      1      0      0      0.7  = 2.7

        expected_result_dir = self.find_data("functional_test_ordered_1_loner_exch_in_clust")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-loner-exch " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")


    def test_ordered_1_clusters_3_loners_w_exchangeable(self):
        # genetic organization of test_7.fasta
        #
        # gene       omf    mfp    abc    mfp    abc    gspd   omf    omf    gspF
        # gene id   01360  01361  01397  01398  01399  01400  01506  01548  02599
        # pos         2      3     11     12     13      14    23     32     46
        # clst      [ML      M ]  [M       M      M      ME]  [ML]    [ML]  [ML]
        # syst                    [abc    mfp    abc    gspd   omf] with equivalent for omf23 [omf32, gspF46]
        # score                     1      1      0      0      .7  = 2.7
        # gspF46 have a score (465.3). omf (90, 111.5, 87) but it's an exhangeable so it canot be the "best loner"

        expected_result_dir = self.find_data("functional_test_ordered_1_cluster_3_loners_w_exchangeable")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-loner-exch " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")


    def test_ordered_multi_loci(self):
        # genetic organization of test_4.fasta
        #
        # gene       abc    mfp    abc    gspd   omf    omf
        # gene id   01397  01398  01399  01400  01548  01562
        # pos        6      7      14      15    26     40
        # clst      [M      A]    [M       ME]
        # 1 syst   [abc6, mfp7    abc14, gspd15]
        # score      1      .5     1      .7  = 1.5 + 1.7 = 3.2 - (1 * 1.5 redundancy penalty) = 2.7
        # in T12SS-simple-exch omf is not a loner

        expected_result_dir = self.find_data("functional_test_ordered_multi_loci")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-simple-exch " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               "--multi-loci functional/T12SS-simple-exch " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"

        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")


    def test_ordered_single_loci(self):
        # genetic organization of test_4.fasta
        #
        # gene       abc    mfp    abc    gspd   omf    omf
        # gene id   01397  01398  01399  01400  01548  01562
        # pos        6      7      14      15    26     40
        # clst      [M      A]     [M      ME]
        # syst    no system

        expected_result_dir = self.find_data("functional_test_ordered_single_loci")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-simple-exch " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"

        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")


    def test_ordered_multi_system(self):
        # genetic organization of test_13.fasta
        #
        # gene      abc    mfp    gspd   omf   gspf     abc    omf    gspd   omf
        # gene id  01397  01398  01400  01360  02599   01399  01506  00409  01562
        # pos       8       9      19     20    21      34     35      36     43
        # clst    [ M       A ]  [ M     M_MS    A ]   [M     M_MS     M]
        # syst 1               [gspd19, omf20, gspf21][abc34, omf35, gspd36]
        # score                     1       1      .5                            = 2.5
        # score                                        1       1       1         = 3.0
        # syst 2   [abc8    mfp9] + [omf20]
        # score      1       .5       .7                                         = 2.2
        # The multi system is in system
        # So it can be used for other clusters to form new occurrence

        expected_result_dir = self.find_data("functional_test_ordered_multi_system")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-multisystem " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"

        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")


    def test_ordered_multi_system_out_system(self):
        # genetic organization of test_12.fasta
        #
        # gene      abc    mfp    gspd   omf    omf    omf
        # gene id  01397  01398  01400  01360  01506  01562
        # pos       8       9      19     20    33      39
        # clst    [ M       A ]  [ M     M_MS ]
        # syst    no system
        # The multi system is not in system
        # So it cannot be used to build new systems

        expected_result_dir = self.find_data("functional_test_ordered_multi_system_out_system")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-multisystem " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"

        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")


    def test_ordered_multi_model(self):
        # genetic organization of test_14.fasta
        #
        # gene      omf    mfp    abc   gspd    pilB    pilW
        # gene id  01506  01398  01399  01400  000980  025680
        # pos        8      9      10    11      12      13
        # clst     [        model C         ]
        #                                [      model D       ]
        # syst     [omf8, mfp9, abc10, gspd11]                   model C score = 2.5
        #                             [gspd11, pilB12, pilW13]   model D score = 2.0
        # gspd is multi_model in C_multi_model and D_multi_model
        # So the 2 systems are in the solution
        expected_result_dir = self.find_data("functional_test_ordered_multi_model")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional  C_multi_model D_multi_model " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"

        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")


    def test_ordered_only_one_multi_model(self):
        # genetic organization of test_14.fasta
        #
        # gene      omf    mfp    abc   gspd    pilB    pilW
        # gene id  01506  01398  01399  01400  000980  025680
        # pos        8      9      10    11      12      13
        # clst     [        model C         ]
        #                                [      model D       ]
        # syst     [omf8, mfp9, abc10, gspd11]                   model C score = 2.5
        #                             [gspd11, pilB12, pilW13]   model D score = 2.0
        # gspd is multi_model in D_multi_model but NOT in C_no_multi_model
        # So the 2 systems are excllusive. msf pick the best system for the best solution (sys c score 2.5)

        expected_result_dir = self.find_data("functional_test_ordered_only_one_multi_model")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional  C_no_multi_model D_multi_model " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"

        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")


    def test_degenerated_systems(self):
        # genetic organization of test_4.fasta
        #
        # gene      mfp   gspd
        # gene id  01398  01400
        # pos        7      15
        # syst     [mfp    gspd]
        # score    [M        A ]   = 1.5
        # inter_gene_max_space="8"
        expected_result_dir = self.find_data("functional_test_degenerated_systems")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional degenerated_systems " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")


    def test_uncomplete_degenerated_systems(self):
        # genetic organization of test_4.fasta
        #
        # gene      mfp(Man)   gspd(acce)
        # gene id    01398        01400
        # pos         7            15
        # syst       mfp
        # score      [1]                   = 1
        # inter_gene_max_space="5"
        expected_result_dir = self.find_data("functional_test_uncomplete_degenerated_systems")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional  uncomplete_degenerated_systems " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")

    def test_2_systems_not_compatible(self):
        # genetic organization of test_9.fasta
        #
        # gene      abc   mfp    gspd   omf    gspf
        # gene id  01397 01398  01400  01506  02599
        # pos       9     10      13    15     17
        # syst A   [abc    mfp   gspd]
        #            A      A     M
        # score      .5    .5     1   = 2.0
        # syst B                [gspd   omf   gspf]
        #                         M      A     A
        # score                   1      .5    .5  = 2.0
        # 2 systems not compatible (share gspd)
        # so 2 solutions

        expected_result_dir = self.find_data("functional_test_2_systems_not_compatible")
        args = "--db-type ordered_replicon " \
               "--replicon-topology linear " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional  A B " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        self._macsyfinder_run(args)

        for file_name in (self.all_systems_tsv,
                          self.all_best_solutions,
                          self.best_solution,
                          self.loners,
                          self.multisystems,
                          self.summary,
                          self.rejected_candidates_tsv):
            with self.subTest(file_name=file_name):
                expected_result = self.find_data(expected_result_dir, file_name)
                get_results = os.path.join(self.out_dir, file_name)
                self.assertTsvEqual(expected_result, get_results, comment="#", tsv_type=file_name)

        self.assertFileEqual(self.find_data(expected_result_dir, self.rejected_candidates_txt),
                             os.path.join(self.out_dir, self.rejected_candidates_txt), comment="#")


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
               f"--index-dir {self._index_dir} " \
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

    def test_unordered_only_forbidden(self):
        # genetic organization of test_10.fasta
        #
        # gene       omf
        # gene id   01506
        # pos        17
        # syst    no Systems
        expected_result_dir = self.find_data("functional_test_unordered_only_forbidden")
        args = "--db-type unordered " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-forbidden " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
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


    def test_unordered_no_hits(self):
        # genetic organization of test_11.fasta
        #
        # gene
        # gene id
        # pos
        # syst    no Systems
        expected_result_dir = self.find_data("functional_test_unordered_no_hits")
        args = "--db-type unordered " \
               f"--models-dir {self.find_data('models')} " \
               "-m functional T12SS-forbidden " \
               "-o {out_dir} " \
               f"--index-dir {self._index_dir} " \
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


    def test_index_dir(self):
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
               f"--index-dir {self._index_dir} " \
               f"--previous-run {expected_result_dir} " \
               "--relative-path"
        sequences_dir = self.find_data('base')
        self._macsyfinder_run(args)

        self.assertTrue(os.path.exists(os.path.join(self._index_dir, "test_3.fasta.idx")))


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
               "--index-dir {out_dir} " \
               "-o {out_dir}"

        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_model_unkwon')
        os.makedirs(self.out_dir)

        args = args.format(out_dir=self.out_dir)
        with self.assertRaises(ValueError) as ctx:
            macsyfinder.main(args=args.split(), loglevel='ERROR')
        self.assertEqual(str(ctx.exception),
                         "Unknown_model does not match with any definitions")


    def test_cfg_n_previous_run(self):
        args = f"--cfg-file foo --previous-run bar " \
               "-o {out_dir}"

        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_cfg_n_previous_run')
        os.makedirs(self.out_dir)

        args = args.format(out_dir=self.out_dir)

        real_exit = sys.exit
        sys.exit = self.fake_exit
        try:
            with self.catch_io(out=True):
                with self.assertRaises(TypeError) as ctx:
                    macsyfinder.main(args=args.split(), loglevel='ERROR')
            self.assertEqual(str(ctx.exception), '2')
        finally:
            sys.exit = real_exit


    def _macsyfinder_run(self, args_tpl):
        try:
            # get the name of the calling function
            test_name = inspect.stack()[1].function
            self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_{}'.format(test_name))
            os.makedirs(self.out_dir)
            args = args_tpl.format(out_dir=self.out_dir)
            # print("\n############################################")
            # print(args)
            # print("##############################################")
            System._id = itertools.count(1)
            RejectedCandidate._id = itertools.count(1)
            macsyfinder.main(args=args.split(),
                             loglevel='ERROR'
                             )
        except Exception as err:
            import traceback
            traceback.print_exc()
            print(err)
            # keep the directory
            i = 0
            new_name = self.out_dir + f'_keep_{err.__class__.__name__}_{i}'
            while os.path.exists(new_name):
                i += 1
                new_name = self.out_dir + f'_keep_{err.__class__.__name__}_{i}'
                if i > 20:
                    break
            shutil.copytree(self.out_dir, new_name)
            shutil.copytree(self._index_dir, os.path.join(new_name, os.path.basename(self._index_dir)))
