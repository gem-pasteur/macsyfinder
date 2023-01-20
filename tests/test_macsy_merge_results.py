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

import tempfile
import shutil
import os
import argparse
import sys
import colorlog
import unittest
import platform
import re

from tests import MacsyTest
from macsypy.scripts import macsy_merge_results


class TestMerge(MacsyTest):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.test_dir = os.path.join(self.tmpdir, 'macsyfinder_test_merge')
        os.mkdir(self.test_dir)

        self.args = argparse.Namespace()

    def tearDown(self):
        try:
            shutil.rmtree(self.tmpdir)
        except:
            pass
        # some function in macsydata script suppress the traceback
        # but without traceback it's hard to debug test :-(
        sys.tracebacklimit = 1000  # the default value

        logger = colorlog.getLogger('macsypy')
        for h in logger.handlers[:]:
            logger.removeHandler(h)

    def test_get_warning(self):

        header = """# parallel_msf 20230113.dev
# merged best_solution.tsv
# Systems found:
"""
        warning ="""# 
# WARNING: The replicon 'GCF_000006765' has been SKIPPED. Cannot be solved before timeout.
#
"""
        first_line = "replicon\thit_id\t...\n"

        path = os.path.join(self.tmpdir, 'fake_result')
        with open(path, 'w') as fake_file:
            print(header + first_line, file=fake_file)
        got_warning = macsy_merge_results.get_warning(path)
        self.assertEqual(got_warning, [])

        with open(path, 'w') as fake_file:
            print(header + warning + first_line, file=fake_file)
        got_warning = macsy_merge_results.get_warning(path)
        self.assertEqual(got_warning, ["# WARNING: The replicon 'GCF_000006765' has been SKIPPED. Cannot be solved before timeout.\n"])


    def test_functional_merge(self):
        # res_1 have hits and systems
        # res_2 have hits and systems
        res_1 = 'results_1'
        src_res_1 = self.find_data('data_set', res_1)
        dest_res_1 = os.path.join(self.test_dir, res_1)
        res_2 = 'results_2'
        src_res_2 = self.find_data('data_set', res_2)
        dest_res_2 = os.path.join(self.test_dir, res_2)
        shutil.copytree(src_res_1, dest_res_1)
        shutil.copytree(src_res_2, dest_res_2)
        merge_dir = os.path.join(self.test_dir, 'merged_results')
        os.mkdir(merge_dir)
        cmd = f"macsymerge -o {merge_dir} --mute {dest_res_1} {dest_res_2}"
        macsy_merge_results.main(args=cmd.split()[1:])

        results_files = ('best_solution.tsv',
                         'best_solution_summary.tsv',
                         'best_solution_loners.tsv',
                         'best_solution_multisystems.tsv',
                         'all_best_solutions.tsv',
                         'all_systems.tsv',
                         'all_systems.txt',
                         'rejected_candidates.txt')

        for res in results_files:
            with self.subTest(res):
                expected = self.find_data('data_set', 'merged_results', f"merged_{res}")
                recieved = os.path.join(merge_dir, f"merged_{res}")
                self.assertFileEqual(expected, recieved, comment='#')


    def test_functional_merge_one_results_empty(self):
        # res_1 has hits
        # res_2 has no hits and no systems
        res_1 = 'results_no_hits'
        src_res_1 = self.find_data('data_set', res_1)
        dest_res_1 = os.path.join(self.test_dir, res_1)
        res_2 = 'results_2'
        src_res_2 = self.find_data('data_set', res_2)
        dest_res_2 = os.path.join(self.test_dir, res_2)
        shutil.copytree(src_res_1, dest_res_1)
        shutil.copytree(src_res_2, dest_res_2)
        merge_dir = os.path.join(self.test_dir, 'merged_results')
        os.mkdir(merge_dir)
        cmd = f"macsymerge -o {merge_dir} --mute -qq {dest_res_1} {dest_res_2}"
        macsy_merge_results.main(args=cmd.split()[1:])

        results_files = ('best_solution.tsv',
                         'best_solution_summary.tsv',
                         'best_solution_loners.tsv',
                         'best_solution_multisystems.tsv',
                         'all_best_solutions.tsv',
                         'all_systems.tsv',
                         'all_systems.txt',
                         'rejected_candidates.txt',
                         'rejected_candidates.tsv')

        for res in results_files:
            with self.subTest(res):
                expected = self.find_data('data_set', 'merged_results_one_res_no_hits', f"merged_{res}")
                recieved = os.path.join(merge_dir, f"merged_{res}")
                self.assertFileEqual(expected, recieved, comment='#')


    def test_functional_merge_results_empty(self):
        # res_1 has no hits and no systems
        # res_2 has no hits and no systems
        res_1 = 'results_1_no_hits'
        src_res_1 = self.find_data('data_set', 'results_no_hits')
        dest_res_1 = os.path.join(self.test_dir, res_1)
        res_2 = 'results_2_no_hits'
        dest_res_2 = os.path.join(self.test_dir, res_2)
        shutil.copytree(src_res_1, dest_res_1)
        shutil.copytree(src_res_1, dest_res_2)
        merge_dir = os.path.join(self.test_dir, 'merged_results')
        os.mkdir(merge_dir)
        cmd = f"macsy_merge_results -o {merge_dir} --mute -qq {dest_res_1} {dest_res_2}"
        macsy_merge_results.main(args=cmd.split()[1:])

        results_files = ('best_solution.tsv',
                         'best_solution_summary.tsv',
                         'best_solution_loners.tsv',
                         'best_solution_multisystems.tsv',
                         'all_best_solutions.tsv',
                         'all_systems.tsv',
                         'all_systems.txt',
                         'rejected_candidates.txt',
                         'rejected_candidates.tsv')

        for res in results_files:
            with self.subTest(res):
                expected = self.find_data('data_set', 'merged_results_no_hits', f"merged_{res}")
                recieved = os.path.join(merge_dir, f"merged_{res}")
                self.assertFileEqual(expected, recieved, comment='#')


    def test_functional_merge_out_error(self):
        res_1 = 'results_1'
        src_res_1 = self.find_data('data_set', res_1)
        dest_res_1 = os.path.join(self.test_dir, res_1)
        res_2 = 'results_no_hits'
        src_res_2 = self.find_data('data_set', res_2)
        dest_res_2 = os.path.join(self.test_dir, res_2)
        merge_dir = os.path.join(self.test_dir, 'merged_results')

        ###########################
        # the -o option is a file
        ###########################
        open(merge_dir, 'w').close()
        cmd = f"macsy_merge_results -o {merge_dir} {dest_res_1} {dest_res_2}"
        with self.catch_io(out=True):
            # we cannot test the log message here
            # because the logger are init when main is called
            # after that the context is establish
            # so when the wrapper is called the handler cannot be substitute by the fake
            # but we can catch stdout
            with self.assertRaises(IOError):
                macsy_merge_results.main(args=cmd.split()[1:], log_level='WARNING')
            stdout = sys.stdout.getvalue().strip()
        self.assertEqual(self.remove_red_ansi_color(stdout),
                         f'{merge_dir} is not a directory')

        os.unlink(merge_dir)


    @unittest.skipIf(platform.system() == 'Windows' or os.getuid() == 0, 'Skip test on Windows or if run as root')
    def test_functional_merge_out_not_writable(self):
        res_1 = 'results_1'
        src_res_1 = self.find_data('data_set', res_1)
        dest_res_1 = os.path.join(self.test_dir, res_1)
        res_2 = 'results_no_hits'
        src_res_2 = self.find_data('data_set', res_2)
        dest_res_2 = os.path.join(self.test_dir, res_2)
        merge_dir = os.path.join(self.test_dir, 'merged_results')
        ######################################
        # the -o option is a dir NOT writable
        #####################################
        os.mkdir(merge_dir, mode=0o444)
        cmd = f"macsy_merge_results -o {merge_dir} {dest_res_1} {dest_res_2}"
        try:
            with self.catch_io(out=True):
                with self.assertRaises(IOError):
                    macsy_merge_results.main(args=cmd.split()[1:], log_level='WARNING')
                stdout = sys.stdout.getvalue().strip()
            self.assertEqual(self.remove_red_ansi_color(stdout),
                             f'{merge_dir} is not writable')
        finally:
            shutil.rmtree(merge_dir)

        ######################################
        # the -o option parent is NOT writable
        ######################################
        os.chmod(self.test_dir, mode=0o444)
        with self.catch_io(out=True):
            with self.assertRaises(IOError):
                macsy_merge_results.main(args=cmd.split()[1:], log_level='WARNING')
            stdout = sys.stdout.getvalue().strip()
        self.assertEqual(self.remove_red_ansi_color(stdout),
                         f"Cannot create {merge_dir} : [Errno 13] Permission denied: '{merge_dir}'")


    def test_merge_and_reindex_invalid_solid(self):
        res1 = self.find_data('data_set', 'results_2', 'all_best_solutions.tsv')
        res2 = self.find_data('data_set', 'results_invalid_output', 'all_best_solutions.tsv')
        merge_dir = os.path.join(self.test_dir, 'merged_results')

        with self.assertRaises(ValueError) as ctx:
            macsy_merge_results.merge_and_reindex([res1, res2], merge_dir, "Systems",
                                                  skip_until=lambda l: l.startswith('sol_id'),
                                                  comment='#')
        self.assertEqual(str(ctx.exception),
                         "Cannot reindex int(GCF_000006845) + 1:"
                         " invalid literal for int() with base 10: 'GCF_000006845'")
