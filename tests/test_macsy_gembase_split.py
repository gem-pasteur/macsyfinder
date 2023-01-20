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
import glob
import unittest
import platform

from tests import MacsyTest
from macsypy.scripts import macsy_gembase_split


class TestSplit(MacsyTest):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.test_dir = os.path.join(self.tmpdir, 'macsyfinder_test_split')
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


    def test_functional_split(self):
        gembase = self.find_data('base', 'gembase_short.fa')
        seq_dir = os.path.join(self.test_dir, 'base')
        cmd = f"macsy_gembase_split -o {seq_dir} --mute -qqq {gembase}"
        with self.catch_io(out=True):
            macsy_gembase_split.main(cmd.split()[1:])
            seq_files = sys.stdout.getvalue().strip().split()

        expected_dir = self.find_data('gembase_split')
        expected_out_files = glob.glob(os.path.join(expected_dir, f'*.fasta'))

        self.assertSetEqual(set([os.path.basename(f) for f in expected_out_files]),
                            set([os.path.basename(f) for f in seq_files])
                            )
        for file in [os.path.basename(f) for f in expected_out_files]:
            with self.subTest(file):
                expected_seq = os.path.join(expected_dir, file)
                split_seq = os.path.join(seq_dir, file)
                self.assertFileEqual(expected_seq, split_seq)


    def test_functional_split_out_error(self):
        gembase = self.find_data('base', 'gembase_short.fa')
        seq_dir = os.path.join(self.test_dir, 'base')

        ###########################
        # the -o option is a file
        ###########################
        open(seq_dir, 'w').close()
        cmd = f"macsy_gembase_split -o {seq_dir} {gembase}"
        with self.catch_io(out=True):
            # we cannot test the log message here
            # because the logger are init when main is called
            # after that the context is establish
            # so when the wrapper is called the handler cannot be substitute by the fake
            # but we can catch stdout
            with self.assertRaises(IOError):
                macsy_gembase_split.main(args=cmd.split()[1:], log_level='WARNING')
            stdout = sys.stdout.getvalue().strip()
        self.assertEqual(self.remove_red_ansi_color(stdout),
                         f'{seq_dir} is not a directory')
        os.unlink(seq_dir)


    @unittest.skipIf(platform.system() == 'Windows' or os.getuid() == 0, 'Skip test on Windows or if run as root')
    def test_functional_split_out_not_writable(self):
        gembase = self.find_data('base', 'gembase_short.fa')
        seq_dir = os.path.join(self.test_dir, 'base')

        ######################################
        # the -o option is a dir NOT writable
        #####################################

        os.mkdir(seq_dir, mode=0o444)
        cmd = f"macsy_gembase_split -o {seq_dir} {gembase}"
        try:
            with self.catch_io(out=True):
                with self.assertRaises(IOError):
                    macsy_gembase_split.main(args=cmd.split()[1:], log_level='WARNING')
                stdout = sys.stdout.getvalue().strip()
            self.assertEqual(self.remove_red_ansi_color(stdout),
                             f'{seq_dir} is not writable')
        finally:
            shutil.rmtree(seq_dir)

        ######################################
        # the -o option parent is NOT writable
        ######################################
        os.chmod(self.test_dir, mode=0o444)
        with self.catch_io(out=True):
            with self.assertRaises(IOError):
                macsy_gembase_split.main(args=cmd.split()[1:], log_level='WARNING')
            stdout = sys.stdout.getvalue().strip()
        self.assertEqual(self.remove_red_ansi_color(stdout),
                         f"Cannot create {seq_dir} : [Errno 13] Permission denied: '{seq_dir}'")
