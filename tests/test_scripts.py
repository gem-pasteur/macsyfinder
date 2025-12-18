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
import shutil
import tempfile
import os.path
import sys

import subprocess
from tests import MacsyTest


class Test_macsyfinder(MacsyTest):

    def test_cmde(self):
        # just test if the command msf_profile is available (there is no error in the pyproject.toml and generated wrapper
        p = subprocess.run("macsyfinder --help", shell=True, check=True, capture_output=True, text=True, encoding='utf8')
        self.assertEqual(p.returncode, 0)


class Test_msf_data(MacsyTest):

    def test_cmde(self):
        # The output of msf_data --help is formatted according to the width of the terminal
        # So it is not so obvious to have a reproducible results
        # just check if the command is available as in bioconda bot
        # the code itself is tested in other modules
        p = subprocess.run("msf_data --help", shell=True, check=True, capture_output=True, text=True, encoding='utf8')
        self.assertEqual(p.returncode, 0)


class Test_msf_profile(MacsyTest):

    def setUp(self):
        self._tmp_dir = tempfile.TemporaryDirectory(prefix='test_msf_scripts_')
        self.tmp_dir = self._tmp_dir.name

    def tearDown(self):
        self._tmp_dir.cleanup()


    def test_cmde(self):
        # just test if the command msf_profile is available (there is no error in the pyproject.toml and generated wrapper
        p = subprocess.run("msf_profile --help", shell=True, check=True, capture_output=True, text=True, encoding='utf8')
        self.assertEqual(p.returncode, 0)


    def test_previous_run(self):
        # test if the command msf_profile on a previous run
        run_name = 'functional_test_ordered_single_loci'
        prev_run = shutil.copytree(self.find_data(run_name),
                    os.path.join(self.tmp_dir, 'functional_test_ordered_single_loci'))

        try:
            p = subprocess.run(f"msf_profile {prev_run}",
                               shell=True,
                               check=True,
                               capture_output=True,
                               text=True,
                               encoding='utf8')
        except subprocess.CalledProcessError as err:
            print(err.stderr, sys.stderr)
            raise err
        self.assertEqual(p.returncode, 0, p.stdout)


class Test_msf_config(MacsyTest):

    def test_cmde(self):
        # just test if the command msf_profile is available (there is no error in the pyproject.toml and generated wrapper
        p = subprocess.run("msf_config --help", shell=True, check=True, capture_output=True, text=True, encoding='utf8')
        self.assertEqual(p.returncode, 0)


class Test_msf_split(MacsyTest):

    def test_cmde(self):
        # just test if the command msf_profile is available (there is no error in the pyproject.toml and generated wrapper
        p = subprocess.run("msf_split --help", shell=True, check=True, capture_output=True, text=True, encoding='utf8')
        self.assertEqual(p.returncode, 0)


class Test_msf_merge(MacsyTest):

    def test_cmde(self):
        # just test if the command msf_profile is available (there is no error in the pyproject.toml and generated wrapper
        p = subprocess.run("msf_merge --help", shell=True, check=True, capture_output=True, text=True, encoding='utf8')
        self.assertEqual(p.returncode, 0)
