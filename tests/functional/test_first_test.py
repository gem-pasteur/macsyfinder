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


import shutil
import tempfile
import platform
import os
from subprocess import Popen
import json

from tests import MacsyTest
from macsypy.utils import which

class Test(MacsyTest):

    def setUp(self):
        if 'MACSY_HOME' in os.environ:
            self.macsy_home = os.environ['MACSY_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.macsy_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..' '..')))
        self.tmp_dir = tempfile.gettempdir()



    def tearDown(self):
        try:
            shutil.rmtree(self.out_dir)
        except:
            pass

    def test_basic_run(self):
        """
        test if returncode of masyfinder is 0 and if json == the expected json
        test the presence of T9SS T3SS T4SS_typeI systems
        with test_aesu.fa sequence db in gembase format
        """
        self.out_dir = os.path.join(self.tmp_dir, 'macsyfinder_test_basic_run')
        os.makedirs(self.out_dir)
        bin = os.path.join(self.macsy_home, 'bin', 'macsyfinder') if self.local_install else which('macsyfinder')

        # bin/macsyfinder
        # --sequence-db /home/bneron/Projects/macsyfinder/data/base/test_aesu.fa
        # --db-type gembase
        # -d /home/bneron/Projects/macsyfinder/data/TXSScan/DEF_TXSS
        # -p /home/bneron/Projects/macsyfinder/data/TXSScan/profiles_TXSS
        #

        command = "{bin} --def={def_dir} --profile-dir={profiles} --out-dir={out_dir} --sequence-db={seq_db} --db-type=gembase {systems}".format(
                    bin=bin,
                    out_dir=self.out_dir,
                    def_dir=os.path.join(self._data_dir, 'data_set_1', 'def'),
                    profiles=os.path.join(self._data_dir, 'data_set_1', 'profiles'),
                    seq_db=os.path.join(self._data_dir, 'base', 'test_aesu.fa'),
                    systems="T9SS T3SS T4SS_typeI",
                    )
        if not bin:
            raise RuntimeError('macsyfinder not found, macsyfinder must be either in your path or MACSY_HOME must be defined')
        # I redirect stdout and stderr in dev null I don't want them on screen
        # I cannot redirect them in output directory as --out-dir expect a non existing directory or an empty one
        # but Popen need to have a file as argument of stdout/err
        try:
            macsy_process = Popen(command,
                                  shell=True,
                                  stdin=None,
                                  stdout=open(os.devnull, 'w'),
                                  stderr=open(os.devnull, 'w'),
                                  close_fds=False
                                  )
        except Exception as err:
            msg = "macsyfinder execution failed: command = {0} : {1}".format(command, err)
            print
            print msg
            raise err

        macsy_process.wait()

        self.assertEqual(macsy_process.returncode, 0, "macsyfinder finished with non zero exit code: {0}".format(macsy_process.returncode))

        expected_result_path = os.path.join(self._data_dir, 'outputs_control', 'basic_run', 'results.macsyfinder.json')
        with open(expected_result_path) as expected_result_file:
            expected_result_json = json.load(expected_result_file)

        test_result_path = os.path.join(self.out_dir, 'results.macsyfinder.json')
        with open(test_result_path) as test_result_file:
            test_result_json = json.load(test_result_file)

        self.assertEqual(expected_result_json, test_result_json)
