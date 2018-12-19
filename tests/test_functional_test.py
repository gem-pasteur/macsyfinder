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
import os
from subprocess import Popen
import json
import unittest

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
        macsy_bin = 'macsyfinder'
        command = "{bin} --out-dir={out_dir} --sequence-db={seq_db} --db-type=gembase --models-dir={models_dir}" \
                  " --models {models}".format(bin=macsy_bin,
                                              out_dir=self.out_dir,
                                              models_dir=os.path.join(self._data_dir, 'data_set_1', 'models'),
                                              seq_db=os.path.join(self._data_dir, 'base', 'test_aesu.fa'),
                                              models="set_1 T9SS T3SS T4SS_typeI",
                                              )

        # I need to prepend the command by setsid because macsyfinder use killpg with group_id to terminated all
        # threads and subprocess when an error occurred in one hmmsearch. It's work fine but when
        # macsyfinder is launched by the tests.py script the kill group kill also the tests.py script
        # so we must run macsyfinder in a new process group
        # but setsid cmd does not exists on darwin so we provide one in utils
        setsid = self.setsid()
        try:
            macsy_process = Popen("{setsid} {cmd}".format(setsid=setsid,
                                                          cmd=command),
                                  shell=True,
                                  stdin=None,
                                  stdout=open(os.devnull, 'w'),
                                  stderr=open(os.devnull, 'w'),
                                  close_fds=False
                                  )
        except Exception as err:
            msg = "macsyfinder execution failed: command = {0} : {1}".format(command, err)
            print(msg)
            raise err

        macsy_process.wait()
        self.assertEqual(macsy_process.returncode, 0,
                         "macsyfinder finished with non zero exit code: {0} command launched=\n{1}".format(
                          macsy_process.returncode,
                          command))

        expected_result_path = os.path.join(self._data_dir, 'data_set_1', 'basic_run_results',
                                            'results.macsyfinder.json')
        with open(expected_result_path) as expected_result_file:
            expected_result_json = json.load(expected_result_file)

        test_result_path = os.path.join(self.out_dir, 'results.macsyfinder.json')
        with open(test_result_path) as test_result_file:
            test_result_json = json.load(test_result_file)

        # it should have only one occurrence of T9SS
        self.assertEqual(len(test_result_json), 1,
                         "different type of systems expected: 1  retrieved: {0}".format(len(test_result_json)))
        expected_result_json = expected_result_json[0]
        test_result_json = test_result_json[0]
        self.assertEqual(expected_result_json['name'],
                         test_result_json['name'],
                         "type of system name expected: {0}   retrieved: {1}".format(expected_result_json['name'],
                                                                                     test_result_json['name']))
        self.assertEqual(expected_result_json['occurrence_number'],
                         test_result_json['occurrence_number'],
                         "occurrence number expected {0}   retrieved: {1}".format(expected_result_json['occurrence_number'],
                                                                                  test_result_json['occurrence_number']))
        self.assertDictEqual(expected_result_json['replicon'],
                             test_result_json['replicon'],
                             "replicon expected {0}   retrieved: {1}".format(expected_result_json['occurrence_number'],
                                                                             test_result_json['occurrence_number']))
        self.assertEqual(expected_result_json['id'],
                         test_result_json['id'],
                         "system occurrence id expected {0}    retrieved: {1}".format(expected_result_json['id'],
                                                                                      test_result_json['id']))
        self.assertDictEqual(expected_result_json['summary']['mandatory'],
                             test_result_json['summary']['mandatory'],
                             "\nmandatory genes expected:  {0}"
                             "\nmandatory genes retrieved: {1}".format(expected_result_json['summary']['mandatory'],
                                                                       test_result_json['summary']['mandatory']))
        self.assertDictEqual(expected_result_json['summary']['accessory'],
                             test_result_json['summary']['accessory'],
                             "\naccessory genes expected:  {0}"
                             "\naccessory genes retrieved: {1}".format(expected_result_json['summary']['accessory'],
                                                                       test_result_json['summary']['accessory']))
        self.assertDictEqual(expected_result_json['summary']['forbidden'],
                             test_result_json['summary']['forbidden'],
                             "\nforbidden genes expected:  {0}"
                             "\nforbidden genes retrieved: {1}".format(expected_result_json['summary']['forbidden'],
                                                                       test_result_json['summary']['forbidden']))
        self.assertListEqual(expected_result_json['genes'], test_result_json['genes'],
                             "\ngenes expected:  {0}"
                             "\ngenes retrieved: {1}".format(expected_result_json['genes'], test_result_json['genes']))


