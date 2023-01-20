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

import os
import shutil
import tempfile
import argparse

from macsypy.registries import ModelRegistry, scan_models_dir
from macsypy.utils import get_def_to_detect, get_replicon_names, threads_available, parse_time
from macsypy.error import MacsypyError

from tests import MacsyTest


class TestUtils(MacsyTest):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()


    def tearDown(self):
        try:
            shutil.rmtree(self.tmp_dir)
        except:
            pass


    def test_get_def_to_detect(self):
        cmd_args = argparse.Namespace()
        cmd_args.models_dir = os.path.join(self._data_dir, 'fake_model_dir')
        cmd_args.models = ('set_1', 'def_1_1', 'def_1_2', 'def_1_3')
        registry = ModelRegistry()
        models_location = scan_models_dir(cmd_args.models_dir)
        for ml in models_location:
            registry.add(ml)

        # case where models are specified on command line
        res, model_family, model_vers = get_def_to_detect(('set_1', ['def_1_1', 'def_1_2', 'def_1_3']), registry)
        model_loc = registry['set_1']
        self.assertEqual(model_family, 'set_1')
        self.assertEqual(model_vers, '0.0b2')
        exp = [model_loc.get_definition(name) for name in ('set_1/def_1_1', 'set_1/def_1_2', 'set_1/def_1_3')]
        self.assertListEqual(res, exp)

        # case we search all models
        res, model_family, model_vers = get_def_to_detect(('set_1', ['all']), registry)
        self.assertEqual(model_family, 'set_1')
        self.assertEqual(model_vers, '0.0b2')
        exp = model_loc.get_all_definitions()
        self.assertListEqual(res, exp)

        # case the models required does not exists
        with self.assertRaises(ValueError):
            get_def_to_detect(('set_1', ['FOO', 'BAR']), registry)


    def test_get_replicon_names_gembase(self):
        replicon_names = get_replicon_names(self.find_data('base', 'gembase.fasta'), 'gembase')
        self.assertListEqual(replicon_names,
                             ['GCF_000005845', 'GCF_000006725', 'GCF_000006745', 'GCF_000006765', 'GCF_000006845',
                              'GCF_000006905', 'GCF_000006925', 'GCF_000006945'])

    def test_get_replicon_names_ordered(self):
        replicon_names = get_replicon_names(self.find_data('base', 'MOBP1_once.prt'), 'ordered_replicon')
        self.assertListEqual(replicon_names,
                             ['MOBP1_once'])

    def test_get_replicon_names_unordered(self):
        replicon_names = get_replicon_names(self.find_data('base', 'MOBP1_once.prt'), 'unordered')
        self.assertListEqual(replicon_names,
                             ['MOBP1_once'])

    def test_get_replicon_names_bad_type(self):
        with self.assertRaises(MacsypyError) as ctx:
            get_replicon_names(self.find_data('base', 'MOBP1_once.prt'), 'bad_dbtype')
        self.assertEqual(str(ctx.exception),
                         'Invalid genome type: bad_dbtype')

    def test_threads_available(self):
        if hasattr(os, "sched_getaffinity"):
            sched_getaffinity_ori = os.sched_getaffinity
        else:
            sched_getaffinity_ori = None
        cpu_count_ori = os.cpu_count

        threads_nb = 7
        cpu_nb = 8

        os.cpu_count = lambda : cpu_nb

        try:
            del os.sched_getaffinity
            self.assertEqual(threads_available(), cpu_nb)
            os.sched_getaffinity = lambda x: [None] * threads_nb
            self.assertEqual(threads_available(), threads_nb)
        finally:
            os.cpu_count = cpu_count_ori
            if sched_getaffinity_ori:
                os.sched_getaffinity = sched_getaffinity_ori
            else:
                del os.sched_getaffinity

    def test_parse_time(self):
        self.assertEqual(parse_time(10), 10)
        self.assertEqual(parse_time('10s'), 10)
        self.assertEqual(parse_time('10m'), 600)
        self.assertEqual(parse_time('1h'), 3600)
        self.assertEqual(parse_time('1d'), 86400)
        self.assertEqual(parse_time(10.5), 10)
        self.assertEqual(parse_time('10m10s1h'), 600 + 10 + 3600)
        self.assertEqual(parse_time('10m 10s 1h'), 600 + 10 + 3600)
        with self.assertRaises(ValueError) as ctx:
            parse_time('10W')
        self.assertEqual(str(ctx.exception),
                         'Not valid time format. Units allowed h/m/s.')
