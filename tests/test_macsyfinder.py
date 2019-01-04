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

import os
import shutil
import tempfile
import argparse

from macsypy.config import Config, MacsyDefaults
from macsypy.registries import ModelRegistry
from macsypy.scripts.macsyfinder import get_models_name_to_detect
from tests import MacsyTest


class Test(MacsyTest):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()


    def tearDown(self):
        try:
            shutil.rmtree(self.tmp_dir)
        except:
            pass

    def test_models_name_to_detect(self):
        cmd_args = argparse.Namespace()
        cmd_args.models_dir = os.path.join(self._data_dir, 'data_set_1', 'models')
        cmd_args.models = [('set_1', 'T9SS', 'T3SS', 'T4SS_typeI')]
        config = Config(MacsyDefaults(models_dir=os.path.join(self._data_dir, 'data_set_1', 'models')),
                        cmd_args)
        registry = ModelRegistry(config)
        res = get_models_name_to_detect([('set_1', ['T9SS', 'T3SS', 'T4SS_typeI'])], registry)
        exp = ['set_1/T9SS', 'set_1/T3SS', 'set_1/T4SS_typeI']
        self.assertListEqual(res, exp)
        with self.assertRaises(ValueError):
            get_models_name_to_detect([('set_1', ['FOO', 'BAR'])], registry)

    def test_get_version_message(self):
        pass

    def test_list_models(self):
        pass

    def parse_args(self):
        pass

