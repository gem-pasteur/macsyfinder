# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Neron                                         #
# Copyright (c) 2014-2019  Institut Pasteur (Paris) and CNRS.                  #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################


import tempfile
import argparse

from macsypy.model import ModelBank
from macsypy.model import Model
from macsypy.config import Config, MacsyDefaults
from tests import MacsyTest


class Test(MacsyTest):
    
    def setUp(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.res_search_dir = tempfile.gettempdir()
        args.log_level = 30
        self.cfg = Config(MacsyDefaults(), args)
        self.system_bank = ModelBank()

    def tearDown(self):
        ModelBank._model_bank = {}

    def test_add_get_system(self):
        model_name = 'foo'
        self.assertRaises(KeyError,  self.system_bank.__getitem__, model_name)
        system_foo = Model(self.cfg, model_name, 10)
        self.system_bank.add_model(system_foo)
        self.assertTrue(isinstance(system_foo, Model))
        self.assertEqual(system_foo,  self.system_bank[model_name])
        with self.assertRaises(KeyError) as ctx:
            self.system_bank.add_model(system_foo)
        self.assertEqual(str(ctx.exception),
                         '"a model named {0} is already registered in the models\' bank"'.format(model_name))

    def test_contains(self):
        system_in = Model(self.cfg, "foo", 10)
        self.system_bank.add_model(system_in)
        self.assertIn(system_in,  self.system_bank)
        system_out = Model(self.cfg, "bar", 10)
        self.assertNotIn(system_out,  self.system_bank)

    def test_iter(self):
        systems = [Model(self.cfg, 'foo', 10), Model(self.cfg, 'bar', 10)]
        for s in systems:
            self.system_bank.add_model(s)
        i = 0
        for s in self.system_bank:
            self.assertIn(s, systems)
            i += 1
        self.assertEqual(i, len(systems))

    def test_get_uniq_object(self):
        system_foo = Model(self.cfg, "foo", 10)
        self.system_bank.add_model(system_foo)
        system_1 = self.system_bank[system_foo.name]
        system_2 = self.system_bank[system_foo.name]
        self.assertEqual(system_1, system_2)
