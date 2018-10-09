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
import logging
from macsypy.system import SystemBank
from macsypy.system import System
from macsypy.config import Config
from tests import MacsyTest


class Test(MacsyTest):
    
    def setUp(self):
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        
        # add only one handler to the macsypy logger
        from macsypy.system import _log
        macsy_log = _log.parent
        log_file = os.devnull
        log_handler = logging.FileHandler(log_file)
        macsy_log.addHandler(log_handler)
        
        self.cfg = Config(hmmer_exe="",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type="gembase",
                          e_value_res=1,
                          i_evalue_sel=0.5,
                          models_dir=self.find_data('models'),
                          res_search_dir=tempfile.gettempdir(),
                          res_search_suffix="",
                          profile_suffix=".hmm",
                          res_extract_suffix="",
                          log_level=30,
                          log_file=log_file
                          )
        self.system_bank = SystemBank()


    def tearDown(self):
        # close loggers filehandles, so they don't block file deletion
        # in shutil.rmtree calls in Windows
        logging.shutdown()
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass
        SystemBank._system_bank = {}

    def test_add_get_system(self):
        model_name = 'foo'
        self.assertRaises(KeyError,  self.system_bank.__getitem__, model_name)
        system_foo = System(self.cfg, model_name, 10)
        self.system_bank.add_system(system_foo)
        self.assertTrue(isinstance(system_foo, System))
        self.assertEqual(system_foo,  self.system_bank[model_name])
        with self.assertRaises(KeyError) as ctx:
            self.system_bank.add_system(system_foo)
        self.assertEqual(ctx.exception.message,
                         "a system named {0} is already registered in the systems' bank".format(model_name))

    def test_contains(self):
        system_in = System(self.cfg, "foo", 10)
        self.system_bank.add_system(system_in)
        self.assertIn(system_in,  self.system_bank)
        system_out = System(self.cfg, "bar", 10)
        self.assertNotIn(system_out,  self.system_bank)

    def test_iter(self):
        systems = [System(self.cfg, 'foo', 10), System(self.cfg, 'bar', 10)]
        for s in systems:
            self.system_bank.add_system(s)
        i = 0
        for s in self.system_bank:
            self.assertIn(s, systems)
            i += 1
        self.assertEqual(i, len(systems))

    def test_get_uniq_object(self):
        system_foo = System(self.cfg, "foo", 10)
        self.system_bank.add_system(system_foo)
        system_1 = self.system_bank[system_foo.name]
        system_2 = self.system_bank[system_foo.name]
        self.assertEqual(system_1, system_2)
