# -*- coding: utf-8 -*-

#===============================================================================
# Created on Jan 15, 2013
# 
# @author: bneron
# @contact: user_email
# @organization: organization_name
# @license: license
#===============================================================================


import os
import unittest
import shutil
from txsscanlib.system import system_bank
from txsscanlib.system import System
from txsscanlib.config import Config


class Test(unittest.TestCase):


    _data_dir = os.path.join(os.path.dirname(__file__), "datatest")
    
    def setUp(self):
        self.cfg = Config( sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type = "gembase",
                           hmmer_exe = "",
                           e_value_res = 1,
                           i_evalue_sel = 0.5,
                           def_dir = os.path.join(self._data_dir, 'DEF'),
                           res_search_dir = "/tmp",
                           res_search_suffix = "",
                           profile_dir = os.path.join(self._data_dir, 'profiles'),
                           profile_suffix = ".fasta-aln_edit.hmm",
                           res_extract_suffix = "",
                           log_level = 30,
                           log_file = '/dev/null'
                           )

    def tearDown(self):
        shutil.rmtree(self.cfg.working_dir)

    def test_add_get_system(self):
        system_name= 'foo'
        self.assertRaises(KeyError, system_bank.__getitem__, system_name)
        system_foo = System(self.cfg, system_name, 10)
        system_bank.add_system(system_foo)
        self.assertTrue( isinstance( system_foo, System ))
        self.assertEqual( system_foo, system_bank[system_name] )

    def test_contains(self):
        system_in = System(self.cfg, "foo", 10)
        system_bank.add_system(system_in)
        self.assertIn(system_in, system_bank)
        system_out = System(self.cfg, "bar", 10)
        self.assertNotIn( system_out, system_bank)

    def test_iter(self):
        systems = [System(self.cfg, 'foo', 10), System(self.cfg, 'bar', 10)]
        for s in systems:
            system_bank.add_system(s)
        i = 0
        for s in system_bank:
            self.assertIn(s, systems)
            i = i + 1
        self.assertEqual(i, len(systems))

    def test_get_uniq_object(self):
        system_foo = System(self.cfg, "foo", 10)
        system_bank.add_system(system_foo)
        system_1 = system_bank[system_foo.name]
        system_2 = system_bank[system_foo.name]
        self.assertEqual( system_1, system_2 )
        
