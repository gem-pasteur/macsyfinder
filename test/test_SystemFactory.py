# -*- coding: utf-8 -*-

#===============================================================================
# Created on Jan 15, 2013
# 
# @author: bneron
# @contact: user_email
# @organization: organization_name
# @license: license
#===============================================================================


import sys
import os
TXSSCAN_HOME = os.path.abspath('..')
if not TXSSCAN_HOME in sys.path:
    sys.path.append(os.path.abspath('..') )

import unittest
import shutil
from txsscanlib.system import system_factory
from txsscanlib.system import System
from txsscanlib.config import Config


class Test(unittest.TestCase):


    def setUp(self):
        self.cfg = Config( sequence_db = ".",
                           db_type = "gembase",
                           hmmer_exe = "",
                           e_value_res = 1,
                           i_evalue_sel = 0.5,
                           def_dir = "../data/DEF",
                           res_search_dir = ".",
                           res_search_suffix = "",
                           profile_dir = "../data/profiles",
                           profile_suffix = ".fasta-aln_edit.hmm",
                           res_extract_suffix = "",
                           log_level = 30,
                           log_file = '/dev/null'
                           )

    def tearDown(self):
        shutil.rmtree(self.cfg.working_dir)

    def test_get_system(self):
        system_foo = system_factory.get_system("foo", 20, self.cfg)
        self.assertTrue( isinstance( system_foo, System ))
        self.assertEqual( system_foo.name, "foo" )
        self.assertEqual( system_foo.inter_gene_max_space, 20 )
    
    def test_get_uniq_object(self):
        system_1 = system_factory.get_system("foo", 20, self.cfg)
        system_2 = system_factory.get_system("foo", 20, self.cfg)
        self.assertEqual( system_1, system_2 )
        
        
                 
if __name__ == "__main__":
    unittest.main()