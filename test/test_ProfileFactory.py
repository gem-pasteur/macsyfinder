'''
Created on Nov 30, 2012

@author: bneron
'''

import os
import unittest
import shutil
from txsscanlib.gene import profile_factory
from txsscanlib.gene import Profile
from txsscanlib.gene import Gene
from txsscanlib.system import System
from txsscanlib.config import Config


class Test(unittest.TestCase):


    def setUp(self):
        self.cfg = Config( sequence_db = ".",
                           db_type = "gembase",
                           hmmer_exe = "",
                           e_value_res = 1,
                           i_evalue_sel = 0.5,
                           def_dir = os.path.join(os.path.dirname(__file__), ".." , "data", "DEF"),
                           res_search_dir = "/tmp",
                           res_search_suffix = "",
                           profile_dir = os.path.join(os.path.dirname(__file__), "..", "data", "profiles"),
                           profile_suffix = ".fasta-aln_edit.hmm",
                           res_extract_suffix = "",
                           log_level = 30,
                           log_file = '/dev/null'
                           )

    def tearDown(self):
        profile_factory._profiles = {}
        shutil.rmtree(self.cfg.working_dir)

    def test_get_profile(self):
        system_foo = System(self.cfg, "foo", 10)
        gene_name = 'sctJ_FLG'
        gene = Gene(self.cfg, gene_name, system_foo)
        profile = profile_factory.get_profile(gene,self.cfg )
        self.assertTrue( isinstance( profile, Profile ))
        self.assertEqual( profile.gene.name, gene_name )

    def test_get_uniq_object(self):
        system_foo = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo)
        profile1 = profile_factory.get_profile(gene,self.cfg )
        profile2 = profile_factory.get_profile(gene,self.cfg )
        self.assertEqual( profile1, profile2 )


