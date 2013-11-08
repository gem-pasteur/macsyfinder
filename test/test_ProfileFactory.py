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
from txsscanlib.registries import ProfilesRegistry
from txsscanlib.txsscan_error import TxsscanError


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
        self.profile_registry = ProfilesRegistry(self.cfg)

    def tearDown(self):
        profile_factory._profiles = {}
        shutil.rmtree(self.cfg.working_dir)

    def test_get_profile(self):
        system_foo = System(self.cfg, "foo", 10)
        gene_name = 'sctJ_FLG'
        gene = Gene(self.cfg, gene_name, system_foo, self.profile_registry)
        profile = profile_factory.get_profile(gene, self.cfg, self.profile_registry )
        self.assertTrue( isinstance( profile, Profile ))
        self.assertEqual( profile.gene.name, gene_name )

    def test_get_uniq_object(self):
        system_foo = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.profile_registry)
        path = self.profile_registry.get('sctJ_FLG')
        profile1 = profile_factory.get_profile(gene, self.cfg, path )
        profile2 = profile_factory.get_profile(gene, self.cfg, path)
        self.assertEqual( profile1, profile2 )

    def test_unknow_profile(self):
        system_foo = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.profile_registry)
        gene.name = "foo"
        self.assertRaises(TxsscanError, profile_factory.get_profile, gene, self.cfg, self.profile_registry)

