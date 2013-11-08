# -*- coding: utf-8 -*-

#===============================================================================
# Created on Nov 30, 2012
# 
# @author: bneron
# @contact: user_email
# @organization: organization_name
# @license: license
#===============================================================================


import os
import unittest
import shutil
from txsscanlib.gene import Profile
from txsscanlib.gene import Gene
from txsscanlib.system import System
from txsscanlib.config import Config
from txsscanlib.registries import ProfilesRegistry
from test import which

class Test(unittest.TestCase):

    _data_dir = os.path.join(os.path.dirname(__file__), "datatest")

    def setUp(self):
        self.cfg = Config( hmmer_exe = "hmmsearch",
                           sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type = "gembase",
                           e_value_res = 1,
                           i_evalue_sel = 0.5,
                           def_dir = os.path.join(self._data_dir, 'DEF'),
                           res_search_dir = '/tmp',
                           res_search_suffix = ".search_hmm.out",
                           profile_dir = os.path.join(self._data_dir, 'profiles'),
                           profile_suffix = ".fasta-aln_edit.hmm",
                           res_extract_suffix = "",
                           log_level = 30,
                           log_file = '/dev/null'
                           )
        self.profile_registry = ProfilesRegistry(self.cfg)


    def tearDown(self):
        shutil.rmtree(self.cfg.working_dir)

    def test_len(self):
        system = System(self.cfg, "T2SS", 10)
        gene = Gene(self.cfg, "abc", system, self.profile_registry)
        path = self.profile_registry.get("abc")
        profile = Profile(gene, self.cfg, path)
        self.assertEqual(len(profile), 501)


    def test_str(self):
        system = System(self.cfg, "T2SS", 10)
        gene = Gene(self.cfg, "abc", system, self.profile_registry)
        path = self.profile_registry.get("abc")
        profile = Profile(gene, self.cfg, path)
        s = "%s : %s" % (gene.name, path)
        self.assertEqual(str(profile), s)

    @unittest.skipIf( not which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_execute(self):
        system = System(self.cfg, "T2SS", 10)
        gene = Gene(self.cfg, "abc", system, self.profile_registry)
        path = self.profile_registry.get("abc")
        profile = Profile(gene, self.cfg, path)
        report = profile.execute()
        hmmer_raw_out = profile.hmm_raw_output
        with open(hmmer_raw_out, 'r') as hmmer_raw_out_file:
            first_l = hmmer_raw_out_file.readline()
            #a hmmsearch output file has been produced
            self.assertTrue(first_l.startswith("# hmmsearch :: search profile(s) against a sequence database"))
            for i in range(5):
                #skip 4 lines
                l  =  hmmer_raw_out_file.readline()
            #a hmmsearch used the abc profile line should become with: "# query HMM file:"
            path = os.path.join(self.cfg.profile_dir, gene.name + self.cfg.profile_suffix)
            self.assertTrue(l.find(path) != -1)

    def test_execute_unknown_binary(self):
        self.cfg.options['hmmer_exe'] = "Nimportnaoik"
        system = System(self.cfg, "T2SS", 10)
        gene = Gene(self.cfg, "abc", system, self.profile_registry)
        path = self.profile_registry.get("abc")
        profile = Profile(gene, self.cfg, path)
        self.assertRaises(RuntimeError, profile.execute)

