# -*- coding: utf-8 -*-

#===============================================================================
# Created on Nov 30, 2012
# 
# @author: bneron
# @contact: user_email
# @organization: organization_name
# @license: license
#===============================================================================


import sys
import os
import shutil


TXSSCAN_HOME = os.path.abspath('..')
if not TXSSCAN_HOME in sys.path: 
    sys.path.append(os.path.abspath('..') )

import unittest
import shutil
from txsscanlib.gene import Profile
from txsscanlib.gene import Gene
from txsscanlib.system import System
from txsscanlib.config import Config

class Test(unittest.TestCase):


    def setUp(self):
        self.cfg = Config( hmmer_exe = "hmmsearch",
                           sequence_db = "./datatest/prru_psae.001.c01.fasta",
                           ordered_db = True,
                           e_value_res = 1,
                           i_evalue_sel = 0.5,
                           def_dir = "../data/DEF",
                           res_search_dir = '.',
                           res_search_suffix = ".search_hmm.out",
                           profile_dir = "../data/profiles",
                           profile_suffix = ".fasta-aln_edit.hmm",
                           res_extract_suffix = "",
                           log_level = 30,
                           log_file = '/dev/null'
                           )

    def tearDown(self):
        shutil.rmtree(self.cfg.working_dir)

    def test_len(self):
        system = System("T2SS", self.cfg)
        gene = Gene(self.cfg, "abc", system)
        profile = Profile(gene, self.cfg)
        self.assertEqual(len(profile), 501)

    def test_unknow_profile(self):
        system = System("T2SS", self.cfg)
        gene = Gene(self.cfg, "abc", system)
        gene.name = "foo"
        self.assertRaises(IOError, Profile, gene, self.cfg)
        
    def test_str(self):
        system = System("T2SS", self.cfg)
        gene = Gene(self.cfg, "abc", system)
        profile = Profile(gene, self.cfg)
        s = "%s : %s" % (gene.name, os.path.join(self.cfg.profile_dir, gene.name + self.cfg.profile_suffix))
        self.assertEqual(str(profile), s)
        
    def test_execute(self):
        system = System("T2SS", self.cfg)
        gene = Gene(self.cfg, "abc", system)
        profile = Profile(gene, self.cfg)
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
        system = System("T2SS", self.cfg)
        gene = Gene(self.cfg, "abc", system)
        profile = Profile(gene, self.cfg)
        self.assertRaises(RuntimeError, profile.execute)
        

            
if __name__ == "__main__":
    unittest.main()