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
from macsypy.gene import Homolog
from macsypy.gene import Gene
from macsypy.system import System
from macsypy.config import Config
from macsypy.registries import ProfilesRegistry

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
        shutil.rmtree(self.cfg.working_dir)


    def test_gene_ref(self):
        system = System(self.cfg, "T2SS", 10)
        gene_ref = Gene(self.cfg, 'sctJ_FLG', system, self.profile_registry)
        gene = Gene(self.cfg, 'sctJ', system, self.profile_registry)
        homolog_1 = Homolog(gene, gene_ref)
        self.assertEqual( homolog_1.gene_ref , gene_ref)
 
    def test_is_aligned(self):
        system = System(self.cfg, "T2SS", 10)
        gene_ref = Gene(self.cfg, 'sctJ_FLG', system, self.profile_registry)
        gene = Gene(self.cfg, 'sctJ', system, self.profile_registry)
        homolog = Homolog( gene, gene_ref)
        self.assertFalse( homolog.is_aligned() )
        homolog = Homolog(gene, gene_ref, aligned = True  )
        self.assertTrue( homolog.is_aligned() )

    def test_delegation(self):
        system = System(self.cfg, "T2SS", 10)
        gene_ref = Gene(self.cfg, 'sctJ_FLG', system, self.profile_registry)
        gene = Gene(self.cfg, 'sctJ', system, self.profile_registry)
        homolog = Homolog( gene, gene_ref)
        self.assertEqual( homolog.system , system )

