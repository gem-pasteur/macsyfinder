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
from txsscanlib.report import Hit
from txsscanlib.config import Config
from txsscanlib.gene import Gene
from txsscanlib.system import System
from txsscanlib.registries import ProfilesRegistry

class Test(unittest.TestCase):

    _data_dir = os.path.join(os.path.dirname(__file__), "datatest")

    def setUp(self):
        self.cfg = Config( hmmer_exe = "",
                           sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type = "gembase",
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

    def test_cmp(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, "gspD",system, self.profile_registry)
        h0 = Hit(gene, system, "PSAE001c01_006940", 803,"PSAE001c01", 3450, float(1.2e-234), float(779.2), float(1.000000), (741.0 - 104.0 + 1)/ 803, 104, 741)
        h1 = Hit(gene, system, "PSAE001c01_013980", 759,"PSAE001c01", 4146, float(3.7e-76), float(255.8), float(1.000000), (736.0 - 105.0 + 1)/ 759, 105, 736)
        self.assertGreater(h1, h0)
        self.assertLess(h0, h1)

    def test_eq(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, "gspD", system, self.profile_registry)
        h0 = Hit(gene, system, "PSAE001c01_006940", 803,"PSAE001c01", 3450, float(1.2e-234), float(779.2), float(1.000000), (741.0 - 104.0 + 1)/ 803, 104, 741)
        h1 = Hit(gene, system, "PSAE001c01_006940", 803,"PSAE001c01", 3450, float(1.2e-234), float(779.2), float(1.000000), (741.0 - 104.0 + 1)/ 803, 104, 741)
        h2 = Hit(gene, system, "PSAE001c01_013980", 759,"PSAE001c01", 4146, float(3.7e-76), float(255.8), float(1.000000), (736.0 - 105.0 + 1)/ 759, 105, 736)
        self.assertEqual(h0, h1)
        self.assertNotEqual(h0, h2)
        
    def test_str(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, "gspD", system, self.profile_registry)
        hit_prop={'id' : "PSAE001c01_006940",
                  'hit_seq_len': 803,
                  'replicon_name' : "PSAE001c01",
                  'position' : 694,
                  'i_eval' : float(1.2e-234),
                  'score' : float(779.2),
                  'gene_name' : gene.name,
                  'system_name' : system.name, 
                  'profil_coverage' : float(1.0),
                  'sequence_coverage' : float(638.000000),
                  'begin' : 104,
                  'end' : 741
                  }
        
        hit = Hit( gene, system, hit_prop['id'], hit_prop['hit_seq_len'], hit_prop['replicon_name'], hit_prop['position'] , hit_prop['i_eval'], hit_prop['score'], 
                   hit_prop['profil_coverage'], hit_prop['sequence_coverage'],hit_prop['begin'],hit_prop['end'])
        s = "%(id)s\t%(replicon_name)s\t%(position)d\t%(hit_seq_len)d\t%(gene_name)s\t%(system_name)s\t%(i_eval)s\t%(score)s\t%(profil_coverage)f\t%(sequence_coverage)f\t%(begin)d\t%(end)d\n" % hit_prop
        self.assertEqual(s,str(hit))
        
