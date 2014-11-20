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
import unittest
import shutil
import tempfile
import platform
from macsypy.report import Hit
from macsypy.config import Config
from macsypy.gene import Gene
from macsypy.system import System
from macsypy.registries import ProfilesRegistry

class Test(unittest.TestCase):

    _data_dir = os.path.join(os.path.dirname(__file__), "datatest")

    def setUp(self):
        self.cfg = Config( hmmer_exe = "",
                           sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type = "gembase",
                           e_value_res = 1,
                           i_evalue_sel = 0.5,
                           def_dir = os.path.join(self._data_dir, 'DEF'),
                           res_search_dir = tempfile.gettempdir(),
                           res_search_suffix = "",
                           profile_dir = os.path.join(self._data_dir, 'profiles'),
                           profile_suffix = ".hmm",
                           res_extract_suffix = "",
                           log_level = 30,
                           log_file = 'NUL' if platform.system() == 'Windows' else '/dev/null'
                           )
        self.profile_registry = ProfilesRegistry(self.cfg)
        
  
    def tearDown(self):
        # close loggers filehandles, so they don't block file deletion
        # in shutil.rmtree calls in Windows
        handlers = self.cfg.options['logger'].handlers[:]
        for handler in handlers:
            handler.close()
            self.cfg.options['logger'].removeHandler(handler)

        handlers = self.cfg.options['out_logger'].handlers[:]
        for handler in handlers:
            handler.close()
            self.cfg.options['out_logger'].removeHandler(handler)

        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass

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
        
