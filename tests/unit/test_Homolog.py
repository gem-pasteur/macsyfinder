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
import logging
from macsypy.gene import Homolog
from macsypy.gene import Gene
from macsypy.system import System
from macsypy.config import Config
from macsypy.registries import ProfilesRegistry

class Test(unittest.TestCase):

    _data_dir = os.path.join(os.path.dirname(__file__), "..", "datatest")

    def setUp(self):
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        
        #add only one handler to the macsypy logger
        from macsypy.gene import _log
        macsy_log = _log.parent
        log_file = 'NUL' if platform.system() == 'Windows' else '/dev/null'
        log_handler = logging.FileHandler(log_file)
        macsy_log.addHandler(log_handler)
        
        self.cfg = Config( sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type = "gembase",
                           hmmer_exe = "",
                           e_value_res = 1,
                           i_evalue_sel = 0.5,
                           def_dir = os.path.join(self._data_dir, 'DEF'),
                           res_search_dir = tempfile.gettempdir(),
                           res_search_suffix = "",
                           profile_dir = os.path.join(self._data_dir, 'profiles'),
                           profile_suffix = ".hmm",
                           res_extract_suffix = "",
                           log_level = 30,
                           log_file = log_file
                           )
        self.profile_registry = ProfilesRegistry(self.cfg)
        

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

