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
from macsypy.gene import Gene
from macsypy.gene import Homolog
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
        logging.shutdown()
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass


    def test_add_homolog(self):
        system_foo = System(self.cfg, "foo", 10)
        system_bar = System(self.cfg, "bar", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.profile_registry)
        gene_ref = Gene(self.cfg, 'sctJ', system_bar, self.profile_registry)
        homolog = Homolog(self.cfg, gene, gene_ref)
        gene.add_homolog( homolog )
        self.assertEqual(len( gene.homologs), 1)
        self.assertEqual(gene.homologs[0], homolog)
    
    
    def test_get_homologs(self):
        system_foo = System(self.cfg, "foo", 10)
        system_bar = System(self.cfg, "bar", 10)
        gene = Gene(self.cfg, 'sctN', system_foo, self.profile_registry)
        sctJ_FLG = Gene(self.cfg, 'sctJ_FLG', system_foo, self.profile_registry)
        sctJ = Gene(self.cfg, 'sctJ', system_bar, self.profile_registry)
        homolog_1 = Homolog(sctJ_FLG, gene)
        gene.add_homolog(homolog_1)
        homolog_2 = Homolog(sctJ, gene)
        gene.add_homolog(homolog_2)
        self.assertEqual(gene.get_homologs(), [homolog_1, homolog_2] )
    
        
    def test_system(self):
        """
        test getter/setter for system property
        """
        system_foo = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.profile_registry)
        self.assertEqual(gene.system, system_foo)
    
    
    def test_loner(self):
        """
        test getter for loner property
        """
        system_foo = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.profile_registry)
        self.assertFalse(gene.loner)
        gene = Gene(self.cfg, 'sctJ', system_foo, self.profile_registry, loner = True)
        self.assertTrue(gene.loner)


    def test_exchangeable(self):
        """
        test getter for exchangeable property
        """
        system_foo = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.profile_registry)
        self.assertFalse(gene.exchangeable)
        gene = Gene(self.cfg, 'sctJ', system_foo, self.profile_registry, exchangeable = True)
        self.assertTrue(gene.exchangeable)
 
    def test_multi_system(self):
        """
        test getter for multi_system property
        """
        system_foo = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.profile_registry)
        self.assertFalse(gene.multi_system)
        gene = Gene(self.cfg, 'sctJ', system_foo, self.profile_registry, multi_system = True)
        self.assertTrue(gene.multi_system)


    def test_inter_gene_max_space(self):
        """
        test getter for inter_gene_max_space property
        """
        system_inter_gene_max_space = 40
        gene_inter_gene_max_space = 50
        system_foo = System(self.cfg, "foo", system_inter_gene_max_space)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.profile_registry)
        self.assertEqual(gene.inter_gene_max_space, system_inter_gene_max_space)
        gene = Gene(self.cfg, 'sctJ', system_foo, self.profile_registry, inter_gene_max_space = gene_inter_gene_max_space)
        self.assertEqual(gene.inter_gene_max_space, gene_inter_gene_max_space)


    def test_str(self):
        """
        """
        system_foo = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.profile_registry)
        system_bar = System(self.cfg, "bar", 20)
        gene_homolog = Gene(self.cfg, 'sctJ', system_bar, self.profile_registry)
        homolog = Homolog( gene_homolog, gene, self.cfg)
        gene.add_homolog( homolog )
        s = """name : sctJ_FLG
inter_gene_max_space: 10
    homologs: sctJ"""
        self.assertEqual( str(gene) , s )
