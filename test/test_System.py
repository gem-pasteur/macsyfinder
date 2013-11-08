# -*- coding: utf-8 -*-

#===============================================================================
# Created on Jan 14, 2013
# 
# @author: bneron
# @contact: user_email
# @organization: organization_name
# @license: license
#===============================================================================

import os
import unittest
import shutil
from txsscanlib.config import Config
from txsscanlib.system import System
from txsscanlib.gene import Gene
from txsscanlib.registries import ProfilesRegistry


class Test(unittest.TestCase):

    _data_dir = os.path.join(os.path.dirname(__file__), "datatest")


    def setUp(self):
        self.cfg = Config(sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
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

    def test_name(self):
        name = 'foo'
        system = System(self.cfg, name, 10)
        self.assertEqual(system.name, name)

    def test_inter_gene_max_space(self):
        name = 'foo'
        inter_gene_max_space = 40
        system = System(self.cfg, name, inter_gene_max_space)
        self.assertEqual(system.inter_gene_max_space, inter_gene_max_space )

    def test_min_genes_required(self):
        name = 'foo'
        min_genes_required = 40
        system = System(self.cfg, name, 10, min_genes_required = min_genes_required)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.profile_registry)
        system.add_mandatory_gene( gene )
        self.assertEqual(system.min_genes_required, min_genes_required)
        #see https://projets.pasteur.fr/issues/1850
        system = System(self.cfg, name, 10)
        self.assertEqual(system.min_genes_required , len(system.mandatory_genes))
        
    def test_min_mandatory_genes_required(self):
        name = 'foo'
        min_mandatory_genes_required = 40
        system = System(self.cfg, name, 10, min_mandatory_genes_required = min_mandatory_genes_required)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.profile_registry)
        system.add_mandatory_gene( gene )
        self.assertEqual( system.min_mandatory_genes_required , min_mandatory_genes_required )    
        #see https://projets.pasteur.fr/issues/1850
        system = System(self.cfg, name, 10)
        self.assertEqual(system.min_mandatory_genes_required , len(system.mandatory_genes))
        
    def test_add_mandatory_gene(self):
        system = System("foo", self.cfg, 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.profile_registry)
        system.add_mandatory_gene( gene )
        self.assertEqual(system._mandatory_genes, [gene])
        self.assertEqual(system._allowed_genes, [])
        self.assertEqual(system._forbidden_genes, [])

    def test_add_allowed_gene(self):
        system = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.profile_registry)
        system.add_allowed_gene( gene )
        self.assertEqual( system._allowed_genes, [gene])
        self.assertEqual( system._mandatory_genes, [])
        self.assertEqual( system._forbidden_genes, [])

    def test_add_forbidden_gene(self):
        system = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.profile_registry)
        system.add_forbidden_gene( gene )
        self.assertEqual( system._forbidden_genes, [gene])
        self.assertEqual( system._allowed_genes, [])
        self.assertEqual( system._mandatory_genes, [])

    def test_mandatory_genes(self):
        system = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.profile_registry)
        system.add_mandatory_gene( gene )
        self.assertEqual( system.mandatory_genes, [gene])

    def test_allowed_genes(self):
        system = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.profile_registry)
        system.add_allowed_gene( gene )
        self.assertEqual( system.allowed_genes, [gene])

    def test_forbidden_genes(self):
        system = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.profile_registry)
        system.add_forbidden_gene( gene )
        self.assertEqual( system.forbidden_genes, [gene])


