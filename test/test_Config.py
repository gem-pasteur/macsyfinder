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
import time
from txsscanlib.config import Config

class Test(unittest.TestCase):

    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass

    def test_build_indexes(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = '/tmp'
                          )
        self.assertFalse(self.cfg.build_indexes)
        self.tearDown()
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                  'db_type' : 'gembase',
                  'def_dir' : os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                  'profile_dir' : os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                  'build_indexes' : True
        }
        config = Config(
                        cfg_file = "nimportnaoik",
                        sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                        db_type = 'gembase',
                        def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                        profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                        res_search_dir = '/tmp',
                        build_indexes = True
                        )
        self.assertTrue(config.build_indexes)


    def test_default(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = '/tmp'
                          )
        self.assertEqual(self.cfg.hmmer_exe, 'hmmsearch')


    def test_coverage_profile(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = '/tmp',
                          )
        self.assertEqual( self.cfg.coverage_profile, 0.5 )
        self.tearDown()
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          coverage_profile = 0.6,
                          res_search_dir = '/tmp',
                          )
        self.assertEqual( self.cfg.coverage_profile, 0.6 )
        self.tearDown()
        #coverage_profile must be a float
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                  'db_type' : 'gembase',
                  'def_dir' : os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                  'profile_dir' : os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                  'coverage_profile' : "foo",
                  'res_search_dir' : '/tmp'
        }
        self.assertRaises(ValueError, Config, **kwargs)


    def test_def_dir(self):
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                  'db_type' : 'gembase',
                  'profile_dir' : os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                  'res_search_dir' : '/tmp',
                  }
        real_def_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'DEF'),
        self.assertRaises(ValueError, Config, **kwargs)
        
        def_dir = os.path.join('/tmp', 'txsscan_DEF')
        if not os.path.exists(def_dir):
            os.mkdir(def_dir)
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = def_dir,
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = '/tmp',
                          ) 
        self.assertEqual( def_dir, self.cfg.def_dir)
        shutil.rmtree(def_dir)
        self.tearDown()
        def_dir = os.path.join(os.path.dirname(__file__),'..', 'data',  'DEF')
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = def_dir,
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = '/tmp',
                          )  
        self.assertEqual( def_dir, self.cfg.def_dir)


    def test_e_value_res(self):
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                  'db_type' : 'gembase',
                  'def_dir' : os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                  'profile_dir' : os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                  'e_value_res' : 'foo',
                  'res_search_dir' : '/tmp'
        }
        self.assertRaises(ValueError, Config, **kwargs)
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = '/tmp'
                          )
        self.assertEqual(self.cfg.e_value_res, 1)
        self.tearDown()
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          e_value_res = 0.7,
                          res_search_dir = '/tmp',
                          )
        self.assertEqual(self.cfg.e_value_res, 0.7)
        self.tearDown()
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                  'db_type' : 'gembase',
                  'def_dir' : os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                  'profile_dir' : os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                  'e_value_res' : 0.7,
                  'i_evalue_sel' : 1,
                  'res_search_dir' : '/tmp',
        }
        self.assertRaises(ValueError, Config, **kwargs)
        
        
    def test_hmmer_exe(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = '/tmp'
                          )
        self.assertEqual(self.cfg.hmmer_exe, 'hmmsearch')
        self.tearDown()
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          hmmer_exe = 'truc',
                          res_search_dir = '/tmp',
                          )
        self.assertEqual(self.cfg.hmmer_exe, 'truc')
    
    
    def test_i_value_sel(self):
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                  'db_type' : 'gembase',
                  'def_dir' : os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                  'profile_dir' : os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                  'i_evalue_sel' : 'foo',
                  'res_search_dir' : '/tmp',
        }
        self.assertRaises(ValueError, Config, **kwargs)
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = '/tmp'
                          )
        self.assertEqual(self.cfg.i_evalue_sel, 0.5)
        self.tearDown()
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          i_evalue_sel = 0.7,
                          res_search_dir = '/tmp',
                          )
        self.assertEqual(self.cfg.i_evalue_sel, 0.7)
        self.tearDown()
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                  'db_type' : 'gembase',
                  'def_dir' : os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                  'profile_dir' : os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                  'e_value_res' : 0.7,
                  'i_evalue_sel' : 1,
                  'res_search_dir' : '/tmp',
        }
        self.assertRaises(ValueError, Config, **kwargs)
      
      
    def test_db_type(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = '/tmp',
                          )
        self.assertEqual( self.cfg.db_type, 'gembase')
        self.tearDown()
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                  'db_type' : 'foo',
                  'def_dir' : os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                  'profile_dir' : os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                  'res_search_dir' : '/tmp',
        }
        self.assertRaises(ValueError, Config, **kwargs)
        
        
    def test_previous_run(self):
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                  'db_type' : 'gembase',
                  'def_dir' : os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                  'profile_dir' : os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                  'previous_run' : 'foo',
                  'res_search_dir' : '/tmp',
        }
        self.assertRaises(ValueError, Config, **kwargs)
        try:
            cfg_base = Config(cfg_file = "nimportnaoik",
                              sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                              db_type = 'gembase',
                              def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                              profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                              res_search_dir = '/tmp',
                              )
            self.assertIsNone(cfg_base.previous_run)
            cfg_base.save( cfg_base.working_dir )
            #wait
            time.sleep(1)
            new_cfg = Config(previous_run = cfg_base.working_dir)
            self.assertEqual(new_cfg.previous_run, cfg_base.working_dir)
        finally:
            try:
                shutil.rmtree(cfg_base.working_dir)
            except:
                pass
            try:
                shutil.rmtree(new_cfg.working_dir)
            except:
                pass
            
            
    def test_profile_dir(self):
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                  'db_type' : 'gembase',
                  'def_dir' : os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                  'profile_dir' : 'foo',
                  'res_search_dir' : '/tmp',
        }
        self.assertRaises(ValueError, Config, **kwargs)
        profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles')
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = profile_dir,
                          res_search_dir = '/tmp',
                          )
        self.assertEqual(self.cfg.profile_dir , profile_dir)
         
         
    def test_profile_suffix(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = '/tmp'
                          )
        self.assertEqual(self.cfg.profile_suffix, '.fasta-aln_edit.hmm')
        self.tearDown()
        profile_suffix = 'foo'
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          profile_suffix = profile_suffix,
                          res_search_dir = '/tmp',
                          )
        self.assertEqual(self.cfg.profile_suffix, profile_suffix)
    
    
    def test_replicon_topology(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = '/tmp',
                          )
        self.assertEqual(self.cfg.replicon_topology, 'circular')
        self.tearDown()
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          replicon_topology = 'linear',
                          res_search_dir = '/tmp',
                          )
        self.assertEqual(self.cfg.replicon_topology, 'linear')
        self.tearDown()
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                  'db_type' : 'gembase',
                  'def_dir' : os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                  'profile_dir' : os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                  'replicon_topology' : 'foo'
        }
        self.assertRaises(ValueError, Config, **kwargs)


    def test_inter_gene_max_space(self):
        inter_gene_max_space = (["T2SS", 32], ['Flagellum', 64])
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          inter_gene_max_space = inter_gene_max_space,
                          res_search_dir = '/tmp',
                          )
        self.assertEqual(self.cfg.inter_gene_max_space('T2SS'), 32)
        self.assertEqual(self.cfg.inter_gene_max_space('Flagellum'), 64)
        self.assertIsNone(self.cfg.inter_gene_max_space('Foo'))
    
    
    def test_min_genes_required(self):
        min_genes_required = (["T2SS", 32], ['Flagellum', 64])
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          min_genes_required = min_genes_required,
                          res_search_dir = '/tmp',
                          )
        self.assertEqual(self.cfg.min_genes_required('T2SS'), 32)
        self.assertEqual(self.cfg.min_genes_required('Flagellum'), 64)
        self.assertIsNone(self.cfg.min_genes_required('Foo'))
    
    
    def test_min_mandatory_genes_required(self):
        min_mandatory_genes_required = (["T2SS", 32], ['Flagellum', 64])
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          min_mandatory_genes_required = min_mandatory_genes_required,
                          res_search_dir = '/tmp',
                          )
        self.assertEqual(self.cfg.min_mandatory_genes_required('T2SS'), 32)
        self.assertEqual(self.cfg.min_mandatory_genes_required('Flagellum'), 64)
        self.assertIsNone(self.cfg.min_mandatory_genes_required('Foo'))
        
        
    def test_res_extract_suffix(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = '/tmp',
                          )
        self.assertEqual(self.cfg.res_extract_suffix, '.res_hmm_extract')
        self.tearDown()
        res_extract_suffix = 'foo'
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_extract_suffix = res_extract_suffix,
                          res_search_dir = '/tmp',
                          )
        self.assertEqual(self.cfg.res_extract_suffix, res_extract_suffix)
        
        
    def test_res_search_dir(self):
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                  'db_type' : 'gembase',
                  'def_dir' : os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                  'profile_dir' : os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                  'res_search_dir' : 'foo',
        }
        self.assertRaises(ValueError, Config, **kwargs)
        res_search_dir = os.path.join(os.path.dirname(__file__), 'datatest')
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = res_search_dir
                          )
        self.assertEqual(self.cfg.res_search_dir , res_search_dir)


    def test_res_search_suffix(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data',  'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = '/tmp',
                          )
        self.assertEqual(self.cfg.res_search_suffix, '.search_hmm.out')
        self.tearDown()
        res_search_suffix = 'foo'
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_suffix = res_search_suffix,
                          res_search_dir = '/tmp',
                          )
        self.assertEqual(self.cfg.res_search_suffix, res_search_suffix)
    
#    def test_save(self):
#        pass 
#      
    def test_sequence_db(self):
        kwargs = {'cfg_file' : "nimportnaoik",
                  'db_type' : 'gembase',
                  'def_dir' : os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                  'profile_dir' : os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                  'worker_nb' : '2.3' ,
                  'res_search_dir' : '/tmp',
        }
        self.assertRaises(ValueError, Config, **kwargs)   
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : "foo",
                  'db_type' : 'gembase',
                  'def_dir' : os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                  'profile_dir' : os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                  'worker_nb' : '2.3',
                  'res_search_dir' : '/tmp',
        }
        self.assertRaises(ValueError, Config, **kwargs) 
        sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta")
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = sequence_db,
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = '/tmp',
                          )  
        self.assertEqual(self.cfg.sequence_db, sequence_db)


    def test_worker_nb(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = '/tmp'
                          )
        self.assertEqual(self.cfg.worker_nb, 0)
        self.tearDown()
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                          db_type = 'gembase',
                          def_dir = os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                          profile_dir = os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                          res_search_dir = '/tmp',
                          worker_nb = 2
                          )
        self.assertEqual(self.cfg.worker_nb, 2)
        self.tearDown()
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : os.path.join(os.path.dirname(__file__), "datatest", "prru_psae.001.c01.fasta"),
                  'db_type' : 'gembase',
                  'def_dir' : os.path.join(os.path.dirname(__file__),'..', 'data', 'DEF'),
                  'profile_dir' : os.path.join(os.path.dirname(__file__), '..', 'data', 'profiles'),
                  'res_search_dir' : '/tmp',
                  'worker_nb' : '2.3'
        }
        self.assertRaises(ValueError, Config, **kwargs)   


