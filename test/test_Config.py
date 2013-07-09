# -*- coding: utf-8 -*-

#===============================================================================
# Created on Jan 14, 2013
# 
# @author: bneron
# @contact: user_email
# @organization: organization_name
# @license: license
#===============================================================================

import sys
import os

TXSSCAN_HOME = os.path.abspath('..')
if not TXSSCAN_HOME in sys.path: 
    sys.path.append(os.path.abspath('..') )

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
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles'
                          )
        self.assertFalse(self.cfg.build_indexes)
        self.tearDown()
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : "./datatest/prru_psae.001.c01.fasta",
                  'db_type' : 'gembase',
                  'def_dir' : '../data/DEF',
                  'profile_dir' : '../data/profiles',
                  'build_indexes' : True
        }
        config = Config(
                        cfg_file = "nimportnaoik",
                        sequence_db = "./datatest/prru_psae.001.c01.fasta",
                        db_type = 'gembase',
                        def_dir = '../data/DEF',
                        profile_dir = '../data/profiles',
                        build_indexes = True
                        )
        self.assertTrue(config.build_indexes)
        
    def test_default(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles'
                          )
        self.assertEqual(self.cfg.hmmer_exe, 'hmmsearch')
    
    def test_coverage_profile(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles'
                          )
        self.assertEqual( self.cfg.coverage_profile, 0.5 )
        self.tearDown()
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          coverage_profile = 0.6
                          )
        self.assertEqual( self.cfg.coverage_profile, 0.6 )
        self.tearDown()
        #coverage_profile must be a float
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : "./datatest/prru_psae.001.c01.fasta",
                  'db_type' : 'gembase',
                  'def_dir' : '../data/DEF',
                  'profile_dir' : '../data/profiles',
                  'coverage_profile' : "foo"}
        self.assertRaises(ValueError, Config, **kwargs)
        
        
    def test_def_dir(self):
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : "./datatest/prru_psae.001.c01.fasta",
                  'db_type' : 'gembase',
                  'profile_dir' : '../data/profiles'}
        self.assertRaises(ValueError, Config, **kwargs)
        def_dir = '../data/DEF'
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = def_dir,
                          profile_dir = '../data/profiles'
                          )  
        self.assertEqual( def_dir, self.cfg.def_dir)                
    
    def test_e_value_res(self):
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : "./datatest/prru_psae.001.c01.fasta",
                  'db_type' : 'gembase',
                  'def_dir' : '../data/DEF',
                  'profile_dir' : '../data/profiles',
                  'e_value_res' : 'foo'
        }
        self.assertRaises(ValueError, Config, **kwargs)
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles')
        self.assertEqual(self.cfg.e_value_res, 1)
        self.tearDown()
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          e_value_res = 0.7
                          )
        self.assertEqual(self.cfg.e_value_res, 0.7)
        self.tearDown()
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : "./datatest/prru_psae.001.c01.fasta",
                  'db_type' : 'gembase',
                  'def_dir' : '../data/DEF',
                  'profile_dir' : '../data/profiles',
                  'e_value_res' : 0.7,
                  'i_evalue_sel' : 1
        }
        self.assertRaises(ValueError, Config, **kwargs)
        
    def test_hmmer_exe(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles')
        self.assertEqual(self.cfg.hmmer_exe, 'hmmsearch')
        self.tearDown()
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          hmmer_exe = 'truc'
                          )
        self.assertEqual(self.cfg.hmmer_exe, 'truc')
    
    def test_i_value_sel(self):
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : "./datatest/prru_psae.001.c01.fasta",
                  'db_type' : 'gembase',
                  'def_dir' : '../data/DEF',
                  'profile_dir' : '../data/profiles',
                  'i_evalue_sel' : 'foo'
        }
        self.assertRaises(ValueError, Config, **kwargs)
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles')
        self.assertEqual(self.cfg.i_evalue_sel, 0.5)
        self.tearDown()
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          i_evalue_sel = 0.7
                          )
        self.assertEqual(self.cfg.i_evalue_sel, 0.7)
        self.tearDown()
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : "./datatest/prru_psae.001.c01.fasta",
                  'db_type' : 'gembase',
                  'def_dir' : '../data/DEF',
                  'profile_dir' : '../data/profiles',
                  'e_value_res' : 0.7,
                  'i_evalue_sel' : 1
        }
        self.assertRaises(ValueError, Config, **kwargs)
      
    def test_db_type(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          )
        self.assertEqual( self.cfg.db_type, 'gembase')
        self.tearDown()
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : "./datatest/prru_psae.001.c01.fasta",
                  'db_type' : 'foo',
                  'def_dir' : '../data/DEF',
                  'profile_dir' : '../data/profiles',
        }
        self.assertRaises(ValueError, Config, **kwargs)
        
    def test_previous_run(self):
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : "./datatest/prru_psae.001.c01.fasta",
                  'db_type' : 'gembase',
                  'def_dir' : '../data/DEF',
                  'profile_dir' : '../data/profiles',
                  'previous_run' : 'foo'
        }
        self.assertRaises(ValueError, Config, **kwargs)
        try:
            cfg_base = Config(cfg_file = "nimportnaoik",
                              sequence_db = "./datatest/prru_psae.001.c01.fasta",
                              db_type = 'gembase',
                              def_dir = '../data/DEF',
                              profile_dir = '../data/profiles',
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
                  'sequence_db' : "./datatest/prru_psae.001.c01.fasta",
                  'db_type' : 'gembase',
                  'def_dir' : '../data/DEF',
                  'profile_dir' : 'foo',
        }
        self.assertRaises(ValueError, Config, **kwargs)
        profile_dir = '../data/profiles'
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = profile_dir,
                          )
        self.assertEqual(self.cfg.profile_dir , profile_dir)
         
    def test_profile_suffix(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          )
        self.assertEqual(self.cfg.profile_suffix, '.fasta-aln_edit.hmm')                 
        self.tearDown()
        profile_suffix = 'foo'
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          profile_suffix = profile_suffix
                          )
        self.assertEqual(self.cfg.profile_suffix, profile_suffix)
    
    def test_replicon_topology(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          )
        self.assertEqual(self.cfg.replicon_topology, 'linear')
        self.tearDown()
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          replicon_topology = 'circular'
                          )
        self.assertEqual(self.cfg.replicon_topology, 'circular')
        self.tearDown()
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : "./datatest/prru_psae.001.c01.fasta",
                  'db_type' : 'gembase',
                  'def_dir' : '../data/DEF',
                  'profile_dir' : '../data/profiles',
                  'replicon_topology' : 'foo'
        }
        self.assertRaises(ValueError, Config, **kwargs)

    def test_inter_gene_max_space(self):
        inter_gene_max_space = (["T2SS", 32], ['Flagellum', 64])
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          inter_gene_max_space = inter_gene_max_space
                          )
        self.assertEqual(self.cfg.inter_gene_max_space('T2SS'), 32)
        self.assertEqual(self.cfg.inter_gene_max_space('Flagellum'), 64)
        self.assertIsNone(self.cfg.inter_gene_max_space('Foo'))
    
    def test_min_genes_required(self):
        min_genes_required = (["T2SS", 32], ['Flagellum', 64])
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          min_genes_required = min_genes_required
                          )
        self.assertEqual(self.cfg.min_genes_required('T2SS'), 32)
        self.assertEqual(self.cfg.min_genes_required('Flagellum'), 64)
        self.assertIsNone(self.cfg.min_genes_required('Foo'))
    
    def test_min_mandatory_genes_required(self):
        min_mandatory_genes_required = (["T2SS", 32], ['Flagellum', 64])
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          min_mandatory_genes_required = min_mandatory_genes_required
                          )
        self.assertEqual(self.cfg.min_mandatory_genes_required('T2SS'), 32)
        self.assertEqual(self.cfg.min_mandatory_genes_required('Flagellum'), 64)
        self.assertIsNone(self.cfg.min_mandatory_genes_required('Foo'))    
        
    def test_res_extract_suffix(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          )
        self.assertEqual(self.cfg.res_extract_suffix, '.res_hmm_extract')                 
        self.tearDown()
        res_extract_suffix = 'foo'
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          res_extract_suffix = res_extract_suffix
                          )
        self.assertEqual(self.cfg.res_extract_suffix, res_extract_suffix)
    
    def test_res_search_dir(self):
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : "./datatest/prru_psae.001.c01.fasta",
                  'db_type' : 'gembase',
                  'def_dir' : '../data/DEF',
                  'profile_dir' : '../data/profiles',
                  'res_search_dir' : 'foo'
        }
        self.assertRaises(ValueError, Config, **kwargs)
        res_search_dir = './datatest/'
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          res_search_dir = res_search_dir
                          )
        self.assertEqual(self.cfg.res_search_dir , res_search_dir)
     
    def test_res_search_suffix(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          )
        self.assertEqual(self.cfg.res_search_suffix, '.search_hmm.out')                 
        self.tearDown()
        res_search_suffix = 'foo'
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          res_search_suffix = res_search_suffix
                          )
        self.assertEqual(self.cfg.res_search_suffix, res_search_suffix)
    
#    def test_save(self):
#        pass 
#      
    def test_sequence_db(self):
        kwargs = {'cfg_file' : "nimportnaoik",
                  'db_type' : 'gembase',
                  'def_dir' : '../data/DEF',
                  'profile_dir' : '../data/profiles',
                  'worker_nb' : '2.3'
        }
        self.assertRaises(ValueError, Config, **kwargs)   
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : "foo",
                  'db_type' : 'gembase',
                  'def_dir' : '../data/DEF',
                  'profile_dir' : '../data/profiles',
                  'worker_nb' : '2.3'
        }
        self.assertRaises(ValueError, Config, **kwargs) 
        sequence_db = "./datatest/prru_psae.001.c01.fasta"
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = sequence_db,
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          )  
        self.assertEqual(self.cfg.sequence_db, sequence_db)
    
    def test_worker_nb(self):
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles'
                          )
        self.assertEqual(self.cfg.worker_nb, 0)
        self.tearDown()
        self.cfg = Config(cfg_file = "nimportnaoik",
                          sequence_db = "./datatest/prru_psae.001.c01.fasta",
                          db_type = 'gembase',
                          def_dir = '../data/DEF',
                          profile_dir = '../data/profiles',
                          worker_nb = 2
                          )
        self.assertEqual(self.cfg.worker_nb, 2)
        self.tearDown()
        kwargs = {'cfg_file' : "nimportnaoik",
                  'sequence_db' : "./datatest/prru_psae.001.c01.fasta",
                  'db_type' : 'gembase',
                  'def_dir' : '../data/DEF',
                  'profile_dir' : '../data/profiles',
                  'worker_nb' : '2.3'
        }
        self.assertRaises(ValueError, Config, **kwargs)   
          




if __name__ == "__main__":
    unittest.main()