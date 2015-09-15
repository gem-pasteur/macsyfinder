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
import time
import tempfile
import logging
import platform

from macsypy.config import Config

class Test(unittest.TestCase):

    _data_dir = os.path.join(os.path.dirname(__file__), "datatest")
    
    def setUp(self):
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        if hasattr(self, 'cfg'):
            try:
                shutil.rmtree(self.cfg.working_dir)
            except Exception as err:
                pass
        self.tmp_dir = tempfile.gettempdir()
        
        
    def tearDown(self):
        logging.shutdown()
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        if hasattr(self, 'cfg'):
            try:
                shutil.rmtree(self.cfg.working_dir)
            except Exception as err:
                pass

  
    def test_build_indexes(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir = self.tmp_dir
                          )
        self.assertFalse(self.cfg.build_indexes)
        self.tearDown()

        config = Config(
                        cfg_file="nimportnaoik",
                        sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                        db_type='gembase',
                        res_search_dir = self.tmp_dir,
                        build_indexes = True
                        )
        self.assertTrue(config.build_indexes)


    def test_default(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.hmmer_exe, 'hmmsearch')


    def test_coverage_profile(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir,
                          )
        self.assertEqual(self.cfg.coverage_profile, 0.5)
        self.tearDown()
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          coverage_profile=0.6,
                          res_search_dir=self.tmp_dir,
                          )
        self.assertEqual(self.cfg.coverage_profile, 0.6)
        self.tearDown()
        #coverage_profile must be a float
        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': os.path.join(self._data_dir, "base", "test_base.fa"),
                  'db_type': 'gembase',
                  'coverage_profile': "foo",
                  'res_search_dir': self.tmp_dir
                  }
        self.assertRaises(ValueError, Config, **kwargs)


    def test_def_dir(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertIsNone(self.cfg.def_dir)
        self.tearDown()

        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': os.path.join(self._data_dir, "base", "test_base.fa"),
                  'db_type': 'gembase',
                  'def_dir': 'foo',
                  'res_search_dir': self.tmp_dir
                  }
        self.assertRaises(ValueError, Config, **kwargs)
        self.tearDown()

        def_dir = os.path.join(tempfile.gettempdir(), 'macsyfinder_DEF')
        if not os.path.exists(def_dir):
            os.mkdir(def_dir)
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          def_dir=def_dir,
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(def_dir, self.cfg.def_dir)
        shutil.rmtree(def_dir)
        self.tearDown()


    def test_models_dir(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                       sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                       db_type='gembase',
                       res_search_dir=self.tmp_dir
                       )
        self.assertIsNone(self.cfg.models_dir)
        self.tearDown()

        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': os.path.join(self._data_dir, "base", "test_base.fa"),
                  'db_type': 'gembase',
                  'models_dir': 'foo',
                  'res_search_dir': self.tmp_dir
                  }
        self.assertRaises(ValueError, Config, **kwargs)
        self.tearDown()

        models_dir = os.path.join(tempfile.gettempdir(), 'macsyfinder_models_dir')
        if not os.path.exists(models_dir):
            os.mkdir(models_dir)
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          models_dir=models_dir,
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(models_dir, self.cfg.models_dir)
        shutil.rmtree(models_dir)
        self.tearDown()


    def test_e_value_res(self):
        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': os.path.join(self._data_dir, "base", "test_base.fa"),
                  'db_type': 'gembase',
                  'def_dir': os.path.join(self._data_dir, 'DEF'),
                  'profile_dir': os.path.join(self._data_dir, 'profiles'),
                  'e_value_res': 'foo',
                  'res_search_dir': self.tmp_dir
                  }
        self.assertRaises(ValueError, Config, **kwargs)
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir = self.tmp_dir
                          )
        self.assertEqual(self.cfg.e_value_res, 1)
        self.tearDown()
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          e_value_res = 0.7,
                          res_search_dir = self.tmp_dir
                          )
        self.assertEqual(self.cfg.e_value_res, 0.7)
        self.tearDown()
        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': os.path.join(self._data_dir, "base", "test_base.fa"),
                  'db_type': 'gembase',
                  'e_value_res': 0.7,
                  'i_evalue_sel': 1,
                  'res_search_dir': self.tmp_dir
                  }
        self.assertRaises(ValueError, Config, **kwargs)


    def test_hmmer_exe(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir = self.tmp_dir
                          )
        self.assertEqual(self.cfg.hmmer_exe, 'hmmsearch')
        self.tearDown()
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          hmmer_exe = 'truc',
                          res_search_dir = self.tmp_dir
                          )
        self.assertEqual(self.cfg.hmmer_exe, 'truc')


    def test_index_db_exe(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir = self.tmp_dir
                          )
        self.assertEqual(self.cfg.index_db_exe, 'makeblastdb')
        self.tearDown()
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          index_db_exe = 'truc',
                          res_search_dir = self.tmp_dir
                          )
        self.assertEqual(self.cfg.index_db_exe, 'truc')


    def test_i_value_sel(self):
        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': os.path.join(self._data_dir, "base", "test_base.fa"),
                  'db_type': 'gembase',
                  'i_evalue_sel': 'foo',
                  'res_search_dir': self.tmp_dir
                  }
        self.assertRaises(ValueError, Config, **kwargs)
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir = self.tmp_dir
                          )
        self.assertEqual(self.cfg.i_evalue_sel, 0.001)
        self.tearDown()
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          i_evalue_sel = 0.7,
                          res_search_dir = self.tmp_dir
                          )
        self.assertEqual(self.cfg.i_evalue_sel, 0.7)
        self.tearDown()
        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': os.path.join(self._data_dir, "base", "test_base.fa"),
                  'db_type': 'gembase',
                  'e_value_res': 0.7,
                  'i_evalue_sel': 1,
                  'res_search_dir': self.tmp_dir
                  }
        self.assertRaises(ValueError, Config, **kwargs)


    def test_db_type(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir = self.tmp_dir
                          )
        self.assertEqual( self.cfg.db_type, 'gembase')
        self.tearDown()
        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': os.path.join(self._data_dir, "base", "test_base.fa"),
                  'db_type': 'foo',
                  'res_search_dir': self.tmp_dir
                  }
        self.assertRaises(ValueError, Config, **kwargs)


    def test_previous_run(self):
        out_dir = os.path.join(self.tmp_dir, 'macsy_test_config')
        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': os.path.join(self._data_dir, "base", "test_base.fa"),
                  'db_type': 'gembase',
                  'models_dir': os.path.join(self._data_dir, 'models'),
                  'previous_run': 'foo',
                  'res_search_dir': out_dir
                  }
        self.assertRaises(ValueError, Config, **kwargs)
        try:
            shutil.rmtree(out_dir)
        except:
            pass
        self.tearDown()
        try:
            cfg_base = Config(cfg_file="nimportnaoik",
                              sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                              db_type='gembase',
                              models_dir=os.path.join(self._data_dir, 'models'),
                              res_search_dir = self.tmp_dir
                              )
            self.assertIsNone(cfg_base.previous_run)
            cfg_base.save( cfg_base.working_dir )
            #wait
            time.sleep(1)
            new_cfg = Config(previous_run = cfg_base.working_dir)
            self.assertEqual(new_cfg.previous_run, cfg_base.working_dir)
        finally:
            # close loggers filehandles, so they don't block file deletion
            # in shutil.rmtree calls in Windows
            logging.shutdown()
            try:
                shutil.rmtree(new_cfg.working_dir)
            except:
                pass
            try:
                shutil.rmtree(cfg_base.working_dir)
            except:
                pass


    def test_profile_dir(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertIsNone(self.cfg.def_dir)
        self.tearDown()

        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': os.path.join(self._data_dir, "base", "test_base.fa"),
                  'db_type': 'gembase',
                  'profile_dir': 'foo',
                  'res_search_dir': self.tmp_dir
                  }
        self.assertRaises(ValueError, Config, **kwargs)
        self.tearDown()

        profile_dir = os.path.join(tempfile.gettempdir(), 'macsyfinder_PROFILE')
        if not os.path.exists(profile_dir):
            os.mkdir(profile_dir)
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          profile_dir=profile_dir,
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.profile_dir, profile_dir)
        shutil.rmtree(profile_dir)


    def test_profile_suffix(self):
         self.cfg = Config(cfg_file="nimportnaoik",
                           sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type='gembase',
                           res_search_dir = self.tmp_dir
                           )
         self.assertEqual(self.cfg.profile_suffix, '.hmm')
         self.tearDown()
         profile_suffix = 'foo'
         self.cfg = Config(cfg_file="nimportnaoik",
                           sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type='gembase',
                           profile_suffix = profile_suffix,
                           res_search_dir = self.tmp_dir
                           )
         self.assertEqual(self.cfg.profile_suffix, profile_suffix)


    def test_replicon_topology(self):
         self.cfg = Config(cfg_file="nimportnaoik",
                           sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type='gembase',
                           res_search_dir = self.tmp_dir
                           )
         self.assertEqual(self.cfg.replicon_topology, 'circular')
         self.tearDown()
         self.cfg = Config(cfg_file="nimportnaoik",
                           sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type='gembase',
                           replicon_topology = 'linear',
                           res_search_dir = self.tmp_dir
                           )
         self.assertEqual(self.cfg.replicon_topology, 'linear')
         self.tearDown()
         kwargs = {'cfg_file': "nimportnaoik",
                   'sequence_db': os.path.join(self._data_dir, "base", "test_base.fa"),
                   'db_type': 'gembase',
                   'replicon_topology': 'foo',
                   'res_search_dir': self.tmp_dir
         }
         self.assertRaises(ValueError, Config, **kwargs)


    def test_inter_gene_max_space(self):
         inter_gene_max_space = (["T2SS", 32], ['Flagellum', 64])
         self.cfg = Config(cfg_file="nimportnaoik",
                           sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type='gembase',
                           inter_gene_max_space = inter_gene_max_space,
                           res_search_dir = self.tmp_dir
                           )
         self.assertEqual(self.cfg.inter_gene_max_space('T2SS'), 32)
         self.assertEqual(self.cfg.inter_gene_max_space('Flagellum'), 64)
         self.assertIsNone(self.cfg.inter_gene_max_space('Foo'))


    def test_min_genes_required(self):
         min_genes_required = (["T2SS", 32], ['Flagellum', 64])
         self.cfg = Config(cfg_file="nimportnaoik",
                           sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type='gembase',
                           min_genes_required = min_genes_required,
                           res_search_dir = self.tmp_dir
                           )
         self.assertEqual(self.cfg.min_genes_required('T2SS'), 32)
         self.assertEqual(self.cfg.min_genes_required('Flagellum'), 64)
         self.assertIsNone(self.cfg.min_genes_required('Foo'))

    def test_max_nb_genes(self):
         max_nb_genes = (["T2SS", 32], ['Flagellum', 64])
         self.cfg = Config(cfg_file="nimportnaoik",
                           sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type='gembase',
                           max_nb_genes = max_nb_genes,
                           res_search_dir = self.tmp_dir
                           )
         self.assertEqual(self.cfg.max_nb_genes('T2SS'), 32)
         self.assertEqual(self.cfg.max_nb_genes('Flagellum'), 64)
         self.assertIsNone(self.cfg.max_nb_genes('Foo'))


    def test_min_mandatory_genes_required(self):
         min_mandatory_genes_required = (["T2SS", 32], ['Flagellum', 64])
         self.cfg = Config(cfg_file="nimportnaoik",
                           sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type='gembase',
                           min_mandatory_genes_required = min_mandatory_genes_required,
                           res_search_dir = self.tmp_dir
                           )
         self.assertEqual(self.cfg.min_mandatory_genes_required('T2SS'), 32)
         self.assertEqual(self.cfg.min_mandatory_genes_required('Flagellum'), 64)
         self.assertIsNone(self.cfg.min_mandatory_genes_required('Foo'))


    def test_multi_loci(self):
         multi_loci = "T2SS,Flagellum"
         self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          multi_loci = multi_loci,
                          res_search_dir = self.tmp_dir
                          )
         self.assertTrue(self.cfg.multi_loci('T2SS'))
         self.assertTrue(self.cfg.multi_loci('Flagellum'))
         self.assertFalse(self.cfg.multi_loci('Foo'))


    def test_res_extract_suffix(self):
         self.cfg = Config(cfg_file="nimportnaoik",
                           sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type='gembase',
                           res_search_dir = self.tmp_dir
                           )
         self.assertEqual(self.cfg.res_extract_suffix, '.res_hmm_extract')
         self.tearDown()
         res_extract_suffix = 'foo'
         self.cfg = Config(cfg_file="nimportnaoik",
                           sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type='gembase',
                           res_extract_suffix = res_extract_suffix,
                           res_search_dir = self.tmp_dir
                           )
         self.assertEqual(self.cfg.res_extract_suffix, res_extract_suffix)


    def test_res_search_dir(self):
         kwargs = {'cfg_file': "nimportnaoik",
                   'sequence_db': os.path.join(self._data_dir, "base", "test_base.fa"),
                   'db_type': 'gembase',
                   'models_dir' :os.path.join(self._data_dir, 'models'),
                   'profile_dir': os.path.join(self._data_dir, 'profiles'),
                   'res_search_dir': 'foo',
                   'log_file': 'NUL' if platform.system() == 'Windows' else '/dev/null'
         }
         self.assertRaises(ValueError, Config, **kwargs)
         self.tearDown()
         res_search_dir = os.path.join(os.path.dirname(__file__), 'datatest')
         self.cfg = Config(cfg_file="nimportnaoik",
                           sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type='gembase',
                           res_search_dir=res_search_dir,
                           log_file='NUL' if platform.system() == 'Windows' else '/dev/null'
                           )
         self.assertEqual(self.cfg.res_search_dir, res_search_dir)


    def test_out_dir(self):
        out_dir = os.path.join(self.tmp_dir, 'macsy_test_config')
        try:
            shutil.rmtree(out_dir)
        except:
            pass
        os.makedirs(out_dir, 0775)
        f = open(os.path.join(out_dir, 'fake'), 'w')
        f.close()
        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': os.path.join(self._data_dir, "base", "test_base.fa"),
                  'db_type': 'gembase',
                  'out_dir': out_dir,
                  'log_file': 'NUL' if platform.system() == 'Windows' else '/dev/null'
        }
        self.assertRaises(ValueError, Config, **kwargs)
        shutil.rmtree(out_dir)
        self.tearDown()

        os.makedirs(out_dir, 0775)
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type='gembase',
                          out_dir=out_dir,
                          log_file = 'NUL' if platform.system() == 'Windows' else '/dev/null'
                          )
        self.assertEqual(self.cfg.working_dir, out_dir)



    def test_res_search_suffix(self):
         self.cfg = Config(cfg_file="nimportnaoik",
                           sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type='gembase',
                           res_search_dir = self.tmp_dir
                           )
         self.assertEqual(self.cfg.res_search_suffix, '.search_hmm.out')
         self.tearDown()
         res_search_suffix = 'foo'
         self.cfg = Config(cfg_file="nimportnaoik",
                           sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type='gembase',
                           res_search_suffix = res_search_suffix,
                           res_search_dir = self.tmp_dir
                           )
         self.assertEqual(self.cfg.res_search_suffix, res_search_suffix)


    def test_sequence_db(self):
         kwargs = {'cfg_file': "nimportnaoik",
                    'db_type': 'gembase',
                    'worker_nb': '2.3' ,
                    'res_search_dir': self.tmp_dir,
                    'log_file': 'NUL' if platform.system() == 'Windows' else '/dev/null'
         }

         self.assertRaises(ValueError, Config, **kwargs)
         kwargs = {'cfg_file': "nimportnaoik",
                    'sequence_db': "foo",
                    'db_type': 'gembase',
                    'worker_nb': '2.3',
                    'res_search_dir': self.tmp_dir,
                    'log_file': 'NUL' if platform.system() == 'Windows' else '/dev/null'
         }
         self.assertRaises(ValueError, Config, **kwargs)
         sequence_db = os.path.join(self._data_dir, "base", "test_base.fa")
         self.cfg = Config(cfg_file="nimportnaoik",
                            sequence_db = sequence_db,
                            db_type='gembase',
                            res_search_dir = self.tmp_dir,
                            log_file = 'NUL' if platform.system() == 'Windows' else '/dev/null'
                            )
         self.assertEqual(self.cfg.sequence_db, sequence_db)


    def test_worker_nb(self):
         self.cfg = Config(cfg_file="nimportnaoik",
                           sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type='gembase',
                           res_search_dir = self.tmp_dir
                           )
         self.assertEqual(self.cfg.worker_nb, 1)
         self.tearDown()
         self.cfg = Config(cfg_file="nimportnaoik",
                           sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type='gembase',
                           res_search_dir = self.tmp_dir,
                           worker_nb = 2
                           )
         self.assertEqual(self.cfg.worker_nb, 2)
         self.tearDown()
         kwargs = {'cfg_file': "nimportnaoik",
                   'sequence_db': os.path.join(self._data_dir, "base", "test_base.fa"),
                   'db_type': 'gembase',
                   'res_search_dir': self.tmp_dir,
                   'worker_nb': '2.3'
         }
         self.assertRaises(ValueError, Config, **kwargs)


