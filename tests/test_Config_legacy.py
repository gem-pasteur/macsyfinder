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
import shutil
import time
import tempfile
import logging

from macsypy.config import Config
from tests import MacsyTest


class TestConfig(MacsyTest):

    def setUp(self):
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        if hasattr(self, 'cfg'):
            try:
                shutil.rmtree(self.cfg.working_dir)
            except Exception:
                pass
        self.tmp_dir = tempfile.gettempdir()


    def tearDown(self):
        logging.shutdown()
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        if hasattr(self, 'cfg'):
            try:
                shutil.rmtree(self.cfg.working_dir)
            except Exception:
                pass


    def test_build_indexes(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertFalse(self.cfg.build_indexes)
        self.tearDown()

        config = Config(
            cfg_file="nimportnaoik",
            sequence_db=self.find_data("base", "test_base.fa"),
            db_type='gembase',
            res_search_dir=self.tmp_dir,
            build_indexes=True
        )
        self.assertTrue(config.build_indexes)


    def test_default(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.hmmer_exe, 'hmmsearch')


    def test_coverage_profile(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir,
                          )
        self.assertEqual(self.cfg.coverage_profile, 0.5)
        self.tearDown()
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          coverage_profile=0.6,
                          res_search_dir=self.tmp_dir,
                          )
        self.assertEqual(self.cfg.coverage_profile, 0.6)
        self.tearDown()
        # coverage_profile must be a float
        with self.assertRaises(ValueError) as ctx:
            Config(cfg_file="nimportnaoik",
                   sequence_db=self.find_data("base", "test_base.fa"),
                   db_type='gembase',
                   coverage_profile="foo",
                   res_search_dir=self.tmp_dir
                   )
        self.assertEqual(str(ctx.exception), "Invalid value for hmmer coverage_profile :foo: (float expected)")
        self.tearDown()

        coverage_profile = 0.6
        cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
        try:
            with open(cfg_file, 'w') as f:
                f.write("""
[hmmer]
coverage_profile = {}
""".format(coverage_profile))
            cfg = Config(cfg_file=cfg_file,
                         sequence_db=self.find_data("base", "test_base.fa"),
                         db_type='gembase',
                         res_search_dir=self.tmp_dir)
            self.assertEqual(cfg.coverage_profile, coverage_profile)
        finally:
            try:
                os.unlink(cfg_file)
                self.tearDown()
            except:
                pass

            coverage_profile = 'foo'
            cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
            try:
                with open(cfg_file, 'w') as f:
                    f.write("""
[hmmer]
coverage_profile = {}
""".format(coverage_profile))
                with self.assertRaises(ValueError) as ctx:
                    self.cfg = Config(cfg_file=cfg_file,
                                      sequence_db=self.find_data("base", "test_base.fa"),
                                      db_type='gembase',
                                      res_search_dir=self.tmp_dir)
                self.assertEqual(str(ctx.exception), "Invalid value for hmmer coverage_profile :foo: (float expected)")
            finally:
                try:
                    os.unlink(cfg_file)
                    self.tearDown()
                except Exception:
                    pass


    def test_def_dir(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertIsNone(self.cfg.def_dir)
        self.tearDown()

        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': self.find_data("base", "test_base.fa"),
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
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          def_dir=def_dir,
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(def_dir, self.cfg.def_dir)
        shutil.rmtree(def_dir)
        self.tearDown()


    def test_models_dir(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.models_dir,
                         os.path.normpath(os.path.join(__file__, '..', '..', 'data', 'models')))
        self.tearDown()

        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': self.find_data("base", "test_base.fa"),
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
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          models_dir=models_dir,
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(models_dir, self.cfg.models_dir)
        shutil.rmtree(models_dir)
        self.tearDown()


    def test_e_value_res(self):
        e_value_res = 'foo'
        cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
        try:
            with open(cfg_file, 'w') as f:
                f.write("""
[hmmer]
e_value_res = {}
""".format(e_value_res))
            with self.assertRaises(ValueError) as ctx:
                Config(cfg_file=cfg_file,
                       sequence_db=self.find_data("base", "test_base.fa"),
                       db_type='gembase',
                       res_search_dir=self.tmp_dir)
            self.assertEqual(str(ctx.exception), "Invalid value for hmmer e_value_res :foo: (float expected)")
        finally:
            try:
                os.unlink(cfg_file)
                self.tearDown()
            except:
                pass

        e_value_res = 0.9
        try:
            with open(cfg_file, 'w') as f:
                f.write("""
[hmmer]
e_value_res = {}
""".format(e_value_res))
            self.cfg = Config(cfg_file=cfg_file,
                              sequence_db=self.find_data("base", "test_base.fa"),
                              db_type='gembase',
                              res_search_dir=self.tmp_dir)
            self.assertEqual(self.cfg.e_value_res, e_value_res)
        finally:
            try:
                os.unlink(cfg_file)
                self.tearDown()
            except:
                pass

        with self.assertRaises(ValueError) as ctx:
            Config(cfg_file="nimportnaoik",
                   sequence_db=self.find_data("base", "test_base.fa"),
                   db_type='gembase',
                   e_value_res='foo',
                   res_search_dir=self.tmp_dir
                   )
        self.assertEqual(str(ctx.exception), "Invalid value for hmmer e_value_res :foo: (float expected)")
        self.tearDown()

        with self.assertRaises(ValueError) as ctx:
            Config(cfg_file="nimportnaoik",
                   sequence_db=self.find_data("base", "test_base.fa"),
                   db_type='gembase',
                   def_dir=os.path.join(self._data_dir, 'DEF'),
                   profile_dir=os.path.join(self._data_dir, 'profiles'),
                   e_value_res='foo',
                   res_search_dir=self.tmp_dir
                   )
        self.assertEqual(str(ctx.exception), "Invalid value for hmmer e_value_res :foo: (float expected)")
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.e_value_res, 1)
        self.tearDown()
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          e_value_res=0.7,
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.e_value_res, 0.7)
        self.tearDown()
        e_value_res = 0.7
        i_evalue_sel = 1
        with self.assertRaises(ValueError) as ctx:
            Config(cfg_file="nimportnaoik",
                   sequence_db=self.find_data("base", "test_base.fa"),
                   db_type='gembase',
                   e_value_res=e_value_res,
                   i_evalue_sel=i_evalue_sel,
                   res_search_dir=self.tmp_dir
                   )
        self.assertEqual(str(ctx.exception),
                         "i_evalue_sel ({:f}) must be lower or equal than e_value_res ({:f})".format(i_evalue_sel,
                                                                                                     e_value_res))

    def test_hmmer_exe(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.hmmer_exe, 'hmmsearch')
        self.tearDown()
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          hmmer_exe='truc',
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.hmmer_exe, 'truc')


    def test_index_db_exe(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.index_db_exe, 'makeblastdb')
        self.tearDown()
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          index_db_exe='truc',
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.index_db_exe, 'truc')


    def test_i_value_sel(self):
        i_evalue_sel = 'foo'
        cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
        try:
            with open(cfg_file, 'w') as f:
                f.write("""
[hmmer]
i_evalue_sel = {}
""".format(i_evalue_sel))
            with self.assertRaises(ValueError) as ctx:
                Config(cfg_file=cfg_file,
                       sequence_db=self.find_data("base", "test_base.fa"),
                       db_type='gembase',
                       res_search_dir=self.tmp_dir)
            self.assertEqual(str(ctx.exception), "Invalid value for hmmer i_evalue_sel :foo: (float expected)")
        finally:
            try:
                os.unlink(cfg_file)
                self.tearDown()
            except:
                pass

        i_evalue_sel = 0.9
        try:
            with open(cfg_file, 'w') as f:
                f.write("""
[hmmer]
i_evalue_sel = {}
""".format(i_evalue_sel))
            self.cfg = Config(cfg_file=cfg_file,
                              sequence_db=self.find_data("base", "test_base.fa"),
                              db_type='gembase',
                              res_search_dir=self.tmp_dir)
            self.assertEqual(self.cfg.i_evalue_sel, i_evalue_sel)
        finally:
            try:
                os.unlink(cfg_file)
                self.tearDown()
            except:
                pass

        with self.assertRaises(ValueError) as ctx:
            Config(cfg_file="nimportnaoik",
                   sequence_db=self.find_data("base", "test_base.fa"),
                   db_type='gembase',
                   i_evalue_sel='foo',
                   res_search_dir=self.tmp_dir
                   )
        self.assertEqual(str(ctx.exception), "Invalid value for hmmer i_evalue_sel :foo: (float expected)")
        self.tearDown()

        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          i_evalue_sel=0.7,
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.i_evalue_sel, 0.7)
        self.tearDown()

        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.i_evalue_sel, 0.001)
        self.tearDown()
        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': self.find_data("base", "test_base.fa"),
                  'db_type': 'gembase',
                  'e_value_res': 0.7,
                  'i_evalue_sel': 1,
                  'res_search_dir': self.tmp_dir
                  }
        self.assertRaises(ValueError, Config, **kwargs)


    def test_db_type(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.db_type, 'gembase')
        self.tearDown()

        with self.assertRaises(ValueError) as ctx:
            self.cfg = Config(cfg_file="nimportnaoik",
                              sequence_db=self.find_data("base", "test_base.fa"),
                              res_search_dir=self.tmp_dir
                              )
        self.assertEqual(str(ctx.exception), 'You must specify the type of the input dataset '
                                             '(unordered_replicon, ordered_replicon, gembase, unordered).')
        self.tearDown()

        with self.assertRaises(ValueError) as ctx:
            self.cfg = Config(cfg_file="nimportnaoik",
                             sequence_db=self.find_data("base", "test_base.fa"),
                             db_type='foo',
                             res_search_dir=self.tmp_dir
                             )
        self.assertEqual(str(ctx.exception), 'Allowed values for the input dataset are : '
                                             'unordered_replicon, ordered_replicon, gembase, unordered')
        self.tearDown()

    def test_previous_run(self):
        with self.assertRaises(ValueError) as ctx:
            self.cfg = Config(cfg_file="nimportnaoik",
                              sequence_db=self.find_data("base", "test_base.fa"),
                              db_type='gembase',
                              models_dir=self.find_data('models'),
                              previous_run='foo',
                              res_search_dir=tempfile.gettempdir()
                              )
        self.assertEqual(str(ctx.exception), "No config file found in dir foo")
        self.tearDown()

        try:
            cfg_base = Config(cfg_file="nimportnaoik",
                              sequence_db=self.find_data("base", "test_base.fa"),
                              db_type='gembase',
                              models_dir=self.find_data('models'),
                              res_search_dir=self.tmp_dir
                              )
            self.assertIsNone(cfg_base.previous_run)
            cfg_base.save(cfg_base.working_dir)
            # wait
            time.sleep(1)
            new_cfg = Config(previous_run=cfg_base.working_dir)
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
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertIsNone(self.cfg.def_dir)
        self.tearDown()

        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': self.find_data("base", "test_base.fa"),
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
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          profile_dir=profile_dir,
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.profile_dir, profile_dir)
        shutil.rmtree(profile_dir)


    def test_profile_suffix(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.profile_suffix, '.hmm')
        self.tearDown()
        profile_suffix = 'foo'
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          profile_suffix=profile_suffix,
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.profile_suffix, profile_suffix)


    def test_replicon_topology(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.replicon_topology, 'circular')
        self.tearDown()
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          replicon_topology='linear',
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.replicon_topology, 'linear')
        self.tearDown()
        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': self.find_data("base", "test_base.fa"),
                  'db_type': 'gembase',
                  'replicon_topology': 'foo',
                  'res_search_dir': self.tmp_dir
                  }
        self.assertRaises(ValueError, Config, **kwargs)


    def test_topology_file(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertIsNone(self.cfg.topology_file)
        self.tearDown()

        with self.assertRaises(ValueError) as ctx:
            self.cfg = Config(cfg_file="nimportnaoik",
                              sequence_db=self.find_data("base", "test_base.fa"),
                              db_type='gembase',
                              res_search_dir=self.tmp_dir,
                              topology_file='foo'
                              )
        self.assertEqual(str(ctx.exception), 'topology_file cannot access foo: No such file')
        self.tearDown()

        topo_file = os.path.join(tempfile.gettempdir(), 'test_macsy_topofile')
        open(topo_file, 'w').close()
        try:
            self.cfg = Config(cfg_file="nimportnaoik",
                                  sequence_db=self.find_data("base", "test_base.fa"),
                                  db_type='gembase',
                                  res_search_dir=self.tmp_dir,
                                  topology_file=topo_file
                                  )
            self.assertEqual(self.cfg.topology_file, topo_file)
        finally:
            self.tearDown()
            try:
                os.unlink(topo_file)
            except:
                pass


        topo_file = os.path.join(tempfile.gettempdir(), 'test_macsy_topofile')
        os.mkdir(topo_file)
        try:
            with self.assertRaises(ValueError) as ctx:
                self.cfg = Config(cfg_file="nimportnaoik",
                                  sequence_db=self.find_data("base", "test_base.fa"),
                                  db_type='gembase',
                                  res_search_dir=self.tmp_dir,
                                  topology_file=topo_file
                                  )
            self.assertEqual(str(ctx.exception), "topology_file {} is not a regular file".format(topo_file))
        finally:
            self.tearDown()
            try:
                shutil.rmtree(topo_file)
            except:
                pass

        topo_file = os.path.join(tempfile.gettempdir(), 'test_macsy_topofile')
        open(topo_file, 'w').close()
        cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')

        try:
            with open(cfg_file, 'w') as f:
                f.write("""
[base]
topology_file = {}
""".format(topo_file))
            self.cfg = Config(cfg_file=cfg_file,
                              sequence_db=self.find_data("base", "test_base.fa"),
                              db_type='gembase',
                              res_search_dir=self.tmp_dir,
                              )
            self.assertEqual(self.cfg.topology_file, topo_file)
        finally:
            self.tearDown()
            try:
                os.unlink(topo_file)
            except:
                pass


    def test_inter_gene_max_space(self):
        inter_gene_max_space = (["T2SS", 32], ['Flagellum', 64])
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          inter_gene_max_space=inter_gene_max_space,
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.inter_gene_max_space('T2SS'), 32)
        self.assertEqual(self.cfg.inter_gene_max_space('Flagellum'), 64)
        self.assertIsNone(self.cfg.inter_gene_max_space('Foo'))

        self.tearDown()

        with self.assertRaises(ValueError) as ctx:
            Config(cfg_file="nimportnaoik",
                   sequence_db=self.find_data("base", "test_base.fa"),
                   db_type='gembase',
                   inter_gene_max_space=(["Foo/T2SS", 'blabla'],),
                   res_search_dir=self.tmp_dir
                   )
        self.assertEqual(str(ctx.exception), "The value for 'inter_gene_max_space' option for system {0} must be an integer, "
                                             "but you provided {1} on command line".format('Foo/T2SS', 'blabla'))
        self.tearDown()

        cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
        try:
            with open(cfg_file, 'w') as f:
                f.write("""
[system]
inter_gene_max_space = {}
""".format(' '.join(["{} {}".format(definition, nb) for definition, nb in inter_gene_max_space])))
            self.cfg = Config(cfg_file=cfg_file,
                         sequence_db=self.find_data("base", "test_base.fa"),
                         db_type='gembase',
                         res_search_dir=self.tmp_dir)
            self.assertEqual(self.cfg.inter_gene_max_space(inter_gene_max_space[0][0]), inter_gene_max_space[0][1])
            self.assertEqual(self.cfg.inter_gene_max_space(inter_gene_max_space[1][0]), inter_gene_max_space[1][1])
        finally:
            try:
                os.unlink(cfg_file)
                self.tearDown()
            except:
                pass

            cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
            try:
                with open(cfg_file, 'w') as f:
                    f.write("""
[system]
inter_gene_max_space = Foo/T2SS blabla
""")
                with self.assertRaises(ValueError) as ctx:
                    Config(cfg_file=cfg_file,
                           sequence_db=self.find_data("base", "test_base.fa"),
                           db_type='gembase',
                           res_search_dir=self.tmp_dir)
                self.assertEqual(str(ctx.exception),
                                 "The value for 'inter_gene_max_space' option for system {} must be an integer, "
                                 "but you provided {} in the configuration file".format('Foo/T2SS', 'blabla'))
            finally:
                try:
                    os.unlink(cfg_file)
                    self.tearDown()
                except:
                    pass

            try:
                with open(cfg_file, 'w') as f:
                    f.write("""
[system]
inter_gene_max_space = Foo/T2SS
""")
                with self.assertRaises(ValueError) as ctx:
                    Config(cfg_file=cfg_file,
                           sequence_db=self.find_data("base", "test_base.fa"),
                           db_type='gembase',
                           res_search_dir=self.tmp_dir)
                self.assertEqual(str(ctx.exception),
                                 "Invalid syntax for 'inter_gene_max_space': you must have a list of systems and "
                                 "corresponding 'inter_gene_max_space' separated by spaces")
            finally:
                try:
                    os.unlink(cfg_file)
                    self.tearDown()
                except:
                    pass


    def test_min_genes_required(self):
        min_genes_required = (["T2SS", 32], ['Flagellum', 64])
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          min_genes_required=min_genes_required,
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.min_genes_required('T2SS'), 32)
        self.assertEqual(self.cfg.min_genes_required('Flagellum'), 64)
        self.assertIsNone(self.cfg.min_genes_required('Foo'))

        self.tearDown()

        with self.assertRaises(ValueError) as ctx:
            Config(cfg_file="nimportnaoik",
                   sequence_db=self.find_data("base", "test_base.fa"),
                   db_type='gembase',
                   min_genes_required=(["Foo/T2SS", 'blabla'],),
                   res_search_dir=self.tmp_dir
                   )
        self.assertEqual(str(ctx.exception), "The value for 'min_genes_required' option for system {} must be an integer, "
                                             "but you provided {} on command line".format('Foo/T2SS', 'blabla'))
        self.tearDown()

        cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
        try:
            with open(cfg_file, 'w') as f:
                f.write("""
[system]
min_genes_required = {}
""".format(' '.join(["{} {}".format(definition, nb) for definition, nb in min_genes_required])))
            self.cfg = Config(cfg_file=cfg_file,
                              sequence_db=self.find_data("base", "test_base.fa"),
                              db_type='gembase',
                             res_search_dir=self.tmp_dir)
            self.assertEqual(self.cfg.min_genes_required(min_genes_required[0][0]), min_genes_required[0][1])
            self.assertEqual(self.cfg.min_genes_required(min_genes_required[1][0]), min_genes_required[1][1])
        finally:
            try:
                os.unlink(cfg_file)
                self.tearDown()
            except:
                pass

            cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
            try:
                with open(cfg_file, 'w') as f:
                    f.write("""
[system]
min_genes_required = Foo/T2SS blabla
""")
                with self.assertRaises(ValueError) as ctx:
                    Config(cfg_file=cfg_file,
                           sequence_db=self.find_data("base", "test_base.fa"),
                           db_type='gembase',
                           res_search_dir=self.tmp_dir)
                self.assertEqual(str(ctx.exception),
                                 "The value for 'min_genes_required' option for system {} must be an integer, "
                                 "but you provided {} in the configuration file".format('Foo/T2SS', 'blabla'))
            finally:
                try:
                    os.unlink(cfg_file)
                    self.tearDown()
                except:
                    pass

            cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
            try:
                with open(cfg_file, 'w') as f:
                    f.write("""
[system]
min_genes_required = Foo/T2SS
""")
                with self.assertRaises(ValueError) as ctx:
                    Config(cfg_file=cfg_file,
                           sequence_db=self.find_data("base", "test_base.fa"),
                           db_type='gembase',
                           res_search_dir=self.tmp_dir)
                self.assertEqual(str(ctx.exception),
                                 "Invalid syntax for 'min_genes_required': you must have a list of systems and "
                                 "corresponding 'min_genes_required' separated by spaces")
            finally:
                try:
                    os.unlink(cfg_file)
                    self.tearDown()
                except:
                    pass


    def test_max_nb_genes(self):
        max_nb_genes = (["Foo/T2SS", 32], ['Foo/Flagellum', 64])
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          max_nb_genes=max_nb_genes,
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.max_nb_genes('Foo/T2SS'), 32)
        self.assertEqual(self.cfg.max_nb_genes('Foo/Flagellum'), 64)
        self.assertIsNone(self.cfg.max_nb_genes('nimportnaoik'))

        self.tearDown()

        with self.assertRaises(ValueError) as ctx:
            Config(cfg_file="nimportnaoik",
                   sequence_db=self.find_data("base", "test_base.fa"),
                   db_type='gembase',
                   max_nb_genes=(["Foo/T2SS", 'blabla'],),
                   res_search_dir=self.tmp_dir
                   )
        self.assertEqual(str(ctx.exception), "The value for 'max_nb_genes' option for system {0} must be an integer, "
                                             "but you provided {1} on command line".format('Foo/T2SS', 'blabla'))
        self.tearDown()

        cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
        try:
            with open(cfg_file, 'w') as f:
                f.write("""
[system]
max_nb_genes = {}
""".format(' '.join(["{} {}".format(definition, nb) for definition, nb in max_nb_genes])))
            self.cfg = Config(cfg_file=cfg_file,
                         sequence_db=self.find_data("base", "test_base.fa"),
                         db_type='gembase',
                         res_search_dir=self.tmp_dir)
            self.assertEqual(self.cfg.max_nb_genes(max_nb_genes[0][0]),max_nb_genes[0][1])
            self.assertEqual(self.cfg.max_nb_genes(max_nb_genes[1][0]),max_nb_genes[1][1])
        finally:
            try:
                os.unlink(cfg_file)
                self.tearDown()
            except:
                pass

            cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
            try:
                with open(cfg_file, 'w') as f:
                    f.write("""
[system]
max_nb_genes = Foo/T2SS blabla
""")
                with self.assertRaises(ValueError) as ctx:
                    Config(cfg_file=cfg_file,
                           sequence_db=self.find_data("base", "test_base.fa"),
                           db_type='gembase',
                           res_search_dir=self.tmp_dir)
                self.assertEqual(str(ctx.exception), "The value for 'max_nb_genes' option for system {0} must be an integer, "
                                                     "but you provided {1} in the configuration file".format('Foo/T2SS', 'blabla'))
            finally:
                try:
                    os.unlink(cfg_file)
                    self.tearDown()
                except:
                    pass

            cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
            try:
                with open(cfg_file, 'w') as f:
                    f.write("""
[system]
max_nb_genes = Foo/T2SS
""")
                with self.assertRaises(ValueError) as ctx:
                    Config(cfg_file=cfg_file,
                           sequence_db=self.find_data("base", "test_base.fa"),
                           db_type='gembase',
                           res_search_dir=self.tmp_dir)
                self.assertEqual(str(ctx.exception),
                                 "Invalid syntax for 'max_nb_genes': you must have a list of systems and "
                                 "corresponding 'max_nb_genes' separated by spaces")
            finally:
                try:
                    os.unlink(cfg_file)
                    self.tearDown()
                except:
                    pass


    def test_min_mandatory_genes_required(self):
        min_mandatory_genes_required = (["T2SS", 32], ['Flagellum', 64])
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          min_mandatory_genes_required=min_mandatory_genes_required,
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.min_mandatory_genes_required('T2SS'), 32)
        self.assertEqual(self.cfg.min_mandatory_genes_required('Flagellum'), 64)
        self.assertIsNone(self.cfg.min_mandatory_genes_required('Foo'))

        self.tearDown()

        with self.assertRaises(ValueError) as ctx:
            Config(cfg_file="nimportnaoik",
                   sequence_db=self.find_data("base", "test_base.fa"),
                   db_type='gembase',
                   min_mandatory_genes_required=(["Foo/T2SS", 'blabla'],),
                   res_search_dir=self.tmp_dir
                   )
        self.assertEqual(str(ctx.exception), "The value for 'min_mandatory_genes_required' option for system {0} must be an integer, "
                                             "but you provided {1} on command line".format('Foo/T2SS', 'blabla'))
        self.tearDown()

        cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
        try:
            with open(cfg_file, 'w') as f:
                f.write("""
[system]
min_mandatory_genes_required = {}
""".format(' '.join(["{} {}".format(definition, nb) for definition, nb in min_mandatory_genes_required])))
            cfg = Config(cfg_file=cfg_file,
                         sequence_db=self.find_data("base", "test_base.fa"),
                         db_type='gembase',
                         res_search_dir=self.tmp_dir)
            self.assertEqual(cfg.min_mandatory_genes_required(min_mandatory_genes_required[0][0]), min_mandatory_genes_required[0][1])
            self.assertEqual(cfg.min_mandatory_genes_required(min_mandatory_genes_required[1][0]), min_mandatory_genes_required[1][1])
        finally:
            try:
                os.unlink(cfg_file)
            except:
                pass

            self.tearDown()

            cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
            try:
                with open(cfg_file, 'w') as f:
                    f.write("""
[system]
min_mandatory_genes_required = Foo/T2SS blabla
""")
                with self.assertRaises(ValueError) as ctx:
                    Config(cfg_file=cfg_file,
                           sequence_db=self.find_data("base", "test_base.fa"),
                           db_type='gembase',
                           res_search_dir=self.tmp_dir)
                self.assertEqual(str(ctx.exception),
                                 "The value for 'min_mandatory_genes_required' option for system {} must be an integer, "
                                 "but you provided {} in the configuration file".format('Foo/T2SS', 'blabla'))
            finally:
                try:
                    os.unlink(cfg_file)
                except:
                    pass

            self.tearDown()

            cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
            try:
                with open(cfg_file, 'w') as f:
                    f.write("""
[system]
min_mandatory_genes_required = Foo/T2SS
""")
                with self.assertRaises(ValueError) as ctx:
                    Config(cfg_file=cfg_file,
                           sequence_db=self.find_data("base", "test_base.fa"),
                           db_type='gembase',
                           res_search_dir=self.tmp_dir)
                self.assertEqual(str(ctx.exception),
                                 "Invalid syntax for 'min_mandatory_genes_required': you must have a list of systems and "
                                 "corresponding 'min_mandatory_genes_required' separated by spaces")
            finally:
                try:
                    os.unlink(cfg_file)
                except:
                    pass

        self.tearDown()


    def test_multi_loci(self):
        multi_loci = "Foo/T2SS,Bar/Flagellum"
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          multi_loci=multi_loci,
                          res_search_dir=self.tmp_dir
                          )
        self.assertTrue(self.cfg.multi_loci('Foo/T2SS'))
        self.assertTrue(self.cfg.multi_loci('Bar/Flagellum'))
        self.assertFalse(self.cfg.multi_loci('Foo'))
        self.tearDown()
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertFalse(self.cfg.multi_loci('Foo'))

    def test_res_extract_suffix(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.res_extract_suffix, '.res_hmm_extract')
        self.tearDown()
        res_extract_suffix = 'foo'
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_extract_suffix=res_extract_suffix,
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.res_extract_suffix, res_extract_suffix)


    def test_res_search_dir(self):
        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': self.find_data("base", "test_base.fa"),
                  'db_type': 'gembase',
                  'models_dir': self.find_data('models'),
                  'res_search_dir': 'foo',
                  'log_file': os.devnull
                  }
        self.assertRaises(ValueError, Config, **kwargs)
        self.tearDown()
        res_search_dir = tempfile.gettempdir()
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=res_search_dir,
                          log_file=os.devnull
                          )
        self.assertEqual(self.cfg.res_search_dir, res_search_dir)


    def test_out_dir(self):
        out_dir = os.path.join(self.tmp_dir, 'macsy_test_config')
        try:
            shutil.rmtree(out_dir)
        except:
            pass
        os.makedirs(out_dir, 0o775)
        f = open(os.path.join(out_dir, 'fake'), 'w')
        f.close()
        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': self.find_data("base", "test_base.fa"),
                  'db_type': 'gembase',
                  'out_dir': out_dir,
                  'log_file': os.devnull
                  }
        self.assertRaises(ValueError, Config, **kwargs)
        shutil.rmtree(out_dir)
        self.tearDown()

        os.makedirs(out_dir, 0o775)
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          out_dir=out_dir,
                          log_file=os.devnull
                          )
        self.assertEqual(self.cfg.working_dir, out_dir)


    def test_res_search_suffix(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.res_search_suffix, '.search_hmm.out')
        self.tearDown()
        res_search_suffix = 'foo'
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_suffix=res_search_suffix,
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.res_search_suffix, res_search_suffix)


    def test_sequence_db(self):
        kwargs = {'cfg_file': "nimportnaoik",
                  'db_type': 'gembase',
                  'worker_nb': '2.3',
                  'res_search_dir': self.tmp_dir,
                  'log_file': os.devnull
                  }

        self.assertRaises(ValueError, Config, **kwargs)
        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': "foo",
                  'db_type': 'gembase',
                  'worker_nb': '2.3',
                  'res_search_dir': self.tmp_dir,
                  'log_file': os.devnull
                  }
        self.assertRaises(ValueError, Config, **kwargs)
        sequence_db = self.find_data("base", "test_base.fa")
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=sequence_db,
                          db_type='gembase',
                          res_search_dir=self.tmp_dir,
                          log_file=os.devnull
                          )
        self.assertEqual(self.cfg.sequence_db, sequence_db)


    def test_worker_nb(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertEqual(self.cfg.worker_nb, 1)
        self.tearDown()
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir,
                          worker_nb=2
                          )
        self.assertEqual(self.cfg.worker_nb, 2)
        self.tearDown()
        kwargs = {'cfg_file': "nimportnaoik",
                  'sequence_db': self.find_data("base", "test_base.fa"),
                  'db_type': 'gembase',
                  'res_search_dir': self.tmp_dir,
                  'worker_nb': '2.3'
                  }
        self.assertRaises(ValueError, Config, **kwargs)

    def test_relative_path(self):
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir
                          )
        self.assertFalse(self.cfg.relative_path)
        self.tearDown()
        self.cfg = Config(cfg_file="nimportnaoik",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type='gembase',
                          res_search_dir=self.tmp_dir,
                          relative_path=True
                          )
        self.assertTrue(self.cfg.relative_path)
