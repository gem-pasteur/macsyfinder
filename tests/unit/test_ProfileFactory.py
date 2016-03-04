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
from macsypy.gene import profile_factory
from macsypy.gene import Profile
from macsypy.gene import Gene
from macsypy.system import System
from macsypy.config import Config
from macsypy.registries import ProfilesRegistry
from macsypy.macsypy_error import MacsypyError
from tests import MacsyTest


class Test(MacsyTest):

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
        profile_factory._profiles = {}
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass

    def test_get_profile(self):
        system_foo = System(self.cfg, "foo", 10)
        gene_name = 'sctJ_FLG'
        gene = Gene(self.cfg, gene_name, system_foo, self.profile_registry)
        profile = profile_factory.get_profile(gene, self.cfg, self.profile_registry )
        self.assertTrue( isinstance( profile, Profile ))
        self.assertEqual( profile.gene.name, gene_name )

    def test_get_uniq_object(self):
        system_foo = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.profile_registry)
        path = self.profile_registry.get('sctJ_FLG')
        profile1 = profile_factory.get_profile(gene, self.cfg, path )
        profile2 = profile_factory.get_profile(gene, self.cfg, path)
        self.assertEqual( profile1, profile2 )

    def test_unknow_profile(self):
        system_foo = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system_foo, self.profile_registry)
        gene.name = "foo"
        self.assertRaises(MacsypyError, profile_factory.get_profile, gene, self.cfg, self.profile_registry)

