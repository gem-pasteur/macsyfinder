# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Neron                                         #
# Copyright (c) 2014-2019  Institut Pasteur (Paris) and CNRS.                  #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################


import os
import shutil
import tempfile
import argparse

from macsypy.profile import ProfileFactory, Profile
from macsypy.gene import Gene
from macsypy.model import Model
from macsypy.config import Config, MacsyDefaults
from macsypy.registries import ModelLocation
from macsypy.error import MacsypyError
from tests import MacsyTest


class TestProfileFactory(MacsyTest):

    def setUp(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.res_search_dir = tempfile.gettempdir()
        args.log_level = 30
        self.cfg = Config(MacsyDefaults(), args)

        self.model_name = 'foo'
        self.models_location = ModelLocation(path=os.path.join(args.models_dir, self.model_name))
        self.profile_factory = ProfileFactory(self.cfg)

    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass


    def test_get_profile(self):
        system_foo = Model("foo", 10)
        gene_name = 'sctJ_FLG'
        gene = Gene(self.profile_factory, gene_name, system_foo, self.models_location)
        profile = self.profile_factory.get_profile(gene, self.models_location)
        self.assertTrue(isinstance(profile, Profile))
        self.assertEqual(profile.gene.name, gene_name)


    def test_get_uniq_object(self):
        system_foo = Model("foo", 10)
        gene = Gene(self.profile_factory, 'sctJ_FLG', system_foo, self.models_location)
        profile1 = self.profile_factory.get_profile(gene, self.models_location)
        profile2 = self.profile_factory.get_profile(gene, self.models_location)
        self.assertEqual(profile1, profile2)


    def test_unknow_profile(self):
        system_foo = Model("foo", 10)
        gene = Gene(self.profile_factory, 'sctJ_FLG', system_foo, self.models_location)
        gene.name = "foo"
        self.assertRaises(MacsypyError, self.profile_factory.get_profile, gene, self.models_location)
