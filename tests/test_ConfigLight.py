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

from macsypy.config_legacy import ConfigLight
from tests import MacsyTest


class TestConfigLight(MacsyTest):

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

    def test_old_data_organization(self):
        cfg = ConfigLight()
        self.assertFalse(cfg.old_data_organization())

    def test_models_dir(self):
        cfg = ConfigLight(cfg_file="nimportnaoik")
        self.assertEqual(cfg.models_dir,
                         os.path.normpath(os.path.join(__file__, '..', '..', 'data', 'models')))

        models_dir = 'foo_bar_baz'
        cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
        try:
            with open(cfg_file, 'w') as f:
                f.write("""
[directories]  
models_dir = {}
""".format(models_dir))
            with self.assertRaises(ValueError) as ctx:
                ConfigLight(cfg_file=cfg_file)
            self.assertEqual(str(ctx.exception), "{}: No such models directory".format(models_dir))
        finally:
            try:
                os.unlink(cfg_file)
            except:
                pass

        models_dir = self.find_data('models')

        cfg = ConfigLight(models_dir=models_dir)
        self.assertEqual(cfg.models_dir, models_dir)

        try:
            with open(cfg_file, 'w') as f:
                f.write("""
[directories]  
models_dir = {}
""".format(models_dir))
            cfg = ConfigLight(cfg_file=cfg_file)
            self.assertEqual(cfg.models_dir, models_dir)
        finally:
            try:
                os.unlink(cfg_file)
            except:
                pass


    def test_profile_suffix(self):
        cfg = ConfigLight(cfg_file="nimportnaoik")
        self.assertEqual(cfg.profile_suffix, '.hmm')
        cfg_file = os.path.join(tempfile.gettempdir(), 'test_macsy_config.cfg')
        profile_suffix = '.foo'
        try:
            with open(cfg_file, 'w') as f:
                f.write("""
[directories]  
profile_suffix = {}
""".format(profile_suffix))
            cfg = ConfigLight(cfg_file=cfg_file)
            self.assertEqual(cfg.profile_suffix, profile_suffix)
        finally:
            try:
                os.unlink(cfg_file)
            except:
                pass

        cfg = ConfigLight(profile_suffix=profile_suffix)
        self.assertEqual(cfg.profile_suffix, profile_suffix)

    def test_previous_run(self):
        try:
            previous_run = os.path.join(tempfile.gettempdir(), 'test_macsy_previous_run')
            os.makedirs(previous_run)
            cfg_file = os.path.join(previous_run, 'macsyfinder.conf')
            profile_suffix = '.foo'
            models_dir = self.find_data('models')

            with self.assertRaises(ValueError) as ctx:
                ConfigLight(previous_run=previous_run)
            self.assertEqual(str(ctx.exception), "No config file found in dir {}".format(previous_run))

            with open(cfg_file, 'w') as f:
                pass
            cfg = ConfigLight(previous_run=previous_run, profile_suffix=profile_suffix, models_dir=models_dir)
            self.assertEqual(cfg.profile_suffix, profile_suffix)
            self.assertEqual(cfg.models_dir, models_dir)

            with open(cfg_file, 'w') as f:
                f.write("""
[directories]  
profile_suffix = {}
models_dir = {}
""".format(profile_suffix, models_dir))
            cfg = ConfigLight(previous_run=previous_run)
            self.assertEqual(cfg.profile_suffix, profile_suffix)
            self.assertEqual(cfg.models_dir, models_dir)
        finally:
            try:
                shutil.rmtree(previous_run)
            except:
                pass



