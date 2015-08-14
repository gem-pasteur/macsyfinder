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
from macsypy.config import ConfigLight
from macsypy.registries import SystemDef


class Test(unittest.TestCase):


    _data_dir = os.path.join(os.path.dirname(__file__), "datatest")
    
    
    def setUp(self):
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        
        #add only one handler to the macsypy logger
        from macsypy.gene import _log
        macsy_log = _log.parent
        log_file = 'NUL' if platform.system() == 'Windows' else '/dev/null'
        log_handler = logging.FileHandler(log_file)
        macsy_log.addHandler(log_handler)
        
        self.cfg = ConfigLight()
        self.tmp_dir = tempfile.mkdtemp()
        self.system_dir = os.path.join(self.tmp_dir, 'systems')
        os.mkdir(self.system_dir)


    def _create_fake_simple_systems_tree(self, sys_def):
        simple_dir = os.path.join(self.system_dir, sys_def['name'])
        os.mkdir(simple_dir)

        profiles_dir = os.path.join(simple_dir, 'profiles')
        os.mkdir(profiles_dir)
        for profile in sys_def['profiles']:
            fh = open(os.path.join(profiles_dir, profile), 'w')
            fh.close()
        for filename in sys_def['not_profiles']:
            fh = open(os.path.join(profiles_dir, filename), 'w')
            fh.close()

        models_dir = os.path.join(simple_dir, 'models')
        os.mkdir(models_dir)
        for model in sys_def['models']:
            fh = open(os.path.join(models_dir, model), 'w')
            fh.close()
        for filename in sys_def['not_models']:
            fh = open(os.path.join(models_dir, filename), 'w')
            fh.close()
        return simple_dir


    def _create_fake_complex_systems_tree(self, system_dir):
        profiles = ()
        models = ()


    def tearDown(self):
        # close loggers filehandles, so they don't block file deletion
        # in shutil.rmtree calls in Windows
        logging.shutdown()
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        try:
            shutil.rmtree(self.tmp_dir)
        except:
            pass


    def test_SystemDef(self):
        simple_system = {'name': 'simple',
                         'profiles': ('prof_1.hmm', 'prof_2.hmm'),
                         'not_profiles': ('not_a_profile', ),
                         'models': ('model_1.xml', 'model_2.xml'),
                         'not_models': ('not_a_model', ),
                         }

        simple_dir = self._create_fake_simple_systems_tree(simple_system)
        sys_def = SystemDef(simple_dir, self.cfg)
        self.assertEqual(sys_def.name, simple_system['name'])
        self.assertEqual(sys_def.path, simple_dir)
        self.assertDictEqual(sys_def._profiles,
                             {os.path.splitext(p)[0]: os.path.join(simple_dir, 'profiles', p) for p in simple_system['profiles']})

        self.assertListEqual(sys_def._models.keys(),
                             [os.path.splitext(m)[0] for m in simple_system['models']])



    def test_get_model(self):
        pass

    def test_models(self):
        pass

    def test_profile(self):
        pass

    def test_str(self):
        pass
    
    
