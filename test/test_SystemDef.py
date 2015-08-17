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
from macsypy.registries import SystemDef, ModelDefLocation


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


    def _create_fake_systems_tree(self, sys_def):
        system_dir = os.path.join(self.system_dir, sys_def['name'])
        os.mkdir(system_dir)

        profiles_dir = os.path.join(system_dir, 'profiles')
        os.mkdir(profiles_dir)
        for profile in sys_def['profiles']:
            fh = open(os.path.join(profiles_dir, profile), 'w')
            fh.close()
        for filename in sys_def['not_profiles']:
            fh = open(os.path.join(profiles_dir, filename), 'w')
            fh.close()

        models_dir = os.path.join(system_dir, 'models')
        os.mkdir(models_dir)

        def create_tree(models, path):
            for level_n in models:
                if isinstance(level_n, str):
                    if models[level_n] is None:
                        fh = open(os.path.join(path, level_n), 'w')
                        fh.close()
                    else:
                        submodel_path = os.path.join(path, level_n)
                        if not os.path.exists(submodel_path):
                            os.mkdir(submodel_path)
                        create_tree(models[level_n], submodel_path)
                elif isinstance(level_n, dict):
                    create_tree(models[level_n], path)
                else:
                    assert False, "error in modeling \"models\" {} {} ".format(level_n, type(level_n))
        create_tree(sys_def['models'], models_dir)
        create_tree(sys_def['not_models'], models_dir)

        return system_dir


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
                         'models': {'model_1.xml': None,
                                    'model_2.xml': None
                                    },
                         'not_models': {'not_a_model': None},
                         }

        simple_dir = self._create_fake_systems_tree(simple_system)
        sys_def = SystemDef(simple_dir, self.cfg)
        self.assertEqual(sys_def.name, simple_system['name'])
        self.assertEqual(sys_def.path, simple_dir)
        self.assertDictEqual(sys_def._profiles,
                             {os.path.splitext(p)[0]: os.path.join(simple_dir, 'profiles', p) for p in simple_system['profiles']})

        self.assertSetEqual(set(sys_def._models.keys()),
                            {os.path.splitext(m)[0] for m in simple_system['models']})

        complex_system = {'name': 'complex',
                          'profiles': ('prof_1.hmm', 'prof_2.hmm'),
                          'not_profiles': ('not_a_profile', ),
                          'models': {'submodel_1': {'model_1_1.xml': None,
                                                    'model_1_2.xml': None
                                                    },
                                     'submodel_2': {'model_2_1.xml': None,
                                                    'model_2_2.xml': None
                                                   },
                                     },
                          'not_models': {'submodel_1': {'not_a_model': None},
                                         'submodel_2': {'not_a_model': None}
                                         },
                          }
        complex_dir = self._create_fake_systems_tree(complex_system)
        sys_def = SystemDef(complex_dir, self.cfg)
        self.assertEqual(sys_def.name, complex_system['name'])
        self.assertEqual(sys_def.path, complex_dir)
        self.assertDictEqual(sys_def._profiles,
                             {os.path.splitext(p)[0]: os.path.join(complex_dir, 'profiles', p) for p in complex_system['profiles']})

        self.assertSetEqual({sm for sm in complex_system['models']}, set(sys_def._models.keys()))
        for submodel_name in complex_system['models']:
            submodel = complex_system['models'][submodel_name]

            self.assertSetEqual({ssm.name for ssm in sys_def._models[submodel_name].submodels.values()},
                                {os.path.splitext(ssm)[0] for ssm in submodel})
            self.assertSetEqual({ssm.path for ssm in sys_def._models[submodel_name].submodels.values()},
                                {os.path.join(complex_dir, 'models', submodel_name,ssm) for ssm in submodel})


    def test_get_model(self):
        simple_system = {'name': 'simple',
                         'profiles': ('prof_1.hmm', 'prof_2.hmm'),
                         'not_profiles': ('not_a_profile', ),
                         'models': {'model_1.xml': None,
                                    'model_2.xml': None
                                    },
                         'not_models': {'not_a_model': None},
                         }

        simple_dir = self._create_fake_systems_tree(simple_system)
        sys_def = SystemDef(simple_dir, self.cfg)

        model_name = 'model_1'
        model_expected = ModelDefLocation(name=model_name,
                                          path=os.path.join(simple_dir, 'models', model_name+'.xml'))
        model_received = sys_def.get_model(model_name)
        self.assertEqual(model_expected, model_received)


    def test_models(self):
        complex_system = {'name': 'complex',
                          'profiles': ('prof_1.hmm', 'prof_2.hmm'),
                          'not_profiles': ('not_a_profile', ),
                          'models': {'submodel_1': {'model_1_1.xml': None,
                                                    'model_1_2.xml': None
                                                    },
                                     'submodel_2': {'model_2_1.xml': None,
                                                    'model_2_2.xml': None
                                                    },
                                     },
                          'not_models': {'submodel_1': {'not_a_model': None},
                                         'submodel_2': {'not_a_model': None}
                                         },
                          }
        complex_dir = self._create_fake_systems_tree(complex_system)
        sys_def = SystemDef(complex_dir, self.cfg)

        models_expected = []
        for model_name in complex_system['models']:
            submodels = {os.path.splitext(m)[0]: ModelDefLocation(name=os.path.splitext(m)[0],
                                                                  path=os.path.join(complex_dir, 'models', model_name, m)) \
                         for m in complex_system['models'][model_name]
                         }
            sub_model = ModelDefLocation(name=model_name,
                                         path=os.path.join(complex_dir, 'models', model_name),
                                         submodels=submodels)
            models_expected.append(sub_model)

        models_received = sys_def.models
        self.assertEqual(len(models_received), len(models_expected))
        self.assertEqual(sorted(models_received), sorted(models_expected))


    def test_get_profile(self):
        simple_system = {'name': 'simple',
                         'profiles': ('prof_1.hmm', 'prof_2.hmm'),
                         'not_profiles': ('not_a_profile', ),
                         'models': {'model_1.xml': None,
                                    'model_2.xml': None
                                    },
                         'not_models': {'not_a_model': None},
                         }

        simple_dir = self._create_fake_systems_tree(simple_system)
        sys_def = SystemDef(simple_dir, self.cfg)

        self.assertEqual(sys_def.get_profile(os.path.splitext(simple_system['profiles'][0])[0]),
                         os.path.join(simple_dir, 'profiles', simple_system['profiles'][0]))

    def test_str(self):
        pass
    
    
