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

from macsypy.config import ConfigLight, Config
from macsypy import registries
from macsypy.registries import ModelLocation, DefinitionLocation, ModelRegistry
from macsypy.macsypy_error import MacsypyError




def _create_fake_models_tree(root_models_dir, sys_def):
    models_dir = os.path.join(root_models_dir, sys_def['name'])
    os.mkdir(models_dir)

    profiles_dir = os.path.join(models_dir, 'profiles')
    os.mkdir(profiles_dir)
    for profile in sys_def['profiles']:
        fh = open(os.path.join(profiles_dir, profile), 'w')
        fh.close()
    for filename in sys_def['not_profiles']:
        fh = open(os.path.join(profiles_dir, filename), 'w')
        fh.close()

    def_dir = os.path.join(models_dir, 'definitions')
    os.mkdir(def_dir)

    def create_tree(definitions, path):
        for level_n in definitions:
            if isinstance(level_n, str):
                if definitions[level_n] is None:
                    fh = open(os.path.join(path, level_n), 'w')
                    fh.close()
                else:
                    subdef_path = os.path.join(path, level_n)
                    if not os.path.exists(subdef_path):
                        os.mkdir(subdef_path)
                    create_tree(definitions[level_n], subdef_path)
            elif isinstance(level_n, dict):
                create_tree(definitions[level_n], path)
            else:
                assert False, "error in modeling \"definitions\" {} {} ".format(level_n, type(level_n))
    create_tree(sys_def['definitions'], def_dir)
    create_tree(sys_def['not_definitions'], def_dir)
    return models_dir


class ModelLocationTest(unittest.TestCase):

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
        self.root_models_dir = os.path.join(self.tmp_dir, 'models')
        os.mkdir(self.root_models_dir)

        self.simple_models = {'name': 'simple',
                              'profiles': ('prof_1.hmm', 'prof_2.hmm'),
                              'not_profiles': ('not_a_profile', ),
                              'definitions': {'def_1.xml': None,
                                              'def_2.xml': None
                                             },
                              'not_definitions': {'not_a_def': None},
                              }

        self.complex_models = {'name': 'complex',
                               'profiles': ('prof_1.hmm', 'prof_2.hmm'),
                               'not_profiles': ('not_a_profile', ),
                               'definitions': {'subdef_1': {'def_1_1.xml': None,
                                                            'def_1_2.xml': None
                                                           },
                                               'subdef_2': {'def_2_1.xml': None,
                                                            'def_2_2.xml': None
                                                           },
                                               },
                               'not_definitions': {'subdef_1': {'not_a_def': None},
                                                   'subdef_2': {'not_a_def': None}
                                                  },
                             }



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


    def test_ModelLocation(self):
        with self.assertRaises(MacsypyError) as cm:
            ModelLocation(self.cfg, path='foo', profile_dir='bar')
        self.assertEqual(cm.exception.message, "'path' and 'profile_dir' are incompatible arguments")

        with self.assertRaises(MacsypyError) as cm:
            ModelLocation(self.cfg, path='foo', def_dir='bar')
        self.assertEqual(cm.exception.message, "'path' and 'def_dir' are incompatible arguments")

        with self.assertRaises(MacsypyError) as cm:
            ModelLocation(self.cfg, def_dir='foo')
        self.assertEqual(cm.exception.message, "if 'profile_dir' is specified 'def_dir' must be specified_too and vice versa")

        with self.assertRaises(MacsypyError) as cm:
            ModelLocation(self.cfg, profile_dir='foo')
        self.assertEqual(cm.exception.message, "if 'profile_dir' is specified 'def_dir' must be specified_too and vice versa")

        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models)
        model_loc = ModelLocation(self.cfg, path=simple_dir)
        self.assertEqual(model_loc.name, self.simple_models['name'])
        self.assertEqual(model_loc.path, simple_dir)
        self.assertDictEqual(model_loc._profiles,
                             {os.path.splitext(p)[0]: os.path.join(simple_dir, 'profiles', p) \
                              for p in self.simple_models['profiles']})

        self.assertSetEqual(set(model_loc._definitions.keys()),
                            {os.path.splitext(m)[0] for m in self.simple_models['definitions']})


        complex_dir = _create_fake_models_tree(self.root_models_dir, self.complex_models)
        model_loc = ModelLocation(self.cfg, path=complex_dir)
        self.assertEqual(model_loc.name, self.complex_models['name'])
        self.assertEqual(model_loc.path, complex_dir)
        self.assertDictEqual(model_loc._profiles,
                             {os.path.splitext(p)[0]: os.path.join(complex_dir, 'profiles', p) \
                              for p in self.complex_models['profiles']})

        self.assertSetEqual({sm for sm in self.complex_models['definitions']}, set(model_loc._definitions.keys()))
        for subdef_name in self.complex_models['definitions']:
            subdef = self.complex_models['definitions'][subdef_name]

            self.assertSetEqual({ssm.name for ssm in model_loc._definitions[subdef_name].subdefinitions.values()},
                                {os.path.splitext(ssm)[0] for ssm in subdef})
            self.assertSetEqual({ssm.path for ssm in model_loc._definitions[subdef_name].subdefinitions.values()},
                                {os.path.join(complex_dir, 'definitions', subdef_name, ssm) for ssm in subdef})



    def test_get_definition(self):
        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models)
        model_loc = ModelLocation(self.cfg, path=simple_dir)

        model_name = os.path.splitext(self.simple_models['definitions'].keys()[0])[0]
        defloc_expected = DefinitionLocation(name=model_name,
                                             path=os.path.join(simple_dir, 'definitions', model_name + '.xml'))
        defloc_received = model_loc.get_definition(model_name)
        self.assertEqual(defloc_expected, defloc_received)

        complex_dir = _create_fake_models_tree(self.root_models_dir, self.complex_models)
        sys_def = ModelLocation(self.cfg, path=complex_dir)

        submodel_name = 'subdef_1'
        model_name = 'def_1_1'
        model_expected = DefinitionLocation(name=model_name,
                                            path=os.path.join(complex_dir, 'definitions', submodel_name, model_name + '.xml'))

        model_received = sys_def.get_definition(submodel_name + '/' + model_name)
        self.assertEqual(model_expected, model_received)


    def test_models(self):
        complex_dir = _create_fake_models_tree(self.root_models_dir, self.complex_models)
        model_loc = ModelLocation(self.cfg, path=complex_dir)

        defs_expected = []
        for def_name in self.complex_models['definitions']:
            subdefinitions = {os.path.splitext(m)[0]: DefinitionLocation(name=os.path.splitext(m)[0],
                                                                    path=os.path.join(complex_dir, 'definitions', def_name, m)) \
                         for m in self.complex_models['definitions'][def_name]
                         }
            sub_def = DefinitionLocation(name=def_name,
                                         path=os.path.join(complex_dir, 'definitions', def_name),
                                         subdefinitions=subdefinitions)
            defs_expected.append(sub_def)

        defs_received = model_loc.definitions
        self.assertEqual(len(defs_expected), len(defs_received))
        self.assertEqual(sorted(defs_expected), sorted(defs_received))


    def test_get_profile(self):
        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models)
        sys_def = ModelLocation(self.cfg, path=simple_dir)

        self.assertEqual(sys_def.get_profile(os.path.splitext(self.simple_models['profiles'][0])[0]),
                         os.path.join(simple_dir, 'profiles', self.simple_models['profiles'][0]))

    def test_str(self):
        pass


class ModelDefLocationTest(unittest.TestCase):


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


    def test_ModelDefLocation(self):
        model_name = 'foo'
        model_path = '/path/to/model.xml'
        mdfl = DefinitionLocation(name=model_name,
                                path=model_path)
        self.assertEqual(mdfl.name, model_name)
        self.assertEqual(mdfl.path, model_path)
        self.assertIsNone(mdfl.subdefinitions)

    def test_add_submodel(self):
        model_name = 'foo'
        model_path = '/path/to/systems/Foo_system/models/foo'
        mdfl = DefinitionLocation(name=model_name,
                                path=model_path)
        submodel_name = 'bar'
        submodel_path = '/path/to/systems/Foo_system/models/foo/bar.xml'
        submodel = DefinitionLocation(name=submodel_name,
                                    path=submodel_path)
        mdfl.add_subdefinition(submodel)

        self.assertEqual(len(mdfl.subdefinitions), 1)
        self.assertEqual(mdfl.subdefinitions[submodel_name], submodel)


    def test_str(self):
        model_name = 'foo'
        model_path = '/path/to/model.xml'
        mdfl = DefinitionLocation(name=model_name,
                                path=model_path)
        self.assertEqual(mdfl.name, model_name)


class SystemRegistryTest(unittest.TestCase):


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
        self.cfg = Config(sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type="gembase",
                          hmmer_exe="",
                          #def_dir=os.path.join(self._data_dir, 'DEF'),
                          res_search_dir=tempfile.gettempdir(),
                          #profile_dir=os.path.join(self._data_dir, 'profiles'),
                          res_extract_suffix="",
                          log_level=30,
                          log_file=log_file)
        self.tmp_dir = tempfile.mkdtemp()
        self._prefix_data_ori = registries._prefix_data
        registries._prefix_data = self.tmp_dir
        self.root_models_dir = os.path.join(self.tmp_dir, 'macsyfinder', 'models')
        os.makedirs(self.root_models_dir)
        simple_models = {'name': 'simple',
                         'profiles': ('prof_1.hmm', 'prof_2.hmm'),
                         'not_profiles': ('not_a_profile', ),
                         'definitions': {'def_1.xml': None,
                                         'def_2.xml': None
                                        },
                         'not_definitions': {'not_a_def': None},
                         }

        complex_models = {'name': 'complex',
                          'profiles': ('prof_1.hmm', 'prof_2.hmm'),
                          'not_profiles': ('not_a_profile', ),
                          'definitions': {'subdef_1': {'def_1_1.xml': None,
                                                       'def_1_2.xml': None
                                                      },
                                          'subdef_2': {'def_2_1.xml': None,
                                                       'def_2_2.xml': None
                                                      },
                                          },
                          'not_definitions': {'subdef_1': {'not_a_def': None},
                                              'subdef_2': {'not_a_def': None}
                                             },
                          }

        self.simple_dir = _create_fake_models_tree(self.root_models_dir, simple_models)
        self.complex_dir = _create_fake_models_tree(self.root_models_dir, complex_models)


    def tearDown(self):
        # close loggers filehandles, so they don't block file deletion
        # in shutil.rmtree calls in Windows
        logging.shutdown()
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass
        try:
            shutil.rmtree(self.tmp_dir)
        except:
            pass
        registries._prefix_data = self._prefix_data_ori


    def test_ModelRegistry(self):
        sr = ModelRegistry(self.cfg)
        self.assertEqual(sorted(sr._registery.keys()), sorted(['simple', 'complex']))

    def test_models(self):
        md = ModelRegistry(self.cfg)
        model_complex_expected = ModelLocation(self.cfg, path=self.complex_dir)
        model_simple_expected = ModelLocation(self.cfg, path=self.simple_dir)
        models_received = md.models()
        self.assertEqual(len(models_received), 2)
        self.assertIn(model_complex_expected, models_received)
        self.assertIn(model_complex_expected, models_received)
