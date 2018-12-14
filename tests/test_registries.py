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
import tempfile
import logging
import copy
import unittest
import imp

from macsypy.config import ConfigLight, Config
from macsypy import registries
from macsypy.registries import ModelLocation, DefinitionLocation, ModelRegistry
from macsypy.macsypy_error import MacsypyError
from tests import MacsyTest


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


class RegitriesUtilsTest(MacsyTest):

    def test_split_def_name(self):
        items = ['CRISPR-Cas', 'typing', 'cas']
        def_name = registries._separator.join(items)
        split = registries.split_def_name(def_name)
        self.assertListEqual(split, items)
        def_name = registries._separator.join(items) + registries._separator
        split = registries.split_def_name(def_name)
        self.assertListEqual(split, items)
        def_name = registries._separator + registries._separator.join(items)
        split = registries.split_def_name(def_name)
        self.assertListEqual(split, items)


    def test_join_def_path(self):
        items = ['CRISPR-Cas', 'typing', 'cas']
        self.assertEqual('/'.join(items), registries.join_def_path(*items))

class ModelLocationTest(MacsyTest):

    def setUp(self):
        l = logging.getLogger()
        l.manager.loggerDict.clear()

        # add only one handler to the macsypy logger
        from macsypy.gene import _log
        macsy_log = _log.parent
        log_file = os.devnull
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
        self.assertEqual(str(cm.exception),
                         "'path' and 'profile_dir' are incompatible arguments")

        with self.assertRaises(MacsypyError) as cm:
            ModelLocation(self.cfg, path='foo', def_dir='bar')
        self.assertEqual(str(cm.exception),
                         "'path' and 'def_dir' are incompatible arguments")

        with self.assertRaises(MacsypyError) as cm:
            ModelLocation(self.cfg, def_dir='foo')
        self.assertEqual(str(cm.exception),
                         "if 'profile_dir' is specified 'def_dir' must be specified_too and vice versa")

        with self.assertRaises(MacsypyError) as cm:
            ModelLocation(self.cfg, profile_dir='foo')
        self.assertEqual(str(cm.exception),
                         "if 'profile_dir' is specified 'def_dir' must be specified_too and vice versa")

        # test new way to specify profiles and defitions
        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models)
        model_loc = ModelLocation(self.cfg, path=simple_dir)
        self.assertEqual(model_loc.name, self.simple_models['name'])
        self.assertEqual(model_loc.path, simple_dir)
        self.assertDictEqual(model_loc._profiles,
                             {os.path.splitext(p)[0]: os.path.join(simple_dir, 'profiles', p)
                              for p in self.simple_models['profiles']})

        self.assertSetEqual(set(model_loc._definitions.keys()),
                            {os.path.splitext(m)[0] for m in self.simple_models['definitions']})

        complex_dir = _create_fake_models_tree(self.root_models_dir, self.complex_models)
        model_loc = ModelLocation(self.cfg, path=complex_dir)
        self.assertEqual(model_loc.name, self.complex_models['name'])
        self.assertEqual(model_loc.path, complex_dir)
        self.assertDictEqual(model_loc._profiles,
                             {os.path.splitext(p)[0]: os.path.join(complex_dir, 'profiles', p)
                              for p in self.complex_models['profiles']})

        self.assertSetEqual({sm for sm in self.complex_models['definitions']}, set(model_loc._definitions.keys()))
        for subdef_name in self.complex_models['definitions']:
            subdef = self.complex_models['definitions'][subdef_name]
            self.assertSetEqual({ssm.name for ssm in model_loc._definitions[subdef_name].subdefinitions.values()},
                                {os.path.splitext(ss_name)[0] for ss_name in subdef.keys()})

            self.assertSetEqual({ssm.fqn for ssm in model_loc._definitions[subdef_name].subdefinitions.values()},
                                {"{m_name}{sep}{sub_d_name}{sep}{d_name}".format(m_name=self.complex_models['name'],
                                                                                 sep=registries._separator,
                                                                                 sub_d_name=subdef_name,
                                                                                 d_name=os.path.splitext(ss_name)[0]) for ss_name in subdef.keys()})

            self.assertSetEqual({ssm.path for ssm in model_loc._definitions[subdef_name].subdefinitions.values()},
                                {os.path.join(complex_dir, 'definitions', subdef_name, ssm) for ssm in subdef})

        # test old way to specify profiles and defitions
        model_loc = ModelLocation(self.cfg,
                                  profile_dir=os.path.join(simple_dir, 'profiles'),
                                  def_dir=os.path.join(simple_dir, 'definitions'))

        self.assertDictEqual(model_loc._profiles,
                             {os.path.splitext(p)[0]: os.path.join(simple_dir, 'profiles', p)
                              for p in self.simple_models['profiles']})

        self.assertSetEqual(set(model_loc._definitions.keys()),
                            {os.path.splitext(m)[0] for m in self.simple_models['definitions']})


    def test_get_definition(self):
        # test new way to specify profiles and defitions
        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models)
        model_loc = ModelLocation(self.cfg, path=simple_dir)

        def_fqn = '{}{}{}'.format(model_loc.name,
                                  registries._separator,
                                  os.path.splitext(list(self.simple_models['definitions'].keys())[0])[0])
        defloc_expected_name = os.path.splitext(list(self.simple_models['definitions'].keys())[0])[0]
        defloc_expected = DefinitionLocation(name=defloc_expected_name,
                                             path=os.path.join(simple_dir, 'definitions', defloc_expected_name + '.xml'))
        defloc_expected.fqn = "{}{}{}".format(model_loc.name,
                                              registries._separator,
                                              os.path.splitext(list(self.simple_models['definitions'].keys())[0])[0])

        defloc_received = model_loc.get_definition(def_fqn)
        self.assertEqual(defloc_expected, defloc_received)

        complex_dir = _create_fake_models_tree(self.root_models_dir, self.complex_models)
        model_loc = ModelLocation(self.cfg, path=complex_dir)

        subdef_name = 'subdef_1'
        def_name = 'def_1_1'
        def_fqn = '{model_name}{sep}{subdef_name}{sep}{def_name}'.format(
                                                                    model_name=model_loc.name,
                                                                    sep=registries._separator,
                                                                    subdef_name=subdef_name,
                                                                    def_name=def_name)

        defloc_expected = DefinitionLocation(name=def_name,
                                             path=os.path.join(complex_dir, 'definitions', subdef_name, def_name + '.xml'))
        defloc_expected.fqn = def_fqn
        defloc_received = model_loc.get_definition(def_fqn)
        self.assertEqual(defloc_expected, defloc_received)

        # test old way to specify profiles and defitions
        model_loc = ModelLocation(self.cfg,
                                  profile_dir=os.path.join(simple_dir, 'profiles'),
                                  def_dir=os.path.join(simple_dir, 'definitions'))

        function_orig = self.cfg.old_data_organization
        self.cfg.old_data_organization = lambda : True
        def_fqn = "definitions/{0}".format(os.path.splitext(list(self.simple_models['definitions'].keys())[0])[0])
        def_name = os.path.splitext(list(self.simple_models['definitions'].keys())[0])[0]
        defloc_expected = DefinitionLocation(name=def_name,
                                             path=os.path.join(simple_dir, 'definitions', def_name + '.xml'))
        defloc_received = model_loc.get_definition(def_fqn)

        self.assertEqual(defloc_expected, defloc_received)
        self.cfg.old_data_organization = function_orig


    def test_get_all_definitions(self):
        complex_dir = _create_fake_models_tree(self.root_models_dir, self.complex_models)
        model_loc = ModelLocation(self.cfg, path=complex_dir)

        defs_expected = []
        for def_name in self.complex_models['definitions']:
            for subdef_name in self.complex_models['definitions'][def_name]:
                new_def = DefinitionLocation(name=os.path.splitext(subdef_name)[0],
                                             path=os.path.join(complex_dir, 'definitions', def_name, subdef_name))
                new_def.fqn = '{model_name}{sep}{def_name}{sep}{subdef_name}'.format(
                                                                    model_name=model_loc.name,
                                                                    sep=registries._separator,
                                                                    subdef_name=os.path.splitext(subdef_name)[0],
                                                                    def_name=def_name)
                defs_expected.append(new_def)

        defs_received = model_loc.get_all_definitions()
        self.assertEqual(sorted(defs_expected), sorted(defs_received))

        defs_expected = []
        def_root_name = 'subdef_1'
        for def_name in self.complex_models['definitions'][def_root_name]:
            new_def = DefinitionLocation(name=os.path.splitext(def_name)[0],
                                         path=os.path.join(complex_dir, 'definitions', def_root_name, def_name))
            new_def.fqn = '{model_name}{sep}{def_root_name}{sep}{def_name}'.format(
                                                                    model_name=model_loc.name,
                                                                    sep=registries._separator,
                                                                    def_root_name=def_root_name,
                                                                    def_name=os.path.splitext(def_name)[0])
            defs_expected.append(new_def)

        defs_received = model_loc.get_all_definitions(root_def_name="{model_name}{sep}{def_root_name}".format(
                                                                                           model_name=model_loc.name,
                                                                                           sep=registries._separator,
                                                                                           def_root_name=def_root_name))
        self.assertEqual(sorted(defs_expected), sorted(defs_received))

        # test old way to specify profiles and defitions
        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models)
        model_loc = ModelLocation(self.cfg,
                                  profile_dir=os.path.join(simple_dir, 'profiles'),
                                  def_dir=os.path.join(simple_dir, 'definitions'))
        defs_expected = [DefinitionLocation(name=os.path.splitext(d)[0],
                                            path=os.path.join(simple_dir, 'definitions', d))
                         for d in self.simple_models['definitions']
                         ]
        defs_received = model_loc.get_all_definitions()
        self.assertEqual(sorted(defs_expected), sorted(defs_received))

        with self.assertRaises(ValueError) as ctx:
            model_loc.get_all_definitions(root_def_name='foobar')
        self.assertEqual(str(ctx.exception), "root_def_name foobar does not match with any definitions")


    def test_get_profile(self):
        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models)
        model_loc = ModelLocation(self.cfg, path=simple_dir)

        self.assertEqual(model_loc.get_profile(os.path.splitext(self.simple_models['profiles'][0])[0]),
                         os.path.join(simple_dir, 'profiles', self.simple_models['profiles'][0]))

        model_loc = ModelLocation(self.cfg,
                                  profile_dir=os.path.join(simple_dir, 'profiles'),
                                  def_dir=os.path.join(simple_dir, 'definitions'))
        self.assertEqual(model_loc.get_profile(os.path.splitext(self.simple_models['profiles'][0])[0]),
                         os.path.join(simple_dir, 'profiles', self.simple_models['profiles'][0]))


    def test_str(self):
        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models)
        model_loc = ModelLocation(self.cfg, path=simple_dir)
        model_loc.name = 'foo20'
        self.assertEqual('foo20', str(model_loc))

    def test_get_definitions(self):
        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models)
        model_loc = ModelLocation(self.cfg, path=simple_dir)

        model_loc._definitions = None
        defs = model_loc.get_definitions()
        self.assertEqual({}, defs)

        model_loc._definitions = {'foo': 'bar'}
        defs = model_loc.get_definitions()
        self.assertListEqual(['bar'], defs)


class DefinitionLocationTest(MacsyTest):


    def setUp(self):
        l = logging.getLogger()
        l.manager.loggerDict.clear()

        # add only one handler to the macsypy logger
        from macsypy.gene import _log
        macsy_log = _log.parent
        log_file = os.devnull
        log_handler = logging.FileHandler(log_file)
        macsy_log.addHandler(log_handler)


    def test_DefinitionLocationn(self):
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
        self.assertEqual('foo', str(mdfl))


    def test_all(self):
        model_name = 'foo'
        model_path = '/path/to/systems/Foo_system/models/foo'
        mdfl = DefinitionLocation(name=model_name,
                                  path=model_path)
        submodel_name1 = 'foo/bar'
        submodel_path1 = '/path/to/systems/Foo_system/models/foo/bar.xml'
        submodel1 = DefinitionLocation(name=submodel_name1,
                                       path=submodel_path1)
        mdfl.add_subdefinition(submodel1)
        submodel_name2 = 'foo/baz'
        submodel_path2 = '/path/to/systems/Foo_system/models/foo/baz.xml'
        submodel2 = DefinitionLocation(name=submodel_name2,
                                       path=submodel_path2)
        mdfl.add_subdefinition(submodel2)
        self.assertListEqual(sorted(mdfl.all()), sorted([submodel1, submodel2]))


class ModelRegistryTest(MacsyTest):


    def setUp(self):
        l = logging.getLogger()
        l.manager.loggerDict.clear()

        # add only one handler to the macsypy logger
        from macsypy.gene import _log
        macsy_log = _log.parent
        log_file = os.devnull
        log_handler = logging.FileHandler(log_file)
        macsy_log.addHandler(log_handler)
        self.cfg = Config(sequence_db=self.find_data("base", "test_base.fa"),
                          db_type="gembase",
                          hmmer_exe="",
                          # def_dir=self.find_data('DEF'),
                          res_search_dir=tempfile.gettempdir(),
                          # profile_dir=self.find_data('profiles'),
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

        # new way models organization
        sr = ModelRegistry(self.cfg)
        self.assertEqual(sorted(sr._registry.keys()), sorted(['simple', 'complex']))

        # old way models organization
        profile_dir = os.path.join(self.simple_dir, 'profiles')
        def_dir = os.path.join(self.simple_dir, 'definitions')

        cfg_old = copy.copy(self.cfg)
        cfg_old.old_data_organization = lambda : True
        cfg_old.options['profile_dir'] = profile_dir
        cfg_old.options['def_dir'] = def_dir

        sr = ModelRegistry(cfg_old)
        self.assertListEqual(list(sr._registry.keys()), ['definitions'])

    def test_models(self):
        md = ModelRegistry(self.cfg)
        model_complex_expected = ModelLocation(self.cfg, path=self.complex_dir)
        model_simple_expected = ModelLocation(self.cfg, path=self.simple_dir)
        models_received = md.models()
        self.assertEqual(len(models_received), 2)
        self.assertIn(model_complex_expected, models_received)
        self.assertIn(model_complex_expected, models_received)

    def test_str(self):
        sr = ModelRegistry(self.cfg)
        self.assertEqual(str(sr), self.output_control_str('001'))
