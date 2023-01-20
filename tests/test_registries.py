#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2023  Institut Pasteur (Paris) and CNRS.           #
# See the COPYRIGHT file for details                                    #
#                                                                       #
# This file is part of MacSyFinder package.                             #
#                                                                       #
# MacSyFinder is free software: you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or     #
# (at your option) any later version.                                   #
#                                                                       #
# MacSyFinder is distributed in the hope that it will be useful,        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
# GNU General Public License for more details .                         #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with MacSyFinder (COPYING).                                     #
# If not, see <https://www.gnu.org/licenses/>.                          #
#########################################################################


import os
import shutil
import tempfile
import argparse
import yaml
import colorlog

import macsypy
from macsypy.config import Config, MacsyDefaults
from macsypy import registries
from macsypy.registries import ModelLocation, DefinitionLocation, ModelRegistry, scan_models_dir
from tests import MacsyTest


def _create_fake_models_tree(root_models_dir, sys_def, metadata=True):
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
    create_tree(sys_def['definitions'], def_dir)
    create_tree(sys_def['not_definitions'], def_dir)

    if metadata:
        with open(os.path.join(models_dir, "metadata.yml"), 'w') as file:
            metadata = {'vers': '1.0'}
            yaml.dump(metadata, file)

    return models_dir


class RegitriesUtilsTest(MacsyTest):

    def test_split_def_name(self):
        items = ['CRISPR-Cas', 'typing', 'cas']
        def_name = registries._SEPARATOR.join(items)
        split = registries.split_def_name(def_name)
        self.assertListEqual(split, items)
        def_name = registries._SEPARATOR.join(items) + registries._SEPARATOR
        split = registries.split_def_name(def_name)
        self.assertListEqual(split, items)
        def_name = registries._SEPARATOR + registries._SEPARATOR.join(items)
        split = registries.split_def_name(def_name)
        self.assertListEqual(split, items)


    def test_join_def_path(self):
        items = ['CRISPR-Cas', 'typing', 'cas']
        self.assertEqual('/'.join(items), registries.join_def_path(*items))


class ModelLocationTest(MacsyTest):

    def setUp(self):
        macsypy.init_logger()
        logger = colorlog.getLogger('macsypy.registries')
        registries._log = logger
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
                                                            'def_2_2.xml': None,
                                                            'sub_subdef': {'def_2_3_1.xml': None,
                                                                           'def_2_3_2.xml': None}
                                                            },
                                               },
                               'not_definitions': {'subdef_1': {'not_a_def': None},
                                                   'subdef_2': {'not_a_def': None}
                                                   },
                               }


    def tearDown(self):
        try:
            shutil.rmtree(self.tmp_dir)
        except Exception as err:
            pass
        logger = colorlog.getLogger('macsypy.registries')
        del logger.manager.loggerDict['macsypy.registries']
        del logger.manager.loggerDict['macsypy']

    def test_ModelLocation(self):
        # test new way to specify profiles and definitions
        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models)
        model_loc = ModelLocation(path=simple_dir)
        self.assertEqual(model_loc.name, self.simple_models['name'])
        self.assertEqual(model_loc.path, simple_dir)
        self.assertDictEqual(model_loc._profiles,
                             {os.path.splitext(p)[0]: os.path.join(simple_dir, 'profiles', p)
                              for p in self.simple_models['profiles']})

        self.assertSetEqual(set(model_loc._definitions.keys()),
                            {os.path.splitext(m)[0] for m in self.simple_models['definitions']})

        complex_dir = _create_fake_models_tree(self.root_models_dir, self.complex_models)
        model_loc = ModelLocation(path=complex_dir)
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
                                                                                 sep=registries._SEPARATOR,
                                                                                 sub_d_name=subdef_name,
                                                                                 d_name=os.path.splitext(ss_name)[0]) for ss_name in subdef.keys()})

            self.assertSetEqual({ssm.path for ssm in model_loc._definitions[subdef_name].subdefinitions.values()},
                                {os.path.join(complex_dir, 'definitions', subdef_name, ssm) for ssm in subdef})

    def test_version(self):
        # test new way to specify profiles and definitions
        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models)
        model_loc = ModelLocation(path=simple_dir)
        self.assertEqual(model_loc.version, '1.0')


    def test_version_no_metadata(self):
        # test new way to specify profiles and definitions
        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models, metadata=False)
        with self.catch_log(log_name='macsypy') as log:
            model_loc = ModelLocation(path=simple_dir)
            log_msg = log.get_value().strip()
        self.assertEqual(model_loc.version, None)
        self.assertEqual(log_msg,
                         "The models package 'simple' is not versioned contact the package manager to fix it.")


    def test_get_definition(self):
        # test new way to specify profiles and definitions
        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models)
        model_loc = ModelLocation(path=simple_dir)

        def_fqn = '{}{}{}'.format(model_loc.name,
                                  registries._SEPARATOR,
                                  os.path.splitext(list(self.simple_models['definitions'].keys())[0])[0])
        defloc_expected_name = os.path.splitext(list(self.simple_models['definitions'].keys())[0])[0]
        defloc_expected = DefinitionLocation(name=defloc_expected_name,
                                             path=os.path.join(simple_dir, 'definitions', defloc_expected_name + '.xml'))
        defloc_expected.fqn = "{}{}{}".format(model_loc.name,
                                              registries._SEPARATOR,
                                              os.path.splitext(list(self.simple_models['definitions'].keys())[0])[0])

        defloc_received = model_loc.get_definition(def_fqn)
        self.assertEqual(defloc_expected, defloc_received)

        complex_dir = _create_fake_models_tree(self.root_models_dir, self.complex_models)
        model_loc = ModelLocation(path=complex_dir)

        subdef_name = 'subdef_1'
        def_name = 'def_1_1'
        def_fqn = '{model_name}{sep}{subdef_name}{sep}{def_name}'.format(
                                                                    model_name=model_loc.name,
                                                                    sep=registries._SEPARATOR,
                                                                    subdef_name=subdef_name,
                                                                    def_name=def_name)

        defloc_expected = DefinitionLocation(name=def_name,
                                             path=os.path.join(complex_dir, 'definitions', subdef_name, def_name + '.xml'))
        defloc_expected.fqn = def_fqn
        defloc_received = model_loc.get_definition(def_fqn)
        self.assertEqual(defloc_expected, defloc_received)


    def test_get_all_definitions(self):
        complex_dir = _create_fake_models_tree(self.root_models_dir, self.complex_models)
        model_loc = ModelLocation(path=complex_dir)

        defs_expected = []
        for def_name in self.complex_models['definitions']:
            for subdef_name in self.complex_models['definitions'][def_name]:
                if subdef_name.endswith('.xml'):
                    new_def = DefinitionLocation(name=os.path.splitext(subdef_name)[0],
                                                 fqn=os.path.splitext(def_name)[0],
                                                 path=os.path.join(complex_dir, 'definitions', def_name, subdef_name))
                    new_def.fqn = '{model_name}{sep}{def_name}{sep}{subdef_name}'.format(
                                                                        model_name=model_loc.name,
                                                                        sep=registries._SEPARATOR,
                                                                        subdef_name=os.path.splitext(subdef_name)[0],
                                                                        def_name=def_name)
                    defs_expected.append(new_def)
                else:
                    for sub_subdef_name in self.complex_models['definitions'][def_name][subdef_name]:
                        new_def = DefinitionLocation(name=os.path.splitext(sub_subdef_name)[0],
                                                     fqn=os.path.splitext(def_name)[0],
                                                     path=os.path.join(complex_dir, 'definitions', subdef_name, sub_subdef_name))
                        new_def.fqn = '{model_name}{sep}{def_name}{sep}{subdef_name}{sep}{sub_subdef_name}'.format(
                                                                        model_name=model_loc.name,
                                                                        sep=registries._SEPARATOR,
                                                                        subdef_name=subdef_name,
                                                                        sub_subdef_name=os.path.splitext(sub_subdef_name)[0],
                                                                        def_name=def_name)
                        defs_expected.append(new_def)
        defs_received = model_loc.get_all_definitions()
        self.assertEqual(sorted(defs_expected), sorted(defs_received))

        defs_expected = []
        def_root_name = 'subdef_1'
        for def_name in self.complex_models['definitions'][def_root_name]:
            new_def = DefinitionLocation(name=os.path.splitext(def_name)[0],
                                         fqn=os.path.splitext(def_name)[0],
                                         path=os.path.join(complex_dir, 'definitions', def_root_name, def_name))
            new_def.fqn = '{model_name}{sep}{def_root_name}{sep}{def_name}'.format(
                                                                    model_name=model_loc.name,
                                                                    sep=registries._SEPARATOR,
                                                                    def_root_name=def_root_name,
                                                                    def_name=os.path.splitext(def_name)[0])
            defs_expected.append(new_def)

        defs_received = model_loc.get_all_definitions(root_def_name="{model_name}{sep}{def_root_name}".format(
                                                                                           model_name=model_loc.name,
                                                                                           sep=registries._SEPARATOR,
                                                                                           def_root_name=def_root_name))
        self.assertEqual(sorted(defs_expected), sorted(defs_received))


    def test_get_profile(self):
        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models)
        model_loc = ModelLocation(path=simple_dir)

        self.assertEqual(model_loc.get_profile(os.path.splitext(self.simple_models['profiles'][0])[0]),
                         os.path.join(simple_dir, 'profiles', self.simple_models['profiles'][0]))


    def test_get_profiles_names(self):
        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models)
        model_loc = ModelLocation(path=simple_dir)

        self.assertEqual(model_loc.get_profile(os.path.splitext(self.simple_models['profiles'][0])[0]),
                         os.path.join(simple_dir, 'profiles', self.simple_models['profiles'][0]))


    def test_str(self):
        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models)
        model_loc = ModelLocation(path=simple_dir)
        model_loc.name = 'foo20'
        self.assertEqual('foo20', str(model_loc))


    def test_get_definitions(self):
        simple_dir = _create_fake_models_tree(self.root_models_dir, self.simple_models)
        model_loc = ModelLocation(path=simple_dir)

        model_loc._definitions = None
        defs = model_loc.get_definitions()
        self.assertEqual({}, defs)

        model_loc._definitions = {'foo': 'bar'}
        defs = model_loc.get_definitions()
        self.assertListEqual(['bar'], defs)


class DefinitionLocationTest(MacsyTest):

    def test_separator(self):
        self.assertEqual(DefinitionLocation.separator, DefinitionLocation._SEPARATOR)

    def test_split_fqn(self):
        self.assertListEqual(DefinitionLocation.split_fqn('/foo/bar'), ['foo', 'bar'])

    def test_root_name(self):
        self.assertEqual(DefinitionLocation.root_name('/foo/bar'), 'foo')

    def test_family_name(self):
        model_name = 'foo'
        model_path = '/path/to/foo.xml'
        family_name = 'family'
        model_fqn = f"{family_name}/{model_name}"
        mdfl = DefinitionLocation(name=model_name,
                                  fqn=model_fqn,
                                  path=model_path)
        self.assertEqual(mdfl.family_name, family_name)

    def test_hash(self):
        mdfl_1 = DefinitionLocation(name='/foo/model_1',
                                    path='/path/to/model_1.xml')
        mdfl_2 = DefinitionLocation(name='/foo/model_1',
                                    path='/path/to/model_1.xml')
        mdfl_3 = DefinitionLocation(name='/foo/model_3',
                                    path='/path/to/model_3.xml')
        self.assertTrue(isinstance(hash(mdfl_1), int))
        self.assertEqual(hash(mdfl_1), hash(mdfl_2))
        self.assertNotEqual(hash(mdfl_1), hash(mdfl_3))


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
        def_name = 'foo'
        def_fqn = 'foo'
        def_path = '/path/to/model.xml'
        mdfl = DefinitionLocation(name=def_name,
                                  fqn=def_fqn,
                                  path=def_path)
        self.assertEqual('foo', str(mdfl))


    def test_all(self):
        def_name = 'foo'
        def_fqn = 'foo'
        def_path = '/path/to/systems/Foo_system/definitions/foo'
        mdfl = DefinitionLocation(name=def_name,
                                  fqn=def_fqn,
                                  path=def_path)
        submodel_name1 = 'bar'
        def_fqn1 = 'foo/bar'
        submodel_path1 = '/path/to/systems/Foo_system/models/foo/bar.xml'
        submodel1 = DefinitionLocation(name=submodel_name1,
                                       fqn=def_fqn1,
                                       path=submodel_path1)
        mdfl.add_subdefinition(submodel1)
        submodel_name2 = 'baz'
        def_fqn2 = 'foo/baz'
        submodel_path2 = '/path/to/systems/Foo_system/models/foo/baz.xml'
        submodel2 = DefinitionLocation(name=submodel_name2,
                                       fqn=def_fqn2,
                                       path=submodel_path2)
        mdfl.add_subdefinition(submodel2)
        self.assertListEqual(sorted(mdfl.all()), sorted([submodel1, submodel2]))

    def test_eq(self):
        model_name = 'foo'
        model_fqn = 'Foo_system/models/foo'
        model_path = '/path/to/systems/Foo_system/models/foo'
        mdfl_1 = DefinitionLocation(name=model_name,
                                    fqn=model_fqn,
                                    path=model_path)
        mdfl_2 = DefinitionLocation(name=model_name,
                                    fqn=model_fqn,
                                    path=model_path)
        self.assertEqual(mdfl_1, mdfl_2)

    def test_lesser(self):
        model_name = 'aaa'
        model_fqn = 'Foo_system/models/aaa'
        model_path = '/path/to/systems/Foo_system/models/aaa'
        mdfl_1 = DefinitionLocation(name=model_name,
                                    fqn=model_fqn,
                                    path=model_path)
        model_name = 'zzz'
        model_fqn = 'Foo_system/models/zzz'
        model_path = '/path/to/systems/Foo_system/models/zzz'
        mdfl_2 = DefinitionLocation(name=model_name,
                                    fqn=model_fqn,
                                    path=model_path)
        self.assertLess(mdfl_1, mdfl_2)

    def test_greater(self):
        model_name = 'aaa'
        model_fqn = 'Foo_system/models/aaa'
        model_path = '/path/to/systems/Foo_system/models/aaa'
        mdfl_1 = DefinitionLocation(name=model_name,
                                    fqn=model_fqn,
                                    path=model_path)
        model_name = 'zzz'
        model_fqn = 'Foo_system/models/zzz'
        model_path = '/path/to/systems/Foo_system/models/zzz'
        mdfl_2 = DefinitionLocation(name=model_name,
                                    fqn=model_fqn,
                                    path=model_path)
        self.assertGreater(mdfl_2, mdfl_1)


class ModelRegistryTest(MacsyTest):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        registries._prefix_data = self.tmp_dir
        self.root_models_dir = os.path.join(self.tmp_dir, 'macsyfinder', 'models')
        self.cfg = Config(MacsyDefaults(models_dir=self.root_models_dir),
                          argparse.Namespace())

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
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass
        try:
            shutil.rmtree(self.tmp_dir)
        except:
            pass

    def test_scan_models_dir(self):
        models_location = scan_models_dir(self.cfg.models_dir())
        models_location_expected = [
            ModelLocation(path=self.simple_dir,
                          profile_suffix='.hmm',
                          relative_path=False),
            ModelLocation(path=self.complex_dir,
                          profile_suffix='.hmm',
                          relative_path=False),
        ]
        self.assertListEqual(sorted(models_location_expected),
                             sorted(models_location))

    def test_add_get(self):
        mr = ModelRegistry()
        model_complex_expected = ModelLocation(path=self.complex_dir)
        with self.assertRaises(KeyError) as ctx:
            mr[model_complex_expected.name]

        self.assertEqual(str(ctx.exception),
                         '"No such model definition: \'complex\'"')
        mr.add(model_complex_expected)
        self.assertEqual(model_complex_expected, mr[model_complex_expected.name])

    def test_models(self):
        mr = ModelRegistry()
        model_complex_expected = ModelLocation(path=self.complex_dir)
        model_simple_expected = ModelLocation(path=self.simple_dir)
        mr.add(model_complex_expected)
        mr.add(model_simple_expected)
        models_received = sorted(mr.models())
        models_expected = sorted([model_complex_expected, model_simple_expected])
        self.assertListEqual(models_received, models_expected)


    def test_str(self):
        mr = ModelRegistry()
        model_complex_expected = ModelLocation(path=self.complex_dir)
        model_simple_expected = ModelLocation(path=self.simple_dir)
        mr.add(model_complex_expected)
        mr.add(model_simple_expected)
        expected_output = """complex
        /subdef_1
                 /def_1_1
                 /def_1_2
        /subdef_2
                 /def_2_1
                 /def_2_2
simple
       /def_1
       /def_2
"""
        self.assertEqual(expected_output, str(mr))
