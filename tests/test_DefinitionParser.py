#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2020  Institut Pasteur (Paris) and CNRS.           #
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

from macsypy.config import Config, MacsyDefaults
from macsypy.model import ModelBank
from macsypy.gene import GeneBank
from macsypy.profile import ProfileFactory
from macsypy.registries import ModelRegistry, scan_models_dir
from macsypy.definition_parser import DefinitionParser
from macsypy.error import MacsypyError, ModelInconsistencyError
from tests import MacsyTest


class TestModelParser(MacsyTest):

    def setUp(self):
        defaults = MacsyDefaults()
        self.args = argparse.Namespace()
        self.args.sequence_db = self.find_data("base", "test_base.fa")
        self.args.db_type = 'gembase'
        self.args.models_dir = self.find_data('models')
        self.args.res_search_dir = tempfile.gettempdir()

        self.cfg = Config(defaults, self.args)
        self.model_bank = ModelBank()
        self.gene_bank = GeneBank()
        self.profile_factory = ProfileFactory(self.cfg)
        self.model_registry = ModelRegistry()
        models_location = scan_models_dir(self.args.models_dir)
        for ml in models_location:
            self.model_registry.add(ml)
        self.parser = DefinitionParser(self.cfg, self.model_bank, self.gene_bank,
                                       self.profile_factory, self.model_registry)
        
        
    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir())
        except:
            pass

    def test_defintion_to_parse(self):
        parsed = set()
        models_2_detect = set()
        models_2_detect.add('foo/model_1')
        def_2_parse = self.parser.definition_to_parse(models_2_detect, parsed)
        self.assertSetEqual(def_2_parse, {s for s in ('foo/model_1', 'foo/model_2')})
        parsed = set()
        models_2_detect = set()
        definition_name = 'foo/nimportnaoik'
        models_2_detect.add(definition_name)
        with self.assertRaises(MacsypyError) as context:
            with self.catch_log():
                self.parser.definition_to_parse(models_2_detect, parsed)
        self.assertEqual(str(context.exception), f'{definition_name}: No such definition')

        parsed = set()
        models_2_detect = set()
        model_name = 'bar'
        definition = f'{model_name}/nimportnaoik'
        models_2_detect.add(definition)
        with self.assertRaises(MacsypyError) as context:
            with self.catch_log():
                self.parser.definition_to_parse(models_2_detect, parsed)
        self.assertEqual(str(context.exception), f'{model_name}: No such Models in {self.cfg.models_dir()}')

        parsed = set()
        models_2_detect = set()
        model_name = 'foo'
        def_name = 'not_xml'
        fqn = f'{model_name}/{def_name}'
        models_2_detect.add(fqn)
        with self.assertRaises(MacsypyError) as context:
            with self.catch_log():
                self.parser.definition_to_parse(models_2_detect, parsed)
        self.assertTrue(str(context.exception).startswith(f'unable to parse model definition "{fqn}" :'))

    def test_parse_with_homologs(self):
        def_2_parse = set()
        def_2_parse.add('foo/model_1')
        parsed = set()
        models_2_detect = self.parser.definition_to_parse(def_2_parse, parsed)
        self.parser.parse(models_2_detect)
        self.assertEqual(len(self.model_bank), 2)

        m2 = self.model_bank['foo/model_2']
        self.assertEqual(m2.name, 'model_2')
        self.assertEqual(m2.fqn, 'foo/model_2')

        m1 = self.model_bank['foo/model_1']
        self.assertEqual(m1.name, 'model_1')
        self.assertEqual(m1.fqn, 'foo/model_1')
        self.assertEqual(m1.inter_gene_max_space, 20)
        self.assertEqual(m1.min_mandatory_genes_required, 4)
        self.assertEqual(m1.min_genes_required, 6)
        self.assertTrue(m1.multi_loci)
        self.assertFalse(m2.multi_loci)
        self.assertIsNone(m1.max_nb_genes)
        self.assertEqual(m2.max_nb_genes, 1)
        self.assertEqual(len(m1.mandatory_genes), 5)
        mandatory_genes_name = [g.name for g in m1.mandatory_genes]
        mandatory_genes_name.sort()
        theoric_list = ["sctJ_FLG", "sctN_FLG", "flgB", "flgC", "fliE"]
        theoric_list.sort()
        self.assertListEqual(mandatory_genes_name, theoric_list)
        sctJ_FLG = [g for g in m1.mandatory_genes if g.name == 'sctJ_FLG'][0]
        sctJ_FLG_homologs = sctJ_FLG.get_homologs()
        self.assertEqual(len(sctJ_FLG_homologs), 1)
        self.assertEqual(sctJ_FLG_homologs[0].name, 'sctJ')
        sctJ = m2.get_gene('sctJ')
        self.assertTrue(sctJ.exchangeable)
        tadZ = m2.get_gene('tadZ')
        self.assertFalse(tadZ.exchangeable)
        self.assertEqual(sctJ_FLG_homologs[0].model, self.model_bank['foo/model_2'])
        self.assertEqual(len(m1.accessory_genes), 1)
        self.assertEqual(m1.accessory_genes[0].name, 'tadZ')
        self.assertEqual(m1.accessory_genes[0].model, self.model_bank['foo/model_2'])
        self.assertEqual(len(m1.forbidden_genes), 1)
        self.assertEqual(m1.forbidden_genes[0].name, 'sctC')
        self.assertEqual(m1.forbidden_genes[0].model, self.model_bank['foo/model_2'])


    def test_parse_with_analogs(self):
        def_2_parse = set()
        def_2_parse.add('foo/model_3')
        parsed = set()
        models_2_detect = self.parser.definition_to_parse(def_2_parse, parsed)
        self.parser.parse(models_2_detect)

        self.assertEqual(len(self.model_bank), 2)

        m4 = self.model_bank['foo/model_4']
        self.assertEqual(m4.name, 'model_4')
        self.assertEqual(m4.fqn, 'foo/model_4')

        m3 = self.model_bank['foo/model_3']
        self.assertEqual(m3.name, 'model_3')
        self.assertEqual(m3.fqn, 'foo/model_3')
        self.assertEqual(m3.inter_gene_max_space, 20)
        self.assertEqual(m3.min_mandatory_genes_required, 4)
        self.assertEqual(m3.min_genes_required, 6)
        self.assertTrue(m3.multi_loci)
        self.assertFalse(m4.multi_loci)
        self.assertIsNone(m3.max_nb_genes)
        self.assertEqual(m4.max_nb_genes, 1)
        self.assertEqual(len(m3.mandatory_genes), 5)
        mandatory_genes_name = [g.name for g in m3.mandatory_genes]
        mandatory_genes_name.sort()
        theoric_list = ["sctJ_FLG", "sctN_FLG", "flgB", "flgC", "fliE"]
        theoric_list.sort()
        self.assertListEqual(mandatory_genes_name, theoric_list)
        sctJ_FLG = [g for g in m3.mandatory_genes if g.name == 'sctJ_FLG'][0]
        sctJ_FLG_homologs = sctJ_FLG.get_homologs()
        self.assertListEqual(sctJ_FLG_homologs, [])
        sctJ_FLG_analogs = sctJ_FLG.get_analogs()
        self.assertEqual(len(sctJ_FLG_analogs), 1)
        self.assertEqual(sctJ_FLG_analogs[0].name, 'sctJ')
        self.assertEqual(sctJ_FLG_analogs[0].model, self.model_bank['foo/model_4'])
        self.assertEqual(len(m3.accessory_genes), 1)
        self.assertEqual(m3.accessory_genes[0].name, 'tadZ')
        self.assertEqual(m3.accessory_genes[0].model, self.model_bank['foo/model_4'])
        self.assertEqual(len(m3.forbidden_genes), 1)
        self.assertEqual(m3.forbidden_genes[0].name, 'sctC')
        self.assertEqual(m3.forbidden_genes[0].model, self.model_bank['foo/model_4'])


    def test_wo_presence(self):
        model_2_detect = ['foo/fail_wo_presence']
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid model definition 'fail_wo_presence': gene without presence")


    def test_invalid_presence(self):
        model_2_detect = ['foo/fail_invalid_presence']
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid model 'fail_invalid_presence' definition: presence value must be either "
                         "[mandatory, accessory, forbidden] not foo_bar")

    def test_gene_no_name(self):
        model_2_detect = ['foo/gene_no_name']
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid model definition 'gene_no_name': gene without a name")

    def test_invalid_homolog(self):
        model_2_detect = ['foo/invalid_homolog']
        with self.assertRaises(ModelInconsistencyError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid model definition 'invalid_homolog': The gene 'foo_bar' described as "
                         "homolog of 'gspD' in model 'invalid_homolog' is not in the 'GeneBank' gene factory")

    def test_invalid_homolog_2(self):
        model_2_detect = ['foo/invalid_homolog_2']
        with self.assertRaises(SyntaxError) as ctx:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(ctx.exception), "Invalid model definition: gene without name")

    def test_invalid_analog(self):
        model_2_detect = ['foo/invalid_analog']
        with self.assertRaises(ModelInconsistencyError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid model definition 'invalid_analog': The gene 'foo_bar' described as "
                         "analog of 'gspD' in model 'invalid_analog' is not in the 'GeneBank' gene factory")

    def test_invalid_analog_2(self):
        model_2_detect = ['foo/invalid_analog_2']
        with self.assertRaises(SyntaxError) as ctx:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(ctx.exception), "Invalid model definition 'invalid_analog_2': gene without name")

    def test_bad_homolog_sys_ref(self):
        model_2_detect = ['foo/bad_homolog_sys_ref']
        with self.assertRaises(ModelInconsistencyError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "Inconsistency in models definitions: the gene 'sctJ' described as homolog of 'sctN' "
                         "with model_ref 'model_1' has an other model in bank (model_2)")

    def test_bad_analog_sys_ref(self):
        model_2_detect = ['foo/bad_analog_sys_ref']
        with self.assertRaises(ModelInconsistencyError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "Inconsistency in models definitions: the gene 'sctJ' described as analog of 'sctN' "
                         "with model_ref 'model_3' has an other model in bank (model_4)")

    def test_bad_min_genes_required(self):
        model_2_detect = ['foo/bad_min_genes_required']
        with self.assertRaises(ModelInconsistencyError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         'model \'bad_min_genes_required\' is not consistent: min_genes_required 16 must be lesser '
                         'or equal than the number of "accessory" and "mandatory" components in the model: 6')

    def test_bad_min_genes_required_2(self):
        model_2_detect = ['foo/bad_min_genes_required_2']
        with self.catch_log():
            with self.assertRaisesRegex(SyntaxError, "Invalid model definition (.*): "
                                                     "min_genes_required must be an integer: 16.5"):
                self.parser.parse(model_2_detect)

    def test_bad_min_mandatory_genes_required(self):
        model_2_detect = ['foo/bad_min_mandatory_genes_required']
        with self.catch_log():
            with self.assertRaises(ModelInconsistencyError) as context:
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         'model \'bad_min_mandatory_genes_required\' is not consistent: min_genes_required 16 must '
                         'be lesser or equal than the number of "accessory" and "mandatory" components in the model: 6')

    def test_bad_min_mandatory_genes_required_2(self):
        model_2_detect = ['foo/bad_min_mandatory_genes_required_2']
        with self.assertRaises(ModelInconsistencyError) as context:
            with self.catch_log():
                # error raised by System initialization
                # which occur before check_consistency
                # the last test : not(model.min_mandatory_genes_required <= model.min_genes_required)
                # seems to be useless
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "foo/bad_min_mandatory_genes_required_2: min_genes_required '6' must be greater or equal"
                         " than min_mandatory_genes_required '8'")

    def test_bad_min_mandatory_genes_required_4(self):
        model_2_detect = ['foo/bad_min_mandatory_genes_required_4']
        with self.assertRaisesRegex(SyntaxError, "Invalid model definition (.*): "
                                                  "min_mandatory_genes_required must be an integer: 12.5"):
            with self.catch_log():
                self.parser.parse(model_2_detect)


    def test_min_mandatory_genes_required_lesser_than_mandatory_genes(self):
        model_2_detect = ['foo/bad_min_mandatory_genes_required_3']
        with self.assertRaises(ModelInconsistencyError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "model 'bad_min_mandatory_genes_required_3' is not consistent:"
                         " 'min_mandatory_genes_required': 6 must be lesser or equal than the number of 'mandatory' "
                         "components in the model: 5")


    def test_bad_max_nb_genes(self):
        model_2_detect = 'foo/bad_max_nb_genes'
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse([model_2_detect])
        model_name, def_name = model_2_detect.split('/')
        self.assertEqual(str(context.exception),
                         "Invalid model definition ({0}.xml): max_nb_genes must be an integer: HOHOHO".format(
                             os.path.join(self.cfg.models_dir(),
                                          model_name,
                                          'definitions',
                                          def_name)))


    def test_bad_inter_gene_max_space(self):
        model_2_detect = ['foo/bad_inter_gene_max_space']
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid model definition ({}): inter_gene_max_space must be an integer: 12.5".format(
                             os.path.join(self.cfg.models_dir(), "foo/definitions/bad_inter_gene_max_space.xml")
                         )
                         )

    def test_no_inter_gene_max_space(self):
        model_2_detect = ['foo/no_inter_gene_max_space']
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)

        self.assertEqual(str(context.exception),
                         "Invalid model definition ({}): inter_gene_max_space must be defined".format(
                             os.path.join(self.cfg.models_dir(), "foo/definitions/no_inter_gene_max_space.xml")
                         )
                         )


    def test_loner(self):
        def_2_parse = set()
        model_fqn = 'foo/model_5'
        def_2_parse.add(model_fqn)
        parsed = set()
        models_2_detect = self.parser.definition_to_parse(def_2_parse, parsed)
        self.parser.parse(models_2_detect)

        m = self.model_bank[model_fqn]
        flgC = m.get_gene('flgC')
        self.assertFalse(flgC.loner)
        tadZ = m.get_gene('tadZ')
        self.assertTrue(tadZ.loner)


    def test_multi_system(self):
        def_2_parse = set()
        model_fqn = 'foo/model_5'
        def_2_parse.add(model_fqn)
        parsed = set()
        models_2_detect = self.parser.definition_to_parse(def_2_parse, parsed)
        self.parser.parse(models_2_detect)

        m = self.model_bank[model_fqn]
        flgC = m.get_gene('flgC')
        self.assertFalse(flgC.multi_system)
        fliE = m.get_gene('fliE')
        self.assertTrue(fliE.multi_system)


    def test_gene_inter_gene_max_space(self):
        def_2_parse = set()
        model_fqn = 'foo/model_5'
        def_2_parse.add(model_fqn)
        parsed = set()
        models_2_detect = self.parser.definition_to_parse(def_2_parse, parsed)
        self.parser.parse(models_2_detect)

        m = self.model_bank[model_fqn]
        self.assertEqual(m.name, 'model_5')
        self.assertEqual(m.fqn, model_fqn)
        self.assertEqual(m.inter_gene_max_space, 20)
        self.assertEqual(m.min_mandatory_genes_required, 2)
        self.assertEqual(m.min_genes_required, 3)
        self.assertFalse(m.multi_loci)
        self.assertEqual(len(m.mandatory_genes), 3)
        mandatory_genes_name = [g.name for g in m.mandatory_genes]
        mandatory_genes_name.sort()
        theoric_list = ["flgB", "flgC", "fliE"]
        theoric_list.sort()
        self.assertListEqual(mandatory_genes_name, theoric_list)
        flgC = m.get_gene('flgC')
        self.assertEqual(flgC.inter_gene_max_space, 2)


    def test_inter_gene_max_space_cfg(self):
        # test inter_gene_max_space is specified from configuration
        # so this value must overload the value read from xml
        def_2_parse = set()
        model_fqn = 'foo/model_5'
        def_2_parse.add(model_fqn)
        parsed = set()

        inter_gene_max_space_cfg = [[model_fqn, '222']]
        self.args.inter_gene_max_space = inter_gene_max_space_cfg

        self.cfg = Config(MacsyDefaults(), self.args)
        self.model_bank = ModelBank()
        self.gene_bank = GeneBank()
        self.profile_factory = ProfileFactory(self.cfg)
        self.model_registry = ModelRegistry()
        models_location = scan_models_dir(self.args.models_dir)
        for ml in models_location:
            self.model_registry.add(ml)
        self.parser = DefinitionParser(self.cfg, self.model_bank, self.gene_bank,
                                       self.profile_factory, self.model_registry)

        models_2_detect = self.parser.definition_to_parse(def_2_parse, parsed)
        self.parser.parse(models_2_detect)
        m = self.model_bank[model_fqn]
        self.assertEqual(m.inter_gene_max_space, 222)


    def test_min_mandatory_genes_required_cfg(self):
        # test min_mandatory_genes_required is specified from configuration
        # so this value must overload the value read from xml
        def_2_parse = set()
        model_fqn = 'foo/model_5'
        def_2_parse.add(model_fqn)
        parsed = set()

        min_mandatory_genes_required = [[model_fqn, '3']]
        self.args.min_mandatory_genes_required = min_mandatory_genes_required

        self.cfg = Config(MacsyDefaults(), self.args)
        self.model_bank = ModelBank()
        self.gene_bank = GeneBank()
        self.profile_factory = ProfileFactory(self.cfg)
        self.model_registry = ModelRegistry()
        models_location = scan_models_dir(self.args.models_dir)
        for ml in models_location:
            self.model_registry.add(ml)
        self.parser = DefinitionParser(self.cfg, self.model_bank, self.gene_bank,
                                       self.profile_factory, self.model_registry)

        models_2_detect = self.parser.definition_to_parse(def_2_parse, parsed)
        self.parser.parse(models_2_detect)
        m = self.model_bank[model_fqn]
        self.assertEqual(m.min_mandatory_genes_required, 3)


    def test_min_genes_required_cfg(self):
        # test min_genes_required is specified from configuration
        # so this value must overload the value read from xml
        def_2_parse = set()
        model_fqn = 'foo/model_5'
        def_2_parse.add(model_fqn)
        parsed = set()

        min_genes_required = [[model_fqn, '4']]
        self.args.min_genes_required = min_genes_required

        self.cfg = Config(MacsyDefaults(), self.args)
        self.model_bank = ModelBank()
        self.gene_bank = GeneBank()
        self.profile_factory = ProfileFactory(self.cfg)
        self.model_registry = ModelRegistry()
        models_location = scan_models_dir(self.args.models_dir)
        for ml in models_location:
            self.model_registry.add(ml)
        self.parser = DefinitionParser(self.cfg, self.model_bank, self.gene_bank,
                                       self.profile_factory, self.model_registry)

        models_2_detect = self.parser.definition_to_parse(def_2_parse, parsed)
        self.parser.parse(models_2_detect)
        m = self.model_bank[model_fqn]
        self.assertEqual(m.min_genes_required, 4)


    def test_max_nb_genes_cfg(self):
        # test max_nb_genes is specified from configuration
        # so this value must overload the value read from xml
        def_2_parse = set()
        model_fqn = 'foo/model_5'
        def_2_parse.add(model_fqn)
        parsed = set()

        max_nb_genes = [[model_fqn, '4']]
        self.args.max_nb_genes = max_nb_genes

        self.cfg = Config(MacsyDefaults(), self.args)
        self.model_bank = ModelBank()
        self.gene_bank = GeneBank()
        self.profile_factory = ProfileFactory(self.cfg)
        self.model_registry = ModelRegistry()
        models_location = scan_models_dir(self.args.models_dir)
        for ml in models_location:
            self.model_registry.add(ml)
        self.parser = DefinitionParser(self.cfg, self.model_bank, self.gene_bank,
                                       self.profile_factory, self.model_registry)

        models_2_detect = self.parser.definition_to_parse(def_2_parse, parsed)
        self.parser.parse(models_2_detect)
        m = self.model_bank[model_fqn]
        self.assertEqual(m.max_nb_genes, 4)


    def test_multi_loci_cfg(self):
        # test multi_loci is specified from configuration
        # so this value must overload the value read from xml
        def_2_parse = set()
        model_fqn = 'foo/model_5'
        def_2_parse.add(model_fqn)
        parsed = set()

        self.args.multi_loci = model_fqn

        self.cfg = Config(MacsyDefaults(), self.args)
        self.model_bank = ModelBank()
        self.gene_bank = GeneBank()
        self.profile_factory = ProfileFactory(self.cfg)
        self.model_registry = ModelRegistry()
        models_location = scan_models_dir(self.args.models_dir)
        for ml in models_location:
            self.model_registry.add(ml)
        self.parser = DefinitionParser(self.cfg, self.model_bank, self.gene_bank,
                                       self.profile_factory, self.model_registry)

        models_2_detect = self.parser.definition_to_parse(def_2_parse, parsed)
        self.parser.parse(models_2_detect)
        m = self.model_bank[model_fqn]
        self.assertTrue(m.multi_loci)


    def test_bad_gene_inter_gene_max_space_2(self):
        models_2_detect = ['foo/bad_inter_gene_max_space_2']
        with self.assertRaises(SyntaxError) as ctx:
            with self.catch_log():
                self.parser.parse(models_2_detect)

        self.assertEqual(str(ctx.exception), "Invalid model definition 'bad_inter_gene_max_space_2': "
                                             "inter_gene_max_space must be an integer: 2.5")


    def test_parse_2genes_defines_in_2models(self):
        models_2_detect = ['foo/model_2', 'foo/model_6']
        with self.assertRaises(MacsypyError) as ctx:
            with self.catch_log():
                self.parser.parse(models_2_detect)
        self.assertRegex(str(ctx.exception),
                         "gene 'sctJ' define in 'foo/model_[26]' model is already defined in an another model")


    def test_parse_model_old_syntax(self):
        models_2_detect = ['old/model_4']
        with self.catch_log(log_name='macsypy') as log:
            self.parser.parse(models_2_detect)
            log_msg = log.get_value()
        self.assertEqual(log_msg, """'system' is deprecated as xml root. Migrate model_4.xml with macsydef_1to2 script.
'system_ref' is deprecated. Migrate model_4.xml with macsydef_1to2 script.
'system' is deprecated as xml root. Migrate model_3.xml with macsydef_1to2 script.
'system_ref' is deprecated. Migrate model_3.xml with macsydef_1to2 script.
""")