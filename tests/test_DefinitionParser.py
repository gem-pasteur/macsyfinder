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

from macsypy.config import Config, MacsyDefaults
from macsypy.model import ModelBank
from macsypy.profile import ProfileFactory
from macsypy.gene import GeneBank, CoreGene, ModelGene, Exchangeable
from macsypy.registries import ModelRegistry, scan_models_dir
from macsypy.definition_parser import DefinitionParser
from macsypy.error import MacsypyError, ModelInconsistencyError
from tests import MacsyTest


class TestModelParser(MacsyTest):

    def setUp(self):
        defaults = MacsyDefaults()
        self.args = argparse.Namespace()
        self.args.sequence_db = self.find_data("base", "test_1.fasta")
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
                                       self.model_registry, self.profile_factory)
        
        
    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir())
        except:
            pass


    def test_parse_with_exchangeable(self):
        model_name = 'model_1'
        model_family = 'foo'
        fqn = f"{model_family}/{model_name}"
        #def_2_parse = set()
        #def_2_parse.add(fqn)
        models_2_detect = [self.model_registry['foo'].get_definition(fqn)]
        self.parser.parse(models_2_detect)
        self.assertEqual(len(self.model_bank), 1)

        m1 = self.model_bank[fqn]
        self.assertEqual(m1.name, model_name)
        self.assertEqual(m1.fqn, fqn)
        self.assertEqual(m1.inter_gene_max_space, 20)
        self.assertEqual(m1.min_mandatory_genes_required, 2)
        self.assertEqual(m1.min_genes_required, 4)
        self.assertTrue(m1.multi_loci)

        self.assertEqual(len(m1.mandatory_genes), 2)
        mandatory_genes_name = sorted([g.name for g in m1.mandatory_genes])
        theoric_list = sorted(["sctJ_FLG", "sctN_FLG"])
        self.assertListEqual(mandatory_genes_name, theoric_list)

        self.assertEqual(len(m1.accessory_genes), 2)
        accessory_genes_name = sorted([g.name for g in m1.accessory_genes])
        theoric_list = sorted(["flgB", "flgC"])
        self.assertListEqual(accessory_genes_name, theoric_list)

        self.assertEqual(len(m1.neutral_genes), 2)
        neutral_genes_name = sorted([g.name for g in m1.neutral_genes])
        theoric_list = sorted(["fliE", "tadZ"])
        self.assertListEqual(neutral_genes_name, theoric_list)

        self.assertEqual(len(m1.forbidden_genes), 1)
        forbidden_genes_name = sorted([g.name for g in m1.forbidden_genes])
        theoric_list = sorted(["sctC"])
        self.assertListEqual(forbidden_genes_name, theoric_list)

        sctJ_FLG = m1.get_gene('sctJ_FLG')
        sctJ_FLG_exchangeables = sctJ_FLG.exchangeables
        self.assertEqual(len(sctJ_FLG_exchangeables), 2)
        self.assertEqual(sctJ_FLG_exchangeables[0].name, 'sctJ')
        self.assertEqual(sctJ_FLG_exchangeables[1].name, 'abc')
        self.assertTrue(isinstance(sctJ_FLG_exchangeables[0], Exchangeable))
        self.assertTrue(isinstance(sctJ_FLG_exchangeables[0]._gene, CoreGene))
        self.assertTrue(isinstance(sctJ_FLG_exchangeables[0].alternate_of(), ModelGene))
        self.assertTrue(sctJ_FLG_exchangeables[0].loner)
        self.assertFalse(sctJ_FLG.is_exchangeable)
        sctJ = m1.get_gene('sctJ')
        self.assertTrue(sctJ.is_exchangeable)


    def test_exchangeable_inheritance(self):
        def_2_parse = set()
        def_2_parse.add('foo/model_1')
        models_2_detect = [self.model_registry['foo'].get_definition('foo/model_1')]
        self.parser.parse(models_2_detect)
        m1 = self.model_bank['foo/model_1']

        sctJ = m1.get_gene('sctJ')
        self.assertTrue(sctJ.is_exchangeable)
        self.assertTrue(sctJ.loner)
        self.assertTrue(sctJ.multi_system)
        self.assertFalse(sctJ.multi_model)
        sctJ_FLG = m1.get_gene('sctJ_FLG')
        self.assertTrue(sctJ_FLG.multi_system)
        abc = m1.get_gene('abc')
        self.assertFalse(abc.multi_system)

        sctN = m1.get_gene('sctN')
        sctN_FLG = m1.get_gene('sctN_FLG')

        self.assertFalse(sctN_FLG.loner)
        self.assertTrue(sctN.loner)
        self.assertIsNone(sctN_FLG.inter_gene_max_space)
        self.assertEqual(sctN.inter_gene_max_space, 10)
        self.assertFalse(sctN_FLG.multi_model)
        self.assertFalse(sctN.multi_model)
        gspD = m1.get_gene('gspD')
        self.assertFalse(sctN_FLG.multi_system)
        self.assertTrue(gspD.multi_model)
        self.assertTrue(gspD.multi_system)


    def test_model_w_unkown_attr(self):
        model_2_detect = [self.model_registry['foo'].get_definition('foo/model_w_unknown_attribute')]
        with self.assertRaises(MacsypyError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "unable to parse model definition 'foo/model_w_unknown_attribute' : "
                         "The model definition model_w_unknown_attribute.xml has an unknow attribute 'multi-loci'. "
                         "Please fix the definition.")

    def test_gene_w_unkown_attr(self):
        model_2_detect = [self.model_registry['foo'].get_definition('foo/gene_w_unknown_attribute')]
        with self.assertRaises(MacsypyError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "unable to parse model definition 'foo/gene_w_unknown_attribute' : "
                         "The model definition gene_w_unknown_attribute.xml has an unknown attribute 'multi-system' for a gene."
                         " Please fix the definition.")

    def test_wo_presence(self):
        model_2_detect = [self.model_registry['foo'].get_definition('foo/fail_wo_presence')]
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid model definition 'foo/fail_wo_presence': gene 'sctN_FLG' without presence")


    def test_invalid_presence(self):
        model_2_detect = [self.model_registry['foo'].get_definition('foo/fail_invalid_presence')]
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid model 'foo/fail_invalid_presence' definition: presence value must be either: "
                         "'mandatory', 'accessory', 'neutral', 'forbidden' not foo_bar")

    def test_gene_no_name(self):
        model_2_detect = [self.model_registry['foo'].get_definition('foo/gene_no_name')]
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid model definition 'foo/gene_no_name': gene without name")

    def test_invalid_homolog(self):
        model_2_detect = [self.model_registry['foo'].get_definition('foo/invalid_homolog')]
        with self.assertRaises(MacsypyError) as context:
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "'foo/foo_bar': No such profile")

    def test_invalid_homolog_2(self):
        model_2_detect = [self.model_registry['foo'].get_definition('foo/invalid_homolog_2')]
        with self.assertRaises(SyntaxError) as ctx:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(ctx.exception), "Invalid model definition 'foo/invalid_homolog_2': gene without name")


    def test_bad_min_genes_required(self):
        model_2_detect = [self.model_registry['foo'].get_definition('foo/bad_min_genes_required')]
        with self.assertRaises(ModelInconsistencyError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         'model \'bad_min_genes_required\' is not consistent: min_genes_required 16 must be lesser '
                         'or equal than the number of "accessory" and "mandatory" components in the model: 6')

    def test_bad_min_genes_required_2(self):
        model_2_detect = [self.model_registry['foo'].get_definition('foo/bad_min_genes_required_2')]
        with self.catch_log():
            with self.assertRaisesRegex(SyntaxError, "Invalid model definition (.*): "
                                                     "min_genes_required must be an integer: 16.5"):
                self.parser.parse(model_2_detect)

    def test_bad_min_mandatory_genes_required(self):
        model_2_detect = [self.model_registry['foo'].get_definition('foo/bad_min_mandatory_genes_required')]
        with self.catch_log():
            with self.assertRaises(ModelInconsistencyError) as context:
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         'model \'bad_min_mandatory_genes_required\' is not consistent: min_genes_required 16 must '
                         'be lesser or equal than the number of "accessory" and "mandatory" components in the model: 6')

    def test_bad_min_mandatory_genes_required_2(self):
        model_2_detect = [self.model_registry['foo'].get_definition('foo/bad_min_mandatory_genes_required_2')]
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
        model_2_detect = [self.model_registry['foo'].get_definition('foo/bad_min_mandatory_genes_required_4')]
        with self.assertRaisesRegex(SyntaxError, "Invalid model definition (.*): "
                                                  "min_mandatory_genes_required must be an integer: 12.5"):
            with self.catch_log():
                self.parser.parse(model_2_detect)


    def test_min_mandatory_genes_required_lesser_than_mandatory_genes(self):
        model_2_detect = [self.model_registry['foo'].get_definition('foo/bad_min_mandatory_genes_required_3')]
        with self.assertRaises(ModelInconsistencyError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "model 'bad_min_mandatory_genes_required_3' is not consistent:"
                         " 'min_mandatory_genes_required': 6 must be lesser or equal than the number of 'mandatory' "
                         "components in the model: 5")

    def test_only_one_accessory(self):
        model_2_detect = [self.model_registry['foo'].get_definition('foo/only_one_accessory')]
        with self.assertRaises(ModelInconsistencyError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         f"model 'only_one_accessory' is not consistent: there is only one gene in your model. " \
                         f"So its status should be 'mandatory'.")


    def test_bad_max_nb_genes(self):
        model_2_detect = [self.model_registry['foo'].get_definition('foo/bad_max_nb_genes')]
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        model_name, def_name = model_2_detect[0].split_fqn(model_2_detect[0].fqn)
        self.assertEqual(str(context.exception),
                         "Invalid model definition ({0}.xml): max_nb_genes must be an integer: HOHOHO".format(
                             os.path.join(self.cfg.models_dir()[0],
                                          model_name,
                                          'definitions',
                                          def_name)))


    def test_bad_inter_gene_max_space(self):
        fqn = 'foo/bad_inter_gene_max_space'
        model_family, model_name = fqn.split('/')
        model_2_detect = [self.model_registry['foo'].get_definition(fqn)]
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid model definition ({}): inter_gene_max_space must be an integer: 12.5".format(
                             os.path.join(self.cfg.models_dir()[0],
                                          model_family, 'definitions', model_name + ".xml")
                         )
                         )

    def test_no_inter_gene_max_space(self):
        model_2_detect = [self.model_registry['foo'].get_definition('foo/no_inter_gene_max_space')]
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)

        self.assertEqual(str(context.exception),
                         "Invalid model definition ({}): inter_gene_max_space must be defined".format(
                             os.path.join(self.cfg.models_dir()[0], "foo/definitions/no_inter_gene_max_space.xml")
                         )
                         )


    def test_loner(self):
        model_fqn = 'foo/model_5'
        model_2_detect = [self.model_registry['foo'].get_definition(model_fqn)]
        self.parser.parse(model_2_detect)

        m5 = self.model_bank[model_fqn]
        m5_flgC = m5.get_gene('flgC')
        self.assertFalse(m5_flgC.loner)
        m5_tadZ = m5.get_gene('tadZ')
        self.assertTrue(m5_tadZ.loner)

        model_fqn = 'foo/model_6'
        model_2_detect = [self.model_registry['foo'].get_definition(model_fqn)]
        self.parser.parse(model_2_detect)
        m6 = self.model_bank[model_fqn]
        m6_flgC = m6.get_gene('flgC')
        self.assertFalse(m6_flgC.loner)
        m6_tadZ = m6.get_gene('tadZ')
        self.assertFalse(m6_tadZ.loner)


    def test_multi_system(self):
        model_fqn = 'foo/model_5'
        model_2_detect = [self.model_registry['foo'].get_definition(model_fqn)]
        self.parser.parse(model_2_detect)

        m = self.model_bank[model_fqn]
        flgC = m.get_gene('flgC')
        self.assertFalse(flgC.multi_system)
        fliE = m.get_gene('fliE')
        self.assertTrue(fliE.multi_system)


    def test_multi_model(self):
        model_fqn = 'foo/model_5'
        model_2_detect = [self.model_registry['foo'].get_definition(model_fqn)]
        self.parser.parse(model_2_detect)

        m = self.model_bank[model_fqn]
        flgC = m.get_gene('flgC')
        self.assertFalse(flgC.multi_model)
        abc = m.get_gene('abc')
        self.assertTrue(abc.multi_model)


    def test_gene_inter_gene_max_space(self):
        model_fqn = ['foo/model_5', 'foo/model_6']
        models_2_detect = [self.model_registry['foo'].get_definition(fqn) for fqn in model_fqn]
        self.parser.parse(models_2_detect)

        m5 = self.model_bank['foo/model_5']
        self.assertEqual(m5.name, 'model_5')
        self.assertEqual(m5.fqn, 'foo/model_5')
        self.assertEqual(m5.inter_gene_max_space, 20)
        m5_flgB = m5.get_gene('flgB')
        m5_flgC = m5.get_gene('flgC')
        self.assertIsNone(m5_flgB.inter_gene_max_space)
        self.assertEqual(m5_flgC.inter_gene_max_space, 2)
        m6 = self.model_bank['foo/model_6']
        m6_flgC = m6.get_gene('flgC')
        self.assertEqual(m6_flgC.inter_gene_max_space, 12)


    def test_inter_gene_max_space_cfg(self):
        # test inter_gene_max_space is specified from configuration
        # so this value must overload the value read from xml
        model_fqn = 'foo/model_5'

        inter_gene_max_space_cfg = [[model_fqn, '222']]
        self.args.inter_gene_max_space = inter_gene_max_space_cfg

        self.cfg = Config(MacsyDefaults(), self.args)
        self.model_bank = ModelBank()
        self.gene_bank = GeneBank()
        self.model_registry = ModelRegistry()
        models_location = scan_models_dir(self.args.models_dir)
        for ml in models_location:
            self.model_registry.add(ml)
        self.parser = DefinitionParser(self.cfg, self.model_bank, self.gene_bank,
                                       self.model_registry, self.profile_factory)

        models_2_detect = [self.model_registry['foo'].get_definition(model_fqn)]
        self.parser.parse(models_2_detect)
        m = self.model_bank[model_fqn]
        self.assertEqual(m.inter_gene_max_space, 222)


    def test_min_mandatory_genes_required_cfg(self):
        # test min_mandatory_genes_required is specified from configuration
        # so this value must overload the value read from xml
        model_fqn = 'foo/model_5'

        min_mandatory_genes_required = [[model_fqn, '3']]
        self.args.min_mandatory_genes_required = min_mandatory_genes_required

        self.cfg = Config(MacsyDefaults(), self.args)
        self.model_bank = ModelBank()
        self.gene_bank = GeneBank()
        self.model_registry = ModelRegistry()
        models_location = scan_models_dir(self.args.models_dir)
        for ml in models_location:
            self.model_registry.add(ml)
        self.parser = DefinitionParser(self.cfg, self.model_bank, self.gene_bank,
                                       self.model_registry, self.profile_factory)

        models_2_detect = [self.model_registry['foo'].get_definition(model_fqn)]
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
        self.model_registry = ModelRegistry()
        models_location = scan_models_dir(self.args.models_dir)
        for ml in models_location:
            self.model_registry.add(ml)
        self.parser = DefinitionParser(self.cfg, self.model_bank, self.gene_bank,
                                       self.model_registry, self.profile_factory)

        models_2_detect = [self.model_registry['foo'].get_definition(model_fqn)]
        self.parser.parse(models_2_detect)
        m = self.model_bank[model_fqn]
        self.assertEqual(m.min_genes_required, 4)


    def test_max_nb_genes_cfg(self):
        self.model_bank = ModelBank()
        self.gene_bank = GeneBank()
        self.model_registry = ModelRegistry()
        models_location = scan_models_dir(self.args.models_dir)
        for ml in models_location:
            self.model_registry.add(ml)

        # max_nb_genes is specified in xml
        # no user configuration on this
        self.cfg = Config(MacsyDefaults(), self.args)
        model_fqn = 'foo/model_6'  # 4 genes in this model but xml specify 3
        self.cfg = Config(MacsyDefaults(), self.args)
        self.parser = DefinitionParser(self.cfg, self.model_bank, self.gene_bank,
                                       self.model_registry, self.profile_factory)


        models_2_detect = [self.model_registry['foo'].get_definition(model_fqn)]
        self.parser.parse(models_2_detect)
        m = self.model_bank[model_fqn]
        self.assertEqual(m.max_nb_genes, 3)

        # max_nb_genes is specified from configuration
        # so this value must overload the value read from xml
        model_fqn = 'foo/model_5' # 4 genes in this model
        max_nb_genes = [[model_fqn, '6']]
        self.args.max_nb_genes = max_nb_genes
        self.cfg = Config(MacsyDefaults(), self.args)
        self.parser = DefinitionParser(self.cfg, self.model_bank, self.gene_bank,
                                       self.model_registry, self.profile_factory)

        models_2_detect = [self.model_registry['foo'].get_definition(model_fqn)]
        self.parser.parse(models_2_detect)
        m = self.model_bank[model_fqn]
        self.assertEqual(m.max_nb_genes, 6)


    def test_multi_loci_cfg(self):
        # test multi_loci is specified from configuration
        # so this value must overload the value read from xml
        model_fqn = 'foo/model_5'

        self.args.multi_loci = model_fqn

        self.cfg = Config(MacsyDefaults(), self.args)
        self.model_bank = ModelBank()
        self.gene_bank = GeneBank()
        self.model_registry = ModelRegistry()
        models_location = scan_models_dir(self.args.models_dir)
        for ml in models_location:
            self.model_registry.add(ml)
        self.parser = DefinitionParser(self.cfg, self.model_bank, self.gene_bank,
                                       self.model_registry, self.profile_factory)

        models_2_detect = [self.model_registry['foo'].get_definition(model_fqn)]
        self.parser.parse(models_2_detect)
        m = self.model_bank[model_fqn]
        self.assertTrue(m.multi_loci)


    def test_bad_gene_inter_gene_max_space_2(self):
        model_fqn = 'foo/bad_inter_gene_max_space_2'
        models_2_detect = [self.model_registry['foo'].get_definition(model_fqn)]
        with self.assertRaises(SyntaxError) as ctx:
            with self.catch_log():
                self.parser.parse(models_2_detect)

        self.assertEqual(str(ctx.exception), "Invalid model definition 'foo/bad_inter_gene_max_space_2': "
                                             "inter_gene_max_space must be an integer: 2.5")


    def test_bad_exchangeable_inter_gene_max_space(self):
        fqn = 'foo/bad_exchangeable_inter_gene_max_space'
        model_2_detect = [self.model_registry['foo'].get_definition(fqn)]
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(model_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid model definition 'foo/bad_exchangeable_inter_gene_max_space': "
                         "inter_gene_max_space must be an integer: 1.5")


    def test_parse_model_old_syntax(self):
        # the attribute vers is not set
        model_fqn = 'foo/model_old_1'
        models_2_detect = [self.model_registry['foo'].get_definition(model_fqn)]
        with self.catch_log(log_name='macsypy') as log:
            with self.assertRaises(MacsypyError) as ctx:
                self.parser.parse(models_2_detect)
            log_msg = log.get_value().strip()
        self.assertEqual(log_msg,
                         "unable to parse model definition 'foo/model_old_1' : "
                         "The model definition model_old_1.xml is not versioned. Please update your model.")

        # the root is system instead of mmodel
        model_fqn = 'foo/model_old_2'
        models_2_detect = [self.model_registry['foo'].get_definition(model_fqn)]
        with self.catch_log(log_name='macsypy') as log:
            with self.assertRaises(MacsypyError) as ctx:
                self.parser.parse(models_2_detect)
            log_msg = log.get_value().strip()
        self.assertEqual(log_msg,
                         f"unable to parse model definition '{model_fqn}' : "
                         "The model definition model_old_2.xml is obsolete. Please update your model.")

        # there still system_ref attribute
        model_fqn = 'foo/model_old_3'
        models_2_detect = [self.model_registry['foo'].get_definition(model_fqn)]
        with self.catch_log(log_name='macsypy') as log:
            with self.assertRaises(MacsypyError) as ctx:
                self.parser.parse(models_2_detect)
            log_msg = log.get_value().strip()
        self.assertEqual(log_msg,
                         f"unable to parse model definition '{model_fqn}' : "
                         "The model definition model_old_3.xml is obsolete. Please update your model.")

        # there still homologs tag
        model_fqn = 'foo/model_old_4'
        models_2_detect = [self.model_registry['foo'].get_definition(model_fqn)]
        with self.catch_log(log_name='macsypy') as log:
            with self.assertRaises(MacsypyError) as ctx:
                self.parser.parse(models_2_detect)
            log_msg = log.get_value().strip()
        self.assertEqual(log_msg,
                         f"unable to parse model definition '{model_fqn}' : "
                         "The model definition model_old_4.xml is obsolete. Please update your model.")

        # there still analogs tag
        model_fqn = 'foo/model_old_5'
        models_2_detect = [self.model_registry['foo'].get_definition(model_fqn)]
        with self.catch_log(log_name='macsypy') as log:
            with self.assertRaises(MacsypyError) as ctx:
                self.parser.parse(models_2_detect)
            log_msg = log.get_value().strip()
        self.assertEqual(log_msg,
                         f"unable to parse model definition '{model_fqn}' : "
                         "The model definition model_old_5.xml is obsolete. Please update your model.")
