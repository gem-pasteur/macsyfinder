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
import platform
import logging
import argparse

from macsypy.config import Config, MacsyDefaults
from macsypy.system import SystemBank
from macsypy.gene import GeneBank
from macsypy.system_parser import SystemParser
from macsypy.error import MacsypyError, SystemInconsistencyError
from tests import MacsyTest


class TestSystemParser(MacsyTest):

    def setUp(self):
        defaults = MacsyDefaults()
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.res_search_dir = tempfile.gettempdir()
        args.log_level = 30

        self.cfg = Config(defaults, args)
        self.system_bank = SystemBank()
        self.system_bank._system_bank = {}
        self.gene_bank = GeneBank()
        self.gene_bank._genes_bank = {}
        self.parser = SystemParser(self.cfg, self.system_bank, self.gene_bank)
        
        
    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir())
        except:
            pass

    def test_defintion_to_parse(self):
        parsed = set()
        models_2_detect = set()
        models_2_detect.add('foo/system_1')
        def_2_parse = self.parser.definition_to_parse(models_2_detect, parsed)
        self.assertSetEqual(def_2_parse, {s for s in ('foo/system_1', 'foo/system_2')})
        parsed = set()
        models_2_detect = set()
        definition_name = 'foo/nimportnaoik'
        models_2_detect.add(definition_name)
        with self.assertRaises(MacsypyError) as context:
            with self.catch_log():
                self.parser.definition_to_parse(models_2_detect, parsed)
        self.assertEqual(str(context.exception), '{}: No such definition'.format(definition_name))

        parsed = set()
        models_2_detect = set()
        model_name = 'bar'
        definition = '{}/nimportnaoik'.format(model_name)
        models_2_detect.add(definition)
        with self.assertRaises(MacsypyError) as context:
            with self.catch_log():
                self.parser.definition_to_parse(models_2_detect, parsed)
        self.assertEqual(str(context.exception), '{}: No such Models in {}'.format(model_name, self.cfg.models_dir()))

        parsed = set()
        models_2_detect = set()
        model_name = 'foo'
        def_name = 'not_xml'
        fqn = '{}/{}'.format(model_name, def_name)
        models_2_detect.add(fqn)
        with self.assertRaises(MacsypyError) as context:
            with self.catch_log():
                self.parser.definition_to_parse(models_2_detect, parsed)
        self.assertTrue(str(context.exception).startswith('unable to parse system definition "{}" :'.format(fqn)))

    def test_parse_with_homologs(self):
        def_2_parse = set()
        def_2_parse.add('foo/system_1')
        parsed = set()
        models_2_detect = self.parser.definition_to_parse(def_2_parse, parsed)
        self.parser.parse(models_2_detect)
        self.assertEqual(len(self.system_bank), 2)

        s2 = self.system_bank['foo/system_2']
        self.assertEqual(s2.name, 'system_2')
        self.assertEqual(s2.fqn, 'foo/system_2')

        s1 = self.system_bank['foo/system_1']
        self.assertEqual(s1.name, 'system_1')
        self.assertEqual(s1.fqn, 'foo/system_1')
        self.assertEqual(s1.inter_gene_max_space, 20)
        self.assertEqual(s1.min_mandatory_genes_required, 4)
        self.assertEqual(s1.min_genes_required, 6)
        self.assertTrue(s1.multi_loci)
        self.assertFalse(s2.multi_loci)
        self.assertIsNone(s1.max_nb_genes)
        self.assertEqual(s2.max_nb_genes, 1)
        self.assertEqual(len(s1.mandatory_genes), 5)
        mandatory_genes_name = [g.name for g in s1.mandatory_genes]
        mandatory_genes_name.sort()
        theoric_list = ["sctJ_FLG", "sctN_FLG", "flgB", "flgC", "fliE"]
        theoric_list.sort()
        self.assertListEqual(mandatory_genes_name, theoric_list)
        sctJ_FLG = [g for g in s1.mandatory_genes if g.name == 'sctJ_FLG'][0]
        sctJ_FLG_homologs = sctJ_FLG.get_homologs()
        self.assertEqual(len(sctJ_FLG_homologs), 1)
        self.assertEqual(sctJ_FLG_homologs[0].name, 'sctJ')
        self.assertEqual(sctJ_FLG_homologs[0].system, self.system_bank['foo/system_2'])
        self.assertEqual(len(s1.accessory_genes), 1)
        self.assertEqual(s1.accessory_genes[0].name, 'tadZ')
        self.assertEqual(s1.accessory_genes[0].system, self.system_bank['foo/system_2'])
        self.assertEqual(len(s1.forbidden_genes), 1)
        self.assertEqual(s1.forbidden_genes[0].name, 'sctC')
        self.assertEqual(s1.forbidden_genes[0].system, self.system_bank['foo/system_2'])


    def test_parse_with_analogs(self):
        def_2_parse = set()
        def_2_parse.add('foo/system_3')
        parsed = set()
        models_2_detect = self.parser.definition_to_parse(def_2_parse, parsed)
        self.parser.parse(models_2_detect)

        self.assertEqual(len(self.system_bank), 2)

        s2 = self.system_bank['foo/system_4']
        self.assertEqual(s2.name, 'system_4')
        self.assertEqual(s2.fqn, 'foo/system_4')

        s1 = self.system_bank['foo/system_3']
        self.assertEqual(s1.name, 'system_3')
        self.assertEqual(s1.fqn, 'foo/system_3')
        self.assertEqual(s1.inter_gene_max_space, 20)
        self.assertEqual(s1.min_mandatory_genes_required, 4)
        self.assertEqual(s1.min_genes_required, 6)
        self.assertTrue(s1.multi_loci)
        self.assertFalse(s2.multi_loci)
        self.assertIsNone(s1.max_nb_genes)
        self.assertEqual(s2.max_nb_genes, 1)
        self.assertEqual(len(s1.mandatory_genes), 5)
        mandatory_genes_name = [g.name for g in s1.mandatory_genes]
        mandatory_genes_name.sort()
        theoric_list = ["sctJ_FLG", "sctN_FLG", "flgB", "flgC", "fliE"]
        theoric_list.sort()
        self.assertListEqual(mandatory_genes_name, theoric_list)
        sctJ_FLG = [g for g in s1.mandatory_genes if g.name == 'sctJ_FLG'][0]
        sctJ_FLG_homologs = sctJ_FLG.get_homologs()
        self.assertListEqual(sctJ_FLG_homologs, [])
        sctJ_FLG_analogs = sctJ_FLG.get_analogs()
        self.assertEqual(len(sctJ_FLG_analogs), 1)
        self.assertEqual(sctJ_FLG_analogs[0].name, 'sctJ')
        self.assertEqual(sctJ_FLG_analogs[0].system, self.system_bank['foo/system_4'])
        self.assertEqual(len(s1.accessory_genes), 1)
        self.assertEqual(s1.accessory_genes[0].name, 'tadZ')
        self.assertEqual(s1.accessory_genes[0].system, self.system_bank['foo/system_4'])
        self.assertEqual(len(s1.forbidden_genes), 1)
        self.assertEqual(s1.forbidden_genes[0].name, 'sctC')
        self.assertEqual(s1.forbidden_genes[0].system, self.system_bank['foo/system_4'])


    def test_wo_presence(self):
        system_2_detect = ['foo/fail_wo_presence']
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(system_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid system definition 'fail_wo_presence': gene without presence")


    def test_invalid_presence(self):
        system_2_detect = ['foo/fail_invalid_presence']
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(system_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid system 'fail_invalid_presence' definition: presence value must be either "
                         "[mandatory, accessory, forbidden] not foo_bar")

    def test_gene_no_name(self):
        system_2_detect = ['foo/gene_no_name']
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(system_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid system definition 'gene_no_name': gene without a name")

    def test_invalid_aligned(self):
        system_2_detect = ['foo/invalid_aligned']
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(system_2_detect)
        self.assertEqual(str(context.exception),
                         'Invalid system definition \'invalid_aligned\': invalid value for an attribute of gene '
                         '\'foo_bar\': \'totote\' allowed values are "1", "true", "True", "0", "false", "False"')

    def test_invalid_homolog(self):
        system_2_detect = ['foo/invalid_homolog']
        with self.assertRaises(SystemInconsistencyError) as context:
            with self.catch_log():
                self.parser.parse(system_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid system definition 'invalid_homolog': The gene 'foo_bar' described as "
                         "homolog of 'gspD' in system 'invalid_homolog' is not in the 'GeneBank' gene factory")

    def test_invalid_homolog_2(self):
        system_2_detect = ['foo/invalid_homolog_2']
        with self.assertRaises(SyntaxError) as ctx:
            with self.catch_log():
                self.parser.parse(system_2_detect)
        self.assertEqual(str(ctx.exception), "Invalid system definition: gene without name")

    def test_invalid_analog(self):
        system_2_detect = ['foo/invalid_analog']
        with self.assertRaises(SystemInconsistencyError) as context:
            with self.catch_log():
                self.parser.parse(system_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid system definition 'invalid_analog': The gene 'foo_bar' described as "
                         "analog of 'gspD' in system 'invalid_analog' is not in the 'GeneBank' gene factory")

    def test_invalid_analog_2(self):
        system_2_detect = ['foo/invalid_analog_2']
        with self.assertRaises(SyntaxError) as ctx:
            with self.catch_log():
                self.parser.parse(system_2_detect)
        self.assertEqual(str(ctx.exception), "Invalid system definition 'invalid_analog_2': gene without name")

    def test_bad_homolog_sys_ref(self):
        system_2_detect = ['foo/bad_homolog_sys_ref']
        with self.assertRaises(SystemInconsistencyError) as context:
            with self.catch_log():
                self.parser.parse(system_2_detect)
        self.assertEqual(str(context.exception),
                         "Inconsistency in systems definitions: the gene 'sctJ' described as homolog of 'sctN' "
                         "with system_ref 'system_1' has an other system in bank (system_2)")

    def test_bad_analog_sys_ref(self):
        system_2_detect = ['foo/bad_analog_sys_ref']
        with self.assertRaises(SystemInconsistencyError) as context:
            with self.catch_log():
                self.parser.parse(system_2_detect)
        self.assertEqual(str(context.exception),
                         "Inconsistency in systems definitions: the gene 'sctJ' described as analog of 'sctN' "
                         "with system_ref 'system_3' has an other system in bank (system_4)")

    def test_bad_min_genes_required(self):
        system_2_detect = ['foo/bad_min_genes_required']
        with self.assertRaises(SystemInconsistencyError) as context:
            with self.catch_log():
                self.parser.parse(system_2_detect)
        self.assertEqual(str(context.exception),
                         'system \'bad_min_genes_required\' is not consistent: min_genes_required 16 must be lesser '
                         'or equal than the number of "accessory" and "mandatory" components in the system: 6')

    def test_bad_min_genes_required_2(self):
        system_2_detect = ['foo/bad_min_genes_required_2']
        with self.catch_log():
            with self.assertRaisesRegex(SyntaxError, "Invalid system definition (.*): "
                                                     "min_genes_required must be an integer: 16.5"):
                self.parser.parse(system_2_detect)

    def test_bad_min_mandatory_genes_required(self):
        system_2_detect = ['foo/bad_min_mandatory_genes_required']
        with self.catch_log():
            with self.assertRaises(SystemInconsistencyError) as context:
                self.parser.parse(system_2_detect)
        self.assertEqual(str(context.exception),
                         'system \'bad_min_mandatory_genes_required\' is not consistent: min_genes_required 16 must '
                         'be lesser or equal than the number of "accessory" and "mandatory" components in the system: 6')

    def test_bad_min_mandatory_genes_required_2(self):
        system_2_detect = ['foo/bad_min_mandatory_genes_required_2']
        with self.assertRaises(SystemInconsistencyError) as context:
            with self.catch_log():
                # error raised by System initialization
                # which occur before check_consistency
                # the last test : not(system.min_mandatory_genes_required <= system.min_genes_required)
                # seems to be useless
                self.parser.parse(system_2_detect)
        self.assertEqual(str(context.exception),
                         "min_genes_required must be greater or equal than min_mandatory_genes_required")

    def test_bad_min_mandatory_genes_required_4(self):
        system_2_detect = ['foo/bad_min_mandatory_genes_required_4']
        with self.assertRaisesRegex(SyntaxError, "Invalid system definition (.*): "
                                                  "min_mandatory_genes_required must be an integer: 12.5"):
            with self.catch_log():
                self.parser.parse(system_2_detect)


    def test_min_mandatory_genes_required_lesser_than_mandatory_genes(self):
        system_2_detect = ['foo/bad_min_mandatory_genes_required_3']
        with self.assertRaises(SystemInconsistencyError) as context:
            with self.catch_log():
                self.parser.parse(system_2_detect)
        self.assertEqual(str(context.exception),
                         "system 'bad_min_mandatory_genes_required_3' is not consistent:"
                         " \"min_mandatory_genes_required\": 6 must be lesser or equal than the number of \"mandatory\" "
                         "components in the system: 5")


    def test_bad_max_nb_genes(self):
        system_2_detect = 'foo/bad_max_nb_genes'
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse([system_2_detect])
        model_name, def_name = system_2_detect.split('/')
        self.assertEqual(str(context.exception),
                         "Invalid system definition ({0}.xml): max_nb_genes must be an integer: HOHOHO".format(
                             os.path.join(self.cfg.models_dir(),
                                          model_name,
                                          'definitions',
                                          def_name)))


    def test_bad_inter_gene_max_space(self):
        system_2_detect = ['foo/bad_inter_gene_max_space']
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(system_2_detect)
        self.assertEqual(str(context.exception),
                         "Invalid system definition ({}): inter_gene_max_space must be an integer: 12.5".format(
                             os.path.join(self.cfg.models_dir(), "foo/definitions/bad_inter_gene_max_space.xml")
                         )
                         )

    def test_no_inter_gene_max_space(self):
        system_2_detect = ['foo/no_inter_gene_max_space']
        with self.assertRaises(SyntaxError) as context:
            with self.catch_log():
                self.parser.parse(system_2_detect)

        self.assertEqual(str(context.exception),
                         "Invalid system definition ({}): inter_gene_max_space must be defined".format(
                             os.path.join(self.cfg.models_dir(), "foo/definitions/no_inter_gene_max_space.xml")
                         )
                         )

    def test_gene_inter_gene_max_space(self):
        def_2_parse = set()
        def_2_parse.add('foo/system_5')
        parsed = set()
        models_2_detect = self.parser.definition_to_parse(def_2_parse, parsed)
        self.parser.parse(models_2_detect)

        s1 = self.system_bank['foo/system_5']
        self.assertEqual(s1.name, 'system_5')
        self.assertEqual(s1.fqn, 'foo/system_5')
        self.assertEqual(s1.inter_gene_max_space, 20)
        self.assertEqual(s1.min_mandatory_genes_required, 2)
        self.assertEqual(s1.min_genes_required, 3)
        self.assertFalse(s1.multi_loci)
        self.assertEqual(len(s1.mandatory_genes), 3)
        mandatory_genes_name = [g.name for g in s1.mandatory_genes]
        mandatory_genes_name.sort()
        theoric_list = ["flgB", "flgC", "fliE"]
        theoric_list.sort()
        self.assertListEqual(mandatory_genes_name, theoric_list)
        flgC = s1.get_gene('flgC')
        self.assertEqual(flgC.inter_gene_max_space, 2)

    def test_bad_gene_inter_gene_max_space_2(self):
        models_2_detect = ['foo/bad_inter_gene_max_space_2']
        with self.assertRaises(SyntaxError) as ctx:
            with self.catch_log():
                self.parser.parse(models_2_detect)

        self.assertEqual(str(ctx.exception), "Invalid system definition 'bad_inter_gene_max_space_2': "
                                             "inter_gene_max_space must be an integer: 2.5")
