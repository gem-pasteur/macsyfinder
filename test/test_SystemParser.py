# -*- coding: utf-8 -*-

#===============================================================================
# Created on Jan 14, 2013
# 
# @author: bneron
# @contact: user_email
# @organization: organization_name
# @license: license
#===============================================================================

import os
import unittest
import shutil
from txsscanlib.config import Config
from txsscanlib.system import SystemBank
from txsscanlib.gene import GeneBank
from txsscanlib.system_parser import SystemParser
from txsscanlib.txsscan_error import TxsscanError, SystemInconsistencyError


class Test(unittest.TestCase):

    _data_dir = os.path.join(os.path.dirname(__file__), "datatest")
    
    def setUp(self):
        self.cfg = Config(sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type = "gembase", 
                          hmmer_exe = "",
                           e_value_res = 1,
                           i_evalue_sel = 0.5,
                           def_dir = os.path.join(self._data_dir, 'DEF'),
                           res_search_dir = "/tmp",
                           res_search_suffix = "",
                           profile_dir = os.path.join(self._data_dir, 'profiles'),
                           profile_suffix = ".fasta-aln_edit.hmm",
                           res_extract_suffix = "",
                           log_level = 30,
                           log_file = '/dev/null'
                           )
        self.system_bank = SystemBank()
        self.system_bank._system_bank = {}
        self.gene_bank = GeneBank()
        self.gene_bank._genes_bank = {}
        self.parser = SystemParser(self.cfg, self.system_bank, self.gene_bank)

    def tearDown(self):
        shutil.rmtree(self.cfg.working_dir)

    def test_system_to_parse(self):
        system_2_detect = ['system_1']
        system_2_parse = self.parser.system_to_parse(system_2_detect)
        system_2_parse.sort()
        self.assertListEqual(system_2_parse, ['system_1', 'system_2'])
        system_2_detect = ['nimportnaoik']
        with self.assertRaises(TxsscanError) as context:
            self.parser.system_to_parse(system_2_detect)

    def test_parse(self):
        system_2_detect = ['system_1']
        system_2_parse = self.parser.system_to_parse(system_2_detect)
        self.parser.parse(system_2_detect)
        self.assertEqual(len(self.system_bank), 2)

        s2 = self.system_bank['system_2']
        self.assertEqual(s2.name, 'system_2')

        s1 = self.system_bank['system_1']
        self.assertEqual(s1.name, 'system_1')
        self.assertEqual(s1.inter_gene_max_space, 20)
        self.assertEqual(s1.min_mandatory_genes_required, 4)
        self.assertEqual(s1.min_genes_required, 6)
        self.assertTrue(s1.multi_loci)
        self.assertFalse(s2.multi_loci)
        self.assertIsNone(s1.max_nb_genes)
        self.assertEqual(s2.max_nb_genes, 1)
        self.assertEqual(len(s1.mandatory_genes), 5)
        mandatory_genes_name = [ g.name for g in s1.mandatory_genes ]
        mandatory_genes_name.sort()
        theoric_list = ["sctJ_FLG", "sctN_FLG", "flgB", "flgC", "fliE"]
        theoric_list.sort()
        self.assertListEqual(mandatory_genes_name, theoric_list)
        sctJ_FLG = [ g for g in s1.mandatory_genes if g.name == 'sctJ_FLG'][0]
        sctJ_FLG_homologs = sctJ_FLG.get_homologs()
        self.assertEqual(len(sctJ_FLG_homologs), 1)
        self.assertEqual(sctJ_FLG_homologs[0].name, 'sctJ')
        self.assertEqual(sctJ_FLG_homologs[0].system, self.system_bank['system_2'])
        self.assertEqual(len(s1.accessory_genes), 1)
        self.assertEqual(s1.accessory_genes[0].name, 'tadZ')
        self.assertEqual(s1.accessory_genes[0].system, self.system_bank['system_2'])
        self.assertEqual(len(s1.forbidden_genes), 1)
        self.assertEqual(s1.forbidden_genes[0].name, 'sctC')
        self.assertEqual(s1.forbidden_genes[0].system, self.system_bank['system_2'])

    def test_wo_presence(self):
        system_2_detect = ['fail_wo_presence']
        with self.assertRaises(SyntaxError) as context:
            self.parser.parse(system_2_detect)
        self.assertEqual(context.exception.message, "Invalid system definition: gene without presence")

    def test_invalid_presence(self):
        system_2_detect = ['fail_invalid_presence']
        with self.assertRaises(SyntaxError) as context:
            self.parser.parse(system_2_detect)
        self.assertEqual(context.exception.message, "Invalid system definition: presence value must be either [mandatory, accessory, forbidden] not foo_bar")

    def test_gene_no_name(self):
        system_2_detect = ['gene_no_name']
        with self.assertRaises(SyntaxError) as context:
            self.parser.parse(system_2_detect)
        self.assertEqual(context.exception.message,"Invalid system definition: gene without a name")

    def test_invalid_aligned(self):
        system_2_detect = ['invalid_aligned']
        with self.assertRaises(SyntaxError) as context:
            self.parser.parse(system_2_detect)
        self.assertEqual(context.exception.message,'Invalid system definition: invalid value for an attribute of gene foo_bar: sctJ allowed values are "1", "true", "True", "0", "false", "False"')
 
    def test_invalid_homolg(self):
        system_2_detect = ['invalid_homolog']
        with self.assertRaises(SystemInconsistencyError) as context:
            self.parser.parse(system_2_detect)
        self.assertEqual(context.exception.message,  'The gene foo_bar described as homolog of sctJ in system invalid_homolog is not in the "GeneBank" gene factory')

    def test_bad_sys_ref(self):
        system_2_detect = ['bad_sys_ref']
        with self.assertRaises(SystemInconsistencyError) as context:
            self.parser.parse(system_2_detect)
        self.assertEqual(context.exception.message,'Inconsistency in systems definitions: the gene sctJ described as homolog of sctJ with system_ref system_1 has an other system in bank (bad_sys_ref)')

    def test_bad_min_genes_required(self):
        system_2_detect = ['bad_min_genes_required']
        with self.assertRaises(SystemInconsistencyError) as context:
            self.parser.parse(system_2_detect)
        self.assertEqual(context.exception.message, 'system bad_min_genes_required is not consistent: min_genes_required 16 must be lesser or equal than the number of "accessory" and "mandatory" components in the system: 6')

    def test_bad_min_mandatory_genes_required(self):
        system_2_detect = ['bad_min_mandatory_genes_required']
        with self.assertRaises(SystemInconsistencyError) as context:
            self.parser.parse(system_2_detect)
        self.assertEqual(context.exception.message, 'system bad_min_mandatory_genes_required is not consistent: min_genes_required 16 must be lesser or equal than the number of "accessory" and "mandatory" components in the system: 6')

    def test_bad_min_mandatory_genes_required_2(self):
        system_2_detect = ['bad_min_mandatory_genes_required_2']
        with self.assertRaises(SystemInconsistencyError) as context:
            self.parser.parse(system_2_detect)
        self.assertEqual(context.exception.message, "min_genes_required must be greater or equal than min_mandatory_genes_required")

    def test_bad_max_nb_genes(self):
        system_2_detect = ['bad_max_nb_genes']
        with self.assertRaises(SyntaxError) as context:
            self.parser.parse(system_2_detect)
        self.assertEqual(context.exception.message, 'Invalid system definition: max_nb_genes must be an integer: HOHOHO')
