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
import argparse

from macsypy.config import Config, MacsyDefaults
from macsypy.system import System
from macsypy.gene import Gene
from macsypy.gene import Homolog
from macsypy.gene import Analog
from macsypy.registries import ModelRegistry
from tests import MacsyTest


class Test(MacsyTest):

    def setUp(self):
        self.args = argparse.Namespace()
        self.args.sequence_db = self.find_data("base", "test_base.fa")
        self.args.db_type = 'gembase'
        self.args.models_dir = self.find_data('models')
        self.args.res_search_dir = tempfile.gettempdir()
        self.args.log_level = 30
        self.args.out_dir = os.path.join(self.args.res_search_dir,
                                    'test_macsyfinder_System')
        if os.path.exists(self.args.out_dir):
            shutil.rmtree(self.args.out_dir)
        os.mkdir(self.args.out_dir)

        self.cfg = Config(MacsyDefaults(), self.args)
        models_registry = ModelRegistry(self.cfg)
        self.model_name = 'foo'
        self.models_location = models_registry[self.model_name]

    def tearDown(self):
        self.clean_working_dir()

    def clean_working_dir(self):
        try:
            shutil.rmtree(self.cfg.working_dir())
        except:
            pass

    def test_name(self):
        name = 'foo'
        system = System(self.cfg, name, 10)
        self.assertEqual(system.name, name)


    def test_inter_gene_max_space(self):
        system_fqn = 'foo/bar'
        inter_gene_max_space_xml = 40
        # test inter_gene_max_space from xml
        system = System(self.cfg, system_fqn, inter_gene_max_space_xml)
        self.assertEqual(system.inter_gene_max_space, inter_gene_max_space_xml)

        self.clean_working_dir()

        # test inter_gene_max_space overloaded by the config
        inter_gene_max_space_cfg = [[system_fqn, '22']]

        self.args.inter_gene_max_space = inter_gene_max_space_cfg
        cfg = Config(MacsyDefaults(), self.args)
        system = System(cfg, system_fqn, inter_gene_max_space_xml)
        self.assertEqual(system.inter_gene_max_space, int(inter_gene_max_space_cfg[0][1]))

    def test_min_genes_required(self):
        system_fqn = 'foo/bar'
        min_genes_required_xml = 40
        system = System(self.cfg, system_fqn, 10, min_genes_required=min_genes_required_xml)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        system.add_mandatory_gene(gene)
        self.assertEqual(system.min_genes_required, min_genes_required_xml)
        # see https://projets.pasteur.fr/issues/1850
        system = System(self.cfg, system_fqn, 10)
        self.assertEqual(system.min_genes_required, len(system.mandatory_genes))

        self.clean_working_dir()

        # test inter_gene_max_space overloaded by the config
        min_genes_required_cfg = [[system_fqn, '22']]

        self.args.min_genes_required = min_genes_required_cfg
        cfg = Config(MacsyDefaults(), self.args)
        system = System(cfg, system_fqn, min_genes_required_xml)
        self.assertEqual(system.min_genes_required, int(min_genes_required_cfg[0][1]))


    def test_min_mandatory_genes_required(self):
        system_fqn = 'foo/bar'
        min_mandatory_genes_required_xml = 40
        system = System(self.cfg, system_fqn, 10, min_mandatory_genes_required=min_mandatory_genes_required_xml)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        system.add_mandatory_gene(gene)
        self.assertEqual(system.min_mandatory_genes_required, min_mandatory_genes_required_xml)
        # see https://projets.pasteur.fr/issues/1850
        system = System(self.cfg, system_fqn, 10)
        self.assertEqual(system.min_mandatory_genes_required, len(system.mandatory_genes))

        self.clean_working_dir()

        # test inter_gene_max_space overloaded by the config
        min_mandatory_genes_required_cfg = [[system_fqn, '22']]

        self.args.min_mandatory_genes_required = min_mandatory_genes_required_cfg
        cfg = Config(MacsyDefaults(), self.args)
        system = System(cfg, system_fqn, 10)
        self.assertEqual(system.min_mandatory_genes_required, int(min_mandatory_genes_required_cfg[0][1]))

    def test_max_nb_genes(self):
        system_fqn = 'foo/bar'
        inter_gene_max_space = 40
        max_nb_genes_xml = 10
        system = System(self.cfg, system_fqn, inter_gene_max_space, max_nb_genes=max_nb_genes_xml)
        self.assertEqual(system.max_nb_genes, max_nb_genes_xml)
        system = System(self.cfg, system_fqn, inter_gene_max_space)
        self.assertIsNone(system.max_nb_genes)

        self.clean_working_dir()

        # test inter_gene_max_space overloaded by the config
        max_nb_genes_cfg = [[system_fqn, '22']]

        self.args.max_nb_genes = max_nb_genes_cfg
        cfg = Config(MacsyDefaults(), self.args)
        system = System(cfg, system_fqn, inter_gene_max_space, max_nb_genes=max_nb_genes_xml)
        self.assertEqual(system.max_nb_genes, int(max_nb_genes_cfg[0][1]))


    def test_multi_loci(self):
        system_fqn = 'foo/True'
        inter_gene_max_space = 40
        system = System(self.cfg, system_fqn, inter_gene_max_space, multi_loci=True)
        self.assertTrue(system.multi_loci)
        system_fqn = 'foo/False'
        inter_gene_max_space = 40
        system = System(self.cfg, system_fqn, inter_gene_max_space)
        self.assertFalse(system.multi_loci)

        self.clean_working_dir()

        self.args.multi_loci = 'foo/False'
        cfg = Config(MacsyDefaults(), self.args)

        system_fqn = 'foo/False'
        inter_gene_max_space = 40
        system = System(cfg, system_fqn, inter_gene_max_space, multi_loci=False)
        self.assertTrue(system.multi_loci)


    def test_add_mandatory_gene(self):
        system = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        system.add_mandatory_gene(gene)
        self.assertEqual(system._mandatory_genes, [gene])
        self.assertEqual(system._accessory_genes, [])
        self.assertEqual(system._forbidden_genes, [])


    def test_add_accessory_gene(self):
        system = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        system.add_accessory_gene(gene)
        self.assertEqual(system._accessory_genes, [gene])
        self.assertEqual(system._mandatory_genes, [])
        self.assertEqual(system._forbidden_genes, [])


    def test_add_forbidden_gene(self):
        system = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        system.add_forbidden_gene(gene)
        self.assertEqual(system._forbidden_genes, [gene])
        self.assertEqual(system._accessory_genes, [])
        self.assertEqual(system._mandatory_genes, [])

    def test_mandatory_genes(self):
        system = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        system.add_mandatory_gene(gene)
        self.assertEqual(system.mandatory_genes, [gene])


    def test_accessory_genes(self):
        system = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        system.add_accessory_gene(gene)
        self.assertEqual(system.accessory_genes, [gene])


    def test_forbidden_genes(self):
        system = System(self.cfg, "foo", 10)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        system.add_forbidden_gene(gene)
        self.assertEqual(system.forbidden_genes, [gene])


    def test_get_gene(self):
        system = System(self.cfg, "foo", 10)
        gene_name = 'sctJ_FLG'
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        for meth in (system.add_forbidden_gene, system.add_accessory_gene, system.add_mandatory_gene):
            system._mandatory_genes = []
            system._accessory_genes = []
            system._forbidden_genes = []
            meth(gene)
            self.assertEqual(gene, system.get_gene(gene_name))

        self.assertRaises(KeyError, system.get_gene, 'bar')

        homolog_name = 'sctJ'
        gene_homolog = Gene(self.cfg, homolog_name, system, self.models_location)
        homolog = Homolog(gene_homolog, gene)
        gene.add_homolog(homolog)
        for meth in (system.add_forbidden_gene, system.add_accessory_gene, system.add_mandatory_gene):
            system._mandatory_genes = []
            system._accessory_genes = []
            system._forbidden_genes = []
            meth(gene)
            self.assertEqual(homolog, system.get_gene(homolog_name))

        analog_name = 'sctC'
        gene_analog = Gene(self.cfg, analog_name, system, self.models_location)
        analog = Analog(gene_analog, gene)
        gene.add_analog(analog)
        for meth in (system.add_forbidden_gene, system.add_accessory_gene, system.add_mandatory_gene):
            system._mandatory_genes = []
            system._accessory_genes = []
            system._forbidden_genes = []
            meth(gene)
            self.assertEqual(analog, system.get_gene(analog_name))

    def test_get_gene_ref(self):
        system = System(self.cfg, "foo", 10)
        gene_name = 'sctJ_FLG'
        gene_ref = Gene(self.cfg, gene_name, system, self.models_location)
        homolog_name = 'sctJ'
        gene_homolg = Gene(self.cfg, homolog_name, system, self.models_location)
        homolog = Homolog(gene_homolg, gene_ref)
        gene_ref.add_homolog(homolog)

        for meth in (system.add_forbidden_gene, system.add_accessory_gene, system.add_mandatory_gene):
            system._mandatory_genes = []
            system._accessory_genes = []
            system._forbidden_genes = []
            meth(gene_ref)
            self.assertEqual(gene_ref, system.get_gene_ref(homolog))
        self.assertIsNone(system.get_gene_ref(gene_ref))
        gene_ukn = Gene(self.cfg, 'abc', system, self.models_location)
        self.assertRaises(KeyError, system.get_gene_ref, gene_ukn)

    def test_str(self):
        system_fqn = "foo/bar"
        system = System(self.cfg, system_fqn, 10)
        mandatory_gene = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        system.add_mandatory_gene(mandatory_gene)
        homolog_name = 'sctJ'
        gene_homolg = Gene(self.cfg, homolog_name, system, self.models_location)
        homolog = Homolog(gene_homolg, mandatory_gene)
        mandatory_gene.add_homolog(homolog)

        accessory_gene = Gene(self.cfg, 'sctN_FLG', system, self.models_location)
        system.add_accessory_gene(accessory_gene)
        analog_name = 'sctN'
        gene_analog = Gene(self.cfg, analog_name, system, self.models_location)
        analog = Analog(gene_analog, accessory_gene)
        accessory_gene.add_analog(analog)

        forbidden_gene = Gene(self.cfg, 'sctC', system, self.models_location)
        system.add_forbidden_gene(forbidden_gene)

        exp_str = """name: bar
fqn: foo/bar
==== mandatory genes ====
sctJ_FLG
==== accessory genes ====
sctN_FLG
==== forbidden genes ====
sctC
============== end pprint system ================
"""
        self.assertEqual(str(system), exp_str)

    def test_eq(self):
        aa1 = System(self.cfg, "aaa", 10)
        aa2 = System(self.cfg, "aaa", 10)
        self.assertEqual(aa1, aa2)

    def test_lt(self):
        aaa = System(self.cfg, "aaa", 10)
        zzz = System(self.cfg, "zzz", 10)
        self.assertLess(aaa, zzz)

    def test_gt(self):
        aaa = System(self.cfg, "aaa", 10)
        zzz = System(self.cfg, "zzz", 10)
        self.assertGreater(zzz, aaa)


