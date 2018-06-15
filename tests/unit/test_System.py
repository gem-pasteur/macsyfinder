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
from macsypy.config import Config
from macsypy.system import System
from macsypy.gene import Gene
from macsypy.gene import Homolog
from macsypy.gene import Analog
from macsypy.registries import ModelRegistry
from tests import MacsyTest


class Test(MacsyTest):

    def setUp(self):
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        
        #add only one handler to the macsypy logger
        from macsypy.gene import _log
        macsy_log = _log.parent
        log_file = os.devnull
        log_handler = logging.FileHandler(log_file)
        macsy_log.addHandler(log_handler)
        
        self.cfg = Config(hmmer_exe="",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type="gembase",
                          e_value_res=1,
                          i_evalue_sel=0.5,
                          models_dir=self.find_data('models'),
                          res_search_dir=tempfile.gettempdir(),
                          res_search_suffix="",
                          profile_suffix=".hmm",
                          res_extract_suffix="",
                          log_level=30,
                          log_file=log_file
                          )
        models_registry = ModelRegistry(self.cfg)
        self.model_name = 'foo'
        self.models_location = models_registry[self.model_name]

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


    def test_name(self):
        name = 'foo'
        system = System(self.cfg, name, 10)
        self.assertEqual(system.name, name)


    def test_inter_gene_max_space(self):
        name = 'foo'
        inter_gene_max_space = 40
        system = System(self.cfg, name, inter_gene_max_space)
        self.assertEqual(system.inter_gene_max_space, inter_gene_max_space)


    def test_min_genes_required(self):
        name = 'foo'
        min_genes_required = 40
        system = System(self.cfg, name, 10, min_genes_required=min_genes_required)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        system.add_mandatory_gene( gene )
        self.assertEqual(system.min_genes_required, min_genes_required)
        # see https://projets.pasteur.fr/issues/1850
        system = System(self.cfg, name, 10)
        self.assertEqual(system.min_genes_required, len(system.mandatory_genes))
        
        
    def test_min_mandatory_genes_required(self):
        name = 'foo'
        min_mandatory_genes_required = 40
        system = System(self.cfg, name, 10, min_mandatory_genes_required=min_mandatory_genes_required)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        system.add_mandatory_gene(gene)
        self.assertEqual(system.min_mandatory_genes_required, min_mandatory_genes_required)
        # see https://projets.pasteur.fr/issues/1850
        system = System(self.cfg, name, 10)
        self.assertEqual(system.min_mandatory_genes_required, len(system.mandatory_genes))


    def test_max_nb_genes(self):
        name = 'foo'
        inter_gene_max_space = 40
        max_nb_genes = 10
        system = System(self.cfg, name, inter_gene_max_space, max_nb_genes=max_nb_genes)
        self.assertEqual(system.max_nb_genes, max_nb_genes)
        name = 'bar'
        system = System(self.cfg, name, inter_gene_max_space)
        self.assertIsNone(system.max_nb_genes)


    def test_multi_loci(self):
        name = 'True'
        inter_gene_max_space = 40
        system = System(self.cfg, name, inter_gene_max_space, multi_loci=True)
        self.assertTrue(system.multi_loci)
        name = 'False'
        inter_gene_max_space = 40
        system = System(self.cfg, name, inter_gene_max_space)
        self.assertFalse(system.multi_loci)


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
