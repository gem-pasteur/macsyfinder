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
from macsypy.report import Hit
from macsypy.search_systems import SystemOccurence
from macsypy.registries import ModelRegistry
from tests import MacsyTest


class Test(MacsyTest):

    def setUp(self):
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        
        # add only one handler to the macsypy logger
        from macsypy.report import _log
        macsy_log = _log.parent
        log_file = os.devnull
        log_handler = logging.FileHandler(log_file)
        macsy_log.addHandler(log_handler)

        self.cfg = Config(hmmer_exe="",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type="gembase",
                          e_value_res=1,
                          i_evalue_sel=0.5,
                          res_search_dir=tempfile.gettempdir(),
                          res_search_suffix="",
                          profile_suffix=".hmm",
                          res_extract_suffix="",
                          log_level=30,
                          models_dir=self.find_data('models'),
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
        
    def test_state(self):
        system = System(self.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        state = system_occurence.state
        self.assertEqual(state, 'empty')

    def test_decision_rule(self):

        # test 'empty' state
        system = System(self.cfg, 'foo', 10, min_mandatory_genes_required=20, min_genes_required=40)
        system_occurence = SystemOccurence(system)
        system_occurence.decision_rule()
        self.assertEqual(system_occurence.state, 'empty')

        # test 'single_locus' state
        system = System(self.cfg, 'foo', 10, min_mandatory_genes_required=1, min_genes_required=2)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        system.add_mandatory_gene(gene)
        gene = Gene(self.cfg, 'tadZ', system, self.models_location)
        system.add_accessory_gene(gene)
        system_occurence = SystemOccurence(system)
        system_occurence.mandatory_genes['sctJ_FLG'] = 1 # simulate match
        system_occurence.accessory_genes['tadZ'] = 1 # simulate match
        system_occurence.nb_cluster = 1
        system_occurence.decision_rule()
        self.assertEqual(system_occurence.state, 'single_locus')

        # test 'multi_loci' state
        system = System(self.cfg, 'foo', 10, min_mandatory_genes_required=1, min_genes_required=2)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        system.add_mandatory_gene(gene)
        gene = Gene(self.cfg, 'tadZ', system, self.models_location)
        system.add_accessory_gene(gene)
        system_occurence = SystemOccurence(system)
        system_occurence.mandatory_genes['sctJ_FLG'] = 1 # simulate match
        system_occurence.accessory_genes['tadZ'] = 1 # simulate match
        system_occurence.nb_cluster = 2
        system_occurence.decision_rule()
        self.assertEqual(system_occurence.state, 'multi_loci')

        # test 'uncomplete' state
        system = System(self.cfg, 'foo', 10, min_mandatory_genes_required=2, min_genes_required=2)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        system.add_mandatory_gene(gene)
        gene = Gene(self.cfg, 'tadZ', system, self.models_location)
        system.add_accessory_gene(gene)
        system_occurence = SystemOccurence(system)
        system_occurence.mandatory_genes['sctJ_FLG'] = 1 # simulate match
        system_occurence.accessory_genes['tadZ'] = 1 # simulate match
        system_occurence.nb_cluster = 2
        system_occurence.decision_rule()
        self.assertEqual(system_occurence.state, 'uncomplete')

        # test 'exclude' state
        system = System(self.cfg, 'foo', 10, min_mandatory_genes_required=2, min_genes_required=2)
        gene = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        system.add_mandatory_gene(gene)
        gene = Gene(self.cfg, 'tadZ', system, self.models_location)
        system.add_accessory_gene(gene)
        gene = Gene(self.cfg, 'fliE', system, self.models_location)
        system.add_forbidden_gene(gene)
        system_occurence = SystemOccurence(system)
        system_occurence.mandatory_genes['sctJ_FLG'] = 1 # simulate match
        system_occurence.accessory_genes['tadZ'] = 1 # simulate match
        system_occurence.forbidden_genes['fliE'] = 1 # simulate match
        system_occurence.nb_cluster = 2
        system_occurence.decision_rule()
        self.assertEqual(system_occurence.state, 'exclude')

    def test_fill_with_multi_systems_genes(self):

        def hit_mock(gene_name):
            li = [None] * 12
            hit = Hit(*li)
            hit.gene = Gene(self.cfg, gene_name, system, self.models_location)
            return hit

        def multi_systems_hits_mock():
            multi_systems_hits = []
            multi_systems_hits.append(hit_mock("tadZ"))
            multi_systems_hits.append(hit_mock("fliE"))
            return multi_systems_hits
            
        system = System(self.cfg, 'foo', 10, min_mandatory_genes_required=20, min_genes_required=40)
        system_occurence = SystemOccurence(system)

        multi_systems_hits = multi_systems_hits_mock() # create multi system genes (genes found in other systems)

        system_occurence.multi_syst_genes = { "tadZ":0 } # create one missing multi system gene

        system_occurence.mandatory_genes = {}
        system_occurence.accessory_genes = { "tadZ":0 } # create one accessory gene
        system_occurence.fill_with_multi_systems_genes(multi_systems_hits)
        self.assertEqual(system_occurence.accessory_genes["tadZ"], 1)

        system_occurence.accessory_genes = {}
        system_occurence.mandatory_genes = { "tadZ":0 } # create one mandatory gene
        system_occurence.fill_with_multi_systems_genes(multi_systems_hits)
        self.assertEqual(system_occurence.mandatory_genes["tadZ"], 1)
