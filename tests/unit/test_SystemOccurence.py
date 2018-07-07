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
from macsypy.gene import Gene, Analog
from macsypy.report import Hit
from macsypy.search_systems import SystemOccurence
from macsypy.registries import ModelRegistry
from macsypy.database import RepliconDB, Indexes
from macsypy.macsypy_error import SystemDetectionError
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

        self.cfg = Config(hmmer_exe="hmmsearch",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type="gembase",
                          e_value_res=1,
                          i_evalue_sel=0.5,
                          res_search_dir=tempfile.gettempdir(),
                          res_search_suffix=".search_hmm.out",
                          profile_suffix=".hmm",
                          res_extract_suffix="",
                          log_level=30,
                          models_dir=self.find_data('models'),
                          log_file=log_file
                         )

        shutil.copy(self.cfg.sequence_db, self.cfg.working_dir)
        self.cfg.options['sequence_db'] = os.path.join(self.cfg.working_dir, os.path.basename(self.cfg.sequence_db))

        idx = Indexes(self.cfg)
        idx._build_my_indexes()

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

        system_occurence.multi_syst_genes = {"tadZ":0} # create one missing multi system gene

        system_occurence.mandatory_genes = {}
        system_occurence.accessory_genes = {"tadZ":0} # create one accessory gene
        system_occurence.fill_with_multi_systems_genes(multi_systems_hits)
        self.assertEqual(system_occurence.accessory_genes["tadZ"], 1)

        system_occurence.accessory_genes = {}
        system_occurence.mandatory_genes = {"tadZ":0} # create one mandatory gene
        system_occurence.fill_with_multi_systems_genes(multi_systems_hits)
        self.assertEqual(system_occurence.mandatory_genes["tadZ"], 1)

    def test_get_gene_counter_output(self):
        system = System(self.cfg, 'foo', 10, min_mandatory_genes_required=20, min_genes_required=40)
        system_occurence = SystemOccurence(system)

        system_occurence.accessory_genes = {"tadZ":0} # create one accessory gene
        system_occurence.mandatory_genes = {"fliE":0} # create one mandatory gene
        system_occurence.forbidden_genes = {"gspD":0} # create one forbiden gene

        out = system_occurence.get_gene_counter_output(True)
        expected = "{'fliE': 0}\t{'tadZ': 0}\t{}"
        self.assertEqual(out, expected)

        out = system_occurence.get_gene_counter_output()
        expected = "{'fliE': 0}\t{'tadZ': 0}\t{'gspD': 0}"
        self.assertEqual(out, expected)

    def test_nb_syst_genes(self):
        system = System(self.cfg, 'foo', 10, min_mandatory_genes_required=20, min_genes_required=40)
        system_occurence = SystemOccurence(system)
        system_occurence._nb_syst_genes = 3
        self.assertEqual(system_occurence.nb_syst_genes, 3)

    def test_get_gene_ref(self):

        def create_analog(system):
            gene_ref = Gene(self.cfg, 'tadZ', system, self.models_location)
            gene_analog = Gene(self.cfg, 'sctC', system, self.models_location)
            analog = Analog(gene_analog, gene_ref)
            return analog

        # test case 1

        system = System(self.cfg, 'foo', 10, min_mandatory_genes_required=20, min_genes_required=40)
        analog = create_analog(system)
        gene = Gene(self.cfg, 'fliE', system, self.models_location) # create regular gene
        gene.add_analog(analog) # attach analog to regular gene

        system.add_mandatory_gene(gene)
        system_occurence = SystemOccurence(system)

        gref = system_occurence.get_gene_ref(analog)
        self.assertEqual(gref.name, 'tadZ')

        # test case 2

        system = System(self.cfg, 'foo', 10, min_mandatory_genes_required=20, min_genes_required=40)
        gene = Gene(self.cfg, 'fliE', system, self.models_location) # create regular gene

        system.add_mandatory_gene(gene)
        system_occurence = SystemOccurence(system)

        gref = system_occurence.get_gene_ref(gene)
        self.assertEqual(gref, None)

    def test_compute_nb_syst_genes(self):
        system = System(self.cfg, 'foo', 10, min_mandatory_genes_required=2, min_genes_required=2)
        gene = Gene(self.cfg, 'sctJ', system, self.models_location)
        system.add_mandatory_gene(gene)
        gene = Gene(self.cfg, 'tadZ', system, self.models_location)
        system.add_accessory_gene(gene)

        system_occurence = SystemOccurence(system)

        system_occurence.mandatory_genes['sctJ'] = 1 # simulate match
        system_occurence.accessory_genes['tadZ'] = 4 # simulate match
        nb = system_occurence.compute_nb_syst_genes()
        self.assertEqual(nb, 2)

    def test_str(self):
        system = System(self.cfg, 'foo', 10, min_mandatory_genes_required=2, min_genes_required=2)

        gene = Gene(self.cfg, 'sctJ', system, self.models_location)
        system.add_mandatory_gene(gene)

        gene = Gene(self.cfg, 'tadZ', system, self.models_location)
        system.add_accessory_gene(gene)

        gene = Gene(self.cfg, 'flgC', system, self.models_location)
        system.add_forbidden_gene(gene)

        system_occurence = SystemOccurence(system)

        gene = Gene(self.cfg, 'gspD', system, self.models_location)
        system_occurence.multi_syst_genes[gene.name] = 0

        out = system_occurence.__str__()
        expected = 'sctJ\t0\ntadZ\t0\nflgC\t0\ngspD\t0\n'
        self.assertEqual(out, expected)

    def test_compute_system_length(self):

        system = System(self.cfg, 'foo', 10, min_mandatory_genes_required=2, min_genes_required=2)
        system_occurence = SystemOccurence(system)

        self.cfg.options['topology_file'] = self.cfg.sequence_db + ".topo"
        db_send = {'ESCO030p01':'circular', 'PSAE001c01':'linear'}
        with open(self.cfg.topology_file, 'w') as f:
            for k, v in db_send.items():
                f.write('{0} : {1}\n'.format(k, v))
        db = RepliconDB(self.cfg)

        rep_info = db['PSAE001c01']
        system_occurence.loci_positions = [(5, 10), (10, 20)]
        length = system_occurence.compute_system_length(rep_info)
        self.assertEqual(length, 17)

        rep_info = db['NC_xxxxx_xx']
        system_occurence.loci_positions = [(10, 5), (20, 10)]
        length = system_occurence.compute_system_length(rep_info)
        self.assertEqual(length, 3)

        rep_info = db['PSAE001c01']
        system_occurence.loci_positions = [(10, 5), (20, 10)]
        with self.assertRaises(SystemDetectionError) as context:
            length = system_occurence.compute_system_length(rep_info)
        self.assertEqual(context.exception.message,
                         "Inconsistency in locus positions in the case of a linear replicon. "
                         "The begin position of a locus cannot be higher than the end position. \n"
                         "Problem with locus found with positions begin: 10 end: 5")

    def test_compute_nb_syst_genes_tot(self):
        system = System(self.cfg, 'foo', 10, min_mandatory_genes_required=2, min_genes_required=2)
        gene = Gene(self.cfg, 'sctJ', system, self.models_location)
        system.add_mandatory_gene(gene)
        gene = Gene(self.cfg, 'tadZ', system, self.models_location)
        system.add_accessory_gene(gene)

        system_occurence = SystemOccurence(system)

        system_occurence.mandatory_genes['sctJ'] = 1 # simulate match
        system_occurence.accessory_genes['tadZ'] = 4 # simulate match
        nb = system_occurence.compute_nb_syst_genes_tot()
        self.assertEqual(nb, 5)

    def test_get_system_name_unordered(self):
        system = System(self.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)

        name = system_occurence.get_system_name_unordered()
        self.assertEqual(name, 'foo_putative')

        name = system_occurence.get_system_name_unordered('_bar')
        self.assertEqual(name, 'foo_bar')

    def test_count_genes(self):
        system = System(self.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        genes = {'ACBA007p01':2, 'ZIIN001c01':0}
        total = system_occurence.count_genes(genes)
        self.assertEqual(total, 1)

    def test_count_genes_tot(self):
        system = System(self.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        genes = {'ACBA007p01':2, 'ZIIN001c01':1}
        total = system_occurence.count_genes_tot(genes)
        self.assertEqual(total, 3)
