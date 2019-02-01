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
import argparse
from operator import attrgetter

from macsypy.config import Config, MacsyDefaults
from macsypy.model import Model
from macsypy.gene import Gene, Analog
from macsypy.report import Hit
from macsypy.search_systems import SystemOccurence, build_clusters
from macsypy.database import RepliconDB, Indexes
from macsypy.search_genes import search_genes
from macsypy.definition_parser import DefinitionParser
from macsypy.error import SystemDetectionError

from tests import MacsyTest
from tests.macsy_test_env import MacsyEnvManager


class Test(MacsyTest, MacsyEnvManager):

    def setUp(self):
        self.load_env('env_001', log_out=False)


    def tearDown(self):
        self.unload_env('env_001')


    def test_state(self):
        system = Model(self.macsy_test_env.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        state = system_occurence.state
        self.assertEqual(state, 'empty')


    def test_decision_rule(self):
        # test 'empty' state
        system = Model(self.macsy_test_env.cfg, 'foo', 10, min_mandatory_genes_required=20, min_genes_required=40)
        system_occurence = SystemOccurence(system)
        system_occurence.decision_rule()
        self.assertEqual(system_occurence.state, 'empty')

        # test 'single_locus' state
        system = Model(self.macsy_test_env.cfg, 'foo', 10, min_mandatory_genes_required=1, min_genes_required=2)
        gene = Gene(self.macsy_test_env.cfg, 'sctJ_FLG', system, self.macsy_test_env.models_location)
        system.add_mandatory_gene(gene)
        gene = Gene(self.macsy_test_env.cfg, 'tadZ', system, self.macsy_test_env.models_location)
        system.add_accessory_gene(gene)
        system_occurence = SystemOccurence(system)
        system_occurence.mandatory_genes['sctJ_FLG'] = 1  # simulate match
        system_occurence.accessory_genes['tadZ'] = 1  # simulate match
        system_occurence.nb_cluster = 1
        system_occurence.decision_rule()
        self.assertEqual(system_occurence.state, 'single_locus')

        # test 'multi_loci' state
        system = Model(self.macsy_test_env.cfg, 'foo', 10, min_mandatory_genes_required=1, min_genes_required=2)
        gene = Gene(self.macsy_test_env.cfg, 'sctJ_FLG', system, self.macsy_test_env.models_location)
        system.add_mandatory_gene(gene)
        gene = Gene(self.macsy_test_env.cfg, 'tadZ', system, self.macsy_test_env.models_location)
        system.add_accessory_gene(gene)
        system_occurence = SystemOccurence(system)
        system_occurence.mandatory_genes['sctJ_FLG'] = 1  # simulate match
        system_occurence.accessory_genes['tadZ'] = 1  # simulate match
        system_occurence.nb_cluster = 2
        system_occurence.decision_rule()
        self.assertEqual(system_occurence.state, 'multi_loci')

        # test 'uncomplete' state
        system = Model(self.macsy_test_env.cfg, 'foo', 10, min_mandatory_genes_required=2, min_genes_required=2)
        gene = Gene(self.macsy_test_env.cfg, 'sctJ_FLG', system, self.macsy_test_env.models_location)
        system.add_mandatory_gene(gene)
        gene = Gene(self.macsy_test_env.cfg, 'tadZ', system, self.macsy_test_env.models_location)
        system.add_accessory_gene(gene)
        system_occurence = SystemOccurence(system)
        system_occurence.mandatory_genes['sctJ_FLG'] = 1  # simulate match
        system_occurence.accessory_genes['tadZ'] = 1  # simulate match
        system_occurence.nb_cluster = 2
        system_occurence.decision_rule()
        self.assertEqual(system_occurence.state, 'uncomplete')

        # test 'exclude' state
        system = Model(self.macsy_test_env.cfg, 'foo', 10, min_mandatory_genes_required=2, min_genes_required=2)
        gene = Gene(self.macsy_test_env.cfg, 'sctJ_FLG', system, self.macsy_test_env.models_location)
        system.add_mandatory_gene(gene)
        gene = Gene(self.macsy_test_env.cfg, 'tadZ', system, self.macsy_test_env.models_location)
        system.add_accessory_gene(gene)
        gene = Gene(self.macsy_test_env.cfg, 'fliE', system, self.macsy_test_env.models_location)
        system.add_forbidden_gene(gene)
        system_occurence = SystemOccurence(system)
        system_occurence.mandatory_genes['sctJ_FLG'] = 1  # simulate match
        system_occurence.accessory_genes['tadZ'] = 1  # simulate match
        system_occurence.forbidden_genes['fliE'] = 1  # simulate match
        system_occurence.nb_cluster = 2
        system_occurence.decision_rule()
        self.assertEqual(system_occurence.state, 'exclude')


    def test_fill_with_multi_systems_genes(self):

        def hit_mock(gene_name):
            li = [None] * 12
            hit = Hit(*li)
            hit.gene = Gene(self.macsy_test_env.cfg, gene_name, system, self.macsy_test_env.models_location)
            return hit

        def multi_systems_hits_mock():
            multi_systems_hits = []
            multi_systems_hits.append(hit_mock("tadZ"))
            multi_systems_hits.append(hit_mock("fliE"))
            return multi_systems_hits

        system = Model(self.macsy_test_env.cfg, 'foo', 10, min_mandatory_genes_required=20, min_genes_required=40)
        system_occurence = SystemOccurence(system)

        multi_systems_hits = multi_systems_hits_mock()  # create multi system genes (genes found in other systems)

        system_occurence.multi_syst_genes = {"tadZ": 0}  # create one missing multi system gene

        system_occurence.mandatory_genes = {}
        system_occurence.accessory_genes = {"tadZ": 0}  # create one accessory gene
        with self.catch_log():
            system_occurence.fill_with_multi_systems_genes(multi_systems_hits)
        self.assertEqual(system_occurence.accessory_genes["tadZ"], 1)

        system_occurence.accessory_genes = {}
        system_occurence.mandatory_genes = {"tadZ": 0}  # create one mandatory gene
        with self.catch_log():
            system_occurence.fill_with_multi_systems_genes(multi_systems_hits)
        self.assertEqual(system_occurence.mandatory_genes["tadZ"], 1)


    def test_get_gene_counter_output(self):
        system = Model(self.macsy_test_env.cfg, 'foo', 10, min_mandatory_genes_required=20, min_genes_required=40)
        system_occurence = SystemOccurence(system)

        system_occurence.accessory_genes = {"tadZ": 0}  # create one accessory gene
        system_occurence.mandatory_genes = {"fliE": 0}  # create one mandatory gene
        system_occurence.forbidden_genes = {"gspD": 0}  # create one forbiden gene

        out = system_occurence.get_gene_counter_output(True)
        expected = "{'fliE': 0}\t{'tadZ': 0}\t{}"
        self.assertEqual(out, expected)

        out = system_occurence.get_gene_counter_output()
        expected = "{'fliE': 0}\t{'tadZ': 0}\t{'gspD': 0}"
        self.assertEqual(out, expected)


    def test_nb_syst_genes(self):
        system = Model(self.macsy_test_env.cfg, 'foo', 10, min_mandatory_genes_required=20, min_genes_required=40)
        system_occurence = SystemOccurence(system)
        system_occurence._nb_syst_genes = 3
        self.assertEqual(system_occurence.nb_syst_genes, 3)


    def test_get_gene_ref(self):

        def create_analog(system):
            gene_ref = Gene(self.macsy_test_env.cfg, 'tadZ', system, self.macsy_test_env.models_location)
            gene_analog = Gene(self.macsy_test_env.cfg, 'sctC', system, self.macsy_test_env.models_location)
            analog = Analog(gene_analog, gene_ref)
            return analog

        # test case 1

        system = Model(self.macsy_test_env.cfg, 'foo', 10, min_mandatory_genes_required=20, min_genes_required=40)
        analog = create_analog(system)
        gene = Gene(self.macsy_test_env.cfg, 'fliE', system, self.macsy_test_env.models_location)  # create regular gene
        gene.add_analog(analog)  # attach analog to regular gene

        system.add_mandatory_gene(gene)
        system_occurence = SystemOccurence(system)

        gref = system_occurence.get_gene_ref(analog)
        self.assertEqual(gref.name, 'tadZ')

        # test case 2

        system = Model(self.macsy_test_env.cfg, 'foo', 10, min_mandatory_genes_required=20, min_genes_required=40)
        gene = Gene(self.macsy_test_env.cfg, 'fliE', system, self.macsy_test_env.models_location)  # create regular gene

        system.add_mandatory_gene(gene)
        system_occurence = SystemOccurence(system)

        gref = system_occurence.get_gene_ref(gene)
        self.assertEqual(gref, None)


    def test_compute_nb_syst_genes(self):
        system = Model(self.macsy_test_env.cfg, 'foo', 10, min_mandatory_genes_required=2, min_genes_required=2)
        gene = Gene(self.macsy_test_env.cfg, 'sctJ', system, self.macsy_test_env.models_location)
        system.add_mandatory_gene(gene)
        gene = Gene(self.macsy_test_env.cfg, 'tadZ', system, self.macsy_test_env.models_location)
        system.add_accessory_gene(gene)

        system_occurence = SystemOccurence(system)

        system_occurence.mandatory_genes['sctJ'] = 1  # simulate match
        system_occurence.accessory_genes['tadZ'] = 4  # simulate match
        nb = system_occurence.compute_nb_syst_genes()
        self.assertEqual(nb, 2)


    def test_str(self):
        system = Model(self.macsy_test_env.cfg, 'foo', 10, min_mandatory_genes_required=2, min_genes_required=2)

        gene = Gene(self.macsy_test_env.cfg, 'sctJ', system, self.macsy_test_env.models_location)
        system.add_mandatory_gene(gene)

        gene = Gene(self.macsy_test_env.cfg, 'tadZ', system, self.macsy_test_env.models_location)
        system.add_accessory_gene(gene)

        gene = Gene(self.macsy_test_env.cfg, 'flgC', system, self.macsy_test_env.models_location)
        system.add_forbidden_gene(gene)

        system_occurence = SystemOccurence(system)

        gene = Gene(self.macsy_test_env.cfg, 'gspD', system, self.macsy_test_env.models_location)
        system_occurence.multi_syst_genes[gene.name] = 0

        out = system_occurence.__str__()
        expected = 'sctJ\t0\ntadZ\t0\nflgC\t0\ngspD\t0\n'
        self.assertEqual(out, expected)


    def test_compute_system_length(self):
        out_dir = MacsyTest.get_uniq_tmp_dir_name()
        defaults = MacsyDefaults()

        seq_ori = MacsyTest.find_data("base", 'test_base.fa')

        args = argparse.Namespace()
        args.out_dir = out_dir
        args.db_type = 'gembase'
        args.log_level = 30
        args.previous_run = None
        args.models_dir = MacsyTest.find_data("data_set_3", "models")

        args.sequence_db = os.path.join(args.out_dir, os.path.basename(seq_ori))
        args.topology_file = args.sequence_db + ".topo"

        if os.path.exists(args.out_dir):
            shutil.rmtree(args.out_dir)
        os.mkdir(args.out_dir)
        try:
            # copy sequence file in working_dir
            shutil.copy(seq_ori, args.out_dir)

            # create topology file in working_dir
            db_send = {'ESCO030p01': 'circular', 'PSAE001c01': 'linear'}
            with open(args.topology_file, 'w') as f:
                for k, v in db_send.items():
                    f.write('{} : {}\n'.format(k, v))

            cfg = Config(defaults, args)
            idx = Indexes(cfg)
            idx.build()

            system = Model(cfg, 'foo', 10, min_mandatory_genes_required=2, min_genes_required=2)
            system_occurence = SystemOccurence(system)
            db = RepliconDB(cfg)
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
            with self.catch_log():
                with self.assertRaises(SystemDetectionError) as context:
                    length = system_occurence.compute_system_length(rep_info)
                self.assertEqual(str(context.exception),
                                 "Inconsistency in locus positions in the case of a linear replicon. "
                                 "The begin position of a locus cannot be higher than the end position. \n"
                                 "Problem with locus found with positions begin: 10 end: 5")
        finally:
            shutil.rmtree(out_dir)

    def test_compute_nb_syst_genes_tot(self):
        system = Model(self.macsy_test_env.cfg, 'foo', 10, min_mandatory_genes_required=2, min_genes_required=2)
        gene = Gene(self.macsy_test_env.cfg, 'sctJ', system, self.macsy_test_env.models_location)
        system.add_mandatory_gene(gene)
        gene = Gene(self.macsy_test_env.cfg, 'tadZ', system, self.macsy_test_env.models_location)
        system.add_accessory_gene(gene)

        system_occurence = SystemOccurence(system)

        system_occurence.mandatory_genes['sctJ'] = 1  # simulate match
        system_occurence.accessory_genes['tadZ'] = 4  # simulate match
        nb = system_occurence.compute_nb_syst_genes_tot()
        self.assertEqual(nb, 5)


    def test_get_system_name_unordered(self):
        system = Model(self.macsy_test_env.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)

        name = system_occurence.get_system_name_unordered()
        self.assertEqual(name, 'foo_putative')

        name = system_occurence.get_system_name_unordered('_bar')
        self.assertEqual(name, 'foo_bar')


    def test_count_genes(self):
        system = Model(self.macsy_test_env.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        genes = {'ACBA007p01': 2, 'ZIIN001c01': 0}
        total = system_occurence.count_genes(genes)
        self.assertEqual(total, 1)


    def test_count_genes_tot(self):
        system = Model(self.macsy_test_env.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        genes = {'ACBA007p01': 2, 'ZIIN001c01': 1}
        total = system_occurence.count_genes_tot(genes)
        self.assertEqual(total, 3)


    def test_get_system_unique_name(self):
        system = Model(self.macsy_test_env.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        name = system_occurence.get_system_unique_name('ACBA007p01')
        self.assertEqual(name, 'ACBA007p01_foo_1')


    def test_compute_missing_genes_list(self):
        system = Model(self.macsy_test_env.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        genes = {'ACBA007p01': 2, 'ZIIN001c01': 0}
        missing = system_occurence.compute_missing_genes_list(genes)
        self.assertEqual(missing, ['ZIIN001c01'])


    def test_count_missing_genes(self):
        system = Model(self.macsy_test_env.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        genes = {'ACBA007p01': 2, 'ZIIN001c01': 0}
        nb = system_occurence.count_missing_genes(genes)
        self.assertEqual(nb, 1)


    def test_is_complete(self):
        system = Model(self.macsy_test_env.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)

        system_occurence._state = 'multi_loci'
        self.assertEqual(system_occurence.is_complete(), True)

        system_occurence._state = 'single_locus'
        self.assertEqual(system_occurence.is_complete(), True)

        system_occurence._state = 'empty'
        self.assertEqual(system_occurence.is_complete(), False)


    def test_get_summary_header(self):
        system = Model(self.macsy_test_env.cfg, 'foo', 10)
        system_occurence = SystemOccurence(system)
        expect = "#Replicon_name\tSystem_Id\tReference_system\tSystem_status\t" \
                 "Nb_loci\tNb_Ref_mandatory\tNb_Ref_accessory\tNb_Ref_Genes_detected_NR\t" \
                 "Nb_Genes_with_match\tSystem_length\tNb_Mandatory_NR\tNb_Accessory_NR\t" \
                 "Nb_missing_mandatory\tNb_missing_accessory\tList_missing_mandatory\t" \
                 "List_missing_accessory\tLoci_positions\tOccur_Mandatory\tOccur_Accessory\t" \
                 "Occur_Forbidden"
        out = system_occurence.get_summary_header()
        self.assertEqual(out, expect)


    def test_get_summary(self):
        out_dir = MacsyTest.get_uniq_tmp_dir_name()
        defaults = MacsyDefaults()

        seq_ori = MacsyTest.find_data("base", 'test_base.fa')

        args = argparse.Namespace()
        args.out_dir = out_dir
        args.db_type = 'gembase'
        args.log_level = 30
        args.previous_run = None
        args.models_dir = MacsyTest.find_data("data_set_3", "models")

        args.sequence_db = os.path.join(args.out_dir, os.path.basename(seq_ori))
        args.topology_file = args.sequence_db + ".topo"

        if os.path.exists(args.out_dir):
            shutil.rmtree(args.out_dir)
        os.mkdir(args.out_dir)
        try:
            # copy sequence file in working_dir
            shutil.copy(seq_ori, args.out_dir)

            # create topology file in working_dir
            db_send = {'ESCO030p01': 'circular', 'PSAE001c01': 'linear'}
            with open(args.topology_file, 'w') as f:
                for k, v in db_send.items():
                    f.write('{} : {}\n'.format(k, v))

            cfg = Config(defaults, args)
            idx = Indexes(cfg)
            idx.build()

            system = Model(cfg, 'foo', 10, min_mandatory_genes_required=2, min_genes_required=2)
            gene = Gene(cfg, 'tadZ', system, self.macsy_test_env.models_location)
            system.add_accessory_gene(gene)

            system_occurence = SystemOccurence(system)
            db = RepliconDB(cfg)
            rep_info = db['PSAE001c01']

            out = system_occurence.get_summary('PSAE001c01', rep_info)
            expect = "PSAE001c01	PSAE001c01_foo_1	foo	empty	0	0	1	0	0	0	0	0	0	1	[]	['tadZ']	[]	{}	{'tadZ': 0}	{}"
            self.assertEqual(out, expect)
        finally:
            shutil.rmtree(out_dir)


    def test_get_summary_unordered(self):
        system = Model(self.macsy_test_env.cfg, 'foo', 10, min_mandatory_genes_required=2, min_genes_required=2)
        gene = Gene(self.macsy_test_env.cfg, 'sctJ', system, self.macsy_test_env.models_location)
        system.add_mandatory_gene(gene)

        system_occurence = SystemOccurence(system)

        out = system_occurence.get_summary_unordered('ESCO030p01')
        expect = "ESCO030p01	foo_putative	foo	empty	None	1	0	0	0	None	0	0	1	0	['sctJ']	[]	None	{'sctJ': 0}	{}	{}"
        self.assertEqual(out, expect)


    def test_fill_with_cluster(self):
        # for this test, we use a specific configuration with a dedicated
        # working directory (i.e. we don't use the generic configuration
        # defined in setUp() method).
        out_dir = self.get_uniq_tmp_dir_name()

        args = argparse.Namespace()
        args.db_type = 'gembase'
        args.models_dir = self.find_data(os.path.join('data_set_1', 'models'))
        args.previous_run = self.find_data(os.path.join('data_set_1', 'complete_run_results'))
        args.res_search_dir = tempfile.gettempdir()
        args.log_level = 30
        args.out_dir = out_dir
        config = Config(MacsyDefaults(), args)

        if os.path.exists(config.working_dir()):
            shutil(config.working_dir())
        os.makedirs(config.working_dir())
        try:
            idx = Indexes(config)
            idx._build_my_indexes()

            model_bank = self.macsy_test_env.model_bank
            gene_bank = self.macsy_test_env.gene_bank
            parser = DefinitionParser(config, model_bank, gene_bank)
            parser.parse(['set_1/T9SS'])

            system = model_bank['set_1/T9SS']
            genes = system.mandatory_genes + system.accessory_genes + system.forbidden_genes
            ex_genes = []
            for g in genes:
                if g.exchangeable:
                    h_s = g.get_homologs()
                    ex_genes += h_s
                    a_s = g.get_analogs()
                    ex_genes += a_s
            all_genes = (genes + ex_genes)
            all_reports = search_genes(all_genes, config)
            all_hits = [hit for subl in [report.hits for report in all_reports] for hit in subl]
            all_hits = sorted(all_hits, key=attrgetter('score'), reverse=True)
            all_hits = sorted(all_hits, key=attrgetter('replicon_name', 'position'))

            db = RepliconDB(config)
            rep_info = db['AESU001c01a']
            with self.catch_log():
                clusters, multi_syst_genes = build_clusters(all_hits, [system], rep_info)
            cluster = clusters.clusters[0]

            # case 1
            system_occurence = SystemOccurence(system)
            system_occurence.fill_with_cluster(cluster)
            self.assertEqual(system_occurence.mandatory_genes['T9SS_sprT'], 1)
            self.assertEqual(len(system_occurence.valid_hits), 1)
            self.assertEqual(system_occurence.loci_positions, [(505, 505)])

            # case 2
            gene = system.get_gene('T9SS_sprT')
            system._mandatory_genes = []
            system._accessory_genes = [gene]
            system_occurence = SystemOccurence(system)
            system_occurence.fill_with_cluster(cluster)
            self.assertEqual(system_occurence.accessory_genes['T9SS_sprT'], 1)
            self.assertEqual(len(system_occurence.valid_hits), 1)

            # case 3
            system._accessory_genes = []
            system._forbidden_genes = [gene]
            system_occurence = SystemOccurence(system)
            system_occurence.fill_with_cluster(cluster)
            self.assertEqual(system_occurence.forbidden_genes['T9SS_sprT'], 1)
            self.assertEqual(len(system_occurence.valid_hits), 0)
            self.assertEqual(system_occurence.nb_cluster, 0)
            self.assertEqual(system_occurence.loci_positions, [])

            # case 4
            gene = Gene(self.macsy_test_env.cfg, 'tadZ', system, self.macsy_test_env.models_location)
            system._mandatory_genes = [gene]
            system._accessory_genes = []
            system._forbidden_genes = []
            system_occurence = SystemOccurence(system)
            system_occurence.exmandatory_genes['T9SS_sprT'] = 'tadZ'
            system_occurence.fill_with_cluster(cluster)
            self.assertEqual(system_occurence.mandatory_genes['tadZ'], 1)
            self.assertEqual(len(system_occurence.valid_hits), 1)

            # case 5
            gene = Gene(self.macsy_test_env.cfg, 'tadZ', system, self.macsy_test_env.models_location)
            system._mandatory_genes = []
            system._accessory_genes = [gene]
            system._forbidden_genes = []
            system_occurence = SystemOccurence(system)
            system_occurence.exaccessory_genes['T9SS_sprT'] = 'tadZ'
            system_occurence.fill_with_cluster(cluster)
            self.assertEqual(system_occurence.accessory_genes['tadZ'], 1)
            self.assertEqual(len(system_occurence.valid_hits), 1)

            # case 6
            gene = Gene(self.macsy_test_env.cfg, 'tadZ', system, self.macsy_test_env.models_location)
            system._mandatory_genes = []
            system._accessory_genes = []
            system._forbidden_genes = [gene]
            system_occurence = SystemOccurence(system)
            system_occurence.exforbidden_genes['T9SS_sprT'] = 'tadZ'
            system_occurence.fill_with_cluster(cluster)
            self.assertEqual(system_occurence.forbidden_genes['tadZ'], 1)
            self.assertEqual(len(system_occurence.valid_hits), 1)

            # case 7
            system._mandatory_genes = []
            system._accessory_genes = []
            system._forbidden_genes = []
            system_occurence = SystemOccurence(system)
            with self.catch_log():
                system_occurence.fill_with_cluster(cluster)
            self.assertEqual(len(system_occurence.forbidden_genes), 0)
            self.assertEqual(len(system_occurence.valid_hits), 0)
            self.assertEqual(system_occurence.nb_cluster, 1)
        finally:
            shutil.rmtree(out_dir)


    def test_fill_with_hits(self):

        # for this test, we use a specific configuration with a dedicated
        # working directory (i.e. we don't use the generic configuration
        # defined in setUp() method).
        out_dir = self.get_uniq_tmp_dir_name()

        args = argparse.Namespace()
        args.db_type = 'gembase'
        args.models_dir = self.find_data(os.path.join('data_set_1', 'models'))
        args.previous_run = self.find_data(os.path.join('data_set_1', 'complete_run_results'))
        args.res_search_dir = tempfile.gettempdir()
        args.log_level = 30
        args.out_dir = out_dir
        config = Config(MacsyDefaults(), args)

        if os.path.exists(config.working_dir()):
            shutil.rmtree(config.working_dir())
        os.makedirs(config.working_dir())

        try:
            idx = Indexes(config)
            idx._build_my_indexes()
            model_bank = self.macsy_test_env.model_bank
            gene_bank = self.macsy_test_env.gene_bank
            parser = DefinitionParser(config, model_bank, gene_bank)
            parser.parse(['set_1/T9SS'])
            system = model_bank['set_1/T9SS']
            genes = system.mandatory_genes + system.accessory_genes + system.forbidden_genes
            ex_genes = []
            for g in genes:
                if g.exchangeable:
                    h_s = g.get_homologs()
                    ex_genes += h_s
                    a_s = g.get_analogs()
                    ex_genes += a_s
            all_genes = (genes + ex_genes)
            all_reports = search_genes(all_genes, config)
            all_hits = [hit for subl in [report.hits for report in all_reports] for hit in subl]
            all_hits = sorted(all_hits, key=attrgetter('score'), reverse=True)
            all_hits = sorted(all_hits, key=attrgetter('replicon_name', 'position'))

            # case 1
            system_occurence = SystemOccurence(system)
            system_occurence.fill_with_hits(all_hits[:1], False)
            self.assertEqual(system_occurence.mandatory_genes['T9SS_sprT'], 1)
            self.assertEqual(len(system_occurence.valid_hits), 1)
            self.assertEqual(system_occurence._state, "no_decision")

            # case 2
            gene = system.get_gene('T9SS_sprT')
            system._mandatory_genes = []
            system._accessory_genes = [gene]
            system_occurence = SystemOccurence(system)
            with self.catch_log():
                system_occurence.fill_with_hits(all_hits[:2], False)
            self.assertEqual(system_occurence.accessory_genes['T9SS_sprT'], 1)
            self.assertEqual(len(system_occurence.valid_hits), 1)

            # case 3
            system._mandatory_genes = []
            system._accessory_genes = []
            system._forbidden_genes = [gene]
            system_occurence = SystemOccurence(system)
            with self.catch_log():
                system_occurence.fill_with_hits(all_hits, True)
            self.assertEqual(system_occurence.forbidden_genes['T9SS_sprT'], 1)
            self.assertEqual(len(system_occurence.valid_hits), 1)

            # case 4
            gene = Gene(self.macsy_test_env.cfg, 'tadZ', system, self.macsy_test_env.models_location)
            system._mandatory_genes = [gene]
            system._accessory_genes = []
            system._forbidden_genes = []
            system_occurence = SystemOccurence(system)
            system_occurence.exmandatory_genes['T9SS_sprT'] = 'tadZ'
            with self.catch_log():
                system_occurence.fill_with_hits(all_hits, True)
            self.assertEqual(system_occurence.mandatory_genes['tadZ'], 1)
            self.assertEqual(len(system_occurence.valid_hits), 1)

            # case 5
            gene = Gene(self.macsy_test_env.cfg, 'tadZ', system, self.macsy_test_env.models_location)
            system._mandatory_genes = []
            system._accessory_genes = [gene]
            system._forbidden_genes = []
            system_occurence = SystemOccurence(system)
            system_occurence.exaccessory_genes['T9SS_sprT'] = 'tadZ'
            with self.catch_log():
                system_occurence.fill_with_hits(all_hits, True)
            self.assertEqual(system_occurence.accessory_genes['tadZ'], 1)
            self.assertEqual(len(system_occurence.valid_hits), 1)

            # case 6
            gene = Gene(self.macsy_test_env.cfg, 'tadZ', system, self.macsy_test_env.models_location)
            system._mandatory_genes = []
            system._accessory_genes = []
            system._forbidden_genes = [gene]
            system_occurence = SystemOccurence(system)
            system_occurence.exforbidden_genes['T9SS_sprT'] = 'tadZ'
            with self.catch_log():
                system_occurence.fill_with_hits(all_hits, True)
            self.assertEqual(system_occurence.forbidden_genes['tadZ'], 1)
            self.assertEqual(len(system_occurence.valid_hits), 1)
        finally:
            shutil.rmtree(out_dir)
