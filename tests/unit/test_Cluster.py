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
import unittest
import shutil
import tempfile
import logging
from operator import attrgetter
from macsypy.config import Config
from macsypy.system import System, system_bank
from macsypy.gene import gene_bank
from macsypy.report import Hit
from macsypy.system_parser import SystemParser
from macsypy.search_systems import Cluster
from macsypy.search_genes import search_genes
from macsypy.registries import ModelRegistry
from macsypy.database import Indexes
from macsypy.macsypy_error import SystemDetectionError
from tests import MacsyTest, md5sum
from tests.unit import MacsyTestEnv


class Test(MacsyTest):

    def setUp(self):
        self.macsy_test_env = MacsyTestEnv()

    def tearDown(self):
        self.macsy_test_env = None

    def test_len(self):
        self.macsy_test_env.load("env_001")

        system = System(self.macsy_test_env.cfg, 'foo', 10)
        cluster = Cluster(system)
        li = [None] * 12
        hit = Hit(*li)
        cluster.hits = [hit]
        self.assertEqual(len(cluster), 1)

        self.macsy_test_env.unload("env_001")

    def test_putative_system(self):
        self.macsy_test_env.load("env_001")

        system_name = 'set_1/T9SS'
        system = System(self.macsy_test_env.cfg, system_name, 10)
        cluster = Cluster(system)
        cluster._putative_system = system_name
        self.assertEqual(cluster.putative_system, system_name)

        self.macsy_test_env.unload("env_001")

    def test_compatible_systems(self):
        self.macsy_test_env.load("env_001")

        system = System(self.macsy_test_env.cfg, 'set_1/T9SS', 10)
        cluster = Cluster(system)
        compatible_system_name = 'set_1/T2SS'
        cluster._compatible_systems.append(compatible_system_name)
        self.assertListEqual(cluster.compatible_systems, [compatible_system_name])

        self.macsy_test_env.unload("env_001")

    def test_state(self):
        self.macsy_test_env.load("env_001")

        system = System(self.macsy_test_env.cfg, 'foo', 4)
        cluster = Cluster(system)
        state = cluster.state
        self.assertEqual(state, '')

        self.macsy_test_env.unload("env_001")

    @unittest.skip("")
    def test_add(self):
        out_dir = MacsyTest.get_uniq_tmp_dir_name()

        # for this test, we use a specific configuration with a dedicated
        # working directory (i.e. we don't use the generic configuration
        # defined in setUp() method).
        # the config create a directory corresponding to out_dir option
        config = Config(hmmer_exe="hmmsearch",
                        out_dir=out_dir,
                        db_type="gembase",
                        previous_run="tests/data/data_set_1/complete_run_results",
                        e_value_res=1,
                        i_evalue_sel=0.5,
                        res_search_suffix=".search_hmm.out",
                        profile_suffix=".hmm",
                        res_extract_suffix="",
                        log_level=30,
                        models_dir="tests/data/data_set_1/models",
                        log_file=os.devnull)

        idx = Indexes(config)
        idx._build_my_indexes()

        parser = SystemParser(config, system_bank, gene_bank)
        parser.parse(['set_1/T9SS'])

        system = system_bank['set_1/T9SS']

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

        # debug
        # print [h.gene.name for h in all_hits]

        h1, h2, h3 = all_hits[:3]

        cluster = Cluster([system])

        cluster.add(h1)
        self.assertEqual(cluster.begin, 505)
        self.assertEqual(cluster.end, 505)
        self.assertEqual(cluster.replicon_name, 'AESU001c01a')
        self.assertEqual(len(cluster.hits), 1)

        cluster.add(h2)
        self.assertEqual(cluster.begin, 505)
        self.assertEqual(cluster.end, 773)
        self.assertEqual(cluster.replicon_name, 'AESU001c01a')
        self.assertEqual(len(cluster.hits), 2)

        h3.replicon_name = 'bar'
        with self.assertRaises(SystemDetectionError) as context:
            cluster.add(h3)
        self.assertEqual(context.exception.message,
                         "Attempting to gather in a cluster hits from different replicons ! ")

        shutil.rmtree(out_dir)

    @unittest.skip("")
    def test_save(self):

        out_dir = MacsyTest.get_uniq_tmp_dir_name()

        # for this test, we use a specific configuration with a dedicated
        # working directory (i.e. we don't use the generic configuration
        # defined in setUp() method).
        config = Config(hmmer_exe="hmmsearch",
                        out_dir=out_dir,
                        db_type="gembase",
                        previous_run="tests/data/data_set_1/complete_run_results",
                        e_value_res=1,
                        i_evalue_sel=0.5,
                        res_search_suffix=".search_hmm.out",
                        profile_suffix=".hmm",
                        res_extract_suffix="",
                        log_level=30,
                        models_dir="tests/data/data_set_1/models",
                        log_file=os.devnull)

        idx = Indexes(config)
        idx._build_my_indexes()

        parser = SystemParser(config, system_bank, gene_bank)
        parser.parse(['set_1/T9SS'])

        system_1 = system_bank['set_1/T9SS']

        def get_hits(system):
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

            return all_hits

        all_hits_1 = get_hits(system_1)

        # debug
        # print [h.gene.name for h in all_hits_1]

        # debug
        # print cluster.systems
        # print cluster._putative_system
        # print cluster._compatible_systems
        # print cluster._state

        # debug
        # fqn = 'set_1/T4P'
        # fqn = 'set_1/CONJ'
        # fqn = 'set_1/Tad'
        # fqn = 'set_1/Flagellum'
        # fqn = 'set_1/T2SS'

        # test case 1

        cluster = Cluster([system_1])
        for h in all_hits_1:
            cluster.add(h)
        cluster.save()
        self.assertEqual(cluster.systems, {'set_1/T9SS': 17})
        self.assertEqual(cluster._putative_system, 'set_1/T9SS')
        self.assertEqual(cluster._compatible_systems, ['set_1/T9SS'])
        self.assertEqual(cluster._state, 'clear')

        # test case 2

        cluster = Cluster([system_1])
        cluster.add(all_hits_1[0])
        cluster.hits[0].gene._loner = False
        cluster.save()
        self.assertEqual(cluster.systems, {'set_1/T9SS': 1})
        self.assertEqual(cluster._putative_system, 'set_1/T9SS')
        self.assertEqual(cluster._compatible_systems, [])
        self.assertEqual(cluster._state, 'ineligible')

        # test case 3

        fqn = 'set_1/T3SS'
        parser.parse([fqn])
        system_2 = system_bank[fqn]
        all_hits_2 = get_hits(system_2)
        cluster = Cluster([system_1, system_2])
        for h in all_hits_1 + all_hits_2:
            cluster.add(h)
        cluster.save()
        self.assertEqual(cluster.systems, {'set_1/T9SS': 17, 'set_1/T3SS': 2})
        self.assertEqual(cluster._putative_system, 'set_1/T3SS')
        self.assertEqual(cluster._compatible_systems, [])
        self.assertEqual(cluster._state, 'ambiguous')

        # test case 4

        """
        FIXME

        Warning: some parts of the 'save()' method are not tested.  To test
        those parts, the try_system() sub-func must return "clear", which seems
        not possible to achieve with the test dataset used.

        DRAFT
        cluster = Cluster([system_1, system_2])
        for h in all_hits_1 + all_hits_2:
            cluster.add(h)
        cluster.save()
        """

        try:
            shutil.rmtree(out_dir)
        except:
            pass

    def test_str(self):
        self.macsy_test_env.load("env_002")
        buffer_ = str(self.macsy_test_env.cluster)
        self.assertEqual(md5sum(str_=buffer_), '752f66b832fee7ec71fecdaef52fd820')
        self.macsy_test_env.unload("env_002")
