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
from operator import attrgetter
from time import strftime
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
from tests import MacsyTest


class Test(MacsyTest):

    def setUp(self):
        logger = logging.getLogger()
        logger.manager.loggerDict.clear()

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
                          log_file=log_file)

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
        logger = logging.getLogger()
        logger.manager.loggerDict.clear()
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass

        system_bank._system_bank = {}
        gene_bank._genes_bank = {}

    def test_len(self):
        system = System(self.cfg, 'foo', 10)
        cluster = Cluster(system)
        li = [None] * 12
        hit = Hit(*li)
        cluster.hits = [hit]
        self.assertEqual(len(cluster), 1)

    def test_putative_system(self):
        system_name = 'set_1/T9SS'
        system = System(self.cfg, system_name, 10)
        cluster = Cluster(system)
        cluster._putative_system = system_name
        self.assertEqual(cluster.putative_system, system_name)

    def test_compatible_systems(self):
        system = System(self.cfg, 'set_1/T9SS', 10)
        cluster = Cluster(system)
        compatible_system_name = 'set_1/T2SS'
        cluster._compatible_systems.append(compatible_system_name)
        self.assertEqual(cluster.compatible_systems, [compatible_system_name])

    def test_state(self):
        system = System(self.cfg, 'foo', 4)
        cluster = Cluster(system)
        state = cluster.state
        self.assertEqual(state, '')

    def test_add(self):
        out_dir = "/tmp/macsyfinder-test_fill_with_cluster-" + strftime("%Y%m%d_%H-%M-%S")

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

        system = System(self.cfg, 'foo', 4)
        cluster = Cluster(system)

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

        try:
            shutil.rmtree(out_dir)
        except:
            pass
