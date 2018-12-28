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

from macsypy.report import Hit
from macsypy.config import Config, MacsyDefaults
from macsypy.gene import Gene
import macsypy.gene
from macsypy.system import System
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

        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.res_search_dir = tempfile.gettempdir()
        args.log_level = 30
        args.log_file = log_file
        self.cfg = Config(MacsyDefaults(), args)

        models_registry = ModelRegistry(self.cfg)
        self.model_name = 'foo'
        self.models_location = models_registry[self.model_name]

        # we need to reset the ProfileFactory
        # because it's a like a singleton
        # so other tests are influenced by ProfileFactory and it's configuration
        # for instance search_genes get profile without hmmer_exe
        macsypy.gene.profile_factory = macsypy.gene.ProfileFactory()
        macsypy.gene.profile_factory._profiles = {}

  
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
        macsypy.gene.profile_factory = macsypy.gene.ProfileFactory()
        macsypy.gene.profile_factory._profiles = {}


    def test_cmp(self):
        system = System(self.cfg, "foo/T2SS", 10)
        gene = Gene(self.cfg, "gspD", system, self.models_location)
        # compare hit with different id (comparison based on seq identifier)
        h0 = Hit(gene, system, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234),
                 float(779.2), float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = Hit(gene, system, "PSAE001c01_013980", 759, "PSAE001c01", 4146, float(3.7e-76),
                 float(255.8), float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)
        self.assertGreater(h1, h0)
        self.assertLess(h0, h1)
        # compare hit with different same id (comparison based on score)
        # score = 779.2
        h0 = Hit(gene, system, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234),
                 float(779.2), float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        # score = 255.8
        h1 = Hit(gene, system, "PSAE001c01_006940", 759, "PSAE001c01", 4146, float(3.7e-76),
                 float(255.8), float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)
        self.assertGreater(h0, h1)
        self.assertLess(h1, h0)
        # compare non homolgous genes
        gene_non_h = Gene(self.cfg, "abc", system, self.models_location)
        h2 = Hit(gene_non_h, system, "PSAE001c01_006940", 759, "PSAE001c01", 4146, float(3.7e-76),
                 float(255.8), float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)
        self.assertGreater(h0, h2)
        self.assertLess(h2, h0)

    def test_eq(self):
        system = System(self.cfg, "foo/T2SS", 10)
        gene = Gene(self.cfg, "gspD", system, self.models_location)
        h0 = Hit(gene, system, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                 float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = Hit(gene, system, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                 float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h2 = Hit(gene, system, "PSAE001c01_013980", 759, "PSAE001c01", 4146, float(3.7e-76), float(255.8),
                 float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)
        self.assertEqual(h0, h1)
        self.assertNotEqual(h0, h2)
        
    def test_str(self):
        system = System(self.cfg, "foo/T2SS", 10)
        gene = Gene(self.cfg, "gspD", system, self.models_location)
        hit_prop = {'id': "PSAE001c01_006940",
                    'hit_seq_len': 803,
                    'replicon_name': "PSAE001c01",
                    'position': 694,
                    'i_eval': float(1.2e-234),
                    'score': float(779.2),
                    'gene_name': gene.name,
                    'system_name': system.name,
                    'profil_coverage': float(1.0),
                    'sequence_coverage': float(638.000000),
                    'begin': 104,
                    'end': 741
                    }
        
        hit = Hit(gene, system, hit_prop['id'], hit_prop['hit_seq_len'], hit_prop['replicon_name'],
                  hit_prop['position'] , hit_prop['i_eval'], hit_prop['score'],
                  hit_prop['profil_coverage'], hit_prop['sequence_coverage'], hit_prop['begin'],hit_prop['end'])
        s = "{id}\t{replicon_name}\t{position:d}\t{hit_seq_len:d}\t{gene_name}\t{system_name}\t{i_eval:.3e}" \
            "\t{score:.3f}\t{profil_coverage:.3f}\t{sequence_coverage:.3f}\t{begin:d}\t{end:d}\n".format(**hit_prop)
        self.assertEqual(s, str(hit))

    def test_get_syst_inter_gene_max_space(self):
        system = System(self.cfg, "foo/T2SS", 10)
        gene = Gene(self.cfg, "gspD", system, self.models_location)
        h0 = Hit(gene, system, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                 float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        self.assertEqual(h0.get_syst_inter_gene_max_space(), 10)
