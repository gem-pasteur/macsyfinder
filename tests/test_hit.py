# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Neron                                         #
# Copyright (c) 2014-2019  Institut Pasteur (Paris) and CNRS.                  #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################


import argparse

from macsypy.hit import Hit, ValidHit, get_best_hits, hit_weight
from macsypy.config import Config, MacsyDefaults
from macsypy.gene import ProfileFactory, Gene, GeneStatus
from macsypy.model import Model
from macsypy.registries import ModelRegistry
from macsypy.error import MacsypyError
from tests import MacsyTest


class HitTest(MacsyTest):

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        self.cfg = Config(MacsyDefaults(), args)

        models_registry = ModelRegistry(self.cfg)
        self.model_name = 'foo'
        self.models_location = models_registry[self.model_name]

        # we need to reset the ProfileFactory
        # because it's a like a singleton
        # so other tests are influenced by ProfileFactory and it's configuration
        # for instance search_genes get profile without hmmer_exe
        self.profile_factory = ProfileFactory(self.cfg)


    def test_cmp(self):
        model = Model("foo/T2SS", 10)
        gene = Gene(self.profile_factory, "gspD", model, self.models_location)
        # compare hit with different id (comparison based on seq identifier)
        h0 = Hit(gene, model, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234),
                 float(779.2), float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = Hit(gene, model, "PSAE001c01_013980", 759, "PSAE001c01", 4146, float(3.7e-76),
                 float(255.8), float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)
        self.assertGreater(h1, h0)
        self.assertLess(h0, h1)
        # compare hit with different same id (comparison based on score)
        score = 779.2
        h0 = Hit(gene, model, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234),
                 float(779.2), float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        # score = 255.8
        h1 = Hit(gene, model, "PSAE001c01_006940", 759, "PSAE001c01", 4146, float(3.7e-76),
                 float(255.8), float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)
        self.assertGreater(h0, h1)
        self.assertLess(h1, h0)
        # compare non homolgous genes
        gene_non_h = Gene(self.profile_factory, "abc", model, self.models_location)

        h2 = Hit(gene_non_h, model, "PSAE001c01_006940", 759, "PSAE001c01", 4146, float(3.7e-76),
                 float(255.8), float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)
        with self.catch_log():
            self.assertGreater(h0, h2)
            self.assertLess(h2, h0)

    def test_eq(self):
        model = Model("foo/T2SS", 10)
        gene = Gene(self.profile_factory, "gspD", model, self.models_location)
        h0 = Hit(gene, model, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                 float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = Hit(gene, model, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                 float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h2 = Hit(gene, model, "PSAE001c01_013980", 759, "PSAE001c01", 4146, float(3.7e-76), float(255.8),
                 float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)
        self.assertEqual(h0, h1)
        self.assertNotEqual(h0, h2)
        
    def test_str(self):
        model = Model("foo/T2SS", 10)
        gene = Gene(self.profile_factory, "gspD", model, self.models_location)
        hit_prop = {'id': "PSAE001c01_006940",
                    'hit_seq_len': 803,
                    'replicon_name': "PSAE001c01",
                    'position': 694,
                    'i_eval': float(1.2e-234),
                    'score': float(779.2),
                    'gene_name': gene.name,
                    'model_name': model.name,
                    'profil_coverage': float(1.0),
                    'sequence_coverage': float(638.000000),
                    'begin': 104,
                    'end': 741
                    }
        
        hit = Hit(gene, model, hit_prop['id'], hit_prop['hit_seq_len'], hit_prop['replicon_name'],
                  hit_prop['position'], hit_prop['i_eval'], hit_prop['score'],
                  hit_prop['profil_coverage'], hit_prop['sequence_coverage'], hit_prop['begin'],hit_prop['end'])
        s = "{id}\t{replicon_name}\t{position:d}\t{hit_seq_len:d}\t{gene_name}\t{model_name}\t{i_eval:.3e}" \
            "\t{score:.3f}\t{profil_coverage:.3f}\t{sequence_coverage:.3f}\t{begin:d}\t{end:d}\n".format(**hit_prop)
        self.assertEqual(s, str(hit))

    def test_get_syst_inter_gene_max_space(self):
        model = Model("foo/T2SS", 10)
        gene = Gene(self.profile_factory, "gspD", model, self.models_location)
        h0 = Hit(gene, model, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                 float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        self.assertEqual(h0.get_syst_inter_gene_max_space(), 10)


class ValidHitTest(MacsyTest):

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        cfg = Config(MacsyDefaults(), args)

        models_registry = ModelRegistry(cfg)
        model_name = 'foo'
        models_location = models_registry[model_name]

        model = Model("foo/T2SS", 10)
        profile_factory = ProfileFactory(cfg)
        self.gene_gspd = Gene(profile_factory, "gspD", model, models_location)
        model.add_mandatory_gene(self.gene_gspd)
        self.gene_sctj = Gene(profile_factory, "sctJ", model, models_location)
        model.add_accessory_gene(self.gene_sctj)

        self.hit_1 = Hit(self.gene_gspd, model, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.hit_2 = Hit(self.gene_sctj, model, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)


    def test_init(self):
        v_hit_1 = ValidHit(self.hit_1, self.gene_gspd, GeneStatus.MANDATORY)
        self.assertEqual(v_hit_1.gene_ref, self.gene_gspd)
        self.assertEqual(v_hit_1.status, GeneStatus.MANDATORY)
        v_hit_2 = ValidHit(self.hit_2, self.gene_gspd, GeneStatus.ACCESSORY)
        self.assertEqual(v_hit_2.gene_ref, self.gene_gspd)
        self.assertEqual(v_hit_2.status, GeneStatus.ACCESSORY)

    def test_hash(self):
        self.assertEqual(hash(self.hit_1), id(self.hit_1))


class GetBestHitTest(MacsyTest):

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        self.cfg = Config(MacsyDefaults(), args)

        models_registry = ModelRegistry(self.cfg)
        self.model_name = 'foo'
        self.models_location = models_registry[self.model_name]

        # we need to reset the ProfileFactory
        # because it's a like a singleton
        # so other tests are influenced by ProfileFactory and it's configuration
        # for instance search_genes get profile without hmmer_exe
        self.profile_factory = ProfileFactory(self.cfg)

    def test_get_best_hits(self):
        model = Model("foo/T2SS", 10)
        gene = Gene(self.profile_factory, "gspD", model, self.models_location)


        #        gene, model, id,            hit_seq_len, replicon_name, position, i_eval,
        #        score,      profil_coverage,      sequence_coverage,     begin,end
        ######################
        # based on the score #
        ######################
        h0 = Hit(gene, model, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234),
                 10, float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = Hit(gene, model, "PSAE001c01_013980", 759, "PSAE001c01", 3450, float(3.7e-76),
                 11, float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)

        h = get_best_hits([h0, h1])
        self.assertEqual(h[0], h1)

        #######################
        # based on the i_eval #
        #######################
        h0 = Hit(gene, model, "PSAE001c01_006940", 803, "PSAE001c01", 3450, 10,
                 10, float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = Hit(gene, model, "PSAE001c01_013980", 759, "PSAE001c01", 3450, 11,
                 10, float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)

        h = get_best_hits([h0, h1], key='i_eval')
        self.assertEqual(h[0], h0)

        #################################
        # based on the profile_coverage #
        #################################
        h0 = Hit(gene, model, "PSAE001c01_006940", 803, "PSAE001c01", 3450, 10,
                 10, 10, (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = Hit(gene, model, "PSAE001c01_013980", 759, "PSAE001c01", 3450, 10,
                 10, 11, (736.0 - 105.0 + 1) / 759, 105, 736)

        h = get_best_hits([h0, h1], key='profile_coverage')
        self.assertEqual(h[0], h1)

        # bad criterion
        with self.assertRaises(MacsypyError) as ctx:
            get_best_hits([h0, h1], key='nimportnaoik')
        self.assertEqual('The criterion for Hits comparison nimportnaoik does not exist or is not available.\n'
                         'It must be either "score", "i_eval" or "profile_coverage".', str(ctx.exception))


class HitWeightTest(MacsyTest):

    def test_hit_weight(self):
        self.assertEqual(hit_weight.mandatory, 1)
        self.assertEqual(hit_weight.accessory, 0.5)
        self.assertEqual(hit_weight.hitself, 1)
        self.assertEqual(hit_weight.homolog, 0.75)
        self.assertEqual(hit_weight.analog, 0.75)