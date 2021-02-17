#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2021  Institut Pasteur (Paris) and CNRS.           #
# See the COPYRIGHT file for details                                    #
#                                                                       #
# This file is part of MacSyFinder package.                             #
#                                                                       #
# MacSyFinder is free software: you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or     #
# (at your option) any later version.                                   #
#                                                                       #
# MacSyFinder is distributed in the hope that it will be useful,        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
# GNU General Public License for more details .                         #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with MacSyFinder (COPYING).                                     #
# If not, see <https://www.gnu.org/licenses/>.                          #
#########################################################################

import os
import argparse

from macsypy.hit import Hit, ValidHit, get_best_hits, HitWeight
from macsypy.config import Config, MacsyDefaults
from macsypy.gene import CoreGene, ModelGene, GeneStatus
from macsypy.profile import ProfileFactory
from macsypy.model import Model
from macsypy.registries import ModelLocation
from macsypy.error import MacsypyError
from tests import MacsyTest


class HitTest(MacsyTest):

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        self.cfg = Config(MacsyDefaults(), args)

        self.model_name = 'foo'
        self.model_location = ModelLocation(path=os.path.join(args.models_dir, self.model_name))

        # we need to reset the ProfileFactory
        # because it's a like a singleton
        # so other tests are influenced by ProfileFactory and it's configuration
        # for instance search_genes get profile without hmmer_exe
        self.profile_factory = ProfileFactory(self.cfg)


    def test_cmp(self):
        gene_name = "gspD"
        gene = CoreGene(self.model_location, gene_name, self.profile_factory)

        # compare hit with different id (comparison based on seq identifier)
        h0 = Hit(gene, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234),
                 float(779.2), float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = Hit(gene, "PSAE001c01_013980", 759, "PSAE001c01", 4146, float(3.7e-76),
                 float(255.8), float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)
        self.assertGreater(h1, h0)
        self.assertLess(h0, h1)
        # compare hit with different same id (comparison based on score)
        # score = 779.2
        h0 = Hit(gene, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234),
                 float(779.2), float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        # score = 255.8
        h1 = Hit(gene, "PSAE001c01_006940", 759, "PSAE001c01", 4146, float(3.7e-76),
                 float(255.8), float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)
        self.assertGreater(h0, h1)
        self.assertLess(h1, h0)


    def test_eq(self):
        gene_name = "gspD"
        gene = CoreGene(self.model_location, gene_name, self.profile_factory)

        h0 = Hit(gene, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                 float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = Hit(gene, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                 float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h2 = Hit(gene, "PSAE001c01_013980", 759, "PSAE001c01", 4146, float(3.7e-76), float(255.8),
                 float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)
        self.assertEqual(h0, h1)
        self.assertNotEqual(h0, h2)


    def test_str(self):
        gene_name = "gspD"
        gene = CoreGene(self.model_location, gene_name, self.profile_factory)

        hit_prop = {'id': "PSAE001c01_006940",
                    'hit_seq_len': 803,
                    'replicon_name': "PSAE001c01",
                    'position': 694,
                    'i_eval': float(1.2e-234),
                    'score': float(779.2),
                    'gene_name': gene.name,
                    'profil_coverage': float(1.0),
                    'sequence_coverage': float(638.000000),
                    'begin': 104,
                    'end': 741
                    }

        hit = Hit(gene, hit_prop['id'], hit_prop['hit_seq_len'], hit_prop['replicon_name'],
                  hit_prop['position'], hit_prop['i_eval'], hit_prop['score'],
                  hit_prop['profil_coverage'], hit_prop['sequence_coverage'], hit_prop['begin'],hit_prop['end'])
        s = "{id}\t{replicon_name}\t{position:d}\t{hit_seq_len:d}\t{gene_name}\t{i_eval:.3e}" \
            "\t{score:.3f}\t{profil_coverage:.3f}\t{sequence_coverage:.3f}\t{begin:d}\t{end:d}\n".format(**hit_prop)
        self.assertEqual(s, str(hit))


    def test_get_position(self):
        gene_name = "gspD"
        gene = CoreGene(self.model_location, gene_name, self.profile_factory)

        h0 = Hit(gene, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                 float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        self.assertEqual(h0.get_position(), 3450)


    def test_hash(self):
        gene_name = "gspD"
        gene = CoreGene(self.model_location, gene_name, self.profile_factory)

        h0 = Hit(gene, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                 float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = Hit(gene, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                 float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h2 = Hit(gene, "PSAE001c01_006941", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                 float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        self.assertTrue(isinstance(hash(h0), int))
        self.assertEqual(hash(h0), hash(h1))
        self.assertNotEqual(hash(h0), hash(h2))


class ValidHitTest(MacsyTest):

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        cfg = Config(MacsyDefaults(), args)

        model_name = 'foo'
        models_location = ModelLocation(path=os.path.join(args.models_dir, model_name))

        model = Model("foo/T2SS", 10)
        profile_factory = ProfileFactory(cfg)

        gene_name = "gspD"
        self.c_gene_gspd = CoreGene(models_location, gene_name, profile_factory)
        self.gene_gspd = ModelGene(self.c_gene_gspd, model)

        gene_name = "sctJ"
        self.c_gene_sctj = CoreGene(models_location, gene_name, profile_factory)
        self.gene_sctj = ModelGene(self.c_gene_sctj, model)

        model.add_mandatory_gene(self.gene_gspd)
        model.add_accessory_gene(self.gene_sctj)

        self.hit_1 = Hit(self.c_gene_gspd, "hit_1", 803, "replicon_id", 2, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.hit_2 = Hit(self.c_gene_sctj, "hit_2", 803, "replicon_id", 3, 1.0, 1.0, 1.0, 1.0, 10, 20)


    def test_init(self):
        v_hit_1 = ValidHit(self.hit_1, self.gene_gspd, GeneStatus.MANDATORY)
        self.assertEqual(v_hit_1.gene_ref, self.gene_gspd)
        self.assertEqual(v_hit_1.status, GeneStatus.MANDATORY)
        v_hit_2 = ValidHit(self.hit_2, self.gene_gspd, GeneStatus.ACCESSORY)
        self.assertEqual(v_hit_2.gene_ref, self.gene_gspd)
        self.assertEqual(v_hit_2.status, GeneStatus.ACCESSORY)

        with self.assertRaises(MacsypyError) as ctx:
            ValidHit(self.hit_1, self.c_gene_gspd, GeneStatus.MANDATORY)
        self.assertEqual(str(ctx.exception),
                         "The ValidHit 'gene_ref' argument must be a ModelGene not <class 'macsypy.gene.CoreGene'>.")

    def test_hash(self):
        self.assertEqual(hash(self.hit_1), hash(self.hit_1))
        self.assertNotEqual(hash(self.hit_1), hash(self.hit_2))

    def test_delegation(self):
        v_hit_1 = ValidHit(self.hit_1, self.gene_gspd, GeneStatus.MANDATORY)
        self.assertEqual(v_hit_1.get_position(), 2)

    def test_eq(self):
        v_hit_0 = ValidHit(self.hit_1, self.gene_sctj, GeneStatus.MANDATORY)
        v_hit_1 = ValidHit(self.hit_1, self.gene_sctj, GeneStatus.MANDATORY)
        v_hit_2 = ValidHit(self.hit_2, self.gene_gspd, GeneStatus.ACCESSORY)

        self.assertEqual(v_hit_0, v_hit_1)
        self.assertNotEqual(v_hit_0, v_hit_2)


class GetBestHitTest(MacsyTest):

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        self.cfg = Config(MacsyDefaults(), args)

        self.model_name = 'foo'
        self.models_location = ModelLocation(path=os.path.join(args.models_dir, self.model_name))

        # we need to reset the ProfileFactory
        # because it's a like a singleton
        # so other tests are influenced by ProfileFactory and it's configuration
        # for instance search_genes get profile without hmmer_exe
        self.profile_factory = ProfileFactory(self.cfg)


    def test_get_best_hits(self):
        model = Model("foo/T2SS", 10)
        gene_name = "gspD"
        c_gene_gspd = CoreGene(self.models_location, gene_name, self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model)

        #        gene, model, id,            hit_seq_len, replicon_name, position, i_eval,
        #        score,      profil_coverage,      sequence_coverage,     begin,end
        ######################
        # based on the score #
        ######################
        h0 = Hit(gene_gspd, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234),
                 10, float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = Hit(gene_gspd, "PSAE001c01_013980", 759, "PSAE001c01", 3450, float(3.7e-76),
                 11, float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)

        h = get_best_hits([h0, h1])
        self.assertEqual(h[0], h1)

        #######################
        # based on the i_eval #
        #######################
        h0 = Hit(gene_gspd, "PSAE001c01_006940", 803, "PSAE001c01", 3450, 10,
                 10, float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = Hit(gene_gspd, "PSAE001c01_013980", 759, "PSAE001c01", 3450, 11,
                 10, float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)

        h = get_best_hits([h0, h1], key='i_eval')
        self.assertEqual(h[0], h0)

        #################################
        # based on the profile_coverage #
        #################################
        h0 = Hit(gene_gspd, "PSAE001c01_006940", 803, "PSAE001c01", 3450, 10,
                 10, 10, (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = Hit(gene_gspd, "PSAE001c01_013980", 759, "PSAE001c01", 3450, 10,
                 10, 11, (736.0 - 105.0 + 1) / 759, 105, 736)

        h = get_best_hits([h0, h1], key='profile_coverage')
        self.assertEqual(h[0], h1)

        # bad criterion
        with self.assertRaises(MacsypyError) as ctx:
            get_best_hits([h0, h1], key='nimportnaoik')
        self.assertEqual('The criterion for Hits comparison nimportnaoik does not exist or is not available.\n'
                         'It must be either "score", "i_eval" or "profile_coverage".', str(ctx.exception))


class HitWeightTest(MacsyTest):

    def test_hit_weight_default(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        cfg = Config(MacsyDefaults(), args)
        hit_weight = HitWeight(**cfg.hit_weights())
        self.assertEqual(hit_weight.mandatory, 1)
        self.assertEqual(hit_weight.accessory, 0.5)
        self.assertEqual(hit_weight.itself, 1)
        self.assertEqual(hit_weight.exchangeable, 0.8)
        self.assertEqual(hit_weight.loner_multi_system, 0.7)


    def test_hit_weight_not_default(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.mandatory_weight = 2.0
        args.accessory_weight = 3.0
        args.neutral_weight = 4.0
        args.exchangeable_weight = 5.0
        args.itself_weight = 6.0
        args.loner_multi_system_weight = 12
        cfg = Config(MacsyDefaults(), args)
        hit_weight = HitWeight(**cfg.hit_weights())
        self.assertEqual(hit_weight.mandatory, 2.0)
        self.assertEqual(hit_weight.accessory, 3.0)
        self.assertEqual(hit_weight.neutral, 4.0)
        self.assertEqual(hit_weight.exchangeable, 5.0)
        self.assertEqual(hit_weight.itself, 6.0)
        self.assertEqual(hit_weight.loner_multi_system, 12.0)
