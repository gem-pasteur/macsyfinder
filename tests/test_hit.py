#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2023  Institut Pasteur (Paris) and CNRS.           #
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

from macsypy.hit import CoreHit, ModelHit, Loner, MultiSystem, LonerMultiSystem, \
    get_best_hits, get_best_hit_4_func, HitWeight, sort_model_hits, compute_best_MSHit
from macsypy.config import Config, MacsyDefaults
from macsypy.gene import CoreGene, ModelGene, Exchangeable, GeneStatus
from macsypy.profile import ProfileFactory
from macsypy.model import Model
from macsypy.registries import ModelLocation
from macsypy.error import MacsypyError
from tests import MacsyTest


class CoreHitTest(MacsyTest):

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
        h0 = CoreHit(gene, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234),
                     float(779.2), float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = CoreHit(gene, "PSAE001c01_013980", 759, "PSAE001c01", 4146, float(3.7e-76),
                     float(255.8), float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)
        self.assertGreater(h1, h0)
        self.assertLess(h0, h1)
        # compare hit with different same id (comparison based on score)
        # score = 779.2
        h0 = CoreHit(gene, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234),
                     float(779.2), float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        # score = 255.8
        h1 = CoreHit(gene, "PSAE001c01_006940", 759, "PSAE001c01", 4146, float(3.7e-76),
                     float(255.8), float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)
        self.assertGreater(h0, h1)
        self.assertLess(h1, h0)


    def test_eq(self):
        gene_name = "gspD"
        gene = CoreGene(self.model_location, gene_name, self.profile_factory)

        h0 = CoreHit(gene, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                     float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = CoreHit(gene, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                     float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h2 = CoreHit(gene, "PSAE001c01_013980", 759, "PSAE001c01", 4146, float(3.7e-76), float(255.8),
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

        hit = CoreHit(gene, hit_prop['id'], hit_prop['hit_seq_len'], hit_prop['replicon_name'],
                      hit_prop['position'], hit_prop['i_eval'], hit_prop['score'],
                      hit_prop['profil_coverage'], hit_prop['sequence_coverage'], hit_prop['begin'], hit_prop['end'])
        s = "{id}\t{replicon_name}\t{position:d}\t{hit_seq_len:d}\t{gene_name}\t{i_eval:.3e}" \
            "\t{score:.3f}\t{profil_coverage:.3f}\t{sequence_coverage:.3f}\t{begin:d}\t{end:d}\n".format(**hit_prop)
        self.assertEqual(s, str(hit))


    def test_get_position(self):
        gene_name = "gspD"
        gene = CoreGene(self.model_location, gene_name, self.profile_factory)

        h0 = CoreHit(gene, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                     float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        self.assertEqual(h0.get_position(), 3450)


    def test_hash(self):
        gene_name = "gspD"
        gene = CoreGene(self.model_location, gene_name, self.profile_factory)

        h0 = CoreHit(gene, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                     float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = CoreHit(gene, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                     float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h2 = CoreHit(gene, "PSAE001c01_006941", 803, "PSAE001c01", 3450, float(1.2e-234), float(779.2),
                     float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        self.assertTrue(isinstance(hash(h0), int))
        self.assertEqual(hash(h0), hash(h1))
        self.assertNotEqual(hash(h0), hash(h2))


class ModelHitTest(MacsyTest):

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
        self.cg_gspd = CoreGene(models_location, gene_name, profile_factory)
        self.mg_gspd = ModelGene(self.cg_gspd, model, multi_system=True)

        gene_name = "sctJ"
        self.cg_sctj = CoreGene(models_location, gene_name, profile_factory)
        self.mg_sctj = ModelGene(self.cg_sctj, model, loner=True)

        gene_name = "abc"
        self.cg_abc = CoreGene(models_location, gene_name, profile_factory)
        self.mg_abc = ModelGene(self.cg_abc, model, multi_model=True)

        model.add_mandatory_gene(self.mg_gspd)
        model.add_accessory_gene(self.mg_sctj)
        model.add_accessory_gene(self.mg_abc)

        self.chit_1 = CoreHit(self.cg_gspd, "hit_1", 803, "replicon_id", 2, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.chit_2 = CoreHit(self.cg_sctj, "hit_2", 803, "replicon_id", 3, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.chit_3 = CoreHit(self.cg_abc, "hit_3", 803, "replicon_id", 4, 1.0, 1.0, 1.0, 1.0, 10, 20)


    def test_init(self):
        mhit_1 = ModelHit(self.chit_1, self.mg_gspd, GeneStatus.MANDATORY)
        self.assertEqual(mhit_1.gene_ref, self.mg_gspd)
        self.assertEqual(mhit_1.status, GeneStatus.MANDATORY)
        mhit_2 = ModelHit(self.chit_2, self.mg_gspd, GeneStatus.ACCESSORY)
        self.assertEqual(mhit_2.gene_ref, self.mg_gspd)
        self.assertEqual(mhit_2.status, GeneStatus.ACCESSORY)

        with self.assertRaises(MacsypyError) as ctx:
            ModelHit(mhit_1, self.cg_gspd, GeneStatus.MANDATORY)
        self.assertEqual(str(ctx.exception),
                         "The ModelHit 'hit' argument must be a CoreHit not <class 'macsypy.hit.ModelHit'>.")
        with self.assertRaises(MacsypyError) as ctx:
            ModelHit(self.chit_1, self.cg_gspd, GeneStatus.MANDATORY)
        self.assertEqual(str(ctx.exception),
                         "The ModelHit 'gene_ref' argument must be a ModelGene not <class 'macsypy.gene.CoreGene'>.")

    def test_hash(self):
        self.assertEqual(hash(self.chit_1), hash(self.chit_1))
        self.assertNotEqual(hash(self.chit_1), hash(self.chit_2))

    def test_delegation(self):
        mhit_1 = ModelHit(self.chit_1, self.mg_gspd, GeneStatus.MANDATORY)
        self.assertEqual(mhit_1.get_position(), 2)
        with self.assertRaises(AttributeError) as ctx:
            mhit_1.nimportnaoik()
        self.assertEqual(str(ctx.exception),
                         "'ModelHit' object has no attribute 'nimportnaoik'")

    def test_eq(self):
        mhit_0 = ModelHit(self.chit_1, self.mg_sctj, GeneStatus.MANDATORY)
        mhit_1 = ModelHit(self.chit_1, self.mg_sctj, GeneStatus.MANDATORY)
        mhit_2 = ModelHit(self.chit_2, self.mg_gspd, GeneStatus.ACCESSORY)

        self.assertEqual(mhit_0, mhit_1)
        self.assertNotEqual(mhit_0, mhit_2)

    def test_multi_system(self):
        mhit_1 = ModelHit(self.chit_1, self.mg_sctj, GeneStatus.MANDATORY)
        mhit_2 = ModelHit(self.chit_2, self.mg_gspd, GeneStatus.MANDATORY)
        self.assertTrue(mhit_2.multi_system)
        self.assertFalse(mhit_1.multi_system)

    def test_multi_model(self):
        mhit_1 = ModelHit(self.chit_1, self.mg_sctj, GeneStatus.MANDATORY)
        mhit_2 = ModelHit(self.chit_3, self.mg_abc, GeneStatus.ACCESSORY)
        self.assertFalse(mhit_1.multi_model)
        self.assertTrue(mhit_2.multi_model)


    def test_loner(self):
        mhit_2 = ModelHit(self.chit_2, self.mg_gspd, GeneStatus.MANDATORY)
        self.assertFalse(mhit_2.loner)


class LonerTest(MacsyTest):

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
        self.cg_gspd = CoreGene(models_location, gene_name, profile_factory)
        self.mg_gspd = ModelGene(self.cg_gspd, model, loner=True)

        gene_name = "sctJ"
        self.cg_sctj = CoreGene(models_location, gene_name, profile_factory)
        self.mg_sctj = ModelGene(self.cg_sctj, model)

        model.add_mandatory_gene(self.mg_gspd)
        model.add_accessory_gene(self.mg_sctj)

        self.chit_1 = CoreHit(self.cg_gspd, "hit_1", 803, "replicon_id", 2, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.chit_2 = CoreHit(self.cg_sctj, "hit_2", 803, "replicon_id", 3, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.chit_3 = CoreHit(self.cg_gspd, "hit_3", 803, "replicon_id", 10, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.mhit_1 = ModelHit(self.chit_1, self.mg_gspd, GeneStatus.MANDATORY)
        self.mhit_2 = ModelHit(self.chit_2, self.mg_sctj, GeneStatus.ACCESSORY)
        self.mhit_3 = ModelHit(self.chit_3, self.mg_gspd, GeneStatus.MANDATORY)

    def test_init(self):
        # create loner from a CoreHit
        l1 = Loner(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY)
        self.assertEqual(l1.gene_ref, self.mg_gspd)
        self.assertEqual(l1.status, GeneStatus.MANDATORY)
        # create loner from a CoreHit with counterpart
        _ = Loner(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY, counterpart=[self.mhit_3])

        # try to create Loner from cCoreHit but without gene_ref nor gene_status
        with self.assertRaises(MacsypyError) as ctx:
            Loner(self.chit_1, gene_ref=self.mg_gspd)
        self.assertEqual(str(ctx.exception),
                         "Cannot Create a Loner hit from CoreHit (gspD, 2) "
                         "without specifying 'gene_ref' and 'gene_status'")
        with self.assertRaises(MacsypyError) as ctx:
            Loner(self.chit_1, gene_status=GeneStatus.MANDATORY)
        self.assertEqual(str(ctx.exception),
                         "Cannot Create a Loner hit from CoreHit (gspD, 2) "
                         "without specifying 'gene_ref' and 'gene_status'")

        # create loner from a ModelHit
        mh1 = ModelHit(self.chit_1, self.mg_gspd, GeneStatus.MANDATORY)
        l3 = Loner(mh1)
        self.assertEqual(l3.gene_ref, self.mg_gspd)
        self.assertEqual(l3.status, GeneStatus.MANDATORY)

        with self.catch_log(log_name='macsypy') as log:
            with self.assertRaises(MacsypyError) as ctx:
                Loner(self.chit_2, gene_ref=self.mg_sctj, gene_status=GeneStatus.ACCESSORY)
        self.assertEqual(str(ctx.exception),
                         "hit_2 cannot be a loner gene_ref 'sctJ' not tag as loner")

        with self.assertRaises(MacsypyError) as ctx:
            Loner(self.chit_1, self.cg_gspd, GeneStatus.MANDATORY)
        self.assertEqual(str(ctx.exception),
                         "The Loner 'gene_ref' argument must be a ModelGene not <class 'macsypy.gene.CoreGene'>.")


    def test_counterpart(self):
        l1 = Loner(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY)
        self.assertSetEqual(l1.counterpart, set())
        l2 = Loner(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY, counterpart=[self.mhit_3])
        self.assertSetEqual(l2.counterpart, set([self.mhit_3]))

    def test_len(self):
        l1 = Loner(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY)
        self.assertEqual(len(l1), 1)
        l2 = Loner(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY, counterpart=[self.mhit_3])
        self.assertEqual(len(l2), 2)


    def test_loner(self):
        l1 = Loner(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY)
        self.assertTrue(l1.loner)

    def test_multi_system(self):
        l1 = Loner(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY)
        self.assertFalse(l1.multi_system)

    def test_str(self):
        l2 = Loner(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY, counterpart=[self.mhit_3])
        self.assertEqual(str(l2),
                         '\t'.join(['hit_1', 'replicon_id', '2', '803', 'gspD', '1.000e+00',
                                    '1.000', '1.000', '1.000', '10', '20', 'hit_3'])
                         )


class MultiSystemTest(MacsyTest):

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
        self.cg_gspd = CoreGene(models_location, gene_name, profile_factory)
        self.mg_gspd = ModelGene(self.cg_gspd, model, multi_system=True)

        gene_name = "sctJ"
        self.cg_sctj = CoreGene(models_location, gene_name, profile_factory)
        self.mg_sctj = ModelGene(self.cg_sctj, model)

        model.add_mandatory_gene(self.mg_gspd)
        model.add_accessory_gene(self.mg_sctj)

        #              CoreHit(gene, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #                       profile_coverage, sequence_coverage, begin_match, end_match
        #                                                                   pos     score
        self.chit_1 = CoreHit(self.cg_gspd, "hit_1", 803, "replicon_id", 2, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.chit_2 = CoreHit(self.cg_sctj, "hit_2", 803, "replicon_id", 3, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.chit_3 = CoreHit(self.cg_gspd, "hit_3", 803, "replicon_id", 10, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.mhit_1 = ModelHit(self.chit_1, self.mg_gspd, GeneStatus.MANDATORY)
        self.mhit_2 = ModelHit(self.chit_2, self.mg_sctj, GeneStatus.ACCESSORY)
        self.mhit_3 = ModelHit(self.chit_3, self.mg_gspd, GeneStatus.MANDATORY)

    def test_init(self):
        # create loner from a CoreHit
        ms1 = MultiSystem(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY)
        self.assertEqual(ms1.gene_ref, self.mg_gspd)
        self.assertEqual(ms1.status, GeneStatus.MANDATORY)
        # create loner from a CoreHit with counterpart
        # test the creation, the content of counterpart is test in test_counterpart
        _ = MultiSystem(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY, counterpart=[self.mhit_3])

        # try to create MS from cCoreHit but without gene_ref nor gene_status
        with self.assertRaises(MacsypyError) as ctx:
            MultiSystem(self.chit_1, gene_ref=self.mg_gspd)
        self.assertEqual(str(ctx.exception),
                         "Cannot Create a MultiSystem hit from CoreHit (gspD, 2) "
                         "without specifying 'gene_ref' and 'gene_status'")
        with self.assertRaises(MacsypyError) as ctx:
            MultiSystem(self.chit_1, gene_status=GeneStatus.MANDATORY)
        self.assertEqual(str(ctx.exception),
                         "Cannot Create a MultiSystem hit from CoreHit (gspD, 2) "
                         "without specifying 'gene_ref' and 'gene_status'")

        # create MultiSystem from a ModelHit
        mh1 = ModelHit(self.chit_1, gene_ref=self.mg_gspd,  gene_status=GeneStatus.MANDATORY)
        ms3 = MultiSystem(mh1)
        self.assertEqual(ms3.gene_ref, self.mg_gspd)
        self.assertEqual(ms3.status, GeneStatus.MANDATORY)

        with self.catch_log(log_name='macsypy'):
            with self.assertRaises(MacsypyError) as ctx:
                MultiSystem(self.chit_2, gene_ref=self.mg_sctj, gene_status=GeneStatus.ACCESSORY)
        self.assertEqual(str(ctx.exception),
                         "hit_2 cannot be a multi systems, gene_ref 'sctJ' not tag as multi_system")

        with self.assertRaises(MacsypyError) as ctx:
            MultiSystem(self.chit_1, self.cg_gspd, GeneStatus.MANDATORY)
        self.assertEqual(str(ctx.exception),
                         "The MultiSystem 'gene_ref' argument must be a ModelGene not <class 'macsypy.gene.CoreGene'>.")


    def test_counterpart(self):
        ms1 = MultiSystem(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY)
        self.assertSetEqual(ms1.counterpart, set())
        ms2 = MultiSystem(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY, counterpart=[self.mhit_3])
        self.assertSetEqual(ms2.counterpart, set([self.mhit_3]))

    def test_len(self):
        ms1 = MultiSystem(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY)
        self.assertEqual(len(ms1), 1)
        ms2 = MultiSystem(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY, counterpart=[self.mhit_3])
        self.assertEqual(len(ms2), 2)


    def test_loner(self):
        ms1 = MultiSystem(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY)
        self.assertFalse(ms1.loner)

    def test_multi_system(self):
        ms1 = MultiSystem(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY)
        self.assertTrue(ms1.multi_system)

    def test_str(self):
        ms2 = MultiSystem(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY, counterpart=[self.mhit_3])
        self.assertEqual(str(ms2),
                         '\t'.join(['hit_1', 'replicon_id', '2', '803', 'gspD', '1.000e+00',
                                    '1.000', '1.000', '1.000', '10', '20', 'hit_3'])
                         )


class LonerMultiSystemTest(MacsyTest):

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
        self.cg_gspd = CoreGene(models_location, gene_name, profile_factory)
        self.mg_gspd = ModelGene(self.cg_gspd, model, loner=True, multi_system=True)

        gene_name = "sctJ"
        self.cg_sctj = CoreGene(models_location, gene_name, profile_factory)
        self.mg_sctj = ModelGene(self.cg_sctj, model)

        model.add_mandatory_gene(self.mg_gspd)
        model.add_accessory_gene(self.mg_sctj)

        self.chit_1 = CoreHit(self.cg_gspd, "hit_1", 803, "replicon_id", 2, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.chit_2 = CoreHit(self.cg_sctj, "hit_2", 803, "replicon_id", 3, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.chit_3 = CoreHit(self.cg_gspd, "hit_3", 803, "replicon_id", 10, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.chit_4 = CoreHit(self.cg_gspd, "hit_4", 803, "replicon_id", 20, 1.0, 1.0, 1.0, 1.0, 10, 20)
        self.mhit_1 = ModelHit(self.chit_1, self.mg_gspd, GeneStatus.MANDATORY)
        self.mhit_2 = ModelHit(self.chit_2, self.mg_sctj, GeneStatus.ACCESSORY)
        self.mhit_3 = ModelHit(self.chit_3, self.mg_gspd, GeneStatus.MANDATORY)
        self.mhit_4 = ModelHit(self.chit_4, self.mg_gspd, GeneStatus.MANDATORY)

    def test_init(self):
        # create loner from a CoreHit
        lms1 = LonerMultiSystem(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY)
        self.assertEqual(lms1.gene_ref, self.mg_gspd)
        self.assertEqual(lms1.status, GeneStatus.MANDATORY)
        # create loner from a CoreHit with counterpart
        # test the creation, the content of counterpart is test in test_counterpart
        _ = LonerMultiSystem(self.chit_1, gene_ref=self.mg_gspd,
                             gene_status=GeneStatus.MANDATORY,
                             counterpart=[self.mhit_3])

        # try to create MS from CoreHit but without gene_ref nor gene_status
        with self.assertRaises(MacsypyError) as ctx:
            LonerMultiSystem(self.chit_1, gene_ref=self.mg_gspd)
        self.assertEqual(str(ctx.exception),
                         "Cannot Create a LonerMultiSystem hit from CoreHit (gspD, 2) "
                         "without specifying 'gene_ref' and 'gene_status'")
        with self.assertRaises(MacsypyError) as ctx:
            LonerMultiSystem(self.chit_1, gene_status=GeneStatus.MANDATORY)
        self.assertEqual(str(ctx.exception),
                         "Cannot Create a LonerMultiSystem hit from CoreHit (gspD, 2) "
                         "without specifying 'gene_ref' and 'gene_status'")

        # create MultiSystem from a ModelHit
        mh1 = ModelHit(self.chit_1, gene_ref=self.mg_gspd,  gene_status=GeneStatus.MANDATORY)
        lms = LonerMultiSystem(mh1)
        self.assertEqual(lms.gene_ref, self.mg_gspd)
        self.assertEqual(lms.status, GeneStatus.MANDATORY)

        with self.catch_log(log_name='macsypy'):
            with self.assertRaises(MacsypyError) as ctx:
                LonerMultiSystem(self.chit_2, gene_ref=self.mg_sctj, gene_status=GeneStatus.ACCESSORY)
        self.assertEqual(str(ctx.exception),
                         "hit_2 cannot be a multi systems, gene_ref 'sctJ' not tag as multi_system")

        with self.assertRaises(MacsypyError) as ctx:
            LonerMultiSystem(self.chit_1, self.cg_gspd, GeneStatus.MANDATORY)
        self.assertEqual(str(ctx.exception),
                         "The LonerMultiSystem 'gene_ref' argument must be a ModelGene"
                         " not <class 'macsypy.gene.CoreGene'>.")

        # create from a MultiSystem
        mh1 = MultiSystem(self.chit_1,
                          gene_ref=self.mg_gspd,
                          gene_status=GeneStatus.MANDATORY,
                          counterpart=[self.mhit_3, self.mhit_4])
        lms = LonerMultiSystem(mh1)
        self.assertEqual(lms.gene_ref, self.mg_gspd)
        self.assertEqual(lms.status, GeneStatus.MANDATORY)
        self.assertSetEqual(set(lms.counterpart), set([self.mhit_3, self.mhit_4]))

    def test_loner(self):
        lms1 = LonerMultiSystem(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY)
        self.assertTrue(lms1.loner)

    def test_multi_system(self):
        lms1 = LonerMultiSystem(self.chit_1, gene_ref=self.mg_gspd, gene_status=GeneStatus.MANDATORY)
        self.assertTrue(lms1.multi_system)

    def test_str(self):
        lms2 = LonerMultiSystem(self.chit_1, gene_ref=self.mg_gspd,
                                gene_status=GeneStatus.MANDATORY,
                                counterpart=[self.mhit_3])
        self.assertEqual(str(lms2),
                         '\t'.join(['hit_1', 'replicon_id', '2', '803', 'gspD', '1.000e+00',
                                    '1.000', '1.000', '1.000', '10', '20', 'hit_3'])
                         )


    def test_counterpart(self):
        mh1 = MultiSystem(self.chit_1,
                          gene_ref=self.mg_gspd,
                          gene_status=GeneStatus.MANDATORY)

        mh1.counterpart = [self.mhit_3, self.mhit_4]
        self.assertSetEqual(mh1.counterpart,
                            {self.mhit_3, self.mhit_4})

        with self.assertRaises(MacsypyError) as ctx:
            with self.catch_log() as log:
                mh1.counterpart = [self.mhit_3, self.mhit_2]
        self.assertEqual(str(ctx.exception),
                         "Try to set counterpart for hit 'gspD' with non compatible hits: ['gspD', 'sctJ']")


class GetBestHitTest(MacsyTest):

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        cfg = Config(MacsyDefaults(), args)

        model_name = 'foo'
        self.models_location = ModelLocation(path=os.path.join(args.models_dir, model_name))

        # we need to reset the ProfileFactory
        # because it's a like a singleton
        # so other tests are influenced by ProfileFactory and it's configuration
        # for instance search_genes get profile without hmmer_exe
        profile_factory = ProfileFactory(cfg)

        model = Model(model_name, 10)
        self.profile_factory = ProfileFactory(cfg)

        gene_name = "gspD"
        c_gene_gspd = CoreGene(self.models_location, gene_name, self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model, multi_system=True)

        gene_name = "sctJ"
        c_gene_sctj = CoreGene(self.models_location, gene_name, self.profile_factory)
        gene_sctj = ModelGene(c_gene_sctj, model, multi_system=True)

        gene_name = "sctN"
        c_gene_sctn = CoreGene(self.models_location, gene_name, self.profile_factory)
        gene_sctn = Exchangeable(c_gene_sctn, gene_sctj)
        gene_sctj.add_exchangeable(gene_sctn)

        model.add_mandatory_gene(gene_gspd)
        model.add_accessory_gene(gene_sctj)

        #        CoreHit(gene, hit_id, hit_seq_length, replicon_name, position, i_eval, score,
        #                       profile_coverage, sequence_coverage, begin_match, end_match
        #                                                        pos      score
        chit_1 = CoreHit(c_gene_gspd, "hit_1", 803, "replicon_id", 2, 1.0, 1.0, 1.0, 1.0, 10, 20)
        chit_2 = CoreHit(c_gene_sctj, "hit_2", 803, "replicon_id", 3, 1.0, 1.0, 1.0, 1.0, 10, 20)
        chit_3 = CoreHit(c_gene_gspd, "hit_3", 803, "replicon_id", 10, 1.0, 3.0, 1.0, 1.0, 10, 20)
        chit_4 = CoreHit(c_gene_sctn, "hit_4", 803, "replicon_id", 14, 1.0, 4.0, 1.0, 1.0, 10, 20)
        chit_5 = CoreHit(c_gene_gspd, "hit_5", 803, "replicon_id", 20, 1.0, 2.0, 1.0, 1.0, 10, 20)

        self.mhit_1 = ModelHit(chit_1, gene_gspd, GeneStatus.MANDATORY)
        self.mhit_2 = ModelHit(chit_2, gene_sctj, GeneStatus.ACCESSORY)
        self.mhit_3 = ModelHit(chit_3, gene_gspd, GeneStatus.MANDATORY)
        self.mhit_4 = ModelHit(chit_4, gene_sctn, GeneStatus.ACCESSORY)
        self.mhit_5 = ModelHit(chit_5, gene_gspd, GeneStatus.MANDATORY)

        self.ms_1 = MultiSystem(chit_1, gene_ref=gene_gspd, gene_status=GeneStatus.MANDATORY)
        self.ms_2 = MultiSystem(chit_2, gene_ref=gene_sctj, gene_status=GeneStatus.ACCESSORY)
        self.ms_3 = MultiSystem(chit_3, gene_ref=gene_gspd, gene_status=GeneStatus.MANDATORY)
        self.ms_4 = MultiSystem(chit_4, gene_ref=gene_sctn, gene_status=GeneStatus.ACCESSORY)
        self.ms_5 = MultiSystem(chit_5, gene_ref=gene_gspd, gene_status=GeneStatus.MANDATORY)

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
        h0 = CoreHit(gene_gspd, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234),
                     10, float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = CoreHit(gene_gspd, "PSAE001c01_013980", 759, "PSAE001c01", 3450, float(3.7e-76),
                     11, float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)

        h = get_best_hits([h0, h1])
        self.assertEqual(h[0], h1)

        #######################
        # based on the i_eval #
        #######################
        h0 = CoreHit(gene_gspd, "PSAE001c01_006940", 803, "PSAE001c01", 3450, 10,
                     10, float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = CoreHit(gene_gspd, "PSAE001c01_013980", 759, "PSAE001c01", 3450, 11,
                     10, float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)

        h = get_best_hits([h0, h1], key='i_eval')
        self.assertEqual(h[0], h0)

        #################################
        # based on the profile_coverage #
        #################################
        h0 = CoreHit(gene_gspd, "PSAE001c01_006940", 803, "PSAE001c01", 3450, 10,
                     10, 10, (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = CoreHit(gene_gspd, "PSAE001c01_013980", 759, "PSAE001c01", 3450, 10,
                     10, 11, (736.0 - 105.0 + 1) / 759, 105, 736)

        h = get_best_hits([h0, h1], key='profile_coverage')
        self.assertEqual(h[0], h1)

        # bad criterion
        with self.assertRaises(MacsypyError) as ctx:
            get_best_hits([h0, h1], key='nimportnaoik')
        self.assertEqual('The criterion for Hits comparison nimportnaoik does not exist or is not available.\n'
                         'It must be either "score", "i_eval" or "profile_coverage".', str(ctx.exception))


    def test_get_best_hits_4_func(self):
        model = Model("foo/T2SS", 10)
        gene_name = "gspD"
        c_gene_gspd = CoreGene(self.models_location, gene_name, self.profile_factory)
        gene_gspd = ModelGene(c_gene_gspd, model, loner=True)

        #        gene, model, id,            hit_seq_len, replicon_name, position, i_eval,
        #        score,      profil_coverage,      sequence_coverage,     begin,end
        ######################
        # based on the score #
        ######################
        h0 = CoreHit(gene_gspd, "PSAE001c01_006940", 803, "PSAE001c01", 3450, float(1.2e-234),
                     10, float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = CoreHit(gene_gspd, "PSAE001c01_013980", 759, "PSAE001c01", 3450, float(3.7e-76),
                     11, float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)
        m0 = ModelHit(h0, gene_gspd, GeneStatus.ACCESSORY)
        m1 = ModelHit(h1, gene_gspd, GeneStatus.ACCESSORY)
        l0 = Loner(h0, gene_ref=gene_gspd, gene_status=GeneStatus.ACCESSORY, counterpart=[m1])
        l1 = Loner(h1, gene_ref=gene_gspd, gene_status=GeneStatus.ACCESSORY, counterpart=[m0])
        l = get_best_hit_4_func(gene_name, [l0, l1])
        self.assertEqual(l, l1)

        #######################
        # based on the i_eval #
        #######################
        h0 = CoreHit(gene_gspd, "PSAE001c01_006940", 803, "PSAE001c01", 3450, 10,
                     10, float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = CoreHit(gene_gspd, "PSAE001c01_013980", 759, "PSAE001c01", 3450, 11,
                     10, float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736)
        m0 = ModelHit(h0, gene_gspd, GeneStatus.ACCESSORY)
        m1 = ModelHit(h1, gene_gspd, GeneStatus.ACCESSORY)
        l0 = Loner(h0, gene_ref=gene_gspd, gene_status=GeneStatus.ACCESSORY, counterpart=[m1])
        l1 = Loner(h1, gene_ref=gene_gspd, gene_status=GeneStatus.ACCESSORY, counterpart=[m0])

        l = get_best_hit_4_func(gene_name, [l0, l1], key='i_eval')
        self.assertEqual(l, l0)

        #################################
        # based on the profile_coverage #
        #################################
        h0 = CoreHit(gene_gspd, "PSAE001c01_006940", 803, "PSAE001c01", 3450, 10,
                     10, 10, (741.0 - 104.0 + 1) / 803, 104, 741)
        h1 = CoreHit(gene_gspd, "PSAE001c01_013980", 759, "PSAE001c01", 3450, 10,
                     10, 11, (736.0 - 105.0 + 1) / 759, 105, 736)
        m0 = ModelHit(h0, gene_gspd, GeneStatus.ACCESSORY)
        m1 = ModelHit(h1, gene_gspd, GeneStatus.ACCESSORY)
        l0 = Loner(h0, gene_ref=gene_gspd, gene_status=GeneStatus.ACCESSORY, counterpart=[m1])
        l1 = Loner(h1, gene_ref=gene_gspd, gene_status=GeneStatus.ACCESSORY, counterpart=[m0])

        l = get_best_hit_4_func(gene_name, [l0, l1], key='profile_coverage')
        self.assertEqual(l, l1)

        # bad criterion
        with self.assertRaises(MacsypyError) as ctx:
            get_best_hits([l0, l1], key='nimportnaoik')
        self.assertEqual('The criterion for Hits comparison nimportnaoik does not exist or is not available.\n'
                         'It must be either "score", "i_eval" or "profile_coverage".', str(ctx.exception))


    def test_sort_multisystem_hits(self):
        # scores                        gspd 1.0     sctj 1.0     gspd 3.0  sctn 4.0 Ex sctj gspd 2.0
        sorted_hits = sort_model_hits([self.mhit_1, self.mhit_2, self.mhit_3, self.mhit_4, self.mhit_5])
        expected = {'gspD': [self.mhit_1, self.mhit_3, self.mhit_5],
                    'sctJ': [self.mhit_2, self.mhit_4]}
        self.assertDictEqual(sorted_hits,
                             expected)

        sorted_hits = sort_model_hits([])
        self.assertDictEqual(sorted_hits,
                             {})


    def test_compute_best_MSHit(self):
        # scores                       gspd 1.0    sctj 1.0   gspd 3.0 sctn 4.0 Ex sctj  gspd 2.0
        sorted_hits = sort_model_hits([self.ms_1, self.ms_2, self.ms_3, self.ms_4, self.ms_5])
        best_hits = compute_best_MSHit(sorted_hits)
        # The order does not matter
        self.assertSetEqual({self.ms_3, self.ms_2}, set(best_hits))
        self.assertSetEqual(self.ms_3.counterpart,
                            {self.ms_1, self.ms_5})
        self.assertSetEqual(self.ms_2.counterpart,
                            {self.ms_4})


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
        self.assertEqual(hit_weight.out_of_cluster, 0.7)


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
        args.out_of_cluster_weight = 12
        cfg = Config(MacsyDefaults(), args)
        hit_weight = HitWeight(**cfg.hit_weights())
        self.assertEqual(hit_weight.mandatory, 2.0)
        self.assertEqual(hit_weight.accessory, 3.0)
        self.assertEqual(hit_weight.neutral, 4.0)
        self.assertEqual(hit_weight.exchangeable, 5.0)
        self.assertEqual(hit_weight.itself, 6.0)
        self.assertEqual(hit_weight.out_of_cluster, 12.0)
