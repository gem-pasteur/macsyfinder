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


import argparse

from macsypy.hit import Hit, HitRegistry, ValidHit
from macsypy.config import Config, MacsyDefaults
from macsypy.gene import Gene, Homolog, Analog, ProfileFactory, GeneStatus
from macsypy.model import Model
from macsypy.registries import ModelRegistry
from macsypy.cluster import Cluster, RejectedClusters
from macsypy.system import PutativeSystem, match

from tests import MacsyTest


class SystemTest(MacsyTest):

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
        self.profile_factory = ProfileFactory()
        self.hit_registry = HitRegistry()

    def test_init(self):
        model = Model(self.cfg, "foo/T2SS", 10)
        # test if id is well incremented
        gene_gspd = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location)
        model.add_mandatory_gene(gene_gspd)
        gene_sctj = Gene(self.cfg, self.profile_factory, "sctJ", model, self.models_location)
        model.add_accessory_gene(gene_sctj)

        hit_1 = Hit(gene_gspd, model, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(gene_sctj, model, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        system_1 = PutativeSystem(model, [Cluster([v_hit_1, v_hit_2], model)])
        self.assertTrue(system_1.id.startswith('replicon_id_T2SS_'))

        system_2 = PutativeSystem(model, [Cluster([v_hit_1, v_hit_2], model)])
        self.assertEqual(int(system_2.id.split('_')[-1]), int(system_1.id.split('_')[-1]) + 1)

    def test_hits(self):
        model = Model(self.cfg, "foo/T2SS", 10)
        gene_gspd = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location)
        model.add_mandatory_gene(gene_gspd)
        gene_sctj = Gene(self.cfg, self.profile_factory, "sctJ", model, self.models_location)
        model.add_accessory_gene(gene_sctj)
        gene_sctn = Gene(self.cfg, self.profile_factory, "sctN", model, self.models_location)
        model.add_accessory_gene(gene_sctn)

        hit_1 = Hit(gene_gspd, model, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(gene_sctj, model, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(gene_sctn, model, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        system_1 = PutativeSystem(model, [Cluster([v_hit_1, v_hit_2], model),
                                               Cluster([v_hit_3], model)])

        self.assertEqual(system_1.hits, [v_hit_1, v_hit_2, v_hit_3])


    def test_match(self):
        model = Model(self.cfg, "foo/T2SS", 10)
        gene_sctn = Gene(self.cfg, self.profile_factory, "sctN", model, self.models_location)
        gene_sctn_flg = Homolog(
            Gene(self.cfg, self.profile_factory, "sctN_FLG", model, self.models_location),
            gene_sctn
        )
        gene_sctn.add_homolog(gene_sctn_flg)
        gene_sctj = Gene(self.cfg, self.profile_factory, "sctJ", model, self.models_location)
        gene_sctj_flg = Analog(
            Gene(self.cfg, self.profile_factory, "sctJ_FLG", model, self.models_location),
            gene_sctj
        )
        gene_sctj.add_analog(gene_sctj_flg)
        gene_gspd = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location)
        gene_abc = Gene(self.cfg, self.profile_factory, "abc", model, self.models_location)

        model.add_mandatory_gene(gene_sctn)
        model.add_mandatory_gene(gene_sctj)
        model.add_accessory_gene(gene_gspd)
        model.add_forbidden_gene(gene_abc)

        h_sctj = Hit(gene_sctj, model, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctj_flg = Hit(gene_sctj_flg, model, "hit_sctj_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = Hit(gene_sctn, model, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn_flg = Hit(gene_sctn_flg, model, "hit_sctn_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd = Hit(gene_gspd, model, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_abc = Hit(gene_gspd, model, "hit_abc", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        #####################
        # test single locus #
        #####################
        model._min_mandatory_genes_required = 2
        model._min_genes_required = 2

        # it lack one mandatory gene
        c1 = Cluster([h_sctj, h_gspd], model)
        res, multi_system_genes = match([c1], model, self.hit_registry)
        self.assertIsInstance(res, RejectedClusters)
        self.assertEqual(res.reason,
                         "The quorum of mandatory genes required (2) is not reached: 1\n"
                         "The quorum of genes required (2) is not reached: 1")

        model._min_mandatory_genes_required = 2
        model._min_genes_required = 1
        c1 = Cluster([h_sctj, h_sctn, h_gspd], model)
        res, multi_system_genes = match([c1], model, self.hit_registry)
        self.assertIsInstance(res, PutativeSystem)