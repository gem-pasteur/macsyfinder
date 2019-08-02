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
import json

from macsypy.hit import Hit, HitRegistry, ValidHit
from macsypy.config import Config, MacsyDefaults
from macsypy.gene import Gene, Homolog, Analog, ProfileFactory, GeneStatus
from macsypy.model import Model
from macsypy.registries import ModelRegistry
from macsypy.cluster import Cluster, RejectedClusters
from macsypy.system import System, match

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
        model = Model("foo/T2SS", 10)
        # test if id is well incremented
        gene_gspd = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location)
        model.add_mandatory_gene(gene_gspd)
        gene_sctj = Gene(self.cfg, self.profile_factory, "sctJ", model, self.models_location)
        model.add_accessory_gene(gene_sctj)

        hit_1 = Hit(gene_gspd, model, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(gene_sctj, model, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        system_1 = System(model, [Cluster([v_hit_1, v_hit_2], model)])
        self.assertTrue(system_1.id.startswith('replicon_id_T2SS_'))

        system_2 = System(model, [Cluster([v_hit_1, v_hit_2], model)])
        self.assertEqual(int(system_2.id.split('_')[-1]), int(system_1.id.split('_')[-1]) + 1)

    def test_hits(self):
        model = Model("foo/T2SS", 10)
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
        system_1 = System(model, [Cluster([v_hit_1, v_hit_2], model),
                                  Cluster([v_hit_3], model)])

        self.assertEqual(system_1.hits, [v_hit_1, v_hit_2, v_hit_3])


    def test_multi_loci(self):
        model = Model("foo/T2SS", 10)
        gene_gspd = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location)
        model.add_mandatory_gene(gene_gspd)
        gene_sctj = Gene(self.cfg, self.profile_factory, "sctJ", model, self.models_location)
        model.add_accessory_gene(gene_sctj)
        gene_sctn = Gene(self.cfg, self.profile_factory, "sctN", model, self.models_location, loner=True)
        model.add_accessory_gene(gene_sctn)

        hit_1 = Hit(gene_gspd, model, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(gene_sctj, model, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(gene_sctn, model, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        c1 = Cluster([v_hit_1, v_hit_2], model)
        c2 = Cluster([v_hit_1, v_hit_3], model)
        sys_single_locus = System(model, [c1])
        self.assertFalse(sys_single_locus.multi_loci)
        sys_multi_loci = System(model, [c1, c2])
        self.assertTrue(sys_multi_loci.multi_loci)
        c1 = Cluster([v_hit_1, v_hit_2], model)
        c3 = Cluster([v_hit_3], model)
        sys_single_locus_plus_loner = System(model, [c1, c3])
        self.assertFalse(sys_single_locus_plus_loner.multi_loci)

    def test_loci(self):
        model = Model("foo/T2SS", 10)
        gene_gspd = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location)
        model.add_mandatory_gene(gene_gspd)
        gene_sctj = Gene(self.cfg, self.profile_factory, "sctJ", model, self.models_location)
        model.add_accessory_gene(gene_sctj)
        gene_sctn = Gene(self.cfg, self.profile_factory, "sctN", model, self.models_location, loner=True)
        model.add_accessory_gene(gene_sctn)

        hit_1 = Hit(gene_gspd, model, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(gene_sctj, model, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(gene_sctn, model, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        c1 = Cluster([v_hit_1, v_hit_2], model)
        c2 = Cluster([v_hit_1, v_hit_3], model)
        sys_single_locus = System(model, [c1])
        self.assertEqual(sys_single_locus.loci, 1)
        sys_multi_loci = System(model, [c1, c2])
        self.assertEqual(sys_multi_loci.loci, 2)
        c1 = Cluster([v_hit_1, v_hit_2], model)
        c3 = Cluster([v_hit_3], model)
        sys_single_locus_plus_loner = System(model, [c1, c3])
        self.assertEqual(sys_single_locus_plus_loner.loci, 1)

    def test_wholeness(self):
        model = Model("foo/T2SS", 10)
        gene_gspd = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location)
        model.add_mandatory_gene(gene_gspd)
        gene_sctj = Gene(self.cfg, self.profile_factory, "sctJ", model, self.models_location)
        model.add_accessory_gene(gene_sctj)
        gene_sctn = Gene(self.cfg, self.profile_factory, "sctN", model, self.models_location, loner=True)
        model.add_accessory_gene(gene_sctn)

        hit_1 = Hit(gene_gspd, model, "hit_1", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(gene_sctj, model, "hit_2", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(gene_sctn, model, "hit_3", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        c1 = Cluster([v_hit_1, v_hit_2], model)
        c2 = Cluster([v_hit_1, v_hit_3], model)
        s = System(model, [c1])
        self.assertEqual(s.wholeness, 2 / 3)
        s = System(model, [c1, c2])
        self.assertEqual(s.wholeness, 3 / 3)


    def test_score(self):
        model = Model("foo/T2SS", 10)
        gene_gspd = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location)
        model.add_mandatory_gene(gene_gspd)
        gene_tadZ = Gene(self.cfg, self.profile_factory, "tadZ", model, self.models_location)
        model.add_mandatory_gene(gene_tadZ)

        gene_sctj = Gene(self.cfg, self.profile_factory, "sctJ", model, self.models_location, exchangeable=True)
        gene_sctJ_FLG = Gene(self.cfg, self.profile_factory, 'sctJ_FLG', model, self.models_location)
        analog = Analog(gene_sctJ_FLG, gene_sctj)
        gene_sctj.add_analog(analog)
        model.add_accessory_gene(gene_sctj)

        gene_sctn = Gene(self.cfg, self.profile_factory, "sctN", model, self.models_location, loner=True)
        model.add_accessory_gene(gene_sctn)

        hit_1 = Hit(gene_gspd, model, "hit_1", 10, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(gene_tadZ, model, "hit_2", 20, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_tadZ, GeneStatus.MANDATORY)

        hit_3 = Hit(gene_sctj, model, "hit_3", 30, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctj, GeneStatus.ACCESSORY)
        hit_3_bis = Hit(gene_sctJ_FLG, model, "hit_3_bis", 30, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3_bis = ValidHit(hit_3_bis, gene_sctj, GeneStatus.ACCESSORY)

        hit_4 = Hit(gene_sctn, model, "hit_4", 40, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_4 = ValidHit(hit_4, gene_sctn, GeneStatus.ACCESSORY)
        # system with 1 cluster 2 mandatory, 1 accessory
        s = System(model, [Cluster([v_hit_1, v_hit_2, v_hit_3, v_hit_4], model)])
        self.assertEqual(s.score, 3)
        # system with 2 clusters 2 mandatory, 2 accessory
        s = System(model, [Cluster([v_hit_1, v_hit_2], model),
                           Cluster([v_hit_3, v_hit_4], model)])
        self.assertEqual(s.score, 3)
        # system with 2 mandatory, 1 accessory , 1 accessory present 3 times
        s = System(model, [Cluster([v_hit_1, v_hit_2, v_hit_3, v_hit_4, v_hit_1, v_hit_4], model)])
        self.assertEqual(s.score, 0)
        # system with 1 cluster 1 mandatory 1 accessory
        s = System(model, [Cluster([v_hit_1, v_hit_3], model)])
        self.assertEqual(s.score, 1.5)
        # system with 1 cluster with 2 mandatory 1 accessory 1 accessory analog
        s = System(model, [Cluster([v_hit_1, v_hit_2, v_hit_3_bis, v_hit_4], model)])
        self.assertEqual(s.score, 2.875)


    def test_to_json(self):
        model = Model("foo/T2SS", 10)
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
        c1 = Cluster([v_hit_1, v_hit_2], model)
        c2 = Cluster([v_hit_3], model)
        sys_multi_loci = System(model, [c1, c2])

        rec_json = sys_multi_loci.to_json()
        exp_json = {'id': sys_multi_loci.id,
                    'model': 'foo/T2SS',
                    'loci_nb': 2,
                    'replicon_name': 'replicon_id',
                    'clusters': [['gspD', 'sctJ'], ['sctN']],
                    'gene_composition':
                        {'mandatory': {'gspD': ['gspD']},
                         'accessory': {'sctJ': ['sctJ'], 'sctN': ['sctN']}
                         }
                    }
        self.assertDictEqual(json.loads(rec_json), exp_json)


    def test_str(self):
        model = Model("foo/T2SS", 10)
        gene_gspd = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location)
        model.add_mandatory_gene(gene_gspd)
        gene_sctj = Gene(self.cfg, self.profile_factory, "sctJ", model, self.models_location)
        model.add_accessory_gene(gene_sctj)
        gene_sctn = Gene(self.cfg, self.profile_factory, "sctN", model, self.models_location)
        model.add_accessory_gene(gene_sctn)

        hit_1 = Hit(gene_gspd, model, "hit_1", 803, "replicon_id", 10, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_1 = ValidHit(hit_1, gene_gspd, GeneStatus.MANDATORY)
        hit_2 = Hit(gene_sctj, model, "hit_2", 803, "replicon_id", 20, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_2 = ValidHit(hit_2, gene_sctj, GeneStatus.ACCESSORY)
        hit_3 = Hit(gene_sctn, model, "hit_3", 803, "replicon_id", 30, 1.0, 1.0, 1.0, 1.0, 10, 20)
        v_hit_3 = ValidHit(hit_3, gene_sctn, GeneStatus.ACCESSORY)
        c1 = Cluster([v_hit_1, v_hit_2], model)
        c2 = Cluster([v_hit_3], model)
        sys_multi_loci = System(model, [c1, c2])
        sys_str = """system id = {}
model = foo/T2SS 
replicon = replicon_id
clusters = [('gspD', 10), ('sctJ', 20)], [('sctN', 30)]
occ = 1
wholeness = 1.000
loci nb = 1
score = 2.000

mandatory genes:
\t- gspD: 1 (gspD)

accessory genes:
\t- sctJ: 1 (sctJ)
\t- sctN: 1 (sctN)
""".format(sys_multi_loci.id)
        self.assertEqual(sys_str, str(sys_multi_loci))


    def test_match(self):
        model = Model("foo/T2SS", 10)
        gene_sctn = Gene(self.cfg, self.profile_factory, "sctN", model, self.models_location, exchangeable=True)
        gene_sctn_flg = Homolog(
            Gene(self.cfg, self.profile_factory, "sctN_FLG", model, self.models_location),
            gene_sctn
        )
        gene_sctn.add_homolog(gene_sctn_flg)
        gene_sctj = Gene(self.cfg, self.profile_factory, "sctJ", model, self.models_location, exchangeable=True)
        gene_sctj_flg = Analog(
            Gene(self.cfg, self.profile_factory, "sctJ_FLG", model, self.models_location),
            gene_sctj
        )
        gene_sctj.add_analog(gene_sctj_flg)
        gene_gspd = Gene(self.cfg, self.profile_factory, "gspD", model, self.models_location, exchangeable=True)
        gene_gspd_an = Analog(
            Gene(self.cfg, self.profile_factory, "flgB", model, self.models_location),
            gene_gspd
        )
        gene_gspd.add_analog(gene_gspd_an)
        gene_abc = Gene(self.cfg, self.profile_factory, "abc", model, self.models_location, exchangeable=True)
        gene_abc_ho = Homolog(
            Gene(self.cfg, self.profile_factory, "tadZ", model, self.models_location),
            gene_abc
        )
        gene_abc.add_homolog(gene_abc_ho)
        model.add_mandatory_gene(gene_sctn)
        model.add_mandatory_gene(gene_sctj)
        model.add_accessory_gene(gene_gspd)
        model.add_forbidden_gene(gene_abc)

        h_sctj = Hit(gene_sctj, model, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctj_flg = Hit(gene_sctj_flg, model, "hit_sctj_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn = Hit(gene_sctn, model, "hit_sctn", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_sctn_flg = Hit(gene_sctn_flg, model, "hit_sctn_flg", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd = Hit(gene_gspd, model, "hit_gspd", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_gspd_an = Hit(gene_gspd_an, model, "hit_gspd_an", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_abc = Hit(gene_abc, model, "hit_abc", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
        h_abc_ho = Hit(gene_abc_ho, model, "hit_abc_ho", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)

        #####################
        # test single locus #
        #####################

        # it lack one mandatory gene
        model._min_mandatory_genes_required = 2
        model._min_genes_required = 3
        c1 = Cluster([h_sctj, h_gspd], model)
        res, multi_system_genes = match([c1], model, self.hit_registry)
        self.assertIsInstance(res, RejectedClusters)
        self.assertEqual(res.reason,
                         "The quorum of mandatory genes required (2) is not reached: 1\n"
                         "The quorum of genes required (3) is not reached: 2")

        # all quorum are reached
        model._min_mandatory_genes_required = 2
        model._min_genes_required = 1
        c1 = Cluster([h_sctj, h_sctn, h_gspd], model)
        res, _ = match([c1], model, self.hit_registry)
        self.assertIsInstance(res, System)

        # with one mandatory analog
        model._min_mandatory_genes_required = 2
        model._min_genes_required = 1
        c1 = Cluster([h_sctj_flg, h_sctn, h_gspd], model)
        res, _ = match([c1], model, self.hit_registry)
        self.assertIsInstance(res, System)

        # with one accessory analog
        model._min_mandatory_genes_required = 2
        model._min_genes_required = 1
        c1 = Cluster([h_sctj, h_sctn, h_gspd_an], model)
        res, _ = match([c1], model, self.hit_registry)
        self.assertIsInstance(res, System)

        # the min_gene_required quorum is not reached
        model._min_mandatory_genes_required = 2
        model._min_genes_required = 4
        c1 = Cluster([h_sctj, h_sctn_flg, h_gspd], model)
        res, _ = match([c1], model, self.hit_registry)
        self.assertIsInstance(res, RejectedClusters)
        self.assertEqual(res.reason,
                         "The quorum of genes required (4) is not reached: 3")

        # the cluster contain a forbidden gene
        model._min_mandatory_genes_required = 2
        model._min_genes_required = 1
        c1 = Cluster([h_sctj, h_sctn, h_gspd, h_abc], model)
        res, _ = match([c1], model, self.hit_registry)
        self.assertIsInstance(res, RejectedClusters)
        self.assertEqual(res.reason, "There is 1 forbidden genes occurrence(s): abc")

        # the cluster contain a forbidden gene homolog
        model._min_mandatory_genes_required = 2
        model._min_genes_required = 1
        c1 = Cluster([h_sctj, h_sctn, h_gspd, h_abc_ho], model)
        res, _ = match([c1], model, self.hit_registry)
        self.assertIsInstance(res, RejectedClusters)
        self.assertEqual(res.reason, "There is 1 forbidden genes occurrence(s): tadZ")

        #####################
        # test multi loci   #
        #####################
        model._min_mandatory_genes_required = 2
        model._min_genes_required = 1
        c1 = Cluster([h_sctj, h_sctn], model)
        c2 = Cluster([h_gspd], model)
        res, _ = match([c1, c2], model, self.hit_registry)
        self.assertIsInstance(res, System)

        # with one analog an one homolog
        model._min_mandatory_genes_required = 2
        model._min_genes_required = 1
        c1 = Cluster([h_sctj_flg, h_sctn_flg], model)
        c2 = Cluster([h_gspd], model)
        res, _ = match([c1, c2], model, self.hit_registry)
        self.assertIsInstance(res, System)

        # with one analog an one homolog and one forbidden in 3 clusters
        model._min_mandatory_genes_required = 2
        model._min_genes_required = 1
        c1 = Cluster([h_sctj_flg, h_sctn_flg], model)
        c2 = Cluster([h_gspd], model)
        c3 = Cluster([h_abc], model)
        res, _ = match([c1, c2, c3], model, self.hit_registry)
        self.assertEqual(res.reason, "There is 1 forbidden genes occurrence(s): abc")
