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
import shutil
import tempfile
import argparse

from macsypy.config import Config, MacsyDefaults
from macsypy.model import Model
from macsypy.gene import CoreGene, ModelGene, Exchangeable
from macsypy.profile import ProfileFactory
from macsypy.hit import Hit
from macsypy.registries import ModelLocation
from tests import MacsyTest


class TestModel(MacsyTest):

    def setUp(self):
        self.args = argparse.Namespace()
        self.args.sequence_db = self.find_data("base", "test_1.fasta")
        self.args.db_type = 'gembase'
        self.args.models_dir = self.find_data('models')
        self.args.res_search_dir = tempfile.gettempdir()
        self.args.log_level = 30
        self.args.out_dir = os.path.join(self.args.res_search_dir,
                                         'test_macsyfinder_Model')
        if os.path.exists(self.args.out_dir):
            shutil.rmtree(self.args.out_dir)
        os.mkdir(self.args.out_dir)

        self.cfg = Config(MacsyDefaults(), self.args)
        self.model_name = 'foo'
        self.model_location = ModelLocation(path=os.path.join(self.args.models_dir, self.model_name))
        self.profile_factory = ProfileFactory(self.cfg)


    def tearDown(self):
        self.clean_working_dir()

    def clean_working_dir(self):
        try:
            shutil.rmtree(self.cfg.working_dir())
        except:
            pass

    def test_fqn(self):
        fqn = 'foo/bla'
        model = Model(fqn, 10)
        self.assertEqual(model.fqn, fqn)

        self.assertEqual(model.name, 'bla')


    def test_inter_gene_max_space(self):
        model_fqn = 'foo/bar'
        inter_gene_max_space_xml = 40
        # test inter_gene_max_space from xml
        model = Model(model_fqn, inter_gene_max_space_xml)
        self.assertEqual(model.inter_gene_max_space, inter_gene_max_space_xml)

        self.clean_working_dir()


    def test_min_genes_required(self):
        model_fqn = 'foo/model_1'
        min_genes_required_xml = 40
        model = Model(model_fqn, 10, min_genes_required=min_genes_required_xml)
        gene_name = 'sctJ_FLG'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene = ModelGene(c_gene, model)
        model.add_mandatory_gene(gene)
        self.assertEqual(model.min_genes_required, min_genes_required_xml)
        model = Model(model_fqn, 10)
        self.assertEqual(model.min_genes_required, len(model.mandatory_genes))

        self.clean_working_dir()


    def test_min_mandatory_genes_required(self):
        model_fqn = 'foo/bar'
        min_mandatory_genes_required_xml = 40
        model = Model(model_fqn, 10, min_mandatory_genes_required=min_mandatory_genes_required_xml)
        gene_name = 'sctJ_FLG'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene = ModelGene(c_gene, model)
        model.add_mandatory_gene(gene)
        self.assertEqual(model.min_mandatory_genes_required, min_mandatory_genes_required_xml)

        system = Model(model_fqn, 10)
        self.assertEqual(system.min_mandatory_genes_required, len(system.mandatory_genes))

        self.clean_working_dir()


    def test_max_nb_genes(self):
        model_fqn = 'foo/bar'
        inter_gene_max_space = 40
        max_nb_genes_xml = 10
        model = Model(model_fqn, inter_gene_max_space, max_nb_genes=max_nb_genes_xml)
        self.assertEqual(model.max_nb_genes, max_nb_genes_xml)

        model = Model(model_fqn, inter_gene_max_space)
        self.assertEqual(model.max_nb_genes, 0)

        c_gene_sctc = CoreGene(self.model_location, "sctC", self.profile_factory)
        gene_sctc = ModelGene(c_gene_sctc, model)

        c_gene_abc = CoreGene(self.model_location, "abc", self.profile_factory)
        gene_abc = ModelGene(c_gene_abc, model)

        model.add_mandatory_gene(gene_sctc)
        model.add_accessory_gene(gene_abc)
        self.assertEqual(model.max_nb_genes, 2)
        self.clean_working_dir()


    def test_multi_loci(self):
        model_fqn = 'foo/True'
        inter_gene_max_space = 40
        model = Model(model_fqn, inter_gene_max_space, multi_loci=True)
        self.assertTrue(model.multi_loci)
        model_fqn = 'foo/False'
        inter_gene_max_space = 40
        model = Model(model_fqn, inter_gene_max_space)
        self.assertFalse(model.multi_loci)

        self.clean_working_dir()

        self.args.multi_loci = 'foo/False'

        model_fqn = 'foo/False'
        inter_gene_max_space = 40
        model = Model(model_fqn, inter_gene_max_space, multi_loci=False)
        self.assertFalse(model.multi_loci)


    def test_accessor_mutator(self):
        model = Model("foo", 10)
        gene_name = 'sctJ_FLG'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene = ModelGene(c_gene, model)
        categories = set(model.gene_category)
        for cat in categories:
            other_cat = categories - {cat}
            getattr(model, f'add_{cat}_gene')(gene)
            self.assertEqual(getattr(model, f'{cat}_genes'), [gene])
            for other in other_cat:
                self.assertEqual(getattr(model, f'{other}_genes'), [])
            # don't forget to reset the model to avoid
            # to accumulate genes
            model = Model("foo", 10)

    def test_get_gene(self):
        model = Model("foo", 10)
        gene_name = 'sctJ_FLG'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene = ModelGene(c_gene, model)
        for meth in [getattr(model, f'add_{cat}_gene') for cat in model.gene_category]:
            for cat in model.gene_category:
                setattr(model, f'_{cat}_genes', [])
            meth(gene)
            self.assertEqual(gene, model.get_gene(gene_name))

        self.assertRaises(KeyError, model.get_gene, 'bar')

        homolog_name = 'sctJ'
        c_gene_homolog = CoreGene(self.model_location, homolog_name, self.profile_factory)
        homolog = Exchangeable(c_gene_homolog, gene)
        gene.add_exchangeable(homolog)
        for meth in [getattr(model, f'add_{cat}_gene') for cat in model.gene_category]:
            for cat in model.gene_category:
                setattr(model, f'_{cat}_genes', [])
            meth(gene)
            self.assertEqual(homolog, model.get_gene(homolog_name))


    def test_str(self):
        model_fqn = "foo/bar"
        model = Model(model_fqn, 10)
        gene_name = 'sctJ_FLG'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        mandatory_gene = ModelGene(c_gene, model)
        model.add_mandatory_gene(mandatory_gene)
        homolog_name = 'sctJ'
        c_gene_homolg = CoreGene(self.model_location, homolog_name, self.profile_factory)
        homolog = Exchangeable(c_gene_homolg, mandatory_gene)
        mandatory_gene.add_exchangeable(homolog)

        gene_name = 'sctN_FLG'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        accessory_gene = ModelGene(c_gene, model)
        model.add_accessory_gene(accessory_gene)
        analog_name = 'sctN'
        c_gene_analog = CoreGene(self.model_location, analog_name, self.profile_factory)
        analog = Exchangeable(c_gene_analog, accessory_gene)
        accessory_gene.add_exchangeable(analog)

        gene_name = 'toto'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        neutral_gene = ModelGene(c_gene, model)
        model.add_neutral_gene(neutral_gene)

        gene_name = 'sctC'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        forbidden_gene = ModelGene(c_gene, model)
        model.add_forbidden_gene(forbidden_gene)

        exp_str = """name: bar
fqn: foo/bar
==== mandatory genes ====
sctJ_FLG
==== accessory genes ====
sctN_FLG
==== neutral genes ====
toto
==== forbidden genes ====
sctC
============== end pprint model ================
"""
        self.assertEqual(str(model), exp_str)

    def test_eq(self):
        aa1 = Model("aaa", 10)
        aa2 = Model("aaa", 10)
        self.assertEqual(aa1, aa2)

    def test_lt(self):
        aaa = Model("aaa", 10)
        zzz = Model("zzz", 10)
        self.assertLess(aaa, zzz)

    def test_gt(self):
        aaa = Model("aaa", 10)
        zzz = Model("zzz", 10)
        self.assertGreater(zzz, aaa)

    def test_filter(self):
        model_fqn = "foo/bar"
        model = Model(model_fqn, 10)
        model_2 = Model("foo/buz", 10)

        gene_name = 'sctJ_FLG'
        sctJ_FLG_core = CoreGene(self.model_location, gene_name, self.profile_factory)
        sctJ_FLG = ModelGene(sctJ_FLG_core, model)
        model.add_mandatory_gene(sctJ_FLG)

        gene_name = 'sctJ'
        sctJ_core = CoreGene(self.model_location, gene_name, self.profile_factory)
        sctj = Exchangeable(sctJ_core, sctJ_FLG)
        sctJ_FLG.add_exchangeable(sctj)

        gene_name = 'sctN_FLG'
        sctN_FLG_core = CoreGene(self.model_location, gene_name, self.profile_factory)
        sctN_FLG = ModelGene(sctN_FLG_core, model)
        model.add_accessory_gene(sctN_FLG)

        gene_name = 'sctN'
        sctN_core = CoreGene(self.model_location, gene_name, self.profile_factory)
        sctn = Exchangeable(sctN_core, sctN_FLG)
        sctN_FLG.add_exchangeable(sctn)

        gene_name = 'sctC'
        sctC_core = CoreGene(self.model_location, gene_name, self.profile_factory)
        sctC = ModelGene(sctC_core, model)
        model.add_forbidden_gene(sctC)

        gene_name = 'toto'
        toto_core = CoreGene(self.model_location, gene_name, self.profile_factory)
        toto = ModelGene(toto_core, model)
        model.add_neutral_gene(toto)

        gene_name = 'totote'
        totote_core = CoreGene(self.model_location, gene_name, self.profile_factory)
        totote = Exchangeable(totote_core, toto)
        toto.add_exchangeable(totote)

        gene_name = 'gspD'
        gspd_core = CoreGene(self.model_location, gene_name, self.profile_factory)
        gspd = ModelGene(gspd_core, model_2)

        gene_name = 'tadZ'
        tadz_core = CoreGene(self.model_location, gene_name, self.profile_factory)
        tadz = Exchangeable(tadz_core, gspd)
        gspd.add_exchangeable(tadz)

        hit_to_keep = []
        for gene in (sctJ_FLG, sctN_FLG, sctC, toto, totote):
            hit_to_keep.append(Hit(gene,
                                   f"PSAE001c01_{gene.name}",
                                   1, "PSAE001c01", 1, 1.0, 1.0, 1.0, 1.0, 1, 2)
                               )
        hit_to_filter_out = []
        for gene in (gspd, tadz):
            hit_to_filter_out.append(Hit(gene,
                                     f"PSAE001c01_{gene.name}",
                                     1, "PSAE001c01", 1, 1.0, 1.0, 1.0, 1.0, 1, 2)
                                     )

        filtered_hits = model.filter(hit_to_keep + hit_to_filter_out)

        self.assertListEqual(sorted(hit_to_keep), sorted(filtered_hits))


    def test_hash(self):
        model_bar = Model('Foo/bar', 10)
        model_bar_bis = Model('Foo/bar', 10)
        model_buz = Model('Foo/buz', 10)
        self.assertTrue(isinstance(hash(model_bar), int))
        self.assertEqual(hash(model_bar), hash(model_bar_bis))
        self.assertNotEqual(hash(model_bar), hash(model_buz))