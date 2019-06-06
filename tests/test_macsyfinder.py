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
from io import StringIO

from macsypy.config import Config, MacsyDefaults
from macsypy.gene import ProfileFactory, Gene, GeneStatus
from macsypy.registries import ModelRegistry
from macsypy.hit import HitRegistry, Hit, ValidHit
from macsypy.model import Model
from macsypy.system import System
from macsypy.cluster import Cluster
from macsypy.scripts.macsyfinder import get_models_name_to_detect, systems_to_file, rejected_clst_to_file
import macsypy
from tests import MacsyTest


class Test(MacsyTest):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()


    def tearDown(self):
        try:
            shutil.rmtree(self.tmp_dir)
        except:
            pass

    def test_models_name_to_detect(self):
        cmd_args = argparse.Namespace()
        cmd_args.models_dir = os.path.join(self._data_dir, 'data_set_1', 'models')
        cmd_args.models = [('set_1', 'T9SS', 'T3SS', 'T4SS_typeI')]
        config = Config(MacsyDefaults(models_dir=os.path.join(self._data_dir, 'data_set_1', 'models')),
                        cmd_args)
        registry = ModelRegistry(config)
        res = get_models_name_to_detect([('set_1', ['T9SS', 'T3SS', 'T4SS_typeI'])], registry)
        exp = ['set_1/T9SS', 'set_1/T3SS', 'set_1/T4SS_typeI']
        self.assertListEqual(res, exp)
        with self.assertRaises(ValueError):
            get_models_name_to_detect([('set_1', ['FOO', 'BAR'])], registry)


    def test_list_models(self):
        cmd_args = argparse.Namespace()
        cmd_args.models_dir = os.path.join(self._data_dir, 'data_set_1', 'models')
        cmd_args.list_models = True
        config = Config(MacsyDefaults(models_dir=os.path.join(self._data_dir, 'data_set_1', 'models')),
                        cmd_args)
        registry = ModelRegistry(config)
        list_models = """set_1
      /CONJ
      /Flagellum
      /T2SS
      /T3SS
      /T4P
      /T4SS_typeF
      /T4SS_typeI
      /T9SS
      /Tad
set_2
      /CONJ
      /Flagellum
      /T2SS
      /T3SS
      /T4P
      /T4SS_typeF
      /T4SS_typeI
      /T9SS
      /Tad
"""
        self.assertEqual(str(registry), list_models)

    def test_systems_to_file(self):
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
        system_1 = System(model, [Cluster([v_hit_1, v_hit_2], model)])

        system_str = """# macsyfinder {}
# tests/run_tests.py -vv tests/test_macsyfinder.py
# Systems found:

system id = replicon_id_T2SS_1
model = foo/T2SS 
loci nb = 1
replicon = replicon_id
clusters = [('gspD', 1), ('sctJ', 1)]
occ = 1

mandatory genes:
\t- gspD: 1 (gspD)

accessory genes:
\t- sctJ: 1 (sctJ)

============================================================
""".format(macsypy.__version__)
        f_out = StringIO()
        systems_to_file([system_1], f_out)
        self.assertMultiLineEqual(system_str, f_out.getvalue())

    def test_rejected_clst_to_file(self):
        pass

    def parse_args(self):
        pass

