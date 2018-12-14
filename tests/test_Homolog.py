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
from macsypy.gene import Homolog
from macsypy.gene import Gene
from macsypy.system import System
from macsypy.config import Config
from macsypy.registries import ModelRegistry
from tests import MacsyTest


class Test(MacsyTest):

    def setUp(self):
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        
        # add only one handler to the macsypy logger
        from macsypy.gene import _log
        macsy_log = _log.parent
        log_file = os.devnull
        log_handler = logging.FileHandler(log_file)
        macsy_log.addHandler(log_handler)
        
        self.cfg = Config(hmmer_exe="",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type="gembase",
                          e_value_res=1,
                          i_evalue_sel=0.5,
                          models_dir=self.find_data('models'),
                          res_search_dir=tempfile.gettempdir(),
                          res_search_suffix="",
                          profile_suffix=".hmm",
                          res_extract_suffix="",
                          log_level=30,
                          log_file=log_file
                          )
        models_registry = ModelRegistry(self.cfg)
        self.model_name = 'foo'
        self.models_location = models_registry[self.model_name]


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

    def test_gene_ref(self):
        system = System(self.cfg, "T2SS", 10)
        gene_ref = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        gene = Gene(self.cfg, 'sctJ', system, self.models_location)
        homolog_1 = Homolog(gene, gene_ref)
        self.assertEqual(homolog_1.gene_ref, gene_ref)
 
    def test_is_aligned(self):
        system = System(self.cfg, "T2SS", 10)
        gene_ref = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        gene = Gene(self.cfg, 'sctJ', system, self.models_location)
        homolog = Homolog(gene, gene_ref)
        self.assertFalse(homolog.is_aligned())
        homolog = Homolog(gene, gene_ref, aligned=True)
        self.assertTrue(homolog.is_aligned())

    def test_delegation(self):
        system = System(self.cfg, "T2SS", 10)
        gene_ref = Gene(self.cfg, 'sctJ_FLG', system, self.models_location)
        gene = Gene(self.cfg, 'sctJ', system, self.models_location)
        homolog = Homolog(gene, gene_ref)
        self.assertEqual(homolog.system, system)

