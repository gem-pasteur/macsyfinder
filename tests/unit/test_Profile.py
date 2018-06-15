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
import unittest
import shutil
import tempfile
import logging
from macsypy.gene import Profile
from macsypy.gene import Gene
from macsypy.system import System
from macsypy.config import Config
from macsypy.registries import ModelRegistry
from macsypy.utils import which
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
        
        self.cfg = Config(hmmer_exe="hmmsearch",
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


    def test_len(self):
        system = System(self.cfg, "foo/T2SS", 10)
        gene = Gene(self.cfg, "abc", system, self.models_location)
        path = self.models_location.get_profile("abc")
        profile = Profile(gene, self.cfg, path)
        self.assertEqual(len(profile), 501)


    def test_str(self):
        system = System(self.cfg, "foo/T2SS", 10)
        gene = Gene(self.cfg, "abc", system, self.models_location)
        path = self.models_location.get_profile("abc")
        profile = Profile(gene, self.cfg, path)
        s = "{0} : {1}".format(gene.name, path)
        self.assertEqual(str(profile), s)


    @unittest.skipIf(not which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_execute(self):
        system = System(self.cfg, "foo/T2SS", 10)
        gene = Gene(self.cfg, "abc", system, self.models_location)
        profile_path = self.models_location.get_profile("abc")
        profile = Profile(gene, self.cfg, profile_path)
        report = profile.execute()
        hmmer_raw_out = profile.hmm_raw_output
        with open(hmmer_raw_out, 'r') as hmmer_raw_out_file:
            first_l = hmmer_raw_out_file.readline()
            # a hmmsearch output file has been produced
            self.assertTrue(first_l.startswith("# hmmsearch :: search profile(s) against a sequence database"))
            for i in range(5):
                # skip 4 lines
                l = hmmer_raw_out_file.readline()
            # a hmmsearch used the abc profile line should become with: "# query HMM file:"
            self.assertTrue(l.find(profile_path) != -1)


    def test_execute_unknown_binary(self):
        self.cfg.options['hmmer_exe'] = "Nimportnaoik"
        system = System(self.cfg, "foo/T2SS", 10)
        gene = Gene(self.cfg, "abc", system, self.models_location)
        path = self.models_location.get_profile("abc")
        profile = Profile(gene, self.cfg, path)
        self.assertRaises(RuntimeError, profile.execute)

