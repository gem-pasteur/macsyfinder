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


import os
import unittest
import shutil
import tempfile
import sysconfig
import argparse

from macsypy.profile import ProfileFactory, Profile
from macsypy.gene import Gene
from macsypy.model import Model
from macsypy.config import Config, MacsyDefaults
from macsypy.registries import ModelLocation
from tests import MacsyTest, which


class TestProfile(MacsyTest):

    def setUp(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.res_search_dir = tempfile.gettempdir()
        args.log_level = 30
        self.cfg = Config(MacsyDefaults(), args)

        if os.path.exists(self.cfg.working_dir()):
            shutil(self.cfg.working_dir())
        os.makedirs(self.cfg.working_dir())

        self.model_name = 'foo'
        self.models_location = ModelLocation(path=os.path.join(args.models_dir, self.model_name))
        self.profile_factory = ProfileFactory(self.cfg)

    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir())
        except:
            pass


    def test_len(self):
        system = Model("foo/T2SS", 10)
        gene = Gene(self.profile_factory, "abc", system, self.models_location)
        path = self.models_location.get_profile("abc", )
        profile = Profile(gene, self.cfg, path)
        self.assertEqual(len(profile), 501)


    def test_str(self):
        system = Model("foo/T2SS", 10)
        gene = Gene(self.profile_factory, "abc", system, self.models_location)
        path = self.models_location.get_profile("abc", )
        profile = Profile(gene, self.cfg, path)
        s = "{0} : {1}".format(gene.name, path)
        self.assertEqual(str(profile), s)


    @unittest.skipIf(not which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_execute(self):
        for db_type in ("gembase", "ordered_replicon", "unordered"):
            self.cfg._set_db_type(db_type)
            system = Model("foo/T2SS", 10)
            gene = Gene(self.profile_factory, "abc", system, self.models_location)
            profile_path = self.models_location.get_profile("abc", )
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

            report_bis = profile.execute()
            self.assertIs(report, report_bis)


    def test_execute_unknown_binary(self):
        self.cfg._options['hmmer'] = "Nimportnaoik"
        system = Model("foo/T2SS", 10)
        gene = Gene(self.profile_factory, "abc", system, self.models_location)
        path = self.models_location.get_profile("abc", )
        profile = Profile(gene, self.cfg, path)
        with self.catch_log():
            with self.assertRaises(RuntimeError):
                profile.execute()


    def test_execute_hmmer_failed(self):
        fake_hmmer = os.path.join(tempfile.gettempdir(), 'hmmer_failed')
        with open(fake_hmmer, 'w') as hmmer:
            hmmer.write("""#! {}
import sys
sys.exit(127)
""".format(sysconfig.sys.executable))
        try:
            os.chmod(hmmer.name, 0o755)
            self.cfg._options['hmmer'] = hmmer.name
            system = Model("foo/T2SS", 10)
            gene = Gene(self.profile_factory, "abc", system, self.models_location)
            path = self.models_location.get_profile("abc", )
            profile = Profile(gene, self.cfg, path)
            with self.catch_log():
                with self.assertRaisesRegex(RuntimeError,
                                         "an error occurred during Hmmer "
                                         "execution: command = .* : return code = 127 .*") as ctx:
                    profile.execute()
        finally:
            try:
                os.unlink(fake_hmmer)
            except Exception:
                pass





