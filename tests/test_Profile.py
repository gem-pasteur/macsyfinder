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
import unittest
import shutil
import tempfile
import sysconfig
import argparse

from macsypy.profile import Profile
from macsypy.gene import CoreGene, ModelGene
from macsypy.model import Model
from macsypy.profile import ProfileFactory
from macsypy.config import Config, MacsyDefaults
from macsypy.registries import ModelLocation
from tests import MacsyTest, which


class TestProfile(MacsyTest):

    def setUp(self):
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.res_search_dir = tempfile.gettempdir()
        args.log_level = 0
        self.cfg = Config(MacsyDefaults(), args)

        if os.path.exists(self.cfg.working_dir()):
            shutil(self.cfg.working_dir())
        os.makedirs(self.cfg.working_dir())

        self.model_name = 'foo'
        self.model_location = ModelLocation(path=os.path.join(args.models_dir, self.model_name))
        self.profile_factory = ProfileFactory(self.cfg)


    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir())
        except:
            pass


    def test_len(self):
        model = Model("foo/T2SS", 10)

        gene_name = 'abc'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene = ModelGene(c_gene, model)

        path = self.model_location.get_profile("abc")
        profile = Profile(gene, self.cfg, path)
        self.assertEqual(len(profile), 501)


    def test_ga_threshold(self):
        model = Model("foo/T2SS", 10)
        gene_name = 'abc'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene = ModelGene(c_gene, model)

        path = self.model_location.get_profile("abc")
        profile = Profile(gene, self.cfg, path)
        self.assertFalse(profile.ga_threshold)

        gene_name = 'T5aSS_PF03797'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene = ModelGene(c_gene, model)

        path = self.model_location.get_profile("T5aSS_PF03797")
        profile = Profile(gene, self.cfg, path)
        self.assertTrue(profile.ga_threshold)

    def test_str(self):
        model = Model("foo/T2SS", 10)

        gene_name = 'abc'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene = ModelGene(c_gene, model)

        path = self.model_location.get_profile("abc")
        profile = Profile(gene, self.cfg, path)
        s = "{0} : {1}".format(gene.name, path)
        self.assertEqual(str(profile), s)


    @unittest.skipIf(not which('hmmsearch'), 'hmmsearch not found in PATH')
    def test_execute(self):
        for db_type in ("gembase", "ordered_replicon", "unordered"):
            self.cfg._set_db_type(db_type)
            model = Model("foo/T2SS", 10)

            gene_name = 'T5aSS_PF03797'
            c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
            gene = ModelGene(c_gene, model)

            # case GA threshold in profile
            profile_path = self.model_location.get_profile("T5aSS_PF03797")
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
                # a hmmsearch used the abc profile line should become with: "# query HMM file: {the path tp hmm profile used}"
                self.assertTrue(l.find(profile_path) != -1)
                for i in range(3):
                    # skip 2 lines
                    l = hmmer_raw_out_file.readline()
                self.assertEqual("# model-specific thresholding:     GA cutoffs", l.strip())
            # test if profile is executed only once per run
            report_bis = profile.execute()
            self.assertIs(report, report_bis)

            # case GA threshold in profile but --no-cut-ga is set
            args = argparse.Namespace()
            args.sequence_db = self.find_data("base", "test_1.fasta")
            args.db_type = 'gembase'
            args.models_dir = self.find_data('models')
            args.res_search_dir = tempfile.gettempdir()
            args.log_level = 0
            args.e_value_search = 0.5
            args.no_cut_ga = True
            cfg = Config(MacsyDefaults(), args)

            profile = Profile(gene, cfg, profile_path)
            report = profile.execute()
            hmmer_raw_out = profile.hmm_raw_output
            with open(hmmer_raw_out, 'r') as hmmer_raw_out_file:
                for i in range(9):
                    l = hmmer_raw_out_file.readline()
                self.assertEqual("# sequence reporting threshold:    E-value <= 0.5", l.strip())


            # case cut-ga but no GA threshold in hmmprofile
            gene_name = 'abc'
            c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
            gene = ModelGene(c_gene, model)

            # case -cut-ga and GA threshold in profile
            profile_path = self.model_location.get_profile("abc")
            profile = Profile(gene, self.cfg, profile_path)

            with self.catch_log() as log:
                report = profile.execute()

            hmmer_raw_out = profile.hmm_raw_output
            with open(hmmer_raw_out, 'r') as hmmer_raw_out_file:
                first_l = hmmer_raw_out_file.readline()
                # a hmmsearch output file has been produced
                self.assertTrue(first_l.startswith("# hmmsearch :: search profile(s) against a sequence database"))
                for i in range(5):
                    # skip 4 lines
                    l = hmmer_raw_out_file.readline()
                # a hmmsearch used the abc profile line should become with: "# query HMM file: {the path tp hmm profile used}"
                self.assertTrue(l.find(profile_path) != -1)
                for i in range(3):
                    # skip 2 lines
                    l = hmmer_raw_out_file.readline()
                self.assertEqual('# sequence reporting threshold:    E-value <= 0.1', l.strip())


    def test_execute_unknown_binary(self):
        self.cfg._options['hmmer'] = "Nimportnaoik"
        model = Model("foo/T2SS", 10)

        gene_name = 'abc'
        c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
        gene = ModelGene(c_gene, model)

        path = self.model_location.get_profile("abc", )
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
            model = Model("foo/T2SS", 10)

            gene_name = 'abc'
            c_gene = CoreGene(self.model_location, gene_name, self.profile_factory)
            gene = ModelGene(c_gene, model)

            path = self.model_location.get_profile("abc", )
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
