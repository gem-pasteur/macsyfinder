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

import tempfile
import shutil
import os
import argparse
import sys
import logging
from io import StringIO
from itertools import groupby

import macsypy
from macsypy.config import MacsyDefaults, Config
from macsypy.database import Indexes
from tests import MacsyTest
from macsypy.scripts import macsyprofile


class TestMacsyprofile(MacsyTest):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

        self.args = argparse.Namespace()
        self.previous_run = self.find_data('functional_test_gembase')

    def tearDown(self):
        try:
            shutil.rmtree(self.tmpdir)
        except:
            pass
        # some function in macsyprofile script suppress the traceback
        # but without traceback it's hard to debug test :-(
        sys.tracebacklimit = 1000  # the default value

        # at each call of macsyprofile.main
        # init_logger is call and new handler is add
        # remove them
        logger = macsyprofile._log
        if logger:
            for handler in logger.handlers[:]:
                logger.removeHandler(handler)


    def test_pasre_args(self):
        cmd = f"macsyprofile {self.previous_run}"
        args = macsyprofile.parse_args(cmd.split()[1:])
        self.assertEqual(args.previous_run, self.previous_run)
        self.assertEqual(args.coverage_profile, -1.)
        self.assertEqual(args.i_evalue_sel, 1.0e9)
        self.assertIsNone(args.best_hits)
        self.assertEqual(args.pattern, '*')
        self.assertIsNone(args.out)
        self.assertFalse(args.force)
        self.assertEqual(args.verbosity, 0)
        self.assertFalse(args.mute)

        cmd = f"macsyprofile -vv {self.previous_run}"
        args = macsyprofile.parse_args(cmd.split()[1:])
        self.assertEqual(args.verbosity, 2)

        cmd = f"macsyprofile -f {self.previous_run}"
        args = macsyprofile.parse_args(cmd.split()[1:])
        self.assertTrue(args.force)

        cmd = f"macsyprofile --out toto {self.previous_run}"
        args = macsyprofile.parse_args(cmd.split()[1:])
        self.assertEqual(args.out, 'toto')

        pattern = '?toto'
        cmd = f"macsyprofile --pattern {pattern} {self.previous_run}"
        args = macsyprofile.parse_args(cmd.split()[1:])
        self.assertEqual(args.pattern, pattern)

        cmd = f"macsyprofile --best-hits score {self.previous_run}"
        args = macsyprofile.parse_args(cmd.split()[1:])
        self.assertEqual(args.best_hits, 'score')

        cmd = f"macsyprofile --i-evalue-sel 1.0 {self.previous_run}"
        args = macsyprofile.parse_args(cmd.split()[1:])
        self.assertEqual(args.i_evalue_sel, 1.0)

        cmd = f"macsyprofile --coverage-profile 1.0 {self.previous_run}"
        args = macsyprofile.parse_args(cmd.split()[1:])
        self.assertEqual(args.coverage_profile, 1.0)


    def test_verbosity_to_log_level(self):
        level = macsyprofile.verbosity_to_log_level(1)
        self.assertEqual(level, 10)
        level = macsyprofile.verbosity_to_log_level(5)
        self.assertEqual(level, 1)


    def test_header(self):
        out = "FOO"
        coverage_profile = 0.1
        version = macsypy.__version__
        cmd = f"macsyprofile --coverage-profile {coverage_profile} --out {out} --index-dir {self.tmpdir} {self.previous_run}"
        expected_header = f"""# macsyprofile {version}
# macsyprofile {' '.join(cmd.split()[1:])}
hit_id\treplicon_name\tposition_hit\thit_sequence_length\tgene_name\ti_eval\tscore\tprofile_coverage\tsequence_coverage\tbegin\tend"""
        got_header = macsyprofile.header(cmd.split()[1:])
        self.assertEqual(expected_header, got_header)


    def test_get_gene_name(self):
        expected_gene_name = "Archaeal-T4P_arCOG01822"
        path = f"tests/data/functional_test_gembase/hmmer_results/{expected_gene_name}.res_hmm_extract"
        suffix = ".res_hmm_extract"
        self.assertEqual(macsyprofile.get_gene_name(path, suffix), expected_gene_name)


    def test_getProfile_len(self):
        path = self.find_data('models', 'functional', 'profiles', 'T1SS_abc.hmm')
        self.assertEqual(macsyprofile.get_profile_len(path), 501)


    def test_LightHit(self):
        data = {'id': 'GCF_000006745_003720',
                'replicon_name': 'GCF_000006745',
                'position': 372,
                'seq_length': 196,
                'gene_name': 'T4bP_pilA',
                'i_eval': 2.000e+01,
                'score': 3.400,
                'profile_coverage': 0.219,
                'sequence_coverage': 0.153,
                'begin_match': 147,
                'end_match': 176,
                 }
        l_hit = macsyprofile.LightHit(data['gene_name'],
                                      data['id'],
                                      data['seq_length'],
                                      data['replicon_name'],
                                      data['position'],
                                      data['i_eval'],
                                      data['score'],
                                      data['profile_coverage'],
                                      data['sequence_coverage'],
                                      data['begin_match'],
                                      data['end_match'],
                                      )
        exp_str = "GCF_000006745_003720\tGCF_000006745\t372\t196\tT4bP_pilA\t2.000e+01\t3.400\t0.219\t0.153\t147\t176"
        self.assertEqual(str(l_hit), exp_str)


    def test_build_my_db(self):
        gene_name = "gspD"
        args = argparse.Namespace()
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.log_level = 30
        self.cfg = Config(MacsyDefaults(), args)
        gspD_hmmer_path = self.find_data(os.path.join('hmm', 'gspD.search_hmm.out'))

        hmm_prof = macsyprofile.HmmProfile(gene_name, 596, gspD_hmmer_path, self.cfg)

        db = hmm_prof._build_my_db(gspD_hmmer_path)
        self.assertDictEqual(db, {'PSAE001c01_031420': None,
                                  'PSAE001c01_051090': None,
                                  'PSAE001c01_018920': None,
                                  'PSAE001c01_043580': None,
                                  'PSAE001c01_017350': None,
                                  'PSAE001c01_013980': None,
                                  'PSAE001c01_026600': None,
                                  'NC_xxxxx_xx_056141': None,
                                  'PSAE001c01_006940': None})


    def test_fill_my_db(self):
        gene_name = "gspD"
        args = argparse.Namespace()
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.log_level = 30
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.index_dir = self.tmpdir
        cfg = Config(MacsyDefaults(), args)
        gspD_hmmer_path = self.find_data('hmm', 'gspD.search_hmm.out')

        idx = Indexes(cfg)
        macsyfinder_idx = idx.build()
        hmm_prof = macsyprofile.HmmProfile(gene_name, 596, gspD_hmmer_path, cfg)

        db = hmm_prof._build_my_db(gspD_hmmer_path)
        hmm_prof._fill_my_db(macsyfinder_idx, db)
        self.assertDictEqual(db, {'PSAE001c01_031420': (658, 73),
                                  'PSAE001c01_051090': (714, 75),
                                  'PSAE001c01_018920': (776, 71),
                                  'PSAE001c01_043580': (416, 74),
                                  'PSAE001c01_017350': (600, 70),
                                  'PSAE001c01_013980': (759, 69),
                                  'PSAE001c01_026600': (273, 72),
                                  'NC_xxxxx_xx_056141': (803, 141),
                                  'PSAE001c01_006940': (803, 68)})


    def test_hit_start(self):
        gene_name = "gspD"
        args = argparse.Namespace()
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.log_level = 30
        self.cfg = Config(MacsyDefaults(), args)
        gspD_hmmer_path = self.find_data(os.path.join('hmm', 'gspD.search_hmm.out'))

        hmm_prof = macsyprofile.HmmProfile(gene_name, 596, gspD_hmmer_path, self.cfg)

        self.assertFalse(hmm_prof._hit_start("NOT starting hit"))
        self.assertTrue(hmm_prof._hit_start(">> starting hit"))


    def test_parse_hmm_header(self):
        gene_name = "gspD"
        args = argparse.Namespace()
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.log_level = 30
        cfg = Config(MacsyDefaults(), args)
        gspD_hmmer_path = self.find_data(os.path.join('hmm', 'gspD.search_hmm.out'))

        hmm_prof = macsyprofile.HmmProfile(gene_name, 596, gspD_hmmer_path, cfg)

        hmm_hit = [">> NC_xxxxx_xx_056141  C ATG TAA 6260390 6261757 Valid PA5567 1368 _NP_254254.1_ PA5567 1 "
                   "6260390 6261757 | tRNA modific"]
        hit_id = hmm_prof._parse_hmm_header(hmm_hit)
        self.assertEqual(hit_id, 'NC_xxxxx_xx_056141')


    def test_parse_hmm_body(self):
        def make_hmm_group(hmm_string):
            hmm_file = StringIO(hmm_string)
            hmm_hits = (x[1] for x in groupby(hmm_file, lambda l: l.startswith('>>')))
            header = next(hmm_hits)
            body = next(hmm_hits)
            return body


        gene_name = "gspD"
        args = argparse.Namespace()
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.log_level = 30
        cfg = Config(MacsyDefaults(), args)
        gspD_hmmer_path = self.find_data(os.path.join('hmm', 'gspD.search_hmm.out'))

        hmm_prof = macsyprofile.HmmProfile(gene_name, 596, gspD_hmmer_path, cfg)

        # with one significant hit
        hmm = """>> NC_xxxxx_xx_056141  C ATG TAA 6260390 6261757 Valid PA5567 1368 _NP_254254.1_ PA5567 1 6260390 6261757 | tRNA modific
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  779.2   5.5  1.4e-237    2e-236       1     596 []     104     741 ..     104     741 .. 0.93

  Alignments for each domain:
"""
        body = make_hmm_group(hmm)
        hits = hmm_prof._parse_hmm_body('NC_xxxxx_xx_056141', 596, 803, 0.5, 'NC_xxxxx_xx', 141, 0.5, body)
        expected_hits = [macsyprofile.LightHit(gene_name, "NC_xxxxx_xx_056141", 803, "NC_xxxxx_xx", 141, float(2e-236), float(779.2),
                           float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)]
        self.assertListEqual(hits, expected_hits)
        # with no significant hit
        hmm = """>> PSAE001c01_051090  C ATG TGA 5675714 5677858 Valid pilQ 2145 _PA5040_NP_253727.1_ PA5040 1 5675714 5677858 | type 4 f
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !   27.1   0.2   6.3e-10   6.6e-07       1     120 [.     286     402 ..     286     407 .. 0.86
   2 !  186.2   0.1   4.2e-58   4.3e-55     294     590 ..     405     709 ..     397     712 .. 0.84

  Alignments for each domain:
"""
        body = make_hmm_group(hmm)
        hits = hmm_prof._parse_hmm_body('NC_xxxxx_xx_056141', 596, 803, 0.5, 'NC_xxxxx_xx', 141, 0.5, body)
        expected_hits = []
        self.assertListEqual(hits, expected_hits)

        # with no hit
        hmm = """>> PSAE001c01_051090  C ATG TGA 5675714 5677858 Valid pilQ 2145 _PA5040_NP_253727.1_ PA5040 1 5675714 5677858 | type 4 f
        bla bla
        """
        body = make_hmm_group(hmm)
        hits = hmm_prof._parse_hmm_body('NC_xxxxx_xx_056141', 596, 803, 0.5, 'NC_xxxxx_xx', 141, 0.5, body)
        expected_hits = []
        self.assertListEqual(hits, expected_hits)

        # with invalid hmm
        hmm = """>> NC_xxxxx_xx_056141  C ATG TAA 6260390 6261757 Valid PA5567 1368 _NP_254254.1_ PA5567 1 6260390 6261757 | tRNA modific
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  779.2   5.5  1.4e-237    foo       1     596 []     104     741 ..     104     741 .. 0.93

  Alignments for each domain:
"""
        body = make_hmm_group(hmm)
        with self.assertRaises(ValueError) as ctx:
            hmm_prof._parse_hmm_body('NC_xxxxx_xx_056141', 596, 803, 0.5, 'NC_xxxxx_xx', 141, 0.5, body)
        self.assertEqual(str(ctx.exception), """Invalid line to parse :   1 !  779.2   5.5  1.4e-237    foo       1     596 []     104     741 ..     104     741 .. 0.93
:could not convert string to float: 'foo'""")


    def test_parse(self):
        gene_name = "gspD"
        args = argparse.Namespace()
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.log_level = 30
        args.sequence_db = self.find_data("base", "test_base.fa")
        args.profile_coverage = -1
        args.i_evalue_sel = 10e9
        args.index_dir = self.tmpdir
        cfg = Config(MacsyDefaults(), args)
        gspD_hmmer_path = self.find_data('hmm', 'gspD.search_hmm.out')

        hmm_prof = macsyprofile.HmmProfile(gene_name, 596, gspD_hmmer_path, cfg)

        expected_hits = [
            macsyprofile.LightHit('gspD', 'NC_xxxxx_xx_056141', 803, 'NC_xxxxx_xx', 141, 2.000e-236, 779.200, 1.000,
                                  (741.0 - 104.0 + 1) / 803, 104, 741),
            macsyprofile.LightHit('gspD', 'PSAE001c01_006940', 803, 'PSAE001c01', 68, 1.2e-234, 779.2, 1.0,
                                  (741.0 - 104.0 + 1) / 803, 104, 741),
            macsyprofile.LightHit('gspD', 'PSAE001c01_031420', 658, 'PSAE001c01', 73, 1.8e-210, 699.3, 1.0,
                                  (614.0 - 55.0 + 1) / 658, 55, 614),
            macsyprofile.LightHit('gspD', 'PSAE001c01_018920', 776, 'PSAE001c01', 71, 6.1e-183, 608.4, 1.0,
                                  (606.0 - 48.0 + 1) / 776, 48, 606),
            macsyprofile.LightHit('gspD', 'PSAE001c01_013980', 759, 'PSAE001c01', 69, 3.7e-76, 255.8, 1.0,
                                  (736.0 - 105.0 + 1) / 759, 105, 736),
            macsyprofile.LightHit('gspD', 'PSAE001c01_017350', 600, 'PSAE001c01', 70, 3.2e-27, 94.2, 0.5,
                                  (506.0 - 226.0 + 1) / 600, 226, 506),
        ]

        hits = hmm_prof.parse()
        self.assertListEqual(expected_hits, hits)

    def test_functional(self):
        out = os.path.join(self.tmpdir, 'test_macsyprofile')
        previous_run = self.find_data('functional_test_degenerated_systems')
        expected_result = self.find_data('results_macsyprofile.tsv')
        cmd = f"macsyprofile -o {out} --index-dir {self.tmpdir} {previous_run} "
        macsyprofile.main(cmd.split()[1:], log_level='WARNING')
        self.assertFileEqual(expected_result, out, comment='#')

    def test_functional_pattern(self):
        out = os.path.join(self.tmpdir, 'test_macsyprofile')
        previous_run = self.find_data('functional_test_degenerated_systems')
        expected_result = self.find_data('results_macsyprofile_pattern.tsv')
        # the argument do not need to be protected (we do not use shell)
        # as on real command line so '*mpf' => *mpf
        cmd = f"macsyprofile -o {out} --index-dir {self.tmpdir} -p *mfp {previous_run} "
        macsyprofile.main(cmd.split()[1:], log_level='WARNING')
        self.assertFileEqual(expected_result, out, comment='#')

    def test_functional_coverage(self):
        out = os.path.join(self.tmpdir, 'test_macsyprofile')
        previous_run = self.find_data('functional_test_degenerated_systems')
        expected_result = self.find_data('results_macsyprofile_coverage.tsv')
        cmd = f"macsyprofile -o {out} --index-dir {self.tmpdir} --coverage-profile 0.5 {previous_run} "
        macsyprofile.main(cmd.split()[1:], log_level='WARNING')
        self.assertFileEqual(expected_result, out, comment='#')

    def test_functional_evalue(self):
        out = os.path.join(self.tmpdir, 'test_macsyprofile')
        previous_run = self.find_data('functional_test_degenerated_systems')
        expected_result = self.find_data('results_macsyprofile_evalue.tsv')
        cmd = f"macsyprofile -o {out} --index-dir {self.tmpdir} --i-evalue-sel 1e-3 {previous_run} "
        macsyprofile.main(cmd.split()[1:], log_level='WARNING')
        self.assertFileEqual(expected_result, out, comment='#')

    def test_functional_best_score(self):
        out = os.path.join(self.tmpdir, 'test_macsyprofile')
        previous_run = self.find_data('functional_test_degenerated_systems')
        expected_result = self.find_data('results_macsyprofile_best_score.tsv')
        cmd = f"macsyprofile -o {out} --index-dir {self.tmpdir} --best-hits score {previous_run} "
        macsyprofile.main(cmd.split()[1:], log_level='WARNING')
        self.assertFileEqual(expected_result, out, comment='#')

    def test_functional_no_hits(self):
        out = os.path.join(self.tmpdir, 'test_macsyprofile')
        previous_run = self.find_data('functional_test_degenerated_systems')
        expected_result = self.find_data('results_macsyprofile_no_hits.tsv')
        cmd = f"macsyprofile -o {out} --index-dir {self.tmpdir} --i-evalue-sel 1e-10 --coverage-profile 2.0 {previous_run}"
        macsyprofile.main(cmd.split()[1:], log_level='WARNING')
        self.assertFileEqual(expected_result, out, comment='#')


    def test_functional_out_exists(self):
        out = os.path.join(self.tmpdir, 'test_macsyprofile')
        open(out, 'w').close()
        previous_run = self.find_data('functional_test_degenerated_systems')
        cmd = f"macsyprofile -o {self.tmpdir} --mute --best-hits score {previous_run} "

        # we cannot test the log message here
        # because the logger are init when main is called
        # after that the context is establish
        # so when the wrapper is called the handler cannot be substitute by the fake
        with self.assertRaises(ValueError):
            macsyprofile.main(cmd.split()[1:])


    def test_functional_no_previous_run(self):
        cmd = "macsyprofile --mute nimportnaoik "

        # we cannot test the log message here
        # because the logger are init when main is called
        # after that the context is establish
        # so when the wrapper is called the handler cannot be substitute by the fake
        with self.assertRaises(FileNotFoundError):
            macsyprofile.main(cmd.split()[1:])

        bad_previous_run = os.path.join(self.tmpdir, 'test_macsyprofile')
        open(bad_previous_run, 'w').close()
        cmd = f"macsyprofile  {bad_previous_run} "

        with self.assertRaises(ValueError):
            macsyprofile.main(cmd.split()[1:], log_level=logging.CRITICAL + 1)

    def test_functional(self):
        old = self.find_data('conf_files', 'macsyfinder-old.conf')
        shutil.copyfile(old, os.path.join(self.tmpdir, 'macsyfinder.conf'))

        previous_run = self.tmpdir
        cmd = f"macsyprofile --index-dir {self.tmpdir} {previous_run}"
        with self.assertRaises(ValueError):
            macsyprofile.main(cmd.split()[1:], log_level=logging.CRITICAL + 1)
