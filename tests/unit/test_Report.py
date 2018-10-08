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
import platform
import logging
from StringIO import StringIO
from itertools import groupby

from macsypy.report import HMMReport, GembaseHMMReport, OrderedHMMReport, Hit
from macsypy.gene import Gene
from macsypy.system import System
from macsypy.config import Config
from macsypy.database import Indexes, RepliconDB
from macsypy.registries import ModelRegistry
from tests import MacsyTest


class TestReport(MacsyTest):


    def setUp(self):
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        
        # add only one handler to the macsypy logger
        from macsypy.gene import _log
        macsy_log = _log.parent
        log_file = 'NUL' if platform.system() == 'Windows' else '/dev/null'
        log_handler = logging.FileHandler(log_file)
        macsy_log.addHandler(log_handler)
        
        self.cfg = Config(hmmer_exe="hmmsearch",
                          sequence_db=self.find_data("base", "test_base.fa"),
                          db_type="gembase",
                          e_value_res=1,
                          i_evalue_sel=0.5,
                          models_dir=self.find_data('models'),
                          res_search_dir=tempfile.gettempdir(),
                          res_search_suffix=".search_hmm.out",
                          profile_suffix=".hmm",
                          res_extract_suffix=".extract",
                          log_level=30,
                          log_file=log_file
                          )
        models_registry = ModelRegistry(self.cfg)
        self.model_name = 'foo'
        self.models_location = models_registry[self.model_name]
        shutil.copy(self.cfg.sequence_db, self.cfg.working_dir)
        self.cfg.options['sequence_db'] = os.path.join(self.cfg.working_dir, os.path.basename(self.cfg.sequence_db))
        idx = Indexes(self.cfg)
        idx._build_my_indexes()


    def tearDown(self):
        # close loggers filehandles, so they don't block file deletion
        # in shutil.rmtree calls in Windows
        logging.shutdown()
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        try:
            shutil.rmtree(self.cfg.working_dir)
            pass
        except Exception:
            pass


class TestHMMReport(TestReport):

    def test_HMMReport(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        shutil.copy(self.find_data("hmm", gene_name + self.cfg.res_search_suffix),
                    self.cfg.working_dir)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        self.assertRaises(TypeError, HMMReport, gene, report_path, self.cfg)

    def test_str(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        shutil.copy(self.find_data("hmm", gene_name + self.cfg.res_search_suffix),
                    self.cfg.working_dir)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        report = GembaseHMMReport(gene, report_path, self.cfg)
        report.extract()

        hits = [Hit(gene, system, "NC_xxxxx_xx_056141", 803, "NC_xxxxx_xx", 141, float(2e-236), float(779.2),
                    float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741),
                Hit(gene, system, "PSAE001c01_006940", 803, "PSAE001c01", 68, float(1.2e-234),
                    float(779.2), float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741),
                Hit(gene, system, "PSAE001c01_013980", 759, "PSAE001c01", 69, float(3.7e-76), float(255.8),
                    float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736),
                Hit(gene, system, "PSAE001c01_017350", 600, "PSAE001c01", 70, float(3.2e-27), float(94.2),
                    float(0.500000), (506.0 - 226.0 + 1) / 600, 226, 506),
                Hit(gene, system, "PSAE001c01_018920", 776, "PSAE001c01", 71, float(6.1e-183), float(608.4),
                    float(1.000000), (606.0 - 48.0 + 1) / 776, 48, 606),
                Hit(gene, system, "PSAE001c01_031420", 658, "PSAE001c01", 73, float(1.8e-210), float(699.3),
                    float(1.000000), (614.0 - 55.0 + 1) / 658, 55, 614)
                ]

        s = "# gene: {0} extract from {1} hmm output\n".format(gene.name, report_path)
        s += "# profile length= {0:d}\n".format(len(gene.profile))
        s += "# i_evalue threshold= {0:.3f}\n".format(self.cfg.i_evalue_sel)
        s += "# coverage threshold= {0:.3f}\n".format(self.cfg.coverage_profile)
        s += "# hit_id replicon_name position_hit hit_sequence_length gene_name gene_system i_eval score " \
             "profile_coverage sequence_coverage begin end\n"
        for h in hits:
            s += str(h)
        self.assertMultiLineEqual(str(report), s)

    def test_save_extract(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        shutil.copy(self.find_data("hmm", gene_name + self.cfg.res_search_suffix),
                    self.cfg.working_dir)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        report = GembaseHMMReport(gene, report_path, self.cfg)
        report.extract()
        report.save_extract()
        extract_filename = gene_name + self.cfg.res_extract_suffix
        extract_path = os.path.join(self.cfg.working_dir, self.cfg.hmmer_dir, extract_filename)
        self.assertTrue(os.path.exists(extract_path))
        self.assertTrue(os.path.isfile(extract_path))

        hits = [Hit(gene, system, "NC_xxxxx_xx_056141", 803, "NC_xxxxx_xx", 141, float(2e-236), float(779.2),
                    float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741),
                Hit(gene, system, "PSAE001c01_006940", 803, "PSAE001c01", 68, float(1.2e-234),
                    float(779.2), float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741),
                Hit(gene, system, "PSAE001c01_013980", 759, "PSAE001c01", 69, float(3.7e-76), float(255.8),
                    float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736),
                Hit(gene, system, "PSAE001c01_017350", 600, "PSAE001c01", 70, float(3.2e-27), float(94.2),
                    float(0.500000), (506.0 - 226.0 + 1) / 600, 226, 506),
                Hit(gene, system, "PSAE001c01_018920", 776, "PSAE001c01", 71, float(6.1e-183), float(608.4),
                    float(1.000000), (606.0 - 48.0 + 1) / 776, 48, 606),
                Hit(gene, system, "PSAE001c01_031420", 658, "PSAE001c01", 73, float(1.8e-210), float(699.3),
                    float(1.000000), (614.0 - 55.0 + 1) / 658, 55, 614)
                ]

        expected_extract_path = os.path.join(self.cfg.working_dir, 'expected_extract')
        with open(expected_extract_path, 'w') as expected_extract:
            extract = """# gene: {name} extract from {path} hmm output
# profile length= {len_profile:d}
# i_evalue threshold= {i_evalue:.3f}
# coverage threshold= {cov:.3f}
# hit_id replicon_name position_hit hit_sequence_length gene_name gene_system i_eval score profile_coverage sequence_coverage begin end
""".format(name=gene.name, path=report_path, len_profile=len(gene.profile),
           i_evalue=self.cfg.i_evalue_sel, cov=self.cfg.coverage_profile)
            expected_extract.write(extract)
            for h in hits:
                expected_extract.write(str(h))

        self.assertFileEqual(extract_path, expected_extract_path)

    def test_best_hit(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        shutil.copy(self.find_data("hmm", gene_name + self.cfg.res_search_suffix),
                    self.cfg.working_dir)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        report = GembaseHMMReport(gene, report_path, self.cfg)
        self.assertIsNone(report.best_hit())
        report.extract()
        best_hit = report.best_hit()
        hit_expected = Hit(gene, system, "NC_xxxxx_xx_056141", 803, "NC_xxxxx_xx", 141, float(2e-236), float(779.2),
                           float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741)
        self.assertEqual(hit_expected, best_hit)


    def test_hit_start(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        shutil.copy(self.find_data("hmm", gene_name + self.cfg.res_search_suffix),
                    self.cfg.working_dir)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        report = GembaseHMMReport(gene, report_path, self.cfg)

        self.assertFalse(report._hit_start("NOT starting hit"))
        self.assertTrue(report._hit_start(">> starting hit"))


    def test_build_my_db(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        report = GembaseHMMReport(gene, report_path, self.cfg)
        gspD_hmmer_path = self.find_data(os.path.join('hmm', 'gspD.search_hmm.out'))

        db = report._build_my_db(gspD_hmmer_path)
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
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        report = GembaseHMMReport(gene, report_path, self.cfg)
        idx = Indexes(self.cfg)
        macsyfinder_idx = idx.find_my_indexes()
        gspD_hmmer_path = self.find_data(os.path.join('hmm', 'gspD.search_hmm.out'))
        db = report._build_my_db(gspD_hmmer_path)
        report._fill_my_db(macsyfinder_idx, db)
        self.assertDictEqual(db, {'PSAE001c01_031420': (658, 73),
                                  'PSAE001c01_051090': (714, 75),
                                  'PSAE001c01_018920': (776, 71),
                                  'PSAE001c01_043580': (416, 74),
                                  'PSAE001c01_017350': (600, 70),
                                  'PSAE001c01_013980': (759, 69),
                                  'PSAE001c01_026600': (273, 72),
                                  'NC_xxxxx_xx_056141': (803, 141),
                                  'PSAE001c01_006940': (803, 68)})


    def test_parse_hmm_header(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        report = GembaseHMMReport(gene, report_path, self.cfg)

        hmm_hit = [">> NC_xxxxx_xx_056141  C ATG TAA 6260390 6261757 Valid PA5567 1368 _NP_254254.1_ PA5567 1 6260390 6261757 | tRNA modific"]
        hit_id = report._parse_hmm_header(hmm_hit)
        self.assertEqual(hit_id, 'NC_xxxxx_xx_056141')

    def test_parse_hmm_body(self):
        def make_hmm_group(hmm_string):
            hmm_file = StringIO(hmm_string)
            hmm_hits = (x[1] for x in groupby(hmm_file, lambda l: l.startswith('>>')))
            header = hmm_hits.next()
            body = hmm_hits.next()
            return body

        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        report = GembaseHMMReport(gene, report_path, self.cfg)

        # with one significant hit
        hmm = """>> NC_xxxxx_xx_056141  C ATG TAA 6260390 6261757 Valid PA5567 1368 _NP_254254.1_ PA5567 1 6260390 6261757 | tRNA modific
   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
   1 !  779.2   5.5  1.4e-237    2e-236       1     596 []     104     741 ..     104     741 .. 0.93

  Alignments for each domain:
"""
        body = make_hmm_group(hmm)
        hits = report._parse_hmm_body('NC_xxxxx_xx_056141', 596, 803, 0.5, 'NC_xxxxx_xx', 141, 0.5, body)
        expected_hits = [Hit(gene, system, "NC_xxxxx_xx_056141", 803, "NC_xxxxx_xx", 141, float(2e-236), float(779.2),
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
        hits = report._parse_hmm_body('NC_xxxxx_xx_056141', 596, 803, 0.5, 'NC_xxxxx_xx', 141, 0.5, body)
        expected_hits = []
        self.assertListEqual(hits, expected_hits)

        # with no hit
        hmm = """>> PSAE001c01_051090  C ATG TGA 5675714 5677858 Valid pilQ 2145 _PA5040_NP_253727.1_ PA5040 1 5675714 5677858 | type 4 f
        bla bla
        """
        body = make_hmm_group(hmm)
        hits = report._parse_hmm_body('NC_xxxxx_xx_056141', 596, 803, 0.5, 'NC_xxxxx_xx', 141, 0.5, body)
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
            hits = report._parse_hmm_body('NC_xxxxx_xx_056141', 596, 803, 0.5, 'NC_xxxxx_xx', 141, 0.5, body)
        self.assertEqual(str(ctx.exception), """Invalid line to parse :   1 !  779.2   5.5  1.4e-237    foo       1     596 []     104     741 ..     104     741 .. 0.93
:could not convert string to float: foo""")


class TestGembaseHMMReport(TestReport):


    def test_extract(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        shutil.copy(self.find_data("hmm", gene_name + self.cfg.res_search_suffix),
                    self.cfg.working_dir)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        report = GembaseHMMReport(gene, report_path, self.cfg)
        report.extract()
        self.assertEqual(len(report.hits), 6)
        #           gene, system,     hit_id,         hit_seq_ length   replicon_name, pos_hit, i_eval,
        #           score,       profile_coverage, sequence_coverage, begin_match, end_match
        hits = [Hit(gene, system, "NC_xxxxx_xx_056141", 803, "NC_xxxxx_xx", 141, float(2e-236), float(779.2),
                    float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741),
                Hit(gene, system, "PSAE001c01_006940", 803, "PSAE001c01", 68, float(1.2e-234), float(779.2),
                    float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741),
                Hit(gene, system, "PSAE001c01_013980", 759, "PSAE001c01", 69, float(3.7e-76), float(255.8),
                    float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736),
                Hit(gene, system, "PSAE001c01_017350", 600, "PSAE001c01", 70, float(3.2e-27), float(94.2),
                    float(0.500000), (506.0 - 226.0 + 1) / 600,  226, 506),
                Hit(gene, system, "PSAE001c01_018920", 776, "PSAE001c01", 71, float(6.1e-183), float(608.4),
                    float(1.000000), (606.0 - 48.0 + 1) / 776, 48, 606),
                Hit(gene, system, "PSAE001c01_031420", 658, "PSAE001c01", 73, float(1.8e-210), float(699.3),
                    float(1.000000), (614.0 - 55.0 + 1) / 658, 55, 614)
        ]
        self.assertListEqual(hits, report.hits)

        report = GembaseHMMReport(gene, report_path, self.cfg)
        report.hits = hits
        self.assertIsNone(report.extract())


    def test_extract_concurent(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        shutil.copy(self.find_data("hmm", gene_name + self.cfg.res_search_suffix),
                    self.cfg.working_dir)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        reports = []
        for i in range(5):
            report = GembaseHMMReport(gene, report_path, self.cfg)
            reports.append(report)

        import threading

        def worker(report):
            report.extract()

        for report in reports:
            t = threading.Thread(target=worker, args=(report,))
            t.start()
        main_thread = threading.currentThread()
        for t in threading.enumerate():
            if t is main_thread:
                continue
        t.join()

        #          gene, system,     hit_id,        hit_seq_length replicon_name, pos_hit, i_eval,  score,
        #          profile_coverage, sequence_coverage, begin_match, end_match
        hits = [Hit(gene, system, "NC_xxxxx_xx_056141", 803, "NC_xxxxx_xx", 141, float(2e-236), float(779.2),
                    float(1.000000), (741.0 - 104.0 + 1)/803, 104, 741),
                Hit(gene, system, "PSAE001c01_006940", 803, "PSAE001c01", 68, float(1.2e-234), float(779.2),
                    float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741),
                Hit(gene, system, "PSAE001c01_013980", 759, "PSAE001c01", 69, float(3.7e-76), float(255.8),
                    float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736),
                Hit(gene, system, "PSAE001c01_017350", 600, "PSAE001c01", 70, float(3.2e-27), float(94.2),
                    float(0.500000), (506.0 - 226.0 + 1) / 600,  226, 506),
                Hit(gene, system, "PSAE001c01_018920", 776, "PSAE001c01", 71, float(6.1e-183), float(608.4),
                    float(1.000000), (606.0 - 48.0 + 1) / 776, 48, 606),
                Hit(gene, system, "PSAE001c01_031420", 658, "PSAE001c01", 73, float(1.8e-210), float(699.3),
                    float(1.000000), (614.0 - 55.0 + 1) / 658, 55, 614)
                ]
        for report in reports:
            report.save_extract()
            self.assertEqual(len(report.hits), len(hits))
            self.assertListEqual(report.hits, hits)

class TestOrderedHMMReport(TestReport):


    def test_extract(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        shutil.copy(self.find_data("hmm", gene_name + self.cfg.res_search_suffix),
                    self.cfg.working_dir)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        report = OrderedHMMReport(gene, report_path, self.cfg)
        report.extract()
        self.assertEqual(len(report.hits), 6)
        #           gene, system,     hit_id,         hit_seq_ length   replicon_name, pos_hit, i_eval,
        #           score,       profile_coverage, sequence_coverage, begin_match, end_match
        hits = [Hit(gene, system, "NC_xxxxx_xx_056141", 803, RepliconDB.ordered_replicon_name, 141, float(2e-236), float(779.2),
                    float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741),
                Hit(gene, system, "PSAE001c01_006940", 803, RepliconDB.ordered_replicon_name, 68, float(1.2e-234), float(779.2),
                    float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741),
                Hit(gene, system, "PSAE001c01_013980", 759, RepliconDB.ordered_replicon_name, 69, float(3.7e-76), float(255.8),
                    float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736),
                Hit(gene, system, "PSAE001c01_017350", 600, RepliconDB.ordered_replicon_name, 70, float(3.2e-27), float(94.2),
                    float(0.500000), (506.0 - 226.0 + 1) / 600,  226, 506),
                Hit(gene, system, "PSAE001c01_018920", 776, RepliconDB.ordered_replicon_name, 71, float(6.1e-183), float(608.4),
                    float(1.000000), (606.0 - 48.0 + 1) / 776, 48, 606),
                Hit(gene, system, "PSAE001c01_031420", 658, RepliconDB.ordered_replicon_name, 73, float(1.8e-210), float(699.3),
                    float(1.000000), (614.0 - 55.0 + 1) / 658, 55, 614)
        ]
        self.assertListEqual(hits, report.hits)

        report = OrderedHMMReport(gene, report_path, self.cfg)
        report.hits = hits
        self.assertIsNone(report.extract())


    def test_extract_concurent(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        shutil.copy(self.find_data("hmm", gene_name + self.cfg.res_search_suffix),
                    self.cfg.working_dir)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        reports = []
        for i in range(5):
            report = OrderedHMMReport(gene, report_path, self.cfg)
            reports.append(report)

        import threading

        def worker(report):
            report.extract()

        for report in reports:
            t = threading.Thread(target=worker, args=(report,))
            t.start()
        main_thread = threading.currentThread()
        for t in threading.enumerate():
            if t is main_thread:
                continue
        t.join()

        #          gene, system,     hit_id,        hit_seq_length replicon_name, pos_hit, i_eval,  score,
        #          profile_coverage, sequence_coverage, begin_match, end_match
        hits = [Hit(gene, system, "NC_xxxxx_xx_056141", 803, RepliconDB.ordered_replicon_name, 141, float(2e-236), float(779.2),
                    float(1.000000), (741.0 - 104.0 + 1)/803, 104, 741),
                Hit(gene, system, "PSAE001c01_006940", 803, RepliconDB.ordered_replicon_name, 68, float(1.2e-234), float(779.2),
                    float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741),
                Hit(gene, system, "PSAE001c01_013980", 759, RepliconDB.ordered_replicon_name, 69, float(3.7e-76), float(255.8),
                    float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736),
                Hit(gene, system, "PSAE001c01_017350", 600, RepliconDB.ordered_replicon_name, 70, float(3.2e-27), float(94.2),
                    float(0.500000), (506.0 - 226.0 + 1) / 600,  226, 506),
                Hit(gene, system, "PSAE001c01_018920", 776, RepliconDB.ordered_replicon_name, 71, float(6.1e-183), float(608.4),
                    float(1.000000), (606.0 - 48.0 + 1) / 776, 48, 606),
                Hit(gene, system, "PSAE001c01_031420", 658, RepliconDB.ordered_replicon_name, 73, float(1.8e-210), float(699.3),
                    float(1.000000), (614.0 - 55.0 + 1) / 658, 55, 614)
                ]
        for report in reports:
            report.save_extract()
            self.assertEqual(len(report.hits), len(hits))
            self.assertListEqual(report.hits, hits)
