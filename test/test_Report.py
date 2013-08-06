# -*- coding: utf-8 -*-

#===============================================================================
# Created on Nov 30, 2012
# 
# @author: bneron
# @contact: user_email
# @organization: organization_name
# @license: license
#===============================================================================


import sys
import os

TXSSCAN_HOME = os.path.abspath('..')
if not TXSSCAN_HOME in sys.path: 
    sys.path.append(os.path.abspath('..') )

import unittest
import shutil
from txsscanlib.report import HMMReport, GeneralHMMReport, GembaseHMMReport, Hit
from txsscanlib.gene import Gene
from txsscanlib.system import System
from txsscanlib.config import Config
from txsscanlib.database import Indexes


class Test(unittest.TestCase):

    _data_dir = "./datatest/res_search" 
    
    def setUp(self):
        self.cfg = Config( hmmer_exe = "hmmsearch",
                           sequence_db = "./datatest/prru_psae.001.c01.fasta",
                           db_type = "gembase",
                           e_value_res = 1,
                           i_evalue_sel = 0.5,
                           def_dir = "../data/DEF",
                           res_search_dir = '.',
                           res_search_suffix = ".search_hmm.out",
                           profile_dir = "../data/profiles",
                           profile_suffix = ".fasta-aln_edit.hmm",
                           res_extract_suffix = "",
                           log_level = 30,
                           log_file = '/dev/null'
                           )
        idx = Indexes(self.cfg)
        txsscan_idx = idx.find_my_indexes()
        if not txsscan_idx:
            idx._build_my_indexes()
        
    def tearDown(self):
        shutil.rmtree(self.cfg.working_dir)
        pass
    
    def test_HMMReport(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system)
        shutil.copy(os.path.join(self._data_dir, gene_name + self.cfg.res_search_suffix), self.cfg.working_dir)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        self.assertRaises(TypeError, HMMReport, gene, report_path, self.cfg)

    def test_GembaseHMMReport_extract(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system)
        shutil.copy(os.path.join(self._data_dir, gene_name + self.cfg.res_search_suffix), self.cfg.working_dir)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        report = GembaseHMMReport(gene, report_path, self.cfg)
        report.extract()
        self.assertEqual(len(report.hits), 5)
                #gene, system,     hit_id,        hit_seq_length replicon_name, pos_hit, i_eval,          score,       profile_coverage, sequence_coverage, begin_match, end_match
        hits=[ Hit(gene, system, "PSAE001c01_006940", 803,"PSAE001c01", 3450, float(1.2e-234), float(779.2), float(1.000000), (741.0 - 104.0 + 1)/ 803, 104, 741),
               Hit(gene, system, "PSAE001c01_013980", 759,"PSAE001c01", 4146, float(3.7e-76), float(255.8), float(1.000000), (736.0 - 105.0 + 1)/ 759, 105, 736),
               Hit(gene, system, "PSAE001c01_017350", 600,"PSAE001c01", 4482, float(3.2e-27), float(94.2), float(0.500000), (506.0 - 226.0 + 1)/ 600,  226, 506),
               Hit(gene, system, "PSAE001c01_018920", 776,"PSAE001c01", 4634, float(6.1e-183), float(608.4), float(1.000000), (606.0 - 48.0 + 1)/ 776, 48, 606),
               Hit(gene, system, "PSAE001c01_031420", 658,"PSAE001c01", 5870, float(1.8e-210), float(699.3), float(1.000000), (614.0 - 55.0 + 1)/ 658, 55, 614)
        ]
        self.assertListEqual(hits, report.hits)

    def test_GembaseHMMReport_extract_concurent(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system)
        shutil.copy(os.path.join(self._data_dir, gene_name + self.cfg.res_search_suffix), self.cfg.working_dir)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        reports = []
        for i in range(5):
            report = GembaseHMMReport(gene, report_path, self.cfg)
            reports.append(report)

        import threading

        def worker(report):
            report.extract()

        for report in reports:
            t = threading.Thread(target = worker, args = (report,))
            t.start()
        main_thread = threading.currentThread()
        for t in threading.enumerate():
            if t is main_thread:
                continue
        t.join()

                        #gene, system,     hit_id,        hit_seq_length replicon_name, pos_hit, i_eval,          score,       profile_coverage, sequence_coverage, begin_match, end_match
        hits=[ Hit(gene, system, "PSAE001c01_006940", 803,"PSAE001c01", 3450, float(1.2e-234), float(779.2), float(1.000000), (741.0 - 104.0 + 1)/ 803, 104, 741),
               Hit(gene, system, "PSAE001c01_013980", 759,"PSAE001c01", 4146, float(3.7e-76), float(255.8), float(1.000000), (736.0 - 105.0 + 1)/ 759, 105, 736),
               Hit(gene, system, "PSAE001c01_017350", 600,"PSAE001c01", 4482, float(3.2e-27), float(94.2), float(0.500000), (506.0 - 226.0 + 1)/ 600,  226, 506),
               Hit(gene, system, "PSAE001c01_018920", 776,"PSAE001c01", 4634, float(6.1e-183), float(608.4), float(1.000000), (606.0 - 48.0 + 1)/ 776, 48, 606),
               Hit(gene, system, "PSAE001c01_031420", 658,"PSAE001c01", 5870, float(1.8e-210), float(699.3), float(1.000000), (614.0 - 55.0 + 1)/ 658, 55, 614)
        ]
        for report in reports:
            report.save_extract()
            self.assertEqual(len(report.hits), len(hits))
            self.assertListEqual(hits, report.hits)

    def test_str(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system)
        shutil.copy(os.path.join(self._data_dir, gene_name + self.cfg.res_search_suffix), self.cfg.working_dir)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        report = GembaseHMMReport(gene, report_path, self.cfg)
        report.extract()
        hits=[ Hit(gene, system, "PSAE001c01_006940", 803,"PSAE001c01", 3450, float(1.2e-234), float(779.2), float(1.000000), (741.0 - 104.0 + 1)/ 803, 104, 741),
               Hit(gene, system, "PSAE001c01_013980", 759,"PSAE001c01", 4146, float(3.7e-76), float(255.8), float(1.000000), (736.0 - 105.0 + 1)/ 759, 105, 736),
               Hit(gene, system, "PSAE001c01_017350", 600,"PSAE001c01", 4482, float(3.2e-27), float(94.2), float(0.500000), (506.0 - 226.0 + 1)/ 600,  226, 506),
               Hit(gene, system, "PSAE001c01_018920", 776,"PSAE001c01", 4634, float(6.1e-183), float(608.4), float(1.000000), (606.0 - 48.0 + 1)/ 776, 48, 606),
               Hit(gene, system, "PSAE001c01_031420", 658,"PSAE001c01", 5870, float(1.8e-210), float(699.3), float(1.000000), (614.0 - 55.0 + 1)/ 658, 55, 614)
        ]
        s = ""
        s = "# gene: %s extract from %s hmm output\n" % (gene.name, report_path)
        s += "# profile length= %d\n" % len(gene.profile)
        s += "# i_evalue threshold= %f\n" % self.cfg.i_evalue_sel
        s += "# coverage threshold= %f\n" % self.cfg.coverage_profile
        s += "# hit_id replicon_name position_hit hit_sequence_length gene_name gene_system i_eval score profile_coverage sequence_coverage begin end\n"
        for h in hits:
            s += str(h)
        self.assertEqual(str(report), s)

if __name__ == "__main__":
    unittest.main()