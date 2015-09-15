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
import platform
import logging
from macsypy.report import HMMReport, GeneralHMMReport, GembaseHMMReport, Hit
from macsypy.gene import Gene
from macsypy.system import System
from macsypy.config import Config
from macsypy.database import Indexes
from macsypy.registries import ModelRegistry


class Test(unittest.TestCase):

    _data_dir = os.path.join(os.path.dirname(__file__), "datatest")

    def setUp(self):
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        
        #add only one handler to the macsypy logger
        from macsypy.gene import _log
        macsy_log = _log.parent
        log_file = 'NUL' if platform.system() == 'Windows' else '/dev/null'
        log_handler = logging.FileHandler(log_file)
        macsy_log.addHandler(log_handler)
        
        self.cfg = Config(hmmer_exe="hmmsearch",
                          sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type="gembase",
                          e_value_res=1,
                          i_evalue_sel=0.5,
                          models_dir=os.path.join(self._data_dir, 'models'),
                          res_search_dir=tempfile.gettempdir(),
                          res_search_suffix=".search_hmm.out",
                          profile_suffix=".hmm",
                          res_extract_suffix="",
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
        except:
            pass

    def test_HMMReport(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        shutil.copy(os.path.join(self._data_dir, "hmm", gene_name + self.cfg.res_search_suffix),
                    self.cfg.working_dir)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        self.assertRaises(TypeError, HMMReport, gene, report_path, self.cfg)


    def test_GembaseHMMReport_extract(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        shutil.copy(os.path.join(self._data_dir, "hmm", gene_name + self.cfg.res_search_suffix),
                    self.cfg.working_dir)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        report = GembaseHMMReport(gene, report_path, self.cfg)
        report.extract()
        self.assertEqual(len(report.hits), 5)
                #gene, system,     hit_id,        hit_seq_length replicon_name, pos_hit, i_eval,          score,       profile_coverage, sequence_coverage, begin_match, end_match
        hits=[ Hit(gene, system, "PSAE001c01_006940", 803,"PSAE001c01", 68, float(1.2e-234), float(779.2), float(1.000000), (741.0 - 104.0 + 1)/ 803, 104, 741),
               Hit(gene, system, "PSAE001c01_013980", 759,"PSAE001c01", 69, float(3.7e-76), float(255.8), float(1.000000), (736.0 - 105.0 + 1)/ 759, 105, 736),
               Hit(gene, system, "PSAE001c01_017350", 600,"PSAE001c01", 70, float(3.2e-27), float(94.2), float(0.500000), (506.0 - 226.0 + 1)/ 600,  226, 506),
               Hit(gene, system, "PSAE001c01_018920", 776,"PSAE001c01", 71, float(6.1e-183), float(608.4), float(1.000000), (606.0 - 48.0 + 1)/ 776, 48, 606),
               Hit(gene, system, "PSAE001c01_031420", 658,"PSAE001c01", 73, float(1.8e-210), float(699.3), float(1.000000), (614.0 - 55.0 + 1)/ 658, 55, 614)
        ]
        self.assertListEqual(hits, report.hits)


    def test_GembaseHMMReport_extract_concurent(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        shutil.copy(os.path.join(self._data_dir, "hmm", gene_name + self.cfg.res_search_suffix),
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

                        #gene, system,     hit_id,        hit_seq_length replicon_name, pos_hit, i_eval,          score,       profile_coverage, sequence_coverage, begin_match, end_match
        hits = [Hit(gene, system, "PSAE001c01_006940", 803,"PSAE001c01", 68, float(1.2e-234), float(779.2), float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741),
                Hit(gene, system, "PSAE001c01_013980", 759,"PSAE001c01", 69, float(3.7e-76), float(255.8), float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736),
                Hit(gene, system, "PSAE001c01_017350", 600,"PSAE001c01", 70, float(3.2e-27), float(94.2), float(0.500000), (506.0 - 226.0 + 1) / 600,  226, 506),
                Hit(gene, system, "PSAE001c01_018920", 776,"PSAE001c01", 71, float(6.1e-183), float(608.4), float(1.000000), (606.0 - 48.0 + 1) / 776, 48, 606),
                Hit(gene, system, "PSAE001c01_031420", 658,"PSAE001c01", 73, float(1.8e-210), float(699.3), float(1.000000), (614.0 - 55.0 + 1) / 658, 55, 614)
               ]
        for report in reports:
            report.save_extract()
            self.assertEqual(len(report.hits), len(hits))
            self.assertListEqual(hits, report.hits)


    def test_str(self):
        system = System(self.cfg, "T2SS", 10)
        gene_name = "gspD"
        gene = Gene(self.cfg, gene_name, system, self.models_location)
        shutil.copy(os.path.join(self._data_dir, "hmm", gene_name + self.cfg.res_search_suffix),
                    self.cfg.working_dir)
        report_path = os.path.join(self.cfg.working_dir, gene_name + self.cfg.res_search_suffix)
        report = GembaseHMMReport(gene, report_path, self.cfg)
        report.extract()
        hits = [Hit(gene, system, "PSAE001c01_006940", 803,"PSAE001c01", 68, float(1.2e-234), float(779.2), float(1.000000), (741.0 - 104.0 + 1) / 803, 104, 741),
                Hit(gene, system, "PSAE001c01_013980", 759,"PSAE001c01", 69, float(3.7e-76), float(255.8), float(1.000000), (736.0 - 105.0 + 1) / 759, 105, 736),
                Hit(gene, system, "PSAE001c01_017350", 600,"PSAE001c01", 70, float(3.2e-27), float(94.2), float(0.500000), (506.0 - 226.0 + 1) / 600,  226, 506),
                Hit(gene, system, "PSAE001c01_018920", 776,"PSAE001c01", 71, float(6.1e-183), float(608.4), float(1.000000), (606.0 - 48.0 + 1) / 776, 48, 606),
                Hit(gene, system, "PSAE001c01_031420", 658,"PSAE001c01", 73, float(1.8e-210), float(699.3), float(1.000000), (614.0 - 55.0 + 1) / 658, 55, 614)
               ]
        s = "# gene: {0} extract from {1} hmm output\n".format(gene.name, report_path)
        s += "# profile length= {0:d}\n".format(len(gene.profile))
        s += "# i_evalue threshold= {0:.3f}\n".format(self.cfg.i_evalue_sel)
        s += "# coverage threshold= {0:.3f}\n".format(self.cfg.coverage_profile)
        s += "# hit_id replicon_name position_hit hit_sequence_length gene_name gene_system i_eval score profile_coverage sequence_coverage begin end\n"
        for h in hits:
            s += str(h)
        self.assertEqual(str(report), s)


if __name__ == "__main__":
    unittest.main()
