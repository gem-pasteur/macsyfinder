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
import platform
import unittest
import shutil
import tempfile
import logging
from macsypy.config import Config
from macsypy.database import Indexes
from tests import MacsyTest


class Test(MacsyTest):

    def __init__(self, methodName = 'runTest'):
        super(Test, self).__init__(methodName)

        def fake_init(obj, cfg):
            obj.cfg = cfg
            obj._fasta_path = cfg.sequence_db
            obj.name = os.path.basename(cfg.sequence_db)
        self.fake_init = fake_init
        self.real_init = Indexes.__init__

    def setUp(self):
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        
        # add only one handler to the macsypy logger
        from macsypy.database import _log
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

        shutil.copy(self.cfg.sequence_db, self.cfg.working_dir)
        self.cfg.options['sequence_db'] = os.path.join(self.cfg.working_dir, os.path.basename(self.cfg.sequence_db))


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

    def test_find_my_indexes(self):
        idx = Indexes(self.cfg)
        self.assertIsNone(idx.find_my_indexes())
        new_idx = os.path.join(os.path.dirname(self.cfg.sequence_db), idx.name + ".idx")
        open(new_idx, 'w')
        self.assertEqual(idx.find_my_indexes(), new_idx)

    def test_build_no_idx(self):
        idx = Indexes(self.cfg)
        idx.build()
        my_idx = idx.find_my_indexes()
        self.assertEqual(my_idx, os.path.join(os.path.dirname(self.cfg.sequence_db), idx.name + ".idx"))


    def test_build_with_idx(self):
        idx = Indexes(self.cfg)
        open(os.path.join(os.path.dirname(self.cfg.sequence_db), idx.name + ".idx"), 'w').close()
        idx.build()
        my_idx = idx.find_my_indexes()
        self.assertEqual(os.path.getsize(my_idx), 0)


    def test_build_force(self):
        idx = Indexes(self.cfg)
        idx.build(force=True)
        my_idx = idx.find_my_indexes()
        self.assertNotEqual(os.path.getsize(my_idx), 0)


    @unittest.skipIf(platform.system() == 'Windows' or os.getuid() == 0, 'Skip test on Windows or if run as root')
    def test_build_not_writable(self):
        # Skip test on Windows, since setting the folder permissions is not affecting files inside
        # in Singularity container tess are run as root and this test as non sense
        idx = Indexes(self.cfg)
        idx_dir = os.path.join(os.path.dirname(self.cfg.sequence_db))
        os.chmod(idx_dir, 0000)
        self.assertRaises(IOError, idx.build)
        os.chmod(idx_dir, 0777)

