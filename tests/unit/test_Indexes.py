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
from macsypy.config import Config
from macsypy.database import Indexes
from macsypy.utils import which
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
        log_file = 'NUL' if platform.system() == 'Windows' else '/dev/null'
        log_handler = logging.FileHandler(log_file)
        macsy_log.addHandler(log_handler)
        
        self.cfg = Config(hmmer_exe="hmmsearch",
                          sequence_db=os.path.join(self._data_dir, "base", "test_base.fa"),
                          db_type="gembase",
                          e_value_res=1,
                          i_evalue_sel=0.5,
                          def_dir=os.path.join(self._data_dir, "DEF"),
                          res_search_dir=tempfile.gettempdir(),
                          res_search_suffix=".search_hmm.out",
                          profile_dir=os.path.join(self._data_dir, "profiles"),
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

    def test_find_hmmer_indexes_no_files(self):
        idx = Indexes(self.cfg)
        # tester pas de fichier
        hmmer_idx = idx.find_hmmer_indexes()
        self.assertListEqual(hmmer_idx, [])

    def test_find_hmmer_indexes_all_files(self):
        idx = Indexes(self.cfg)
        suffixes = ('.phr', '.pin', '.psd', '.psi', '.psq')
        files_2_find = []
        for s in suffixes:
            new_idx = os.path.join(self.cfg.sequence_db + s)
            open(new_idx, 'w')
            files_2_find.append(new_idx)
        hmmer_idx = idx.find_hmmer_indexes()
        self.assertListEqual(hmmer_idx, files_2_find)


    def test_find_hmmer_indexes_all_files_and_pal(self):
        idx = Indexes(self.cfg)
        # tester tous les fichiers + pal
        suffixes = ('.phr', '.pin', '.psd', '.psi', '.psq', '.pal')
        for s in suffixes:
            new_idx = os.path.join(self.cfg.sequence_db + s)
            open(new_idx, 'w')
        self.assertRaises(RuntimeError, idx.find_hmmer_indexes)


    def test_find_hmmer_indexes_some_files(self):
        idx = Indexes(self.cfg)
        # tester pas tous les fichiers
        suffixes = ('.phr', '.pin', '.psd', '.psi')
        for s in suffixes:
            new_idx = os.path.join(self.cfg.sequence_db + s)
            open(new_idx, 'w')
        self.assertRaises(RuntimeError, idx.find_hmmer_indexes)


    def test_find_hmmer_indexes_lack_pal(self):
        idx = Indexes(self.cfg)
        # tester plusieurs index pas de pal
        suffixes = ('.phr', '.pin', '.psd', '.psi', '.psq')
        for s in suffixes:
            for i in range(2):
                new_idx = os.path.join(self.cfg.sequence_db + str(i) + s)
                open(new_idx, 'w')
        self.assertRaises(RuntimeError, idx.find_hmmer_indexes)


    def test_find_hmmer_indexes_all_files_and_2virtual(self):
        idx = Indexes(self.cfg)
        # tester 1 fichier index + pal
        suffixes = ('.phr', '.pin', '.psd', '.psi', '.psq', '.pal')
        files_2_find = []
        for s in suffixes:
            for i in range(2):
                new_idx = os.path.join(self.cfg.sequence_db + str(i) + s)
                open(new_idx, 'w')
                files_2_find.append(new_idx)
        self.assertRaises(RuntimeError, idx.find_hmmer_indexes)


    def test_find_hmmer_indexes_all_files_and_virtual(self):
        idx = Indexes(self.cfg)
        # tester index + pal
        suffixes = ('.phr', '.pin', '.psd', '.psi', '.psq')
        files_2_find = []
        for s in suffixes:
            for i in range(2):
                new_idx = os.path.join("{0}.{1:d}.{2}".format(self.cfg.sequence_db, i, s))
                open(new_idx, 'w')
                files_2_find.append(new_idx)
        new_idx = os.path.join(self.cfg.sequence_db + '.pal')
        open(new_idx, 'w')
        files_2_find.append(new_idx)
        files_2_find.sort()
        hmmer_idx = idx.find_hmmer_indexes()
        hmmer_idx.sort()
        self.assertListEqual(hmmer_idx, files_2_find)


    def test_find_my_indexes(self):
        idx = Indexes(self.cfg)
        self.assertIsNone(idx.find_my_indexes())
        new_idx = os.path.join(os.path.dirname(self.cfg.sequence_db), idx.name + ".idx")
        open(new_idx, 'w')
        self.assertEqual(idx.find_my_indexes(), new_idx)

    @unittest.skipIf(not (which('makeblastdb') or which('formatdb')), 'neither makeblast nor formatdb found in PATH')
    def test_build_no_idx(self):
        if not which('makeblastdb') and which('formatdb'):
            self.cfg.options['index_db_exe'] = 'formatdb'
        idx = Indexes(self.cfg)
        idx.build()
        my_idx = idx.find_my_indexes()
        hmmer_idx = idx.find_hmmer_indexes()
        self.assertEqual(my_idx, os.path.join(os.path.dirname(self.cfg.sequence_db), idx.name + ".idx"))
        self.assertEqual(hmmer_idx, [self.cfg.sequence_db + suffix for suffix in ('.phr', '.pin', '.psd', '.psi', '.psq')])
        

    @unittest.skipIf(not (which('makeblastdb') or which('formatdb')), 'neither makeblast nor formatdb found in PATH')
    def test_build_with_idx(self):
        if not which('makeblastdb') and which('formatdb'):
            self.cfg.options['index_db_exe'] = 'formatdb'
        # put fake hmmer indexes
        suffixes = ('.phr', '.pin', '.psd', '.psi', '.psq')
        for s in suffixes:
            new_idx = os.path.join(self.cfg.sequence_db + s)
            open(new_idx, 'w')
        idx = Indexes(self.cfg)
        new_idx = open(os.path.join( os.path.dirname(self.cfg.sequence_db), idx.name + ".idx"), 'w')
        idx.build()
        my_idx = idx.find_my_indexes()
        hmmer_idx = idx.find_hmmer_indexes()
        for f in hmmer_idx + [my_idx]:
            self.assertEqual(os.path.getsize(f), 0)

    @unittest.skipIf(not (which('makeblastdb') or which('formatdb')), 'neither makeblast nor formatdb found in PATH')
    def test_build_force(self):
        # put fake hmmer indexes
        if not which('makeblastdb') and which('formatdb'):
            self.cfg.options['index_db_exe'] = 'formatdb'
       
        suffixes = ('.phr', '.pin', '.psd', '.psi', '.psq')
        for s in suffixes:
            new_idx = os.path.join( self.cfg.sequence_db + s)
            open(new_idx, 'w')
        idx = Indexes(self.cfg)
        idx.build(force=True)
        my_idx = idx.find_my_indexes()
        hmmer_idx = idx.find_hmmer_indexes()
        for f in hmmer_idx + [my_idx]:
            self.assertNotEqual(os.path.getsize(f), 0)
            
    
    def test_build_not_writable(self):
        idx = Indexes(self.cfg)
        idx_dir = os.path.join(os.path.dirname(self.cfg.sequence_db))

        # Skip test on Windows, since setting the folder permissions is not affecting files inside
        if platform.system() != 'Windows':
            os.chmod(idx_dir, 0000)
            self.assertRaises(IOError, idx.build)
            os.chmod(idx_dir, 0777)

