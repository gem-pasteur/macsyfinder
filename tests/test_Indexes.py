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
import platform
import unittest
import shutil
import tempfile
import argparse

from macsypy.config import Config, MacsyDefaults
from macsypy.database import Indexes
from macsypy.error import MacsypyError
from tests import MacsyTest


class Test(MacsyTest):


    def setUp(self):
        args = argparse.Namespace()

        args.db_type = 'gembase'
        args.e_value_res = 1
        args.i_evalue_sel = 0.5
        args.models_dir = self.find_data('models')
        args.res_search_suffix = ''
        args.log_level = 30

        args.out_dir = os.path.join(tempfile.gettempdir(), 'test_macsyfinder_indexes')
        if os.path.exists(args.out_dir):
            shutil.rmtree(os.path.join(tempfile.gettempdir(), 'test_macsyfinder_indexes'))
        os.makedirs(args.out_dir)
        seq_db = self.find_data("base", "test_base.fa")
        shutil.copy(seq_db, args.out_dir)

        args.sequence_db = os.path.join(args.out_dir, os.path.basename(seq_db))

        self.cfg = Config(MacsyDefaults(), args)


    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir())
        except:
            pass

    def test_find_my_indexes(self):
        idx = Indexes(self.cfg)
        self.assertIsNone(idx.find_my_indexes())
        new_idx = os.path.join(os.path.dirname(self.cfg.sequence_db()), idx.name + ".idx")
        with open(new_idx, 'w'):
            pass
        self.assertEqual(idx.find_my_indexes(), new_idx)

    def test_build_no_idx(self):
        idx = Indexes(self.cfg)
        idx.build()
        my_idx = idx.find_my_indexes()
        self.assertEqual(my_idx, os.path.join(os.path.dirname(self.cfg.sequence_db()), idx.name + ".idx"))


    def test_build_with_idx(self):
        idx = Indexes(self.cfg)
        open(os.path.join(os.path.dirname(self.cfg.sequence_db()), idx.name + ".idx"), 'w').close()
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
        idx_dir = os.path.join(os.path.dirname(self.cfg.sequence_db()))
        os.chmod(idx_dir, 0000)
        try:
            with self.assertRaises(IOError) as ctx:
                with self.catch_log():
                    idx.build()
            self.assertRegex(str(ctx.exception),
                             "cannot build indexes, \(.+/test_macsyfinder_indexes\) is not writable")
        finally:
            os.chmod(idx_dir, 0o777)


    def test_build_my_indexes(self):
        args = argparse.Namespace()

        args.db_type = 'gembase'
        args.e_value_res = 1
        args.i_evalue_sel = 0.5
        args.models_dir = self.find_data('models')
        args.res_search_suffix = ''
        args.log_level = 30

        args.out_dir = os.path.join(tempfile.gettempdir(), 'test_macsyfinder_indexes')
        if os.path.exists(args.out_dir):
            shutil.rmtree(os.path.join(tempfile.gettempdir(), 'test_macsyfinder_indexes'))
        os.makedirs(args.out_dir)
        seq_db = self.find_data("base", "test_base_with_errors.fa")
        shutil.copy(seq_db, args.out_dir)
        args.sequence_db = os.path.join(args.out_dir, os.path.basename(seq_db))
        cfg = Config(MacsyDefaults(), args)

        idx = Indexes(cfg)
        with self.assertRaises(MacsypyError) as e:
            with self.catch_log():
                idx._build_my_indexes()
        self.assertTrue(str(e.exception).startswith("unable to index the sequence dataset:"))
