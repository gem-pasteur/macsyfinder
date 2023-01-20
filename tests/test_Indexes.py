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


import os
import platform
import time
import unittest
import shutil
import tempfile
import argparse

from macsypy.config import Config, MacsyDefaults
from macsypy.database import Indexes, fasta_iter
from macsypy.error import MacsypyError
from tests import MacsyTest

class TestFastaIter(MacsyTest):

    def setUp(self) -> None:
        self.sequences = [
            ('seq1', 'comment 1', 'ATGCATGC' ),
            ('seq2', 'comment 2', 'GGGGCCCCTT')
        ]
        self.tmpdir = os.path.join(tempfile.gettempdir(), 'test_macsyfinder_fasta_iter')
        if os.path.exists(self.tmpdir):
            shutil.rmtree(self.tmpdir)
        os.makedirs(self.tmpdir)

    def tearDown(self):
        try:
            shutil.rmtree(self.tmpdir)
        except:
            pass

    def test_fasta_iter(self):
        fasta_path = os.path.join(self.tmpdir, "sequence.fa")
        with open(fasta_path, 'w') as fasta:
            for id_, comment, seq in self.sequences:
                fasta.write(f">{id_} {comment}\n{seq}\n")

        with open(fasta_path) as fasta:
            f_iter = fasta_iter(fasta)
            for seq_read, seq_expected in zip(f_iter, self.sequences):
                self.assertEqual(seq_read[0], seq_expected[0])
                self.assertEqual(seq_read[1], seq_expected[1])
                self.assertEqual(seq_read[2], len(seq_expected[2]))

    def test_fasta_iter_bad_fasta(self):
        fasta_path = os.path.join(self.tmpdir, "sequence.fa")
        with open(fasta_path, 'w') as fasta:
            first = True
            for id_, comment, seq in self.sequences:
                if first:
                    fasta.write(f" >{id_} {comment}\n{seq}\n")
                    first = False
                else:
                    fasta.write(f">{id_} {comment}\n{seq}\n")

        with open(fasta_path) as fasta:
            f_iter = fasta_iter(fasta)
            with self.assertRaises(MacsypyError) as ctx:
                with self.catch_log():
                    for seq in f_iter:
                        pass
            self.assertEqual(str(ctx.exception),
                             f"Error during sequence '{fasta_path}' parsing: Check the fasta format.")

class TestIndex(MacsyTest):

    def setUp(self):
        args = argparse.Namespace()
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.out_dir = os.path.join(tempfile.gettempdir(), 'test_macsyfinder_indexes')
        if os.path.exists(args.out_dir):
            shutil.rmtree(args.out_dir)
        os.makedirs(args.out_dir)
        seq_db = self.find_data("base", "test_1.fasta")
        shutil.copy(seq_db, args.out_dir)

        args.index_dir = args.out_dir
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
        my_idx =idx.build()
        self.assertEqual(my_idx, os.path.join(os.path.dirname(self.cfg.sequence_db()), idx.name + ".idx"))

    def test_build_idx_older_than_fasta(self):
        # test if idx.build
        # if index is present and newer than fasta the index are not rebuild
        # if index is present but older than fasta the index is rebuild
        fasta_path = self.cfg.sequence_db()
        idx = Indexes(self.cfg)
        idx_path = idx.build()
        first_build_stamp = os.path.getmtime(idx_path)
        idx_path = idx.build()
        second_build_stamp = os.path.getmtime(idx_path)
        time.sleep(.2)
        f = open(fasta_path, 'a')
        f.write('\n')
        f.close()
        idx_path = idx.build()
        third_build_stamp = os.path.getmtime(idx_path)
        self.assertEqual(first_build_stamp, second_build_stamp)
        self.assertGreater(third_build_stamp, second_build_stamp)

    def test_build_idx_point_wrong_fasta(self):
        # test if idx.build
        # if index is present and newer than fasta the index are not rebuild
        # if index is present but older than fasta the index is rebuild
        fasta_path = self.cfg.sequence_db()
        idx = Indexes(self.cfg)
        idx_path = idx.build()
        first_build_stamp = os.path.getmtime(idx_path)
        with open(idx_path, 'r') as f:
            idx_content = f.readlines()
        first_header = idx_content[0]
        bad_header = f"{first_header[:-1]}_fake\n"
        idx_content[0] = bad_header
        time.sleep(.2)
        with open(idx_path, 'w') as f:
            f.write(''.join(idx_content))
        with self.catch_log():
            idx_path = idx.build()
        second_build_stamp = os.path.getmtime(idx_path)
        self.assertGreater(second_build_stamp, first_build_stamp)
        with open(idx_path, 'r') as f:
            new_header = f.readline()
        self.assertEqual(new_header, first_header)


    def test_build_with_idx(self):
        idx = Indexes(self.cfg)
        # case new style idx
        with open(os.path.join(os.path.dirname(self.cfg.sequence_db()), idx.name + ".idx"), 'w') as idx_file:
            idx_content_new = f"{self.cfg.sequence_db()}\nVICH001.B.00001.C001_01359{idx._field_separator}200{idx._field_separator}1\n"
            idx_file.write(idx_content_new)
        my_idx = idx.build()
        self.assertEqual(os.path.getsize(idx_file.name), len(idx_content_new))

        # case old style no path as first line
        idx_path = os.path.join(os.path.dirname(self.cfg.sequence_db()), idx.name + ".idx")
        with open(idx_path, 'w') as idx_file:
            idx_content_old = "VICH001.B.00001.C001_01359;200;1\n"
            idx_file.write(idx_content_old)
        with self.catch_log(log_name='macsypy') as log:
            _ = idx.build()
            log_msg = log.get_value().strip()
        self.assertEqual(log_msg,
                         f"The '{idx_path}' index file is in old format. Force index building.")

        # case old style bad separator
        idx_path = os.path.join(os.path.dirname(self.cfg.sequence_db()), idx.name + ".idx")
        with open(idx_path, 'w') as idx_file:
            idx_content_old = f"{self.cfg.sequence_db()}\nVICH001.B.00001.C001_01359;200;1\n"
            idx_file.write(idx_content_old)
        with self.catch_log(log_name='macsypy') as log:
            _ = idx.build()
            log_msg = log.get_value().strip()
        self.assertEqual(log_msg,
                         f"The '{idx_path}' index file is in old format. Force index building.")

        # case idx seems valid read it
        with open(os.path.join(os.path.dirname(self.cfg.sequence_db()), idx.name + ".idx")) as idx_file_test:
            data = idx_file_test.read()

        new_content = f"""{self.cfg.sequence_db()}
VICH001.B.00001.C001_01359{idx._field_separator}200{idx._field_separator}1
VICH001.B.00001.C001_01360{idx._field_separator}484{idx._field_separator}2
VICH001.B.00001.C001_01361{idx._field_separator}406{idx._field_separator}3
VICH001.B.00001.C001_01390{idx._field_separator}326{idx._field_separator}4
VICH001.B.00001.C001_01391{idx._field_separator}54{idx._field_separator}5
VICH001.B.00001.C001_01392{idx._field_separator}206{idx._field_separator}6
VICH001.B.00001.C001_01393{idx._field_separator}477{idx._field_separator}7
VICH001.B.00001.C001_01394{idx._field_separator}126{idx._field_separator}8
VICH001.B.00001.C001_01395{idx._field_separator}405{idx._field_separator}9
VICH001.B.00001.C001_01396{idx._field_separator}572{idx._field_separator}10
VICH001.B.00001.C001_01397{idx._field_separator}721{idx._field_separator}11
VICH001.B.00001.C001_01398{idx._field_separator}467{idx._field_separator}12
VICH001.B.00001.C001_01399{idx._field_separator}720{idx._field_separator}13
VICH001.B.00001.C001_01400{idx._field_separator}559{idx._field_separator}14
VICH001.B.00001.C001_01401{idx._field_separator}153{idx._field_separator}15
VICH001.B.00001.C001_01402{idx._field_separator}4558{idx._field_separator}16
VICH001.B.00001.C001_01500{idx._field_separator}120{idx._field_separator}17
VICH001.B.00001.C001_01501{idx._field_separator}344{idx._field_separator}18
VICH001.B.00001.C001_01502{idx._field_separator}478{idx._field_separator}19
VICH001.B.00001.C001_01503{idx._field_separator}724{idx._field_separator}20
VICH001.B.00001.C001_01504{idx._field_separator}309{idx._field_separator}21
VICH001.B.00001.C001_01505{idx._field_separator}390{idx._field_separator}22
VICH001.B.00001.C001_01506{idx._field_separator}419{idx._field_separator}23
VICH001.B.00001.C001_01540{idx._field_separator}353{idx._field_separator}24
VICH001.B.00001.C001_01541{idx._field_separator}229{idx._field_separator}25
VICH001.B.00001.C001_01542{idx._field_separator}267{idx._field_separator}26
VICH001.B.00001.C001_01543{idx._field_separator}328{idx._field_separator}27
VICH001.B.00001.C001_01544{idx._field_separator}258{idx._field_separator}28
VICH001.B.00001.C001_01545{idx._field_separator}228{idx._field_separator}29
VICH001.B.00001.C001_01546{idx._field_separator}538{idx._field_separator}30
VICH001.B.00001.C001_01547{idx._field_separator}77{idx._field_separator}31
VICH001.B.00001.C001_01548{idx._field_separator}476{idx._field_separator}32
VICH001.B.00001.C001_01549{idx._field_separator}324{idx._field_separator}33
VICH001.B.00001.C001_01550{idx._field_separator}387{idx._field_separator}34
VICH001.B.00001.C001_01551{idx._field_separator}382{idx._field_separator}35
VICH001.B.00001.C001_01552{idx._field_separator}149{idx._field_separator}36
VICH001.B.00001.C001_01553{idx._field_separator}319{idx._field_separator}37
VICH001.B.00001.C001_01554{idx._field_separator}237{idx._field_separator}38
VICH001.B.00001.C001_01555{idx._field_separator}74{idx._field_separator}39
VICH001.B.00001.C001_01556{idx._field_separator}362{idx._field_separator}40
VICH001.B.00001.C001_01557{idx._field_separator}170{idx._field_separator}41
VICH001.B.00001.C001_01558{idx._field_separator}77{idx._field_separator}42
VICH001.B.00001.C001_01559{idx._field_separator}296{idx._field_separator}43
VICH001.B.00001.C001_01560{idx._field_separator}405{idx._field_separator}44
VICH001.B.00001.C001_01561{idx._field_separator}182{idx._field_separator}45
VICH001.B.00001.C001_01562{idx._field_separator}445{idx._field_separator}46
VICH001.B.00001.C001_01563{idx._field_separator}212{idx._field_separator}47
VICH001.B.00001.C001_01564{idx._field_separator}387{idx._field_separator}48
VICH001.B.00001.C001_01565{idx._field_separator}414{idx._field_separator}49
"""

        self.assertEqual(data, new_content)




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
            self.assertEqual(f"The '{idx_dir}' dir is not writable.",
                             str(ctx.exception))
        finally:
            os.chmod(idx_dir, 0o777)

    @unittest.skipIf(platform.system() == 'Windows' or os.getuid() == 0, 'Skip test on Windows or if run as root')
    def test_index_dir(self):
        # case index-dir is not specify sequence-db dir is not writable
        args = argparse.Namespace()
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.out_dir = os.path.join(tempfile.gettempdir(), 'test_macsyfinder_indexes')
        args.sequence_db = os.path.join(args.out_dir, os.path.basename(self.cfg.sequence_db()))
        cfg = Config(MacsyDefaults(), args)
        idx = Indexes(cfg)
        index_dir = idx._index_dir(build=False)
        expc_idx_dir = os.path.dirname(cfg.sequence_db())
        self.assertEqual(index_dir, expc_idx_dir)
        try:
            os.chmod(index_dir, 0000)
            with self.assertRaises(ValueError) as ctx:
                _ = idx._index_dir(build=True)
            self.assertEqual(f"The '{index_dir}' dir is not writable. Change rights or specify --index-dir.",
                             str(ctx.exception))
        finally:
            os.chmod(index_dir, 0o777)


        args = argparse.Namespace()
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        args.out_dir = os.path.join(tempfile.gettempdir(), 'test_macsyfinder_indexes')
        args.index_dir = os.path.join(args.out_dir, 'index_dir')
        args.sequence_db = os.path.join(args.out_dir, os.path.basename(self.cfg.sequence_db()))
        cfg = Config(MacsyDefaults(), args)
        idx = Indexes(cfg)

        # case --index-dir does not exists
        with self.assertRaises(ValueError) as ctx:
            _ = idx._index_dir(build=False)
        self.assertEqual(str(ctx.exception),
                         f"No such directory: {args.index_dir}")

        # case --index-dir is not writable
        os.makedirs(args.index_dir)
        os.chmod(args.index_dir, 0000)
        try:
            # but I do not care I only read
            index_dir = idx._index_dir(build=False)
            self.assertEqual(index_dir, args.index_dir)

            # it's important to build indexes
            with self.assertRaises(ValueError) as ctx:
                _ = idx._index_dir(build=True)
            self.assertEqual(str(ctx.exception),
                             f"The '{index_dir}' dir is not writable.")
        finally:
            os.chmod(args.index_dir, 0o777)

        # case the sequence_db value is just a filename not a path
        current_dir = os.getcwd()
        try:
            os.chdir(args.out_dir)
            args = argparse.Namespace()
            args.db_type = 'gembase'
            args.models_dir = self.find_data('models')
            args.out_dir = os.path.join(tempfile.gettempdir(), 'test_macsyfinder_indexes')
            args.sequence_db = os.path.basename(self.cfg.sequence_db())
            cfg = Config(MacsyDefaults(), args)
            idx = Indexes(cfg)
            idx_dir = idx._index_dir(build=True)
            self.assertEqual(idx_dir, os.getcwd())
        finally:
            os.chdir(current_dir)

    def test_build_my_indexes(self):
        args = argparse.Namespace()
        args.db_type = 'gembase'

        args.out_dir = os.path.join(tempfile.gettempdir(), 'test_macsyfinder_indexes')
        if os.path.exists(args.out_dir):
            shutil.rmtree(os.path.join(tempfile.gettempdir(), 'test_macsyfinder_indexes'))
        os.makedirs(args.out_dir)
        seq_db = self.find_data("base", "test_base_with_errors.fa")
        shutil.copy(seq_db, args.out_dir)
        args.sequence_db = os.path.join(args.out_dir, os.path.basename(seq_db))
        self.cfg = Config(MacsyDefaults(), args)

        idx = Indexes(self.cfg)
        with self.assertRaises(MacsypyError) as e:
            # the directory for index exist and is writable but
            # the sequence file is corrupted and cannot be read correctly
            with self.catch_log():
                idx._build_my_indexes(args.out_dir)
        self.assertTrue(str(e.exception).startswith("unable to index the sequence dataset:"))


    def test_iter(self):
        idx = Indexes(self.cfg)
        with self.assertRaises(MacsypyError) as ctx:
            next(iter(idx))

        self.assertEqual(str(ctx.exception),
                         'Build index before to use it.')

        idx.build()
        expected_idx = [('VICH001.B.00001.C001_01359', 200, 1), ('VICH001.B.00001.C001_01360', 484, 2),
                        ('VICH001.B.00001.C001_01361', 406, 3), ('VICH001.B.00001.C001_01390', 326, 4),
                        ('VICH001.B.00001.C001_01391', 54, 5), ('VICH001.B.00001.C001_01392', 206, 6),
                        ('VICH001.B.00001.C001_01393', 477, 7), ('VICH001.B.00001.C001_01394', 126, 8),
                        ('VICH001.B.00001.C001_01395', 405, 9), ('VICH001.B.00001.C001_01396', 572, 10),
                        ('VICH001.B.00001.C001_01397', 721, 11), ('VICH001.B.00001.C001_01398', 467, 12),
                        ('VICH001.B.00001.C001_01399', 720, 13), ('VICH001.B.00001.C001_01400', 559, 14),
                        ('VICH001.B.00001.C001_01401', 153, 15), ('VICH001.B.00001.C001_01402', 4558, 16),
                        ('VICH001.B.00001.C001_01500', 120, 17), ('VICH001.B.00001.C001_01501', 344, 18),
                        ('VICH001.B.00001.C001_01502', 478, 19), ('VICH001.B.00001.C001_01503', 724, 20),
                        ('VICH001.B.00001.C001_01504', 309, 21), ('VICH001.B.00001.C001_01505', 390, 22),
                        ('VICH001.B.00001.C001_01506', 419, 23), ('VICH001.B.00001.C001_01540', 353, 24),
                        ('VICH001.B.00001.C001_01541', 229, 25), ('VICH001.B.00001.C001_01542', 267, 26),
                        ('VICH001.B.00001.C001_01543', 328, 27), ('VICH001.B.00001.C001_01544', 258, 28),
                        ('VICH001.B.00001.C001_01545', 228, 29), ('VICH001.B.00001.C001_01546', 538, 30),
                        ('VICH001.B.00001.C001_01547', 77, 31), ('VICH001.B.00001.C001_01548', 476, 32),
                        ('VICH001.B.00001.C001_01549', 324, 33), ('VICH001.B.00001.C001_01550', 387, 34),
                        ('VICH001.B.00001.C001_01551', 382, 35), ('VICH001.B.00001.C001_01552', 149, 36),
                        ('VICH001.B.00001.C001_01553', 319, 37), ('VICH001.B.00001.C001_01554', 237, 38),
                        ('VICH001.B.00001.C001_01555', 74, 39), ('VICH001.B.00001.C001_01556', 362, 40),
                        ('VICH001.B.00001.C001_01557', 170, 41), ('VICH001.B.00001.C001_01558', 77, 42),
                        ('VICH001.B.00001.C001_01559', 296, 43), ('VICH001.B.00001.C001_01560', 405, 44),
                        ('VICH001.B.00001.C001_01561', 182, 45), ('VICH001.B.00001.C001_01562', 445, 46),
                        ('VICH001.B.00001.C001_01563', 212, 47), ('VICH001.B.00001.C001_01564', 387, 48),
                        ('VICH001.B.00001.C001_01565', 414, 49)]
        self.assertListEqual(list(iter(idx)), expected_idx)

