# -*- coding: utf-8 -*-

#===============================================================================
# Created on Nov 30, 2012
# 
# @author: bneron
# @contact: user_email
# @organization: organization_name
# @license: license
#===============================================================================


import os
import unittest
import shutil
from txsscanlib.config import Config
from txsscanlib.database import RepliconDB, Indexes, RepliconInfo


class Test(unittest.TestCase):

    _data_dir = os.path.join(os.path.dirname(__file__), "datatest")

    def __init__(self, methodName = 'runTest'):
        super(Test, self).__init__(methodName)
        def fake_init(obj, cfg):
            obj.cfg = cfg
            idx = Indexes(self.cfg)
            obj.sequence_idx = idx.find_my_indexes()
            obj.topology_file = self.cfg.topology_file
            obj._DB = {}
        self.fake_init = fake_init
        self.real_init = RepliconDB.__init__


    def setUp(self):
        self.cfg = Config( hmmer_exe = "hmmsearch",
                           sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type = "gembase",
                           e_value_res = 1,
                           i_evalue_sel = 0.5,
                           def_dir = os.path.join(self._data_dir, 'DEF'),
                           res_search_dir = '/tmp',
                           res_search_suffix = ".search_hmm.out",
                           profile_dir = os.path.join(self._data_dir, 'profiles'),
                           profile_suffix = ".fasta-aln_edit.hmm",
                           res_extract_suffix = "",
                           log_level = 30,
                           log_file = '/dev/null'
                           )

        shutil.copy(self.cfg.sequence_db, self.cfg.working_dir)
        self.cfg.options['sequence_db'] = os.path.join(self.cfg.working_dir, os.path.basename(self.cfg.sequence_db))
        
        self.ESCO030p01_genes =['000010', '000020', '000030', '000040', '000050', '000060', '000070', '000080', '000090',
                                '000100', '000110', '000120', '000130', '000140', '000150', '000160', '000170', '000180',
                                '000190', '000200', '000210', '000220', '000230', '000240', '000250', '000260', '000270',
                                '000280', '000290', '000300', '000310', '000320', '000330', '000340', '000350', '000360',
                                '000370', '000380', '000390', '000400', '000410', '000420', '000430', '000440', '000450',
                                '000460', '000470', '000480', '000490', '000500', '000510', '000520', '000530', '000540',
                                '000550', '000560', '000570', '000580', '000590', '000600', '000610', '000620', '000630',
                                '000640', '000650', '000660', '000670', '000670']
        self.PSAE001c01_genes = ['006940', '013980', '017350', '018920', '026600', '031420', '043580', '051090', '055870', 
                                 '055880', '055890', '055900', '055910', '055920', '055930', '055940', '055950', '055960', 
                                 '055970', '055980', '055990', '056000', '056010', '056020', '056030', '056040', '056050', 
                                 '056060', '056070', '056080', '056090', '056100', '056110', '056120', '056130', '056140', 
                                 '056150', '056160', '056170', '056180', '056190', '056200', '056210', '056220', '056230', 
                                 '056240', '056250', '056260', '056270', '056280', '056290', '056300', '056310', '056320', 
                                 '056330', '056340', '056350', '056360', '056370', '056380', '056390', '056400', '056410', 
                                 '056420', '056430', '056440', '056440']
        
        idx = Indexes(self.cfg)
        idx._build_my_indexes()

    def tearDown(self):
        try:
            shutil.rmtree(self.cfg.working_dir)
        except:
            pass
        RepliconDB.__init__ = self.real_init


    def test_fill_topology(self):
        self.cfg.options['topology_file'] = self.cfg.sequence_db + ".topo"
        db_send = {'ESCO030p01' : 'circular',
                   'PSAE001c01' : 'linear'
                   }
        with open(self.cfg.topology_file , 'w') as f:
            for k, v in db_send.items():
                f.write('%s : %s\n' % (k,v))
        RepliconDB.__init__ = self.fake_init
        db = RepliconDB(self.cfg)
        rcv_topo = db._fill_topology()
        self.assertDictEqual(db_send, rcv_topo)
    

    def test_fill_ordered_replicon_min_max(self):
        self.tearDown()
        self.cfg = Config( hmmer_exe = "hmmsearch",
                           sequence_db = os.path.join(self._data_dir, "base", "ordered_replicon_base"),
                           db_type = "ordered_replicon",
                           e_value_res = 1,
                           i_evalue_sel = 0.5,
                           def_dir = os.path.join(self._data_dir, 'DEF'),
                           res_search_dir = '/tmp',
                           res_search_suffix = ".search_hmm.out",
                           profile_dir = os.path.join(self._data_dir, 'profiles'),
                           profile_suffix = ".fasta-aln_edit.hmm",
                           res_extract_suffix = "",
                           log_level = 30,
                           log_file = '/dev/null'
                           )

        shutil.copy(self.cfg.sequence_db, self.cfg.working_dir)
        self.cfg.options['sequence_db'] = os.path.join(self.cfg.working_dir, os.path.basename(self.cfg.sequence_db))

        idx = Indexes(self.cfg)
        idx._build_my_indexes()
        RepliconDB.__init__ = self.fake_init
        db = RepliconDB(self.cfg)
        db._fill_ordered_min_max(self.cfg.replicon_topology)

        self.assertEqual(len(db._DB), 1)
        rep = db[RepliconDB.ordered_replicon_name]
        self.assertEqual(rep.topology, self.cfg.replicon_topology)
        self.assertEqual(rep.min, 1)
        self.assertEqual(rep.max, 52)


    def test_fill_gembase_min_max_default_topology(self):
        RepliconDB.__init__ = self.fake_init
        db = RepliconDB(self.cfg)
        db._fill_gembase_min_max({}, self.cfg.replicon_topology)
        self.assertEqual(len(db._DB), 2)
        PRRU001c01 = db['ESCO030p01']
        self.assertEqual(PRRU001c01.topology, 'circular')
        self.assertEqual(PRRU001c01.min, 1)
        self.assertEqual(PRRU001c01.max, 67)
        PSAE001c01 = db['PSAE001c01']
        self.assertEqual(PSAE001c01.topology, 'circular')
        self.assertEqual(PSAE001c01.min, 68)
        self.assertEqual(PSAE001c01.max, 133)
 
 
    def test_fill_gembase_min_max_with_topology(self):
        self.cfg.options['topology_file'] = self.cfg.sequence_db + ".topo"
        with open(self.cfg.topology_file , 'w') as f:
            f.write('ESCO030p01 : circular\nPSAE001c01 : linear\n')
        RepliconDB.__init__ = self.fake_init
        db = RepliconDB(self.cfg)
        topo_dict = db._fill_topology()
        db._fill_gembase_min_max(topo_dict, 'circular')
        self.assertEqual(len(db._DB), 2)
        ESCO030p01 = db['ESCO030p01']
        self.assertEqual(ESCO030p01.topology, 'circular')
        self.assertEqual(ESCO030p01.min, 1)
        self.assertEqual(ESCO030p01.max, 67)
        PSAE001c01 = db['PSAE001c01']
        self.assertEqual(PSAE001c01.topology, 'linear')
        self.assertEqual(PSAE001c01.min, 68)
        self.assertEqual(PSAE001c01.max, 133)
         
 
    def test_in(self):
        db = RepliconDB(self.cfg)
        self.assertIn('ESCO030p01', db)
        self.assertIn('PSAE001c01', db)
        self.assertNotIn('toto', db)
 
    def test_getitem(self):
        db = RepliconDB(self.cfg)
        ESCO030p01 = RepliconInfo(self.cfg.replicon_topology, 1, 67, self.ESCO030p01_genes)
        PSAE001c01 = RepliconInfo(self.cfg.replicon_topology, 68, 133, self.PSAE001c01_genes)
        self.assertEqual(ESCO030p01, db['ESCO030p01'])
        self.assertEqual(PSAE001c01, db['PSAE001c01'])
        self.assertRaises(KeyError, db.__getitem__, 'foo')
 
    def test_get(self):
        db = RepliconDB(self.cfg)
        ESCO030p01 = RepliconInfo(self.cfg.replicon_topology, 1, 67, self.ESCO030p01_genes)
        PSAE001c01 = RepliconInfo(self.cfg.replicon_topology, 68, 133, self.PSAE001c01_genes)
        self.assertEqual(ESCO030p01, db.get('ESCO030p01'))
        self.assertEqual(PSAE001c01, db.get('PSAE001c01'))
        self.assertIsNone(db.get('foo'))
        self.assertEqual('bar', db.get('foo', 'bar'))

    def test_items(self):
        db = RepliconDB(self.cfg)
        ESCO030p01 = RepliconInfo(self.cfg.replicon_topology, 1, 67, self.ESCO030p01_genes)
        PSAE001c01 = RepliconInfo(self.cfg.replicon_topology, 68, 133, self.PSAE001c01_genes)
        self.assertItemsEqual(db.items(), [('ESCO030p01',ESCO030p01),('PSAE001c01',PSAE001c01)])
        db = RepliconDB(self.cfg)
        PRRU001c01 = RepliconInfo(self.cfg.replicon_topology, 1, 67, self.ESCO030p01_genes)
        PSAE001c01 = RepliconInfo(self.cfg.replicon_topology, 68, 133, self.PSAE001c01_genes)
        
