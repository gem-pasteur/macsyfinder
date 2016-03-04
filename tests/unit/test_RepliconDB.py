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
from macsypy.database import RepliconDB, Indexes, RepliconInfo


class Test(unittest.TestCase):

    _data_dir = os.path.join(os.path.dirname(__file__), "..", "data")

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
        l = logging.getLogger()
        l.manager.loggerDict.clear()
        
        #add only one handler to the macsypy logger
        from macsypy.database import _log
        macsy_log = _log.parent
        log_file = 'NUL' if platform.system() == 'Windows' else '/dev/null'
        log_handler = logging.FileHandler(log_file)
        macsy_log.addHandler(log_handler)
        
        self.cfg = Config( hmmer_exe = "hmmsearch",
                           sequence_db = os.path.join(self._data_dir, "base", "test_base.fa"),
                           db_type = "gembase",
                           e_value_res = 1,
                           i_evalue_sel = 0.5,
                           def_dir = os.path.join(self._data_dir, 'DEF'),
                           res_search_dir = tempfile.gettempdir(),
                           res_search_suffix = ".search_hmm.out",
                           profile_dir = os.path.join(self._data_dir, 'profiles'),
                           profile_suffix = ".hmm",
                           res_extract_suffix = "",
                           log_level = 30,
                           log_file = log_file
                           )

        shutil.copy(self.cfg.sequence_db, self.cfg.working_dir)
        self.cfg.options['sequence_db'] = os.path.join(self.cfg.working_dir, os.path.basename(self.cfg.sequence_db))
        
        self.ESCO030p01_genes = [('000010', '886'), ('000020', '291'), ('000030', '656'), ('000040', '500'), ('000050', '407'),
                                 ('000060', '144'), ('000070', '183'), ('000080', '121'), ('000090', '199'), ('000100', '325'),
                                 ('000110', '425'), ('000120', '171'), ('000130', '277'), ('000140', '133'), ('000150', '108'),
                                 ('000160', '295'), ('000170', '273'), ('000180', '367'), ('000190', '573'), ('000200', '343'),
                                 ('000210', '295'), ('000220', '108'), ('000230', '117'), ('000240', '153'), ('000250', '479'),
                                 ('000260', '706'), ('000270', '998'), ('000280', '171'), ('000290', '108'), ('000300', '295'),
                                 ('000310', '165'), ('000320', '243'), ('000330', '295'), ('000340', '108'), ('000350', '1755'),
                                 ('000360', '248'), ('000370', '286'), ('000380', '186'), ('000390', '83'), ('000400', '153'),
                                 ('000410', '69'), ('000420', '295'), ('000430', '108'), ('000440', '145'), ('000450', '59'),
                                 ('000460', '124'), ('000470', '246'), ('000480', '325'), ('000490', '54'), ('000500', '95'),
                                 ('000510', '83'), ('000520', '56'), ('000530', '401'), ('000540', '320'), ('000550', '256'),
                                 ('000560', '73'), ('000570', '144'), ('000580', '258'), ('000590', '133'), ('000600', '140'),
                                 ('000610', '63'), ('000620', '138'), ('000630', '68'), ('000640', '169'), ('000650', '127'),
                                 ('000660', '295'), ('000670', '108'), ('000670', '108')]
                                
        self.PSAE001c01_genes = [('006940', '803'), ('013980', '759'), ('017350', '600'), ('018920', '776'), ('026600', '273'), 
                                 ('031420', '658'), ('043580', '416'), ('051090', '714'), ('055870', '449'), ('055880', '447'), 
                                 ('055890', '588'), ('055900', '292'), ('055910', '262'), ('055920', '166'), ('055930', '288'), 
                                 ('055940', '194'), ('055950', '567'), ('055960', '188'), ('055970', '247'), ('055980', '252'), 
                                 ('055990', '455'), ('056000', '450'), ('056010', '260'), ('056020', '246'), ('056030', '70'), 
                                 ('056040', '133'), ('056050', '284'), ('056060', '585'), ('056070', '435'), ('056080', '342'), 
                                 ('056090', '252'), ('056100', '122'), ('056110', '213'), ('056120', '400'), ('056130', '134'), 
                                 ('056140', '138'), ('056150', '397'), ('056160', '298'), ('056170', '186'), ('056180', '445'), 
                                 ('056190', '414'), ('056200', '132'), ('056210', '674'), ('056220', '319'), ('056230', '394'), 
                                 ('056240', '207'), ('056250', '401'), ('056260', '611'), ('056270', '257'), ('056280', '169'), 
                                 ('056290', '454'), ('056300', '141'), ('056310', '458'), ('056320', '286'), ('056330', '514'), 
                                 ('056340', '178'), ('056350', '156'), ('056360', '85'), ('056370', '289'), ('056380', '126'), 
                                 ('056390', '290'), ('056400', '262'), ('056410', '214'), ('056420', '630'), ('056430', '127'), 
                                 ('056440', '455'), ('056440', '455')]
        self.NCDB_genes = [('056134', '289'), ('056135', '126'), ('056136', '290'),
                           ('056137', '262'), ('056138', '214'), ('056139', '630'),
                           ('056140', '127'), ('056141', '803'), ('056141', '803')]

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
        RepliconDB.__init__ = self.real_init


    def test_fill_topology(self):
        self.cfg.options['topology_file'] = self.cfg.sequence_db + ".topo"
        db_send = {'ESCO030p01' : 'circular',
                   'PSAE001c01' : 'linear'
                   }
        with open(self.cfg.topology_file , 'w') as f:
            for k, v in db_send.items():
                f.write('{0} : {1}\n'.format(k,v))
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
                           res_search_dir = tempfile.gettempdir(),
                           res_search_suffix = ".search_hmm.out",
                           profile_dir = os.path.join(self._data_dir, 'profiles'),
                           profile_suffix = ".hmm",
                           res_extract_suffix = "",
                           log_level = 30,
                           log_file = 'NUL' if platform.system() == 'Windows' else '/dev/null'
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
        self.assertEqual(len(db._DB), 3)
        self.assertEqual(set(db._DB.keys()), set(['ESCO030p01', 'PSAE001c01', 'NC_xxxxx_xx']))
        PRRU001c01 = db['ESCO030p01']
        self.assertEqual(PRRU001c01.topology, 'circular')
        self.assertEqual(PRRU001c01.min, 1)
        self.assertEqual(PRRU001c01.max, 67)
        self.assertEqual(PRRU001c01.genes, self.ESCO030p01_genes)
        PSAE001c01 = db['PSAE001c01']
        self.assertEqual(PSAE001c01.topology, 'circular')
        self.assertEqual(PSAE001c01.min, 68)
        self.assertEqual(PSAE001c01.max, 133)
        self.assertEqual(PSAE001c01.genes, self.PSAE001c01_genes)
        DBNC = db['NC_xxxxx_xx']
        self.assertEqual(DBNC.topology, 'circular')
        self.assertEqual(DBNC.min, 134)
        self.assertEqual(DBNC.max, 141)
        self.assertEqual(DBNC.genes, self.NCDB_genes)

 
 
    def test_fill_gembase_min_max_with_topology(self):
        self.cfg.options['topology_file'] = self.cfg.sequence_db + ".topo"
        with open(self.cfg.topology_file , 'w') as f:
            f.write('# topology file\nESCO030p01 : circular\nPSAE001c01 : linear\n')
        RepliconDB.__init__ = self.fake_init
        db = RepliconDB(self.cfg)
        topo_dict = db._fill_topology()
        db._fill_gembase_min_max(topo_dict, 'circular')
        self.assertEqual(len(db._DB), 3)
        self.assertEqual(set(db._DB.keys()), set(['ESCO030p01', 'PSAE001c01', 'NC_xxxxx_xx']))
        ESCO030p01 = db['ESCO030p01']
        self.assertEqual(ESCO030p01.topology, 'circular')
        self.assertEqual(ESCO030p01.min, 1)
        self.assertEqual(ESCO030p01.max, 67)
        self.assertEqual(ESCO030p01.genes, self.ESCO030p01_genes)
        PSAE001c01 = db['PSAE001c01']
        self.assertEqual(PSAE001c01.topology, 'linear')
        self.assertEqual(PSAE001c01.min, 68)
        self.assertEqual(PSAE001c01.max, 133)
        self.assertEqual(PSAE001c01.genes, self.PSAE001c01_genes)
        DBNC = db['NC_xxxxx_xx']
        self.assertEqual(DBNC.topology, 'circular')
        self.assertEqual(DBNC.min, 134)
        self.assertEqual(DBNC.max, 141)
        self.assertEqual(DBNC.genes, self.NCDB_genes)
         
 
    def test_in(self):
        db = RepliconDB(self.cfg)
        self.assertIn('ESCO030p01', db)
        self.assertIn('PSAE001c01', db)
        self.assertIn('NC_xxxxx_xx', db)
        self.assertNotIn('toto', db)
 
    def test_getitem(self):
        db = RepliconDB(self.cfg)
        ESCO030p01 = RepliconInfo(self.cfg.replicon_topology, 1, 67, self.ESCO030p01_genes)
        PSAE001c01 = RepliconInfo(self.cfg.replicon_topology, 68, 133, self.PSAE001c01_genes)
        NCXX = RepliconInfo("circular", 134, 141, self.NCDB_genes)
        self.assertEqual(ESCO030p01, db['ESCO030p01'])
        self.assertEqual(PSAE001c01, db['PSAE001c01'])
        self.assertEqual(NCXX, db['NC_xxxxx_xx'])
        self.assertRaises(KeyError, db.__getitem__, 'foo')

    def test_get(self):
        db = RepliconDB(self.cfg)
        ESCO030p01 = RepliconInfo(self.cfg.replicon_topology, 1, 67, self.ESCO030p01_genes)
        PSAE001c01 = RepliconInfo(self.cfg.replicon_topology, 68, 133, self.PSAE001c01_genes)
        NCXX = RepliconInfo("circular", 134, 141, self.NCDB_genes)
        self.assertEqual(ESCO030p01, db.get('ESCO030p01'))
        self.assertEqual(PSAE001c01, db.get('PSAE001c01'))
        self.assertEqual(NCXX, db.get('NC_xxxxx_xx', 'foo'))
        self.assertIsNone(db.get('foo'))
        self.assertEqual('bar', db.get('foo', 'bar'))

    def test_items(self):
        db = RepliconDB(self.cfg)
        ESCO030p01 = RepliconInfo(self.cfg.replicon_topology, 1, 67, self.ESCO030p01_genes)
        PSAE001c01 = RepliconInfo(self.cfg.replicon_topology, 68, 133, self.PSAE001c01_genes)
        NCXX = RepliconInfo("circular", 134, 141, self.NCDB_genes)
        self.assertItemsEqual(db.items(), [('ESCO030p01',ESCO030p01), ('NC_xxxxx_xx',NCXX),
                                           ('PSAE001c01',PSAE001c01)])