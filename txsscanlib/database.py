# -*- coding: utf-8 -*-

#===============================================================================
# Created on Feb 20, 2013
# 
# @author: bneron
# @contact: bneron@pasteur.fr
# @organization: Institut Pasteur
# @license: license
#===============================================================================



from itertools import groupby
from collections import namedtuple
from glob import glob
import os.path
import logging
_log = logging.getLogger('txsscan.' + __name__)
from subprocess import Popen
from bsddb3 import db

""" 
SequenceInfo contains 3 fields 
 
* id : the identifier of the sequence.
* length : the length of the sequence.
* position : it's rank in the sequences base in fasta format.

"""
SequenceInfo = namedtuple('SequenceInfo', 'id, length, position')


def fasta_iter(fasta_file):
    """
    :param fasta_file: the file containing all sequences in fasta format.
    :type fasta_file: file object
    :author: http://biostar.stackexchange.com/users/36/brentp
    :return: given a fasta file. return an iterator which yield tuples
             (string id, string comment, int sequence length)
    :rtype: iterator
    """
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fasta_file , lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        header = header.split()
        id = header[0]
        comment = ' '.join(header[1:])
        seq = ''.join(s.strip() for s in faiter.next())
        length = len(seq)
        yield id, comment, length


class Database(object):
    """
    Handle the indexes for txsscan
    
    # find the indexes for hmmer, or build them using formatdb external tool
    # find the indexes need by txsscan to compute some scores, or build them. These indexes are store in a Berkeley DB
      
    """


    def __init__(self, cfg):
        """
        Constructor retrieve file of indexes, if they are not present
        or the user ask for build indexes (--idx) launch the indexes building. 

        :param cfg: the configuration 
        :type cfg: :class:`txsscanlib.config.Config` object
        """
        self.cfg = cfg
        self._fasta_path = cfg.sequence_db
        self.name = os.path.basename(cfg.sequence_db)
        self._hmmer_indexes = None
        self._my_indexes = None
        self._my_open_db = None

    def build(self, force = False):
        """
        build the indexes from the sequences base in fasta format
        """
        hmmer_indexes = self._find_hmmer_indexes()
        my_indexes = self._find_my_indexes()

        ###########################
        # build indexes if needed #
        ###########################
        index_dir = os.path.dirname(self.cfg.sequence_db)
        if force or not hmmer_indexes or not my_indexes:
            if not os.access(index_dir, os.W_OK):
                msg = "cannot build indexes, (%s) is not writable" % index_dir
                _log.critical(msg)
                raise IOError(msg)
        
        if force or not hmmer_indexes:
            #self._build_hmmer_indexes() is asynchron
            hmmer_indexes_proc = self._build_hmmer_indexes()
        if force or not my_indexes:
            #self._build_my_indexes() is synchron
            self._build_my_indexes()

        ################################# 
        # synchronization point between #
        # hmmer_indexes and my_indexes  #
        #################################
        if force or not hmmer_indexes:
            hmmer_indexes_proc.wait()
            if hmmer_indexes_proc.returncode != 0:
                msg = "an error occurred during databases indexation see formatdb.log f"
                _log.error( msg, exc_info = True )
                raise RuntimeError(msg)
        self._hmmer_indexes = self._find_hmmer_indexes()
        self._my_indexes = self._find_my_indexes()
        assert self._hmmer_indexes , "cannot create hmmer indexes"
        assert self._my_indexes, "cannot create txsscan indexes"

    def _find_hmmer_indexes(self):
        """
        :return: the files wich belongs to the hmmer indexes. 
                 If indexes are inconsistent (lack file) a Runtime Error is raised
        :rtype: list of string 
        """
        suffixes = ('.phr', '.pin', '.psd', '.psi', '.psq', '.pal')
        idx = []
        file_nb = 0
        for suffix in suffixes:
            index_files = glob( "%s*%s" % (self._fasta_path, suffix))
            nb_of_index = len(index_files)
            if suffix != '.pal':
                if file_nb and file_nb != nb_of_index:
                    msg = "some indexes lack. remove indexes (*.phr, *.pin, *.psd, *.psi, *.psq, *.pal) and try to rebuild them."
                    _log.critical(msg)
                    raise  RuntimeError(msg)
            else:
                if nb_of_index > 1:
                    msg = "too many .pal file . remove indexes (*.phr, *.pin, *.psd, *.psi, *.psq, *.pal) and try to rebuild them."
                    _log.critical(msg)
                    raise  RuntimeError(msg)    
                elif file_nb > 1 and nb_of_index == 0:
                    msg = "some indexes lack. remove indexes (*.phr, *.pin, *.psd, *.psi, *.psq, *.pal) and try to rebuild them."
                    _log.critical(msg)
                    raise  RuntimeError(msg)
                elif file_nb == 1 and nb_of_index == 1:
                    msg = "a virtual index is detected (.pal) but there is only one file per index type. remove indexes (*.phr, *.pin, *.psd, *.psi, *.psq, *.pal) and try to rebuild them."
                    _log.critical(msg)
                    raise  RuntimeError(msg)
            idx.extend(index_files)
            file_nb = nb_of_index
        return idx
    

    def _find_my_indexes(self):
        """
        :return: the file of txsscan if exits, None otherwise. 
        :rtype: string
        """ 
        path = os.path.join( os.path.dirname(self.cfg.sequence_db), self.name + ".idx")
        if os.path.exists(path):
            return path


 
    def _build_hmmer_indexes(self):
        """
        build the indexes for hmmer using formatdb tool
        """
        #formatdb create indexes in the same directory as the sequence_db
        #so it must be writable
        #if the directory is not writable, formatdb do a Segmentation fault
        index_dir = os.path.dirname(self.cfg.sequence_db)
        command = "formatdb -t %s -i %s -p T -o T -s T" % ( self.name,
                                                            self.cfg.sequence_db
                                                          )

        err_path = os.path.join(index_dir, "formatdb.err")
        with  open(err_path, 'w') as err_file:
            try:
                formatdb = Popen( command ,
                                  shell = True ,
                                  stdout = None ,
                                  stdin  = None ,
                                  stderr = err_file ,
                                  close_fds = False ,
                                  )
            except Exception, err:
                msg = "unable to format the sequence base : %s : %s" % ( command , err)
                _log.critical( msg, exc_info = True )
                raise err
            return formatdb


    def _build_my_indexes(self):
        """
        build txsscan indexes. This index is stored in a berkeley DB
        """
        my_index = db.DB()
        my_index.open(os.path.join( os.path.dirname(self.cfg.sequence_db), self.name + ".idx"),
                      dbname = self.name,
                      dbtype = db.DB_HASH,
                      flags = db.DB_CREATE)
        try:
            with open(self._fasta_path, 'r') as fasta_file:
                f_iter = fasta_iter(fasta_file)
                seq_nb = 0
                for seqid, comment, length in f_iter:
                    seq_nb += 1
                    my_index.put(seqid, "%d;%d" % (length, seq_nb))
        except Exception, err:
            msg = "unable to index the sequence base: %s : %s" % (self.cfg.sequence_db, err )
            _log.critical( msg, exc_info = True )
            raise err
        finally:
            my_index.close()
    
    def open(self):
        """
        open a connection. This method must be called before to get items or use the with statement::
        
          db = Database(self.cfg)
          db.build()
          db.open()
          db[ my_key ]
          db.close()
        
        or::
        
          with db.open():
             db[ my_key ] 
        
        """
        if not self._my_indexes:
            msg = "no index available, the database must be build before to open it"
            _log.critical(msg)
            raise IOError(msg)
        mydb = db.DB()
        mydb.open(self._my_indexes,
                      dbname = self.name,
                      dbtype = db.DB_HASH,
                      flags = db.DB_RDONLY)
        self._my_open_db = mydb
        return self

    def close(self):
        """
        close the connection. This method must be called when the search in the db is ended.
        or use the with statement::
        
          db = Database(self.cfg)
          db.build()
          db.open()
          db[ my_key ]
          db.close()
        
        or::
        
          with db.open():
             db[ my_key ] 
        
        """
        if self._my_open_db:
            self._my_open_db.close()
        self._my_open_db = None
        
    def __enter__(self):
        return self.open()
    
    def __exit__(self, exctype, exc, traceback):
        if traceback is None: #no traceback means no error
            self.close()
            success = True
        else:
            self.close()
            success = False
        return success    
    
    
    def __getitem__(self, seq_id):
        """
        allow to use the following notation to retrieve sequence information ::
         
          db = Database()
          with db.open():
             seq_info = db[ 'my_id' ] 
          seq_info.lenght
         
        :param seq_id: the sequence identifier
        :type seq_id: string
        :return: the SequenceInfo corresponding to the seq_id.
                 if the seq_id does not exist in the database a KeyError is raised.
        :rtype:  :class:`txsscanlib.database.SeqInfo` object
        """
        if not self._my_open_db:
            msg = msg = "the database must be open before to get items"
            _log.critical(msg)
            raise IOError(msg)
        data = self._my_open_db.get(seq_id)
        if data:
            length, seq_nb = data.split(';')
            return SequenceInfo(seq_id, int(length), int(seq_nb))
        else:
            raise KeyError(str(seq_id))


    def get(self, seq_id, default = None):
        """
        Return the :class:`txsscanlib.database.SeqInfo` object for given seq_id if seq_id is in the dictionary, else default. 
        If default is not given, it defaults to None, so that this method never raises a KeyError.::
        
          db = Database()
          with db.open():
             seq_info = db.get( 'my_id', None ) 
          seq_info.lenght

        :param seq_id: the sequence identifier
        :type seq_id: string
        :param default: the value return if the seq_id is not in the database
        :type default: any
        """
        if not self._my_open_db:
            msg = "the database must be open before to get items"
            _log.critical(msg)
            raise IOError(msg)
        data = self._my_open_db.get(seq_id)
        if data:
            length, seq_nb = data.split(';')
            return SequenceInfo(seq_id, int(length), int(seq_nb))
        else:
            return default
