#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2021  Institut Pasteur (Paris) and CNRS.           #
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


from itertools import groupby
from collections import namedtuple
import os.path
import logging
from macsypy.error import MacsypyError
_log = logging.getLogger(__name__)


def fasta_iter(fasta_file):
    """
    :param fasta_file: the file containing all input sequences in fasta format.
    :type fasta_file: file object
    :author: http://biostar.stackexchange.com/users/36/brentp
    :return: for a given fasta file, it returns an iterator which yields tuples
             (string id, string comment, int sequence length)
    :rtype: iterator
    """
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fasta_file, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = next(header)[1:].strip()
        header = header.split()
        _id = header[0]
        comment = ' '.join(header[1:])
        seq = ''.join(s.strip() for s in next(faiter))
        length = len(seq)
        yield _id, comment, length


class Indexes:
    """
    Handle the indexes for macsyfinder:

     - find the indexes required by macsyfinder to compute some scores, or build them.
    """


    def __init__(self, cfg):
        """
        The constructor retrieves the file of indexes in the case they are not present
        or the user asked for build indexes (--idx) 
        Launch the indexes building. 

        :param cfg: the configuration 
        :type cfg: :class:`macsypy.config.Config` object
        """
        self.cfg = cfg
        self._fasta_path = cfg.sequence_db()
        self.name = os.path.basename(self._fasta_path)
        self._my_indexes = None  # path

    def build(self, force=False):
        """
        Build the indexes from the sequence data set in fasta format

        :param force: If True, force the index building even
                      if the index files are present in the sequence data set folder
        :type force: boolean
        :return: the path to the index
        :rtype: str
        """
        my_indexes = self.find_my_indexes()

        ###########################
        # build indexes if needed #
        ###########################
        index_dir = os.path.abspath(os.path.dirname(self.cfg.sequence_db()))

        if force or not my_indexes:
            # formatdb create indexes in the same directory as the sequence_db
            # so it must be writable
            # if the directory is not writable, formatdb do a Segmentation fault
            if not os.access(index_dir, os.R_OK | os.W_OK):
                msg = f"cannot build indexes, ({index_dir}) is not writable"
                _log.critical(msg)
                raise IOError(msg)

            self._build_my_indexes()

        self._my_indexes = self.find_my_indexes()
        assert self._my_indexes, "failed create macsyfinder indexes"
        return self._my_indexes


    def find_my_indexes(self):
        """
        :return: the file of macsyfinder indexes if it exists in the dataset folder, None otherwise. 
        :rtype: string
        """ 
        path = os.path.join(os.path.dirname(self.cfg.sequence_db()), self.name + ".idx")
        if os.path.exists(path):
            return path

    def _build_my_indexes(self):
        """
        Build macsyfinder indexes. These indexes are stored in a file.

        The file format is the following:
         - one entry per line, with each line having this format:
         - sequence id;sequence length;sequence rank

        """
        try:
            with open(self._fasta_path, 'r') as fasta_file:
                with open(os.path.join(os.path.dirname(self.cfg.sequence_db()), self.name + ".idx"), 'w') as my_base:
                    f_iter = fasta_iter(fasta_file)
                    seq_nb = 0
                    for seq_id, comment, length in f_iter:
                        seq_nb += 1
                        my_base.write(f"{seq_id};{length:d};{seq_nb:d}\n")
        except Exception as err:
            msg = f"unable to index the sequence dataset: {self.cfg.sequence_db()} : {err}"
            _log.critical(msg, exc_info=True)
            raise MacsypyError(msg)


"""handle name, topology type, and min/max positions in the sequence dataset for a replicon and list of genes.
each genes is representing by a tuple (seq_id, length)"""
RepliconInfo = namedtuple('RepliconInfo', 'topology, min, max, genes')


class RepliconDB:
    """
    Stores information (topology, min, max, [genes]) for all replicons in the sequence_db
    the Replicon object must be instantiated only for sequence_db of type 'gembase' or 'ordered_replicon'
    """

    ordered_replicon_name = 'UserReplicon'

    def __init__(self, cfg):
        """
        :param cfg: The configuration object
        :type cfg: :class:`macsypy.config.Config` object

        .. note ::
            This class can be instanciated only if the db_type is 'gembase' or 'ordered_replicon' 
        """
        self.cfg = cfg
        assert self.cfg.db_type() in ('gembase', 'ordered_replicon')
        idx = Indexes(self.cfg)
        self.sequence_idx = idx.find_my_indexes()
        self.topology_file = self.cfg.topology_file()
        self._DB = {}
        if self.topology_file:
            topo_dict = self._fill_topology()
        else:
            topo_dict = {}
        if self.cfg.db_type() == 'gembase':
            self._fill_gembase_min_max(topo_dict, default_topology=self.cfg.replicon_topology())
        else:
            self._fill_ordered_min_max(self.cfg.replicon_topology())


    def _fill_topology(self):
        """
        Fill the internal dictionary with min and max positions for each replicon_name of the sequence_db
        """
        topo_dict = {}
        with open(self.topology_file) as topo_f:
            for line in topo_f:
                if line.startswith('#'):
                    continue
                replicon_name, topo = line.split(':')
                replicon_name = replicon_name.strip()
                topo = topo.strip().lower()
                topo_dict[replicon_name] = topo
        return topo_dict

    def _fill_ordered_min_max(self, default_topology=None):
        """
        For the replicon_name of the ordered_replicon sequence base, fill the internal dict with RepliconInfo

        :param default_topology: the topology provided by config.replicon_topology 
        :type default_topology: string
        """
        _min = 1
        # self.sequence_idx is a file with the following structure seq_id;seq_length;seq_rank\n
        with open(self.sequence_idx) as idx_f:
            _max = 0
            genes = []
            for line in idx_f:
                seq_id, length, _rank = line.split(";")
                genes.append((seq_id, length))
                _max += 1
            self._DB[self.ordered_replicon_name] = RepliconInfo(default_topology, _min, _max, genes)


    def _fill_gembase_min_max(self, topology, default_topology):
        """
        For each replicon_name of a gembase dataset, it fills the internal dictionary with a namedtuple RepliconInfo

        :param topology: the topologies for each replicon 
                         (parsed from the file specified with the option --topology-file)
        :type topology: dict
        :param default_topology: the topology provided by the config.replicon_topology 
        :type default_topology: string
        """
        def grp_replicon(line):
            """
            in gembase the identifier of fasta sequence follows the following schema: 
            <replicon-name>_<seq-name> with eventually '_' inside the <replicon_name>
            but not in the <seq-name>.
            so grp_replicon allow to group sequences belonging to the same replicon.
            """
            return "_".join(line.split('_')[: -1])

        def parse_entry(entry):
            """
            parse an entry in the index file (.idx)
            an entry have the following format sequence_id;sequence length;sequence rank in replicon
            """
            entry = entry.rstrip()
            seq_id, length, rank = entry.split(';')
            replicon_name = "_".join(seq_id.split('_')[: -1])
            seq_name = seq_id.split('_')[-1]
            return replicon_name, seq_name, length, int(rank)

        with open(self.sequence_idx) as idx_f:
            replicons = (x[1] for x in groupby(idx_f, grp_replicon))
            for replicon in replicons:
                genes = []
                entry = next(replicon)
                replicon_name, seq_name, seq_length, _min = parse_entry(entry)
                genes.append((seq_name, seq_length))
                for entry in replicon:
                    # pass all sequence of the replicon until the last one
                    _, seq_name, seq_length, _ = parse_entry(entry)
                    genes.append((seq_name, seq_length))
                _, seq_name, seq_length, _max = parse_entry(entry)
                genes.append((seq_name, seq_length))
                if replicon_name in topology:
                    self._DB[replicon_name] = RepliconInfo(topology[replicon_name], _min, _max, genes)
                else:
                    self._DB[replicon_name] = RepliconInfo(default_topology, _min, _max, genes)


    def __contains__(self, replicon_name):
        """
        :param replicon_name: the name of the replicon
        :type replicon_name: string
        :returns: True if replicon_name is in the repliconDB, false otherwise.
        :rtype: boolean
        """
        return replicon_name in self._DB


    def __getitem__(self, replicon_name):
        """
        :param replicon_name: the name of the replicon to get information on
        :type replicon_name: string
        :returns: the RepliconInfo for the provided replicon_name 
        :rtype: :class:`RepliconInfo` object
        :raise: KeyError if replicon_name is not in repliconDB
        """
        return self._DB[replicon_name]


    def get(self, replicon_name, default=None):
        """
        :param replicon_name: the name of the replicon to get informations
        :type replicon_name: string
        :param default: the value to return if the replicon_name is not in the RepliconDB
        :type default: any
        :returns: the RepliconInfo for replicon_name if replicon_name is in the repliconDB, else default.
                  If default is not given, it is set to None, so that this method never raises a KeyError.
        :rtype: :class:`RepliconInfo` object
        """
        return self._DB.get(replicon_name, default)


    def items(self):
        """
        :return: a copy of the RepliconDB as a list of (replicon_name, RepliconInfo) pairs
        """
        return list(self._DB.items())


    def iteritems(self):
        """
        :return: an iterator over the RepliconDB as a list (replicon_name, RepliconInfo) pairs
        """
        return iter(self._DB.items())


    def replicon_names(self):
        """
        :return: a copy of the RepliconDB as a list of replicon_names
        """
        return list(self._DB.keys())


    def replicon_infos(self):
        """
        :return: a copy of the RepliconDB as list of replicons info
        :rtype: RepliconInfo instance
        """
        return list(self._DB.values())
