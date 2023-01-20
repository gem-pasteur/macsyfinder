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

"""
Module to handle sequences and their indexes
"""

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
        try:
            seq = ''.join(s.strip() for s in next(faiter))
        except StopIteration:
            # the sequence was not start by '>'
            # bad fasta format
            msg = f"Error during sequence '{fasta_file.name}' parsing: Check the fasta format."
            _log.critical(msg)
            raise MacsypyError(msg)
        length = len(seq)
        yield _id, comment, length


class Indexes:
    """
    Handle the indexes for macsyfinder:

     - find the indexes required by macsyfinder to compute some scores, or build them.
    """

    _field_separator = "^^"

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


    def build(self, force=False):
        """
        Build the indexes from the sequence data set in fasta format,

        :param force: If True, force the index building even
                      if the index files are present in the sequence data set folder
        :type force: boolean
        :return: the path to the index
        :rtype: str
        """
        my_indexes = self.find_my_indexes()  # check read

        ###########################
        # build indexes if needed #
        ###########################
        if my_indexes and not force:
            with open(my_indexes) as idx:
                seq_path = next(idx).strip()
                try:
                    first_item = next(idx).strip()
                except StopIteration:
                    # there is only one line in file
                    first_item = None
            if seq_path.count(';') == 2:
                # there is no path in idx, it's an old index
                _log.warning(f"The '{my_indexes}' index file is in old format. Force index building.")
                force = True
            elif seq_path != self._fasta_path:
                _log.warning(f"The '{my_indexes}' index file does not point to '{self._fasta_path}'. Force building")
                force = True
            if not force and first_item:
                # the first line of idx is a valid path
                if first_item.count(self._field_separator) == 0:
                    # the separator is different than the actual separator
                    _log.warning(f"The '{my_indexes}' index file is in old format. Force index building.")
                    force = True
            # if fasta file is newer than idx
            stamp_fasta = os.path.getmtime(self._fasta_path)
            stamp_idx = os.path.getmtime(my_indexes)
            if stamp_idx < stamp_fasta:
                _log.debug("the sequence index is older than sequence file: rebuild the index.")
                force = True

        if force or not my_indexes:
            try:
                index_dir = self._index_dir(build=True)  # check build
            except ValueError as err:
                msg = str(err)
                _log.critical(msg)
                raise IOError(msg) from None

            my_indexes = self._build_my_indexes(index_dir)
        return my_indexes


    def find_my_indexes(self):
        """
        :return: the file of macsyfinder indexes if it exists in the dataset folder, None otherwise.
        :rtype: string
        """
        index_dir = self._index_dir(build=False)
        path = os.path.join(index_dir, self.name + ".idx")
        if os.path.exists(path):
            return path


    def _index_dir(self, build=False):
        """
        search where to store(build=True) read indexes

        :param bool build: if check the index-dir permissions to write
        :return: The directory where read or write the indexes
        :rtype: str
        :raise ValueError: if the directory specify by --index-dir option does not exists
                           or if build = True index-dir is not writable
        """
        index_dir = self.cfg.index_dir()
        if index_dir:
            if not os.path.exists(index_dir):
                raise ValueError(f"No such directory: {index_dir}")
            elif build and not os.access(index_dir, os.W_OK):
                raise ValueError(f"The '{index_dir}' dir is not writable.")
            else:
                return index_dir
        else:
            # we need abspath because if user provide filename not path for sequence_db
            # for instance my_seq.faste instead of ./my_seq.fasta
            # then index_dir is empty string
            # and os.access return False
            index_dir = os.path.dirname(os.path.abspath(self.cfg.sequence_db()))
            if build and not os.access(index_dir, os.W_OK):
                raise ValueError(f"The '{index_dir}' dir is not writable. Change rights or specify --index-dir.")
            else:
                return index_dir


    def _build_my_indexes(self, index_dir):
        """
        Build macsyfinder indexes. These indexes are stored in a file.

        The file format is the following:
         - the first line is the path of the sequence-db indexed
         - one entry per line, with each line having this format:
         - sequence id;sequence length;sequence rank

        """
        index_file = os.path.join(index_dir, self.name + ".idx")
        try:
            with open(self._fasta_path, 'r') as fasta_file:
                with open(index_file, 'w') as my_base:
                    my_base.write(self._fasta_path + '\n')
                    f_iter = fasta_iter(fasta_file)
                    seq_nb = 0
                    for seq_id, comment, length in f_iter:
                        seq_nb += 1
                        my_base.write(f"{seq_id}{self._field_separator}{length:d}{self._field_separator}{seq_nb:d}\n")
                    my_base.flush()
        except Exception as err:
            msg = f"unable to index the sequence dataset: {self.cfg.sequence_db()} : {err}"
            _log.critical(msg, exc_info=True)
            raise MacsypyError(msg) from err
        return index_file


    def __iter__(self):
        """
        :raise MacsypyError: if the indexes are not buid
        :return: an iterator on the indexes

        To use it the index must be build.
        """
        path = self.find_my_indexes()
        if path is None:
            raise MacsypyError("Build index before to use it.")
        with open(path) as idx_file:
            # The first line of index is the path to the data
            # It is not an index
            _ = next(idx_file)
            for line in idx_file:
                try:
                    seq_id, length, _rank = line.split(self._field_separator)
                except Exception as err:
                    raise MacsypyError(f"fail to parse database index {path} at line: {line}."
                                       f"Try to rebuild index with --idx option or remove file."
                                       f"If error persist feel free to submit an issue at"
                                       f"https://github.com/gem-pasteur/macsyfinder/issues/new?assignees=&labels=bug&template=bug_report.md&title=%5BBUG%5D ", err) from err
                length = int(length)
                _rank = int(_rank)
                yield seq_id, length, _rank


RepliconInfo = namedtuple('RepliconInfo', ('topology', 'min', 'max', 'genes'))
"""
handle information about a replicon

.. py:attribute:: topology
    :noindex:

    The type of replicon topology 'linear or 'circular'

.. py:attribute:: min
    :noindex:

    The position of the last gene of the replicon in the sequence dataset.

.. py:attribute:: max
    :noindex:

    The position of the last gene of the replicon in the sequence dataset.

.. py:attribute:: genes
    :noindex:

    A list of genes beloging to the replicon. Each genes is representing by a tuple (str seq_id, int length)
"""


class RepliconDB:
    """
    Stores information (topology, min, max, [genes]) for all replicons in the sequence_db
    the Replicon object must be instantiated only for sequence_db of type 'gembase' or 'ordered_replicon'
    """


    def __init__(self, cfg):
        """
        :param cfg: The configuration object
        :type cfg: :class:`macsypy.config.Config` object

        .. note ::
            This class can be instanciated only if the db_type is 'gembase' or 'ordered_replicon'
        """
        self.cfg = cfg
        assert self.cfg.db_type() in ('gembase', 'ordered_replicon')
        self._idx = Indexes(self.cfg)
        self.topology_file = self.cfg.topology_file()
        self._DB = {}
        if self.topology_file:
            topo_dict = self._fill_topology()
        else:
            topo_dict = {}
        if self.cfg.db_type() == 'gembase':
            self._fill_gembase_min_max(topo_dict, default_topology=self.cfg.replicon_topology())
        else:
            self._ordered_replicon_name = os.path.splitext(os.path.basename(self.cfg.sequence_db()))[0]
            self._fill_ordered_min_max(self.cfg.replicon_topology())

    @property
    def ordered_replicon_name(self):
        return self._ordered_replicon_name

    def guess_if_really_gembase(self):
        """
        Count the number of replicon with only on sequence
        if this number is above a threshold may be it's not gembase. for instance the folowing sequence
        have id compliant with the gembase id syntax but it's not it only contains one replicon ('ordered replicon')

        | >1E10S0A0cP00_0010 D GTG TGA 483 2027 Valid dnaA 1545 _PA0001_NP_064721.1_ PA0001 1 483 2027
        | MSVELWQQCVDLLRDELPSQQFNTWIRPLQVEAEGDELRVYAPNRFVLDW
        | >0200S001A0c_0P1E0 D ATG TAA 2056 3159 Valid dnaN 1104 _PA0002_NP_064722.1_ PA0002 1 2056 3159
        | MHFTIQREALLKPLQLVAGVVERRQTLPVLSNVLLVVEGQQLSLTGTDLE
        | >0000310E00S0c_1PA D ATG TGA 3169 4278 Valid recF 1110 _PA0003_NP_064723.1_ PA0003 1 3169 4278
        | MSLTRVSVTAVRNLHPVTLSPSPRINILYGDNGSGKTSVLEAIHLLGLAR
        | >c_01000A0PS00014E D ATG TGA 4275 6695 Valid gyrB 2421 _PA0004_NP_064724.1_ PA0004 1 4275 6695
        | MSENNTYDSSSIKVLKGLDAVRKRPGMYIGDTDDGTGLHHMVFEVVDNSI
        | >07700ES100A0cP01_ C ATG TGA 91521 94826 Valid icmF1 3306 _PA0077_NP_248767.1_ PA0077 1 91521 94826
        | MQSLAEVSAPDAASVAT

        :return: False if most of replicon contains only one seaquence, True otherwise
        :rtype: bool
        """
        all_len = [rep.max - rep.min for rep in self._DB.values()]
        replicon_with_one_seq = all_len.count(0)
        if replicon_with_one_seq > len(all_len) * 0.8:
            return False
        else:
            return True

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
        _max = 0
        genes = []
        for seq_id, length, _rank in self._idx:
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
        def grp_replicon(entry):
            """
            in gembase the identifier of fasta sequence follows the following schema:
            <replicon-name>_<seq-name> with eventually '_' inside the <replicon_name>
            but not in the <seq-name>.
            so grp_replicon allow to group sequences belonging to the same replicon.
            """
            return "_".join(entry[0].split('_')[: -1])

        def parse_seq_id(seq_id):
            """
            parse a gemabse sequence id (.idx)
            seq_id has the following format <replicon-name>_<seq-name> with eventually '_' inside the <replicon_name>
            but not in the <seq-name>.
            """
            *replicon_name, seq_name = seq_id.split('_')
            replicon_name = "_".join(replicon_name)
            return replicon_name, seq_name

        replicons = (x[1] for x in groupby(self._idx, grp_replicon))
        for replicon in replicons:
            genes = []
            seq_id, seq_length, _min = next(replicon)
            
            replicon_name, seq_name = parse_seq_id(seq_id)
            genes.append((seq_name, seq_length))
            for seq_id, seq_length, rank in replicon:
                # pass all sequence of the replicon until the last one
                _, seq_name = parse_seq_id(seq_id)
                genes.append((seq_name, seq_length))
            _, seq_name = parse_seq_id(seq_id)
            try:
                _max = rank
            except UnboundLocalError:
                # there is only one sequence for this replicon
                _max = _min
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
