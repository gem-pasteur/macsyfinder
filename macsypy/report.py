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
Extract informations from the results of hmmsearch
"""

import os
import logging
_log = logging.getLogger(__name__)

import abc
from threading import Lock
from itertools import groupby

from .database import Indexes, RepliconDB
from .hit import CoreHit
from .error import MacsypyError


class HMMReport(metaclass=abc.ABCMeta):
    """
    Handle the results from the HMM search. Extract a synthetic report from the raw hmmer output,
    after having applied a hit filtering.
    This class is an **abstract class**. There are two implementations of this abstract class
    depending on whether the input sequence dataset is "ordered" ("gembase" or "ordered_replicon" db_type)
    or not ("unordered" db_type).
    """

    def __init__(self, gene, hmmer_output, cfg):
        """
        :param gene: the gene corresponding to the profile search reported here
        :type gene: :class:`macsypy.gene.CoreGene` object
        :param hmmer_output: The path to the raw Hmmer output file
        :type hmmer_output: string
        :param cfg: the configuration object
        :type cfg: :class:`macsypy.config.Config` object
        """
        self.gene = gene
        self._hmmer_raw_out = hmmer_output
        self._extract_out = None
        self.hits = []
        self.cfg = cfg
        self._lock = Lock()

    @abc.abstractmethod
    def _get_replicon_name(self, hit_id):
        """
        This method is used by extract method and must be implemented by concrete class

        :param str hit_id: the id of the current hit extract from hmm output.
        :return: The name of the replicon
        """
        rep_name = os.path.splitext(os.path.basename(self.cfg.sequence_db()))[0]
        return rep_name


    def extract(self):
        """
        Parse the output file of hmmer compute from an unordered genes base
        and produced a new synthetic report file.
        """

        with self._lock:
            # so the extract of a given HMM output is executed only once per run
            # if this method is called several times the first call induce the parsing of HMM out
            # the other calls do nothing
            if self.hits:
                return

            idx = Indexes(self.cfg)
            my_db = self._build_my_db(self._hmmer_raw_out)
            self._fill_my_db(my_db)

            with open(self._hmmer_raw_out, 'r') as hmm_out:
                i_evalue_sel = self.cfg.i_evalue_sel()
                coverage_threshold = self.cfg.coverage_profile()
                gene_profile_lg = len(self.gene.profile)
                hmm_hits = (x[1] for x in groupby(hmm_out, self._hit_start))
                # drop summary
                next(hmm_hits)
                for hmm_hit in hmm_hits:
                    hit_id = self._parse_hmm_header(hmm_hit)
                    try:
                        seq_lg, position_hit = my_db[hit_id]
                    except TypeError as err:
                        if my_db[hit_id] is None:
                            msg = f"hit id '{hit_id}' was not indexed, rebuild sequence '{idx.name}' index"
                            _log.critical(msg)
                            raise MacsypyError(msg) from err
                    replicon_name = self._get_replicon_name(hit_id)

                    body = next(hmm_hits)
                    c_hit = self._parse_hmm_body(hit_id, gene_profile_lg, seq_lg, coverage_threshold,
                                                 replicon_name, position_hit, i_evalue_sel, body)
                    self.hits += c_hit
                self.hits.sort()
                return self.hits


    def __str__(self):
        """
        :return: string representation of this report
        :rtype: str
        """
        rep = f"""# gene: {self.gene.name} extract from {self._hmmer_raw_out} hmm output
# profile length= {len(self.gene.profile):d}
# i_evalue threshold= {self.cfg.i_evalue_sel():.3f}
# coverage threshold= {self.cfg.coverage_profile():.3f}
# hit_id replicon_name position_hit hit_sequence_length gene_name gene_system i_eval score profile_coverage sequence_coverage begin end
"""

        for c_hit in self.hits:
            rep += str(c_hit)
        return rep


    def save_extract(self):
        """
        Write the string representation of the extract report in a file.
        The name of this file is the concatenation of the gene name and of the "res_extract_suffix"
        from the config object
        """
        with self._lock:
            extract_out_name = self.gene.name + self.cfg.res_extract_suffix()
            self._extract_out = os.path.join(self.cfg.working_dir(), self.cfg.hmmer_dir(), extract_out_name)
            with open(self._extract_out, 'w') as _file:
                _file.write(str(self))


    def best_hit(self):
        """
        Return the best hit among multiple hits
        """
        try:
            return self.hits[0]
        except IndexError:
            return None


    def _hit_start(self, line):
        """
        :param line: the line to parse
        :type line: string
        :return: True if it's the beginning of a new hit in Hmmer raw output files.
         False otherwise
        :rtype: boolean.
        """
        return line.startswith(">>")


    def _build_my_db(self, hmm_output):
        """
        Build the keys of a dictionary object to store sequence identifiers of hits.

        :param hmm_output: the path to the hmmsearch output to parse.
        :type hmm_output: string
        :return: a dictionary containing a key for each sequence id of the hits
        :rtype: dict
        """
        db = {}
        with open(hmm_output) as hmm_file:
            hits = (x[1] for x in groupby(hmm_file, self._hit_start) if x[0])
            for hit in hits:
                db[self._parse_hmm_header(hit)] = None
        return db


    def _fill_my_db(self, db):
        """
        Fill the dictionary with information on the matched sequences

        :param db: the database containing all sequence id of the hits.
        :type db: dict
        """
        idx = Indexes(self.cfg)
        # the indexes are already build
        # just use them
        for seqid, length, rank in idx:
            if seqid in db:
                db[seqid] = (int(length), int(rank))


    def _parse_hmm_header(self, h_grp):
        """
        :param h_grp: the sequence of string return by groupby function representing the header of a hit
        :type h_grp: sequence of string (<itertools._grouper object at 0x7ff9912e3b50>)
        :returns: the sequence identifier from a set of lines that corresponds to a single hit
        :rtype: string
        """
        for line in h_grp:
            hit_id = line.split()[1]
        return hit_id


    def _parse_hmm_body(self, hit_id, gene_profile_lg, seq_lg, coverage_threshold, replicon_name,
                        position_hit, i_evalue_sel, b_grp):
        """
        Parse the raw Hmmer output to extract the hits, and filter them with threshold criteria selected
        ("coverage_profile" and "i_evalue_select" command-line parameters)

        :param str hit_id: the sequence identifier
        :param int gene_profile_lg: the length of the profile matched
        :paramint  seq_lg: the length of the sequence
        :param float coverage_threshold: the minimal coverage of the profile to be reached in the Hmmer alignment
                                        for hit selection.
        :param str replicon_name: the identifier of the replicon
        :param int position_hit: the rank of the sequence matched in the input dataset file
        :param float i_evalue_sel: the maximal i-evalue (independent evalue) for hit selection
        :param b_grp: the Hmmer output lines to deal with (grouped by hit)
        :type b_grp: list of list of strings
        :returns: a sequence of hits
        :rtype: list of :class:`macsypy.report.CoreHit` objects

        """
        first_line = next(b_grp)
        if not first_line.startswith('   #    score'):
            return []
        else:
            hits = []
            for line in b_grp:
                if line[0] == '\n':
                    return hits
                elif line.startswith(" ---   ------ ----- --------"):
                    pass
                else:
                    fields = line.split()
                    try:
                        # fields[2] = score
                        # fields[5] = i_evalue
                        # fields[6] = hmmfrom
                        # fields[7] = hmm to
                        # fields[9] = alifrom
                        # fields[10] = ali to
                        if len(fields) > 1 and float(fields[5]) <= i_evalue_sel:
                            cov_profile = (float(fields[7]) - float(fields[6]) + 1) / gene_profile_lg
                            begin = int(fields[9])
                            end = int(fields[10])
                            cov_gene = (float(end) - float(begin) + 1) / seq_lg  # To be added in Gene: sequence_length
                            if cov_profile >= coverage_threshold:
                                i_eval = float(fields[5])
                                score = float(fields[2])
                                hits.append(CoreHit(self.gene,
                                                    hit_id,
                                                    seq_lg,
                                                    replicon_name,
                                                    position_hit,
                                                    i_eval,
                                                    score,
                                                    cov_profile,
                                                    cov_gene,
                                                    begin,
                                                    end))
                    except ValueError as err:
                        msg = f"Invalid line to parse :{line}:{err}"
                        _log.debug(msg)
                        raise ValueError(msg) from err


class GeneralHMMReport(HMMReport):
    """
    Handle HMM report. Extract a synthetic report from the raw hmmer output.
    Dedicated to any type of 'unordered' datasets.
    """
    def _get_replicon_name(self, hit_id):
        return super()._get_replicon_name(hit_id)



class OrderedHMMReport(HMMReport):
    """
    Handle HMM report. Extract a synthetic report from the raw hmmer output.
    Dedicated to 'ordered_replicon' datasets.
    """
    def _get_replicon_name(self, hit_id):
        return super()._get_replicon_name(hit_id)

class GembaseHMMReport(HMMReport):
    """
    Handle HMM report. Extract a synthetic report from the raw hmmer output.
    Dedicated to 'gembase' format datasets.
    """

    def _get_replicon_name(self, hit_id):
        replicon_name = "_".join(hit_id.split('_')[:-1])
        return replicon_name
