#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2020  Institut Pasteur (Paris) and CNRS.           #
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

from operator import attrgetter
import logging
from dataclasses import dataclass

from macsypy.error import MacsypyError

_log = logging.getLogger(__name__)


class Hit:
    """
    Handle the hits filtered from the Hmmer search. The hits are instanciated by :py:meth:`HMMReport.extract` method
    """


    def __init__(self, gene, model, hit_id, hit_seq_length, replicon_name,
                 position_hit, i_eval, score, profile_coverage, sequence_coverage, begin_match, end_match):
        """
        :param gene: the gene corresponding to this profile
        :type gene: :class:`macsypy.gene.Gene` object
        :param model: the model to which this gene belongs
        :type model: :class:`macsypy.model.Model` object
        :param hit_id: the identifier of the hit
        :type hit_id: string
        :param hit_seq_length: the length of the hit sequence
        :type hit_seq_length: integer
        :param replicon_name: the name of the replicon
        :type replicon_name: string
        :param position_hit: the rank of the sequence matched in the input dataset file
        :type position_hit: integer
        :param i_eval: the best-domain evalue (i-evalue, "independent evalue")
        :type i_eval: float
        :param score: the score of the hit
        :type score: float
        :param profile_coverage: percentage of the profile that matches the hit sequence
        :type profile_coverage: float
        :param sequence_coverage: percentage of the hit sequence that matches the profile
        :type sequence_coverage: float
        :param begin_match: where the hit with the profile starts in the sequence
        :type begin_match: integer
        :param end_match: where the hit with the profile ends in the sequence
        :type end_match: integer
        """
        self.gene = gene
        self.model = model
        self.id = hit_id
        self.seq_length = hit_seq_length
        self.replicon_name = replicon_name
        self.position = position_hit
        self.i_eval = i_eval
        self.score = score
        self.profile_coverage = profile_coverage
        self.sequence_coverage = sequence_coverage
        self.begin_match = begin_match
        self.end_match = end_match
        self._systems = set()

    def __hash__(self):
        """To be hashable, it's needed to be put in a set or used as dict key"""
        return id(self)


    def __str__(self):
        """
        Print useful information on the Hit: regarding Hmmer statistics, and sequence information
        """
        return f"{self.id}\t{self.replicon_name}\t{self.position:d}\t{self.seq_length:d}\t{self.gene.name}\t" \
               f"{self.model.name}\t{self.i_eval:.3e}\t{self.score:.3f}\t{self.profile_coverage:.3f}\t" \
               f"{self.sequence_coverage:.3f}\t{self.begin_match:d}\t{self.end_match:d}\n"


    def __lt__(self, other):
        """
        compare two Hits. If the sequence identifier is the same, do the comparison on the score.
        Otherwise, do it on alphabetical comparison of the sequence identifier.

        :param other: the hit to compare to the current object
        :type other: :class:`macsypy.report.Hit` object
        :return: True if self is < other, False otherwise
        """
        if self.id == other.id:
            if not self.gene.is_homolog(other.gene):
                _log.warning(f"Non homologs match: {self.gene.name} {self.model.name}) {other.gene.name} "
                             f"({other.model.name}) for {self.id}")
            return self.score < other.score
        else:
            return self.id < other.id


    def __gt__(self, other):
        """
        compare two Hits. If the sequence identifier is the same, do the comparison on the score.
        Otherwise, do it on alphabetical comparison of the sequence identifier.

        :param other: the hit to compare to the current object
        :type other: :class:`macsypy.report.Hit` object
        :return: True if self is > other, False otherwise
        """
        if self.id == other.id:
            if not self.gene.is_homolog(other.gene):
                _log.warning(f"Non homologs match: {self.gene.name} ({self.model.name}) {other.gene.name} "
                             f"({other.model.name}) for {self.id}")
            return self.score > other.score
        else:
            return self.id > other.id


    def __eq__(self, other):
        """
        Return True if two hits are totally equivalent, False otherwise.

        :param other: the hit to compare to the current object
        :type other: :class:`macsypy.report.Hit` object
        :return: the result of the comparison
        :rtype: boolean
        """
        epsilon = 0.001
        return (self.gene.name == other.gene.name and
                self.model.name == other.model.name and
                self.id == other.id and
                self.seq_length == other.seq_length and
                self.replicon_name == other.replicon_name and
                self.position == other.position and
                abs(self.i_eval - other.i_eval) <= epsilon and
                abs(self.score - other.score) <= epsilon and
                abs(self.profile_coverage - other.profile_coverage) <= epsilon and
                abs(self.sequence_coverage - other.sequence_coverage) <= epsilon and
                self.begin_match == other.begin_match and
                self.end_match == other.end_match
                )


    def get_position(self):
        """
        :returns: the position of the hit (rank in the input dataset file)
        :rtype: integer
        """
        return self.position


    def get_syst_inter_gene_max_space(self):
        """
        :returns: the 'inter_gene_max_space' parameter defined for the gene of the hit
        :rtype: integer
        """
        return self.gene.model.inter_gene_max_space


class ValidHit:
    """
    Encapsulates a :class:`macsypy.report.Hit`
    This class stores a Hit that has been attributed to a putative system.
    Thus, it also stores:

    - the system,
    - the status of the gene in this system, ('mandatory', 'accessory', ...
    - the gene in the model for which it's an occurrence
    """

    def __init__(self, hit, gene_ref, gene_status):
        """
        :param hit:
        :type hit: :class:`macsypy.hit.Hit` object
        :param gene_ref:
        :type gene_ref: :class:`macsypy.gene.Gene` object
        :param gene_status:
        :type gene_status: :class:`macsypy.gene.GeneStatus` object
        """
        self.hit = hit
        self.gene_ref = gene_ref
        self.status = gene_status


    def __getattr__(self, item):
        return getattr(self.hit, item)


@dataclass(frozen=True)
class HitWeight:

    hitself: float = 1
    homolog: float = 0.75
    analog: float = 0.75
    mandatory: float = 1
    accessory: float = 0.5


hit_weight = HitWeight()


def get_best_hits(hits, key='score'):
    """
    If several hits match the same protein, keep only the best match based either on

        * score
        * i_evalue
        * profile_coverage

    :param hits: the hits to filter, all hits must match the same protein.
    :type hits: [ :class:`macsypy.hit.Hit` object, ...]
    :param str key: The criterion used to select the best hit 'score', i_evalue', 'profile_coverage'
    :return: the list of the best hits
    :rtype: [ :class:`macsypy.hit.Hit` object, ...]
    """
    hits_register = {}
    for hit in hits:
        register_key = hit.replicon_name, hit.position
        if register_key in hits_register:
            hits_register[register_key].append(hit)
        else:
            hits_register[register_key] = [hit]

    best_hits = []
    for hits_on_same_prot in hits_register.values():
        if key == 'score':
            hits_on_same_prot.sort(key=attrgetter(key), reverse=True)
        elif key == 'i_eval':
            hits_on_same_prot.sort(key=attrgetter(key))
        elif key == 'profile_coverage':
            hits_on_same_prot.sort(key=attrgetter(key), reverse=True)
        else:
            raise MacsypyError(f'The criterion for Hits comparison {key} does not exist or is not available.\n'
                               f'It must be either "score", "i_eval" or "profile_coverage".')
        best_hits.append(hits_on_same_prot[0])
    return best_hits
