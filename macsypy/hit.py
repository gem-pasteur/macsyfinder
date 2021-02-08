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

from operator import attrgetter
import logging
from dataclasses import dataclass

from macsypy.gene import ModelGene
from macsypy.error import MacsypyError

_log = logging.getLogger(__name__)


class Hit:
    """
    Handle the hits filtered from the Hmmer search. The hits are instanciated by :py:meth:`HMMReport.extract` method
    """


    def __init__(self, gene, hit_id, hit_seq_length, replicon_name,
                 position_hit, i_eval, score, profile_coverage, sequence_coverage, begin_match, end_match):
        """
        :param gene: the gene corresponding to this profile
        :type gene: :class:`macsypy.gene.CoreGene` object
        :param str hit_id: the identifier of the hit
        :param int hit_seq_length: the length of the hit sequence
        :param str replicon_name: the name of the replicon
        :param int position_hit: the rank of the sequence matched in the input dataset file
        :param float i_eval: the best-domain evalue (i-evalue, "independent evalue")
        :param float score: the score of the hit
        :param float profile_coverage: percentage of the profile that matches the hit sequence
        :param float sequence_coverage: percentage of the hit sequence that matches the profile
        :param int begin_match: where the hit with the profile starts in the sequence
        :param int end_match: where the hit with the profile ends in the sequence
        """
        self.gene = gene
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
        return hash((self.gene.name, self.id, self.seq_length, self.position, self.i_eval))


    def __str__(self):
        """
        :return: Useful information on the Hit: regarding Hmmer statistics, and sequence information
        :rtype: str
        """
        return f"{self.id}\t{self.replicon_name}\t{self.position:d}\t{self.seq_length:d}\t{self.gene.name}\t" \
               f"{self.i_eval:.3e}\t{self.score:.3f}\t{self.profile_coverage:.3f}\t" \
               f"{self.sequence_coverage:.3f}\t{self.begin_match:d}\t{self.end_match:d}\n"


    def __lt__(self, other):
        """
        Compare two Hits. If the sequence identifier is the same, do the comparison on the score.
        Otherwise, do it on alphabetical comparison of the sequence identifier.

        :param other: the hit to compare to the current object
        :type other: :class:`macsypy.report.Hit` object
        :return: True if self is < other, False otherwise
        """
        if self.id == other.id:
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
        :param hit: a match between a hmm profile and a replicon
        :type hit: :class:`macsypy.hit.Hit` object
        :param gene_ref: The ModelGene link to this hit
                         The ModeleGene have the same name than the CoreGene
                         But one hit can be link to several ModelGene (several Model)
                         To know for what gene this hit play role use the
                         :meth:`macsypy.gene.ModelGene.alternate_of` ::

                            hit.gene_ref.alternate_of()

        :type gene_ref: :class:`macsypy.gene.ModelGene` object
        :param gene_status:
        :type gene_status: :class:`macsypy.gene.GeneStatus` object
        """
        self.hit = hit
        if not isinstance(gene_ref, ModelGene):
            raise MacsypyError(f"The ValidHit 'gene_ref' argument must be a ModelGene not {type(gene_ref)}.")
        self.gene_ref = gene_ref
        self.status = gene_status

    @property
    def multi_system(self):
        return self.gene_ref.multi_system

    def __getattr__(self, item):
        return getattr(self.hit, item)

    def __gt__(self, other):
        return self.hit > other

    def __eq__(self, other):
        return self.hit == other and self.gene_ref.name == self.gene_ref.name

    def __lt__(self, other):
        return self.hit < other


@dataclass(frozen=True)
class HitWeight:
    """
    The weight to compute the cluster and system score
    see user documentation macsyfinder functionning for further details
    by default

        * itself = 1
        * exchangeable = 0.8

        * mandatory = 1
        * accessory = 0.5
        * neutral = 0

        * loner_multi_system = 0.7
    """
    itself: float = 1
    exchangeable: float = 0.8
    mandatory: float = 1
    accessory: float = 0.5
    neutral: float = 0
    loner_multi_system: float = 0.7


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
