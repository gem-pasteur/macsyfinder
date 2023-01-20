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

import abc
from operator import attrgetter
import logging
from dataclasses import dataclass

from macsypy.gene import ModelGene
from macsypy.error import MacsypyError

_log = logging.getLogger(__name__)


class CoreHit:
    """
    Handle the hits filtered from the Hmmer search.
    The hits are instanciated by :py:meth:`HMMReport.extract` method
    In one run of MacSyFinder, there exists only one CoreHit per gene
    These hits are independent of any :class:`macsypy.model.Model` instance.
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
        :return: Useful information on the CoreHit: regarding Hmmer statistics, and sequence information
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
        :type other: :class:`macsypy.report.CoreHit` object
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
        :type other: :class:`macsypy.report.CoreHit` object
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
        :type other: :class:`macsypy.report.CoreHit` object
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


class ModelHit:
    """
    Encapsulates a :class:`macsypy.report.CoreHit`
    This class stores a CoreHit that has been attributed to a putative system.
    Thus, it also stores:

    - the system,
    - the status of the gene in this system, ('mandatory', 'accessory', ...
    - the gene in the model for which it's an occurrence

    for one gene it can exist several ModelHit instance one for each Model containing this gene
    """

    def __init__(self, hit, gene_ref, gene_status):
        """
        :param hit: a match between a hmm profile and a replicon
        :type hit: :class:`macsypy.hit.CoreHit` object
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
        if not isinstance(hit, CoreHit):
            raise MacsypyError(f"The {self.__class__.__name__} 'hit' argument must be a CoreHit not {type(hit)}.")
        self._hit = hit
        if not isinstance(gene_ref, ModelGene):
            raise MacsypyError(f"The {self.__class__.__name__} 'gene_ref' argument must be a ModelGene "
                               f"not {type(gene_ref)}.")
        self.gene_ref = gene_ref
        self.status = gene_status


    def __str__(self):
        return str(self._hit)


    def __hash__(self):
        """To be hashable, it's needed to be put in a set or used as dict key"""
        return hash((hash(self.hit), self.gene_ref.model.fqn))


    @property
    def hit(self):
        """
        :return: The CoreHit below this ModelHit
        :rtype: :class:`macsypy.hit.CoreHit` oject
        """
        return self._hit


    @property
    def multi_system(self):
        """
        :return: True if the hit represent a `multi_system` :class:`macsypy.Gene.ModelGene`, False otherwise.
        :rtype: bool
        """
        return self.gene_ref.multi_system


    @property
    def multi_model(self):
        """
        :return: True if the hit represent a `multi_model` :class:`macsypy.Gene.ModelGene`, False otherwise.
        :rtype: bool
        """
        return self.gene_ref.multi_model


    @property
    def loner(self):
        """
        :return: True if the hit represent a `loner` :class:`macsypy.Gene.ModelGene`, False otherwise.
                 A True Loner is a hit representing a gene with the attribute loner and which does not include in a cluster.

                 - a hit representing a loner gene but include in a cluster is not a true loner
                 - a hit which is not include with other gene in a cluster but does not represent a gene loner is not a
                   True loner (This situation may append when min_genes_required = 1)
        :rtype: bool
        """
        return False


    def __getattr__(self, item):
        try:
            return getattr(self._hit, item)
        except AttributeError:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{item}'") from None

    def __gt__(self, other):
        return self._hit > other

    def __eq__(self, other):
        return self._hit == other and self.gene_ref.name == self.gene_ref.name

    def __lt__(self, other):
        return self._hit < other

    @property
    def counterpart(self):
        # used polymophism in serialization, loop over all ModelHit
        # instead of testing the type of ModelHit to now if there is a counterpart method
        # defined one in the ModelHit which return empty counterpart
        return []


class AbstractCounterpartHit(ModelHit, metaclass=abc.ABCMeta):
    """
    Abstract Class to handle ModelHit wit equivalent for instance Loner or MultiSystem hit
    """

    def __init__(self, hit, gene_ref=None, gene_status=None, counterpart=None):
        if isinstance(hit, CoreHit) and not (gene_ref and gene_status):
            raise MacsypyError(f"Cannot Create a {self.__class__.__name__} hit from "
                               f"CoreHit ({hit.gene.name}, {hit.position}) "
                                "without specifying 'gene_ref' and 'gene_status'")
        elif isinstance(hit, CoreHit):
            super().__init__(hit, gene_ref, gene_status)
        elif isinstance(hit, ModelHit):
            super().__init__(hit.hit, gene_ref=hit.gene_ref, gene_status=hit.gene_ref.status)
        self.counterpart = counterpart


    def __getattr__(self, item):
        return getattr(self._hit, item)


    @property
    def counterpart(self):
        """
        :return: The set of hits that can play the same role
        """
        return set(self._counterpart)


    @counterpart.setter
    def counterpart(self, counterparts):
        """

        :param counterpart:
        """
        if not counterparts:
            self._counterpart = set()
        elif all([hit.gene_ref.alternate_of().name is self.gene_ref.alternate_of().name for hit in counterparts]):
            self._counterpart = set(counterparts)
        else:
            msg = f"Try to set counterpart for hit '{self.gene_ref.name}' with non compatible hits: " \
                  f"{[hit.gene_ref.name for hit in counterparts]}"
            _log.error(msg)
            raise MacsypyError(msg)

    def __str__(self):
        ch_str = str(self._hit)[:-1]

        return ch_str + '\t' + ','.join([h.id for h in self.counterpart])


    def __len__(self):
        return len(self.counterpart) + 1

    @property
    def loner(self):
        return False

    @property
    def multi_system(self):
        return False


class Loner(AbstractCounterpartHit):
    """
    Handle hit which encode for a gene tagged as loner and which not clustering with other hit.
    """

    def __init__(self, hit, gene_ref=None, gene_status=None, counterpart=None):
        """
        hit that is outside a cluster, the gene_ref is a loner

        :param hit: a match between a hmm profile and a replicon
        :type hit: :class:`macsypy.hit.CoreHit` object
        :param gene_ref: The ModelGene link to this hit
                         The ModeleGene have the same name than the CoreGene
                         But one hit can be link to several ModelGene (several Model)
                         To know for what gene this hit play role use the
                         :meth:`macsypy.gene.ModelGene.alternate_of` ::

                            hit.gene_ref.alternate_of()

        :type gene_ref: :class:`macsypy.gene.ModelGene` object
        :param gene_status:
        :type gene_status: :class:`macsypy.gene.GeneStatus` object
        :param counterpart: the other occurence of the gene or exchangeable in the replicon
        :type counterpart: list of :class:`macsypy.hit.CoreHit`
        """
        super().__init__(hit, gene_ref=gene_ref, gene_status=gene_status, counterpart=counterpart)

        if not self.gene_ref.loner:
            msg = f"{hit.id} cannot be a loner gene_ref '{gene_ref.name}' not tag as loner"
            _log.critical(msg)
            raise MacsypyError(msg)


    @property
    def loner(self):
        return True


class MultiSystem(AbstractCounterpartHit):
    """
    Handle hit which encode for a gene tagged as loner and which not clustering with other hit.
    """

    def __init__(self, hit, gene_ref=None, gene_status=None, counterpart=None):
        """
        hit that is outside a cluster, the gene_ref is a loner

        :param hit: a match between a hmm profile and a replicon
        :type hit: :class:`macsypy.hit.CoreHit` object
        :param gene_ref: The ModelGene link to this hit
                         The ModeleGene have the same name than the CoreGene
                         But one hit can be link to several ModelGene (several Model)
                         To know for what gene this hit play role use the
                         :meth:`macsypy.gene.ModelGene.alternate_of` ::

                            hit.gene_ref.alternate_of()

        :type gene_ref: :class:`macsypy.gene.ModelGene` object
        :param gene_status:
        :type gene_status: :class:`macsypy.gene.GeneStatus` object
        :param counterpart: the other occurence of the gene or exchangeable in the replicon
        :type counterpart: list of :class:`macsypy.hit.CoreHit`
        """
        super().__init__(hit, gene_ref=gene_ref, gene_status=gene_status, counterpart=counterpart)

        if not self.gene_ref.multi_system:
            msg = f"{hit.id} cannot be a multi systems, gene_ref '{gene_ref.name}' not tag as multi_system"
            _log.critical(msg)
            raise MacsypyError(msg)


    @property
    def multi_system(self):
        return True


class LonerMultiSystem(Loner, MultiSystem):
    """
    Handle hit which encode for a gene
     * gene tagged as multi-system
     * and gene tagged as loner also
     * and the hit do not clustering with other hits.
    """

    def __init__(self, hit, gene_ref=None, gene_status=None, counterpart=None):
        """
        hit that is outside a cluster, the gene_ref is loner and multi_system

        :param hit: a match between a hmm profile and a replicon
        :type hit: :class:`macsypy.hit.CoreHit` | :class:`macsypy.hit.ModelHit` | :class:`macsypy.hit.MultiSystem` object
        :param gene_ref: The ModelGene link to this hit
                         The ModeleGene have the same name than the CoreGene
                         But one hit can be link to several ModelGene (several Model)
                         To know for what gene this hit play role use the
                         :meth:`macsypy.gene.ModelGene.alternate_of` ::

                            hit.gene_ref.alternate_of()

        :type gene_ref: :class:`macsypy.gene.ModelGene` object
        :param gene_status:
        :type gene_status: :class:`macsypy.gene.GeneStatus` object
        :param counterpart: the other occurence of the gene or exchangeable in the replicon
        :type counterpart: list of :class:`macsypy.hit.CoreHit`
        """
        if isinstance(hit, (Loner, MultiSystem)):
            super().__init__(hit,
                             gene_ref=hit.gene_ref,
                             gene_status=hit.status,
                             counterpart=hit.counterpart)
        else:
            super().__init__(hit, gene_ref=gene_ref, gene_status=gene_status, counterpart=counterpart)


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

        * out_of_cluster = 0.7
    """
    itself: float = 1
    exchangeable: float = 0.8
    mandatory: float = 1
    accessory: float = 0.5
    neutral: float = 0
    out_of_cluster: float = 0.7


def get_best_hit_4_func(function, hits, key='score'):
    """
    select the best Loner among several ones encoding for same function

        * score
        * i_evalue
        * profile_coverage

    :param str function: the name of the function fulfill by the hits (all hits must have same function)
    :param hits: the hits to filter.
    :type hits: sequence of :class:`macsypy.hit.ModelHit` object
    :param str key: The criterion used to select the best hit 'score', i_evalue', 'profile_coverage'
    :return: the best hit
    :rtype: :class:`macsypy.hit.ModelHit` object
    """
    originals = []
    exchangeables = []
    for hit in hits:
        if hit.gene_ref.name == function:
            originals.append(hit)
        else:
            exchangeables.append(hit)
    if originals:
        hits = originals
    else:
        hits = exchangeables
    if key == 'score':
        hits.sort(key=attrgetter(key), reverse=True)
    elif key == 'i_eval':
        hits.sort(key=attrgetter(key))
    elif key == 'profile_coverage':
        hits.sort(key=attrgetter(key), reverse=True)
    else:
        raise MacsypyError(f'The criterion for Loners comparison {key} does not exist or is not available.\n')
    return hits[0]


def sort_model_hits(model_hits):
    """
    Sort :class:`macsypy.hit.ModelHit` per function

    :param model_hits: a sequence of :class:`macsypy.hit.ModelHit`
    :return: dict {str function name: [model_hit, ...] }
    """
    ms_registry = {}
    for hit in model_hits:
        func_name = hit.gene_ref.alternate_of().name
        if func_name in ms_registry:
            ms_registry[func_name].append(hit)
        else:
            ms_registry[func_name] = [hit]

    return ms_registry


def compute_best_MSHit(ms_registry):
    """

    :param ms_registry:
    :return:
    """

    best_multisystem_hits = []
    for func_name in ms_registry:
        equivalent_ms = ms_registry[func_name]
        best_ms = get_best_hit_4_func(func_name, equivalent_ms, key='score')
        equivalent_ms.remove(best_ms)
        best_ms.counterpart = equivalent_ms
        best_multisystem_hits.append(best_ms)
    return best_multisystem_hits


def get_best_hits(hits, key='score'):
    """
    If several hits match the same protein, keep only the best match based either on

        * score
        * i_evalue
        * profile_coverage

    :param hits: the hits to filter, all hits must match the same protein.
    :type hits: [ :class:`macsypy.hit.CoreHit` object, ...]
    :param str key: The criterion used to select the best hit 'score', i_evalue', 'profile_coverage'
    :return: the list of the best hits
    :rtype: [ :class:`macsypy.hit.CoreHit` object, ...]
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
