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
import itertools
import statistics
from itertools import chain
from operator import attrgetter
import logging
_log = logging.getLogger(__name__)

from.gene import GeneStatus
from .cluster import Cluster
from .error import MacsypyError


class MatchMaker(metaclass=abc.ABCMeta):
    """
    Is an abstract class for (Odered/Unordered)MatchMaker
    the `match` class method must be implemented in concrete classes
    """

    def __init__(self, model):
        self._model = model
        # init my structures to count gene occurrences
        self.mandatory_counter = {g.name: 0 for g in model.mandatory_genes}
        self.exchangeable_mandatory = self._create_exchangeable_map(model.mandatory_genes)

        self.accessory_counter = {g.name: 0 for g in model.accessory_genes}
        self.exchangeable_accessory = self._create_exchangeable_map(model.accessory_genes)

        self.forbidden_counter = {g.name: 0 for g in model.forbidden_genes}
        self.exchangeable_forbidden = self._create_exchangeable_map(model.forbidden_genes)

        self.neutral_counter = {g.name: 0 for g in model.neutral_genes}
        self.exchangeable_neutral = self._create_exchangeable_map(model.neutral_genes)


    def _create_exchangeable_map(self, genes):
        """
        create a map between an exchangeable (formly homolog or analog) gene name and it's gene reference

        :param genes: The genes to get the exchangeable genes
        :type genes: list of :class:`macsypy.gene.ModelGene` objects
        :rtype: a dict with keys are the exchangeable gene_name and the value the reference gene
        """
        map = {}
        for gene in genes:
            for ex_gene in gene.exchangeables:
                map[ex_gene.name] = gene
        return map


    def sort_hits_by_status(self, hits):
        """
        sort :class:`macsypy.hit.ModelHit` according the
        the status of the gene the hit code for.

        :param hits: list of :class:`macsypy.hit.ModelHit` object
        :return: the valid hits according their status
        :rtype: a tuple of 4 lists

              * :class:`macsypy.hit.ModelHit` for MANDATORY genes
              * :class:`macsypy.hit.ModelHit` for ACCESSORY genes
              * :class:`macsypy.hit.ModelHit` for NEUTRAL genes
              * :class:`macsypy.hit.ModelHit` for FORBIDDEN genes

        """
        mandatory_hits = []
        accessory_hits = []
        neutral_hits = []
        forbidden_hits = []
        for hit in hits:
            gene_name = hit.gene.name
            # the ModelHit need to be linked to the
            # gene of the model
            if gene_name in self.mandatory_counter:
                self.mandatory_counter[hit.gene.name] += 1
                mandatory_hits.append(hit)
            elif gene_name in self.exchangeable_mandatory:
                gene_ref = self.exchangeable_mandatory[gene_name]
                self.mandatory_counter[gene_ref.name] += 1
                mandatory_hits.append(hit)

            elif gene_name in self.accessory_counter:
                self.accessory_counter[gene_name] += 1
                accessory_hits.append(hit)
            elif gene_name in self.exchangeable_accessory:
                gene_ref = self.exchangeable_accessory[gene_name]
                self.accessory_counter[gene_ref.name] += 1
                accessory_hits.append(hit)

            elif gene_name in self.neutral_counter:
                self.neutral_counter[gene_name] += 1
                neutral_hits.append(hit)
            elif gene_name in self.exchangeable_neutral:
                gene_ref = self.exchangeable_neutral[gene_name]
                self.neutral_counter[gene_ref.name] += 1
                neutral_hits.append(hit)

            elif gene_name in self.forbidden_counter:
                self.forbidden_counter[gene_name] += 1
                forbidden_hits.append(hit)
            elif gene_name in self.exchangeable_forbidden:
                gene_ref = self.exchangeable_forbidden[gene_name]
                self.forbidden_counter[gene_ref.name] += 1
                forbidden_hits.append(hit)
            else:
                model = hit.gene_ref.model
                msg = f"Gene '{gene_name}' not found in model '{model.fqn}'"
                _log.critical(msg)
                raise MacsypyError(msg)

        return mandatory_hits, accessory_hits, neutral_hits, forbidden_hits


    def present_genes(self):
        """
        :return: the lists of genes name in model which are present in the replicon (included exchangeable)
        :rtype: tuple of 4 lists for mandatory, accessory, neutral and forbidden
                ([str gene_name, ...], [str gene_name], [str gene_name], [str gene_name])
        """
        mandatory_genes = [g_name for g_name, occ in self.mandatory_counter.items() if occ > 0]
        accessory_genes = [g_name for g_name, occ in self.accessory_counter.items() if occ > 0]
        neutral_genes = [g_name for g_name, occ in self.neutral_counter.items() if occ > 0]
        forbidden_genes = [g_name for g_name, occ in self.forbidden_counter.items() if occ > 0]
        return mandatory_genes, accessory_genes, neutral_genes, forbidden_genes


    @abc.abstractmethod
    def match(self, clusters):
        pass


class OrderedMatchMaker(MatchMaker):
    """
    check if a set of hits match the quorum for ordered replicons (ordered_replicon or gembase)
    """

    def __init__(self, model, redundancy_penalty):
        self._redundancy_penalty = redundancy_penalty
        super().__init__(model)


    def match(self, clusters):
        """
        Check a set of clusters fill model constraints.
        If yes create a :class:`macsypy.system.System` otherwise create
        a :class:`macsypy.cluster.RejectedCandidate`.

        :param clusters: The list of cluster to check if fit the model
        :type clusters: list of :class:`macsypy.cluster.Cluster` objects
        :return: either a System or a RejectedCandidates
        :rtype: :class:`macsypy.system.System` or :class:`macsypy.system.RejectedCandidate` object
        """
        # count the hits
        # and track for each hit for which gene it counts for
        valid_clusters = []
        forbidden_hits = []
        for cluster in clusters:
            # sort model hits between forbidden and the other
            *one_clst_allowed_hits, one_clst_forbidden_hits = self.sort_hits_by_status(cluster.hits)
            if one_clst_forbidden_hits:
                forbidden_hits.extend(one_clst_forbidden_hits)
            # merge MANDATORY, ACCESSORY, NEUTRAL ModelHit
            one_clst_allowed_hits = [mh for hits in one_clst_allowed_hits for mh in hits]
            one_clst_allowed_hits.sort(key=attrgetter('position'))
            valid_clusters.append(Cluster(one_clst_allowed_hits + one_clst_forbidden_hits,
                                          self._model, cluster._hit_weights))

        mandatory_genes, accessory_genes, neutral_genes, forbidden_genes = self.present_genes()
        # the count is finished
        # check if the quorum is reached
        # count how many different genes are represented in the clusters
        # the neutral genes belong to the cluster
        # but they do not count for the quorum

        _log.debug("#" * 50)
        _log.debug(f"mandatory_genes: {mandatory_genes}")
        _log.debug(f"accessory_genes: {accessory_genes}")
        _log.debug(f"neutral_genes: {neutral_genes}")
        _log.debug(f"forbidden_genes: {forbidden_genes}")

        reasons = []
        is_a_system = True
        if forbidden_genes:
            is_a_system = False
            msg = f"There is {len(forbidden_hits)} forbidden genes occurrence(s):" \
                  f" {', '.join(h.gene.name for h in forbidden_hits)}"
            reasons.append(msg)
            _log.debug(msg)
        if len(mandatory_genes) < self._model.min_mandatory_genes_required:
            is_a_system = False
            msg = f'The quorum of mandatory genes required ({self._model.min_mandatory_genes_required}) ' \
                  f'is not reached: {len(mandatory_genes)}'
            reasons.append(msg)
            _log.debug(msg)
        if len(accessory_genes) + len(mandatory_genes) < self._model.min_genes_required:
            is_a_system = False
            msg = f'The quorum of genes required ({self._model.min_genes_required}) is not reached:' \
                  f' {len(accessory_genes) + len(mandatory_genes)}'
            reasons.append(msg)
            _log.debug(msg)

        if is_a_system:
            res = System(self._model, valid_clusters, self._redundancy_penalty)
            _log.debug("is a system")
        else:
            res = RejectedCandidate(self._model, valid_clusters, reasons)
        _log.debug("#" * 50)
        return res


class UnorderedMatchMaker(MatchMaker):

    def match(self, hits):
        """
        :param hits:
        :return:
        """
        # count the hits
        # and track for each hit for which gene it counts for
        mandatory_hits, accessory_hits, neutral_hits, forbidden_hits = self.sort_hits_by_status(hits)
        # the count is finished
        # check if the quorum is reached
        # count how many different genes are represented in the clusters
        # the neutral genes belong to the cluster
        # but they do not count for the quorum
        mandatory_genes, accessory_genes, neutral_genes, forbidden_genes = self.present_genes()

        _log.debug("#" * 50)
        _log.debug(f"mandatory_genes: {mandatory_genes}")
        _log.debug(f"accessory_genes: {accessory_genes}")
        _log.debug(f"neutral_genes: {neutral_genes}")
        _log.debug(f"forbidden_genes: {forbidden_genes}")

        is_a_potential_system = True
        reasons = []
        if forbidden_genes:
            msg = f"There is {len(forbidden_hits)} forbidden genes occurrence(s):" \
                  f" {', '.join(h.gene.name for h in forbidden_hits)}"
            _log.debug(msg)
        if len(mandatory_genes) < self._model.min_mandatory_genes_required:
            is_a_potential_system = False
            msg = f'The quorum of mandatory genes required ({self._model.min_mandatory_genes_required}) is not reached: ' \
                  f'{len(mandatory_genes)}'
            reasons.append(msg)
            _log.debug(msg)
        if len(accessory_genes) + len(mandatory_genes) < self._model.min_genes_required:
            is_a_potential_system = False
            msg = f'The quorum of genes required ({self._model.min_genes_required}) is not reached:' \
                  f' {len(accessory_genes) + len(mandatory_genes)}'
            reasons.append(msg)
            _log.debug(msg)

        if is_a_potential_system:
            res = LikelySystem(self._model, mandatory_hits, accessory_hits, neutral_hits, forbidden_hits)
            _log.debug("There is a genetic potential for a system")
        elif any((mandatory_hits, accessory_hits, neutral_hits)):
            res = UnlikelySystem(self._model, mandatory_hits, accessory_hits, neutral_hits, forbidden_hits, reasons)
        else:
            res = None
        _log.debug("#" * 50)
        return res


class HitSystemTracker(dict):
    """
    track in which system is implied each hit
    """

    def __init__(self, systems):
        super(HitSystemTracker, self).__init__()
        for system in systems:
            m_hits = system.hits
            for m_hit in m_hits:
                c_hit = m_hit.hit
                if c_hit not in self:
                    self[c_hit] = set()
                self[c_hit].add(system)


class MetaSetOfHits(abc.ABCMeta):
    """
    This metaclass control the AbstractSetOfHits class creation.
    In this metaclass we inject on the fly several attributes and properties
    two private attributes and one public property corresponding to each value
    of _supported_status class attribute defined in the concrete classes.
    for instance for System class

     * the attributes
        * self._mandatory
        * self._mandatory_occ
        * self._accessory
        * self._accessory_occ
        * self._neutral
        * self._neutral_occ
     * and the properties
        * mandatory
        * accessory
        * neutral

    are automatically injected

    The value for attributes `_<status>_occ` are filled by the `count` method
    which is defined in AbstractSetOfHits
    """

    def getter_maker(status):
        """
        Create a property which allow to access to the gene corresponding of the cat of the model

        :param str cat: the type of gene category to which we create the getter
        :return: unbound method
        """
        def getter(self):
            occ = getattr(self, f"_{status}_occ")
            return {k: v for k, v in occ.items()}
        return getter

    def __call__(cls, *args, **kwargs):
        new_system_inst = super().__call__(*args, **kwargs)
        # print()
        # print("##################### MetaSystem __call__ #################")
        # print("### new_system_inst", new_system_inst)
        # print("### new_system_inst._supported_status", new_system_inst._supported_status)
        for status in [str(s) for s in new_system_inst._supported_status]:
            # set the private attribute in the Model instance
            setattr(new_system_inst, f"_{status}_occ", {})
            # set the public property in the Model class
            setattr(cls, f"{status}_occ", property(MetaSetOfHits.getter_maker(status)))
        # count fill _mandatory_occ, accessory_occ, ... attributes
        new_system_inst.count()
        return new_system_inst


class AbstractSetOfHits(metaclass=MetaSetOfHits):
    """
    Is the mother class of  System, RejectedCandidates, LikelySystems UnlikelySystem, ...
    """

    def __init__(self, model):
        self.model = model

    def _sort_hits(self, hits):
        hits = sorted(hits, key=attrgetter('position'))
        return hits

    @property
    @abc.abstractmethod
    def hits(self):
        pass


    @property
    def replicon_name(self):
        """
        :return: The name of the replicon
        :rtype: str
        """
        return self._replicon_name


    @property
    def position(self):
        """
        :return: The position of the first and last hit,
                 excluded the hit coding for loners.
                 If the system is composed only by loners, used loners to compute position
        :rtype: tuple (start: int, end:int)
        """
        # hits are sorted by their positions
        hits = [h.position for h in self.hits if not h.gene_ref.loner]
        if hits:
            hits.sort()
            pos = hits[0], hits[-1]
        else:
            # there are only loners
            # take them
            pos = self.hits[0].position, self.hits[-1].position
        return pos


    def count(self):
        """
        fill structures one for supported status mandatory, accessory, ...
        each structure count how many hit for each gene of the model
        mandatory_occ = { gene_name : [ModelHit, ...]
        :return: None
        """
        for status in [str(s) for s in self._supported_status]:
            setattr(self,
                    f"_{status}_occ",
                    {g.name: [] for g in getattr(self.model, f"{status}_genes")}
                    )
        # all the hits are ModelHit
        for hit in self.hits:
            name = hit.gene_ref.alternate_of().name
            status = str(hit.status)  # transform gene status in lower string
            try:
                getattr(self, f"_{status}_occ")[name].append(hit)
            except AttributeError:
                pass
            except KeyError:
                raise MacsypyError(f"gene '{name}' does not belong to '{status}' genes in model '{self.model.name}'")


    @property
    def wholeness(self):
        """

        :return: a score indicating the genes ratio of the model which have at least one hit
                 by default full system is mandatory + accessory ('neutral' genes do not count)
                 but for special corner case it can be sepcified in model definition (xml)
                 or on the command line
        :rtype: float
        """
        # model completude
        # the neutral hit do not participate to the model completude
        score = sum([1 for hits in chain(self._mandatory_occ.values(), self._accessory_occ.values()) if hits]) / \
                self.model.max_nb_genes
        return score


class AbstractClusterizedHits(AbstractSetOfHits):
    """

    """

    def __init__(self, model, clusters):
        if isinstance(clusters, Cluster):
            self.clusters = [clusters]
        else:
            self.clusters = clusters
        self._replicon_name = clusters[0].replicon_name
        super().__init__(model)


    def fulfilled_function(self, *genes):
        """

        :param gene: The genes which must be tested.
        :type genes: :class:`macsypy.gene.ModelGene` object or string representing the gene name
        :return: the common functions between genes and this system.
        :rtype: set of string
        """
        # we do not filter out neutral from the model
        common_functions = set()
        for cluster in self.clusters:
            common_functions.update(cluster.fulfilled_function(*genes))
        return common_functions


class System(AbstractClusterizedHits):
    """
    Modelize as system. a system is an occurrence of a given model on a replicon.
    """

    _supported_status = (GeneStatus.MANDATORY,
                         GeneStatus.ACCESSORY,
                         GeneStatus.NEUTRAL)

    _id = itertools.count(1)

    def __init__(self, model, clusters, redundancy_penalty=1.5):
        """

        :param model:  The model which has ben used to build this system
        :type model: :class:`macsypy.model.Model` object
        :param clusters: The list of cluster that form this system
        :type clusters: list of :class:`macsypy.cluster.Cluster` objects
        """
        super().__init__(model, clusters)
        self.id = f"{self.replicon_name}_{model.name}_{next(self._id)}"
        self.redundancy_penalty = redundancy_penalty
        self._score = None


    @property
    def score(self):
        """
        :return: a score take in account
            * if a hit match for the gene or it is an exchangeable gene
            * if a hit is duplicated and already present in the system or the cluster
            * if a hit match for mandatory/accessory gene of the model
        :rtype: float
        """
        if self._score is not None:
            # The score of the system is called for each clique
            # So to avoid computation we cached it
            return self._score

        def clst_func(clsts):
            """

            :param clsts: list of clusters
            :return: return the functions which are represented in these clusters.
            :rtype: dict with func str as keys and list of :class:`macsypy.cluster.Cluster` as values.
            """
            func_in_clst = {}
            for clst in clsts:
                for m_hit in clst.hits:
                    func = m_hit.gene_ref.alternate_of().name
                    if func in func_in_clst:
                        func_in_clst[func].append(clst)
                    else:
                        func_in_clst[func] = [clst]
            return func_in_clst

        _log.debug(f"=================== score computation for system {self.id} ===================")
        # split clusters in 2
        # the clusters true loners and multi systems (out of regular cluster)
        # and the others: regular cluster
        regular_clsts = []
        loner_multi_syst_clsts = []
        for clst in self.clusters:
            if clst.loner or clst.multi_system:
                # clst.multi_system is True
                # only if the cluster is composed of only one MultiSystemHit
                # that mean the hit come from a loner or an other cluster
                # So we have to apply a weight
                loner_multi_syst_clsts.append(clst)
            else:
                regular_clsts.append(clst)

        # Compute score of regular clusters
        clst_scores = [clst.score for clst in regular_clsts]
        score = sum(clst_scores)
        _log.debug(f"regular clusters scores sum({clst_scores}) = {score}")
        for gene in self.model.mandatory_genes + self.model.accessory_genes:
            _log.debug("compute penalty redundancy")
            # it cannot be forbidden gene in System instance
            # the neutral genes do not play a role in score (only to build clusters)
            clst_having_hit = sum([1 for clst in regular_clsts if clst.fulfilled_function(gene)])
            _log.debug(f"nb of clusters which fulfill function of {gene.name} = {clst_having_hit}")
            if clst_having_hit:
                clst_penalty = (clst_having_hit - 1) * self.redundancy_penalty
                _log.debug(f"clst_penalty {- clst_penalty}")
                score -= clst_penalty

        # compute score of loners
        # and multi systems out of regular cluster
        # loners_multi_syst_functions = clst_func(loner_multi_syst_clsts)
        # transform keys of dict in set
        regular_functions = set(clst_func(regular_clsts))
        _log.debug("compute score of loner or multi systems")
        for clst in loner_multi_syst_clsts:
            loner_or_ms = clst.hits[0]  # len(clst) == 1
            funct = loner_or_ms.gene_ref.alternate_of().name
            if not funct in regular_functions:
                _log.debug(f"{funct} is not already in regular clusters {regular_functions}")
                # call the cluster score
                # because it's in this method that the out of cluster penalty is applied
                loner_score = clst.score
                _log.debug(f"score for {funct} = {loner_score}")
                score += loner_score
            else:
                # if the biological funct is already encoded by regular clusters
                # we do not increase the score
                _log.debug(f"{funct} is already in regular clusters")
                pass

        self._score = score
        _log.debug(f"score of system {self.id} = {score:.2f}")
        return score


    def occurrence(self):
        """
        sometimes several systems collocates so they form only one cluster
        so macsyfinder build only one system
        the occurrence is an indicator of how many systems are
        it's based on the number of occurrence of each mandatory genes
        The multi_system genes are not take in account.

        :return: a predict number of biologic systems
        """
        genes = {g.name: g for g in self.model.genes()}
        occ_per_gene = [len(hits) for gene_name, hits in self._mandatory_occ.items()
                        if not genes[gene_name].multi_system]
        # if a systems contains 5 gene with occ of 1 and 5 gene with 0 occ
        # the median is 0.5
        # round(0.5) = 0
        # so I fix a floor value at 1
        return max(1, round(statistics.median(occ_per_gene)))


    @property
    def hits(self):
        """
        :return: The list of all hits that compose this system
        :rtype: [:class:`macsypy.hit.ValidHits` , ... ]
        """
        hits = self._sort_hits([h for cluster in self.clusters for h in cluster.hits])
        return hits

    @property
    def loci_num(self):
        """
        :return: the number of the corresponding locus for each cluster
                 the cluster made of only one Loner are not considered as a loci
                 so these clusters have a negative locus_num
        :rtype: list of int
        """
        loci = []
        loci_num = 0
        loners_num = 0
        # we do not take loners in account
        for clst in self.clusters:
            if not clst.loner:
                loci_num += 1
                loci.append(loci_num)
            else:
                loners_num -= 1
                loci.append(loners_num)
        return loci


    @property
    def loci_nb(self):
        """
        :return: The number of loci of this system (loners are not considered)
        :rtype: int >= 0
        """
        loci_nb = len([1 for c in self.clusters if not c.loner])
        return loci_nb


    @property
    def multi_loci(self):
        """
        :return: True if the systems is encoded in multiple loci. False otherwise
        :rtype: bool
        """
        return self.loci_nb > 1


    def is_compatible(self, other):
        """
        :param other: the other systems to test compatibility
        :type other: :class:`macsypy.system.System` object
        :return: True if other system is compatible with this one. False otherwise.
                 Two systems are compatible if they do not share :class:`macsypy.hit.CoreHit`
                 except hit corresponding to a multi_system gene in the model.

                 .. note::
                    This method is used to compute the best combination of systems.
        """
        # self and other are 2 System they can be related to different Model
        # a CoreHit correspond to 1 gene in replicon
        # it can be only one CoreHit by gene
        # but several ModelGene can exist on the same gene
        # So to know if two systems share same genes we have to work on the CoreHit
        # which is hold by the ModelHit hit attribute
        if self.model == other.model:
            # The Multi_systems are allowed to be "shared" by several oocurences of the same Model
            # The loners not but msf cannot decide which occurence is the right, so all occurences
            # are proposed to the user with eventually a warning
            from macsypy.hit import Loner, MultiSystem, LonerMultiSystem
            other_hits = {mh.hit for mh in other.hits if not isinstance(mh, (Loner, MultiSystem, LonerMultiSystem))}
            my_hits = {mh.hit for mh in self.hits if not isinstance(mh, (Loner, MultiSystem, LonerMultiSystem))}
            return not (my_hits & other_hits)
        else:
            # When we conpare 2 systems from 2 diffrent models
            # Only multi_model gene are allowed in same combination
            # to be allowed the hit from self and other must be multi_model
            other_hits = {mh.hit: mh for mh in other.hits}
            my_hits = {mh.hit: mh for mh in self.hits}
            common_hits = set(my_hits.keys()) & set(other_hits.keys())

            for ch in common_hits:
                other_hit = other_hits[ch]
                my_hit = my_hits[ch]
                if not(other_hit.multi_model and my_hit.multi_model):
                    return False
            return True





    def get_loners(self):
        """
        :return: The True Loners (Loner which not colocalize with an other hit) belonging to the systems
        :rtype: set of :class:`macsypy.hit.Loner` object
        """
        # a model hit is a loner only if it's a true loner
        return {mh for mh in self.hits if mh.loner}


    def get_hits_encoding_multisystem(self):
        """

        :return: The hits codding for a gene taged as multi system
        :rtype: set of :class:`macsypy.hit.ModelHit` object
        """
        return {mh for mh in self.hits if mh.gene_ref.multi_system}


    def get_multisystems(self):
        """
        :return: The MultiSystem hit (comming from out system (other cluster or loner) and tag as multisystem)
        :rtype: set of :class:`macsypy.hit.MultiSystem` | :class:`macsypy.hit.LonerMultiSystem` object
        """
        return {mh for mh in self.hits if mh.multi_system}


class RejectedCandidate(AbstractClusterizedHits):
    """
    Handle a set of clusters which has been rejected during the :func:`macsypy.system.match`  step
    This clusters (can be one) does not fill the requirements or contains forbidden genes.
    """
    _supported_status = (GeneStatus.MANDATORY,
                         GeneStatus.ACCESSORY,
                         GeneStatus.NEUTRAL,
                         GeneStatus.FORBIDDEN)

    _id = itertools.count(1)

    def __init__(self, model, clusters, reasons):
        """
        :param model:
        :type model: :class:`macsypy.model.Model` object
        :param clusters: list of clusters. These Clusters should be created with
                         :class:`macsypy.cluster.Cluster` of :class:`macsypy.hit.ModelHit` objects
        :type clusters: list of :class:`macsypy.cluster.Cluster` objects
        :param reasons: the reason why these clusters have been rejected
        :type reasons: list of string
        """
        super().__init__(model, clusters)
        self.id = f"{self.replicon_name}_{model.name}_{next(self._id)}"  # for testing purpose only
        self._reasons = reasons if isinstance(reasons, list) else [reasons]


    def __str__(self):
        """

        :return: a string representation of this RejectedCandidates
        """
        s = ''
        for c in self.clusters:
            s += str(c)
            s += '\n'
        s += "This candidate has been rejected because:\n"
        for r in self.reasons:
            s += f"\t- {r}\n"
        return s


    @property
    def hits(self):
        """

        :return: The list of all hits that compose this system
        :rtype: [:class:`macsypy.hit.ModelHit` , ... ]
        """
        hits = self._sort_hits([h for cluster in self.clusters for h in cluster.hits])
        return hits


    @property
    def reasons(self):
        return self._reasons


class AbstractUnordered(AbstractSetOfHits):
    """
    Technical abstract class to factorize code share between
    LikelySystem and UnlikelySystem
    """

    _id = itertools.count(1)

    def __init__(self, model, mandatory_hits, accessory_hits, neutral_hits, forbidden_hits):
        self._mandatory_hits = mandatory_hits
        self._accessory_hits = accessory_hits
        self._neutral_hits = neutral_hits
        self._forbidden_hits = forbidden_hits
        self._replicon_name = self.allowed_hits[0].replicon_name
        self.id = f"{self.replicon_name}_{model.name}_{next(self._id)}"
        super().__init__(model)

    @property
    def hits(self):
        """
        :return: The list of all hits sorted by their position
        :rtype: list of :class:`macsypy.hit.ModelHit` objects
        """
        return self._sort_hits(self._mandatory_hits + self._accessory_hits + self._neutral_hits + self.forbidden_hits)

    @property
    def mandatory_hits(self):
        """
        :return: The list of mandatory hits
        :rtype: list of :class:`macsypy.hit.ModelHit` objects
        """
        return self._mandatory_hits[:]

    @property
    def accessory_hits(self):
        """
        :return: The list of accesory hits
        :rtype: list of :class:`macsypy.hit.ModelHit` objects
        """
        return self._accessory_hits[:]

    @property
    def neutral_hits(self):
        """
        :return: The list of neutral hits
        :rtype: list of :class:`macsypy.hit.ModelHit` objects
        """
        return self._neutral_hits[:]

    @property
    def forbidden_hits(self):
        """
        :return: The list of forbidden hits
        :rtype: list of :class:`macsypy.hit.ModelHit` objects
        """
        return self._forbidden_hits[:]

    @property
    def allowed_hits(self):
        """
        :return: The list of allowed (mandatory, accessory, neutral) hits
        :rtype: list of :class:`macsypy.hit.ModelHit` objects
        """

        return self._mandatory_hits + self._accessory_hits + self._neutral_hits


class LikelySystem(AbstractUnordered):
    """"
    Handle components that fill the quorum requirements defined in model.
    with no idea about genetic organization (gene cluster)
    so we cannot take in account forbidden genes

    .. note:
        do not forget that this class inherits from MetaSetOfHits
        so the getter to mandatory, accessory, neutral, forbidden is dynamically injected
        by the meta class base on  _supported_status
    """

    _supported_status = (GeneStatus.MANDATORY,
                         GeneStatus.ACCESSORY,
                         GeneStatus.NEUTRAL,
                         GeneStatus.FORBIDDEN)


    def __str__(self):
        """

        :return: a string representation of this LikelySystem
        """
        return ', '.join([f"({h.id}, {h.gene.name}, {h.position})" for h in self.hits])


class UnlikelySystem(AbstractUnordered):
    """
    Handle components that not fill the quorum requirements defined in model.
    """

    _supported_status = (GeneStatus.MANDATORY,
                         GeneStatus.ACCESSORY,
                         GeneStatus.NEUTRAL,
                         GeneStatus.FORBIDDEN)


    def __init__(self, model, mandatory_hits, accessory_hits, neutral_hits, forbidden_hits, reasons):
        """

        :param model:  The model which has ben used to build this system
        :type model: :class:`macsypy.model.Model` object
        :param mandatory_hits: The list of mandatory hits (encode for a gene tagged as mandatory)
        :type mandatory_hits: list of :class:`macsypy.hit.ModelHit` objects
        :param accessory_hits: The list of accessory hits (encode for a gene tagged as accessory)
        :type accessory_hits: list of :class:`macsypy.hit.ModelHit` objects
        :param neutral_hits: The list of neutral hits (encode for a gene tagged as neutral)
        :type neutral_hits: list of :class:`macsypy.hit.ModelHit` objects
        :param forbidden_hits: The list of hits that are forbidden
        :type forbidden_hits: list of :class:`macsypy.hit.ModelHit` objects
        :param reasons: the reasons why this set of hits has been rejected
        :type reasons: List of str
        """
        super().__init__(model, mandatory_hits, accessory_hits, neutral_hits, forbidden_hits)
        self._reasons = reasons if isinstance(reasons, list) else [reasons]


    def __str__(self):
        """

        :return: a string representation of this UnlikelySystem
        """
        s = ', '.join([f"({h.id}, {h.gene.name}, {h.position})" for h in self.hits])
        s += ': These hits does not probably constitute a system because:\n'
        s += '\n'.join(self.reasons)
        return s


    @property
    def reasons(self):
        """
        :return: The reasons why it probably not a system
        :rtype: list of string
        """
        return self._reasons
