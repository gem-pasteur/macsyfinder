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

import abc
import itertools
import statistics
from itertools import chain
from operator import attrgetter
import logging
_log = logging.getLogger(__name__)

from.gene import GeneStatus
from .cluster import Cluster
from .hit import ValidHit


# la liste des clusters a fournir est a generer avant match
# si len(clusters) = 1 single_loci
# si len(clusters) > 1 multi_loci
# il faut genegerer la liste de toutes les combinaisons
# et appeler cette fonction pour chaqu'une entre elles
# from itertools import combinations

# combinations('ABCD', 1) => inutile mais generique => single_loucs
# combinations('ABCD', 2) => multi_locus a ne faire que si model.multi_locus= True
# combinations('ABCD', 3)
# combinations('ABCD', len("ABCD")) => inutile mais generique => recheche parmis tous les clusters


def match(clusters, model):
    """
    Check a set of clusters fill model constraints.
    If yes create a :class:`macsypy.system.System` otherwise create
    a :class:`macsypy.cluster.RejectedClusters`.

    :param clusters: The list of cluster to check if fit the model
    :type clusters: list of :class:`macsypy.cluster.Cluster` objects
    :param model:  The model to consider
    :type model: :class:`macsypy.model.Model` object
    :return: either a System or a RejectedClusters
    :rtype: :class:`macsypy.system.System` or :class:`macsypy.cluster.RejectedClusters` object
    """
    def create_exchangeable_map(genes):
        """
        create a map between an exchangeable (formly homolog or analog) gene name and it's gene reference

        :param genes: The genes to get the exchangeable genes
        :type genes: list of :class:`macsypy.gene.ModelGene` objects
        :rtype: a dict with keys are the exchangeable gene_name and the value the reference gene name
        """
        map = {}
        for gene in genes:
            for ex_gene in gene.exchangeables:
                map[ex_gene.name] = gene
        return map

    # init my structures to count gene occurrences
    mandatory_counter = {g.name: 0 for g in model.mandatory_genes}
    exchangeable_mandatory = create_exchangeable_map(model.mandatory_genes)

    accessory_counter = {g.name: 0 for g in model.accessory_genes}
    exchangeable_accessory = create_exchangeable_map(model.accessory_genes)

    forbidden_counter = {g.name: 0 for g in model.forbidden_genes}
    exchangeable_forbidden = create_exchangeable_map(model.forbidden_genes)

    neutral_counter = {g.name: 0 for g in model.neutral_genes}
    exchangeable_neutral = create_exchangeable_map(model.neutral_genes)

    # count the hits
    # and track for each hit for which gene it counts for
    valid_clusters = []
    forbidden_hits = []
    for cluster in clusters:
        valid_hits = []
        for hit in cluster.hits:
            gene_name = hit.gene.name
            # the ValidHit need to be linked to the
            # gene of the model
            gene = model.get_gene(gene_name)
            if gene_name in mandatory_counter:
                mandatory_counter[hit.gene.name] += 1
                valid_hits.append(ValidHit(hit, gene, GeneStatus.MANDATORY))
            elif gene_name in exchangeable_mandatory:
                gene_ref = exchangeable_mandatory[gene_name]
                mandatory_counter[gene_ref.name] += 1
                valid_hits.append(ValidHit(hit, gene, GeneStatus.MANDATORY))
            elif gene_name in accessory_counter:
                accessory_counter[gene_name] += 1
                valid_hits.append(ValidHit(hit, gene, GeneStatus.ACCESSORY))
            elif gene_name in exchangeable_accessory:
                gene_ref = exchangeable_accessory[gene_name]
                accessory_counter[gene_ref.name] += 1
                valid_hits.append(ValidHit(hit, gene, GeneStatus.ACCESSORY))
            elif gene_name in neutral_counter:
                neutral_counter[gene_name] += 1
                valid_hits.append(ValidHit(hit, gene, GeneStatus.NEUTRAL))
            elif gene_name in exchangeable_neutral:
                gene_ref = exchangeable_neutral[gene_name]
                neutral_counter[gene_ref.name] += 1
                valid_hits.append(ValidHit(hit, gene, GeneStatus.NEUTRAL))
            elif gene_name in forbidden_counter:
                forbidden_counter[gene_name] += 1
                # valid_hits.append(ValidHit(hit, hit.gene, GeneStatus.FORBIDDEN))
                forbidden_hits.append(hit)
            elif gene_name in exchangeable_forbidden:
                gene_ref = exchangeable_forbidden[gene_name]
                forbidden_counter[gene_ref.name] += 1
                # valid_hits.append(ValidHit(hit, hit.gene.ref, GeneStatus.FORBIDDEN))
                forbidden_hits.append(hit)
        if valid_hits:
            valid_clusters.append(Cluster(valid_hits, model))

    # the count is finished
    # check if the quorum is reached
    # count how many different genes are represented in the clusters
    # the neutral genes belong to the cluster
    # but they do not count for the quorum
    mandatory_genes = [g for g, occ in mandatory_counter.items() if occ > 0]
    accessory_genes = [g for g, occ in accessory_counter.items() if occ > 0]
    neutral_genes = [g for g, occ in neutral_counter.items() if occ > 0]
    forbidden_genes = [g for g, occ in forbidden_counter.items() if occ > 0]
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
    if len(mandatory_genes) < model.min_mandatory_genes_required:
        is_a_system = False
        msg = f'The quorum of mandatory genes required ({model.min_mandatory_genes_required}) is not reached: ' \
              f'{len(mandatory_genes)}'
        reasons.append(msg)
        _log.debug(msg)
    if len(accessory_genes) + len(mandatory_genes) < model.min_genes_required:
        is_a_system = False
        msg = f'The quorum of genes required ({model.min_genes_required}) is not reached:' \
              f' {len(accessory_genes) + len(mandatory_genes)}'
        reasons.append(msg)
        _log.debug(msg)

    if is_a_system:
        res = System(model, valid_clusters)
        _log.debug("is a system")
    else:
        reason = '\n'.join(reasons)
        res = RejectedClusters(model, clusters, reason)
    _log.debug("#" * 50)
    return res


def unordered_match(hits, model):
    def create_exchangeable_map(genes):
        """
        create a map between an exchangeable (formly homolog or analog) gene name and it's gene reference

        :param genes: The genes to get the exchangeable genes
        :type genes: list of :class:`macsypy.gene.ModelGene` objects
        :rtype: a dict with keys are the exchangeable gene_name and the value the reference gene name
        """
        map = {}
        for gene in genes:
            for ex_gene in gene.exchangeables:
                map[ex_gene.name] = gene
        return map

    # init my structures to count gene occurrences
    mandatory_counter = {g.name: 0 for g in model.mandatory_genes}
    exchangeable_mandatory = create_exchangeable_map(model.mandatory_genes)

    accessory_counter = {g.name: 0 for g in model.accessory_genes}
    exchangeable_accessory = create_exchangeable_map(model.accessory_genes)

    forbidden_counter = {g.name: 0 for g in model.forbidden_genes}
    exchangeable_forbidden = create_exchangeable_map(model.forbidden_genes)

    neutral_counter = {g.name: 0 for g in model.neutral_genes}
    exchangeable_neutral = create_exchangeable_map(model.neutral_genes)

    # count the hits
    # and track for each hit for which gene it counts for
    forbidden_hits = []
    valid_hits = []
    for hit in hits:
        gene_name = hit.gene.name
        # the ValidHit need to be linked to the
        # gene of the model
        gene = model.get_gene(gene_name)
        if gene_name in mandatory_counter:
            mandatory_counter[hit.gene.name] += 1
            valid_hits.append(ValidHit(hit, gene, GeneStatus.MANDATORY))
        elif gene_name in exchangeable_mandatory:
            gene_ref = exchangeable_mandatory[gene_name]
            mandatory_counter[gene_ref.name] += 1
            valid_hits.append(ValidHit(hit, gene, GeneStatus.MANDATORY))
        elif gene_name in accessory_counter:
            accessory_counter[gene_name] += 1
            valid_hits.append(ValidHit(hit, gene, GeneStatus.ACCESSORY))
        elif gene_name in exchangeable_accessory:
            gene_ref = exchangeable_accessory[gene_name]
            accessory_counter[gene_ref.name] += 1
            valid_hits.append(ValidHit(hit, gene, GeneStatus.ACCESSORY))
        elif gene_name in neutral_counter:
            neutral_counter[gene_name] += 1
            valid_hits.append(ValidHit(hit, gene, GeneStatus.NEUTRAL))
        elif gene_name in exchangeable_neutral:
            gene_ref = exchangeable_neutral[gene_name]
            neutral_counter[gene_ref.name] += 1
            valid_hits.append(ValidHit(hit, gene, GeneStatus.NEUTRAL))
        elif gene_name in forbidden_counter:
            forbidden_counter[gene_name] += 1
            valid_hits.append(ValidHit(hit, hit.gene, GeneStatus.FORBIDDEN))
            # forbidden_hits.append(hit)
        elif gene_name in exchangeable_forbidden:
            gene_ref = exchangeable_forbidden[gene_name]
            forbidden_counter[gene_ref.name] += 1
            valid_hits.append(ValidHit(hit, hit.gene.ref, GeneStatus.FORBIDDEN))
            # forbidden_hits.append(hit)

    # the count is finished
    # check if the quorum is reached
    # count how many different genes are represented in the clusters
    # the neutral genes belong to the cluster
    # but they do not count for the quorum
    mandatory_genes = [g for g, occ in mandatory_counter.items() if occ > 0]
    accessory_genes = [g for g, occ in accessory_counter.items() if occ > 0]
    neutral_genes = [g for g, occ in neutral_counter.items() if occ > 0]
    forbidden_genes = [g for g, occ in forbidden_counter.items() if occ > 0]
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
    if len(mandatory_genes) < model.min_mandatory_genes_required:
        is_a_potential_system = False
        msg = f'The quorum of mandatory genes required ({model.min_mandatory_genes_required}) is not reached: ' \
              f'{len(mandatory_genes)}'
        reasons.append(msg)
        _log.debug(msg)
    if len(accessory_genes) + len(mandatory_genes) < model.min_genes_required:
        is_a_potential_system = False
        msg = f'The quorum of genes required ({model.min_genes_required}) is not reached:' \
              f' {len(accessory_genes) + len(mandatory_genes)}'
        reasons.append(msg)
        _log.debug(msg)

    if is_a_potential_system:
        res = PotentialSystem(model, valid_hits)
        _log.debug("There is a genetic potential for a system")
    else:
        reason = '\n'.join(reasons)
        res = NotPotentialSystem(model, valid_hits, reason)
    _log.debug("#" * 50)
    return res


class HitSystemTracker(dict):
    """
    track in which system is implied each hit
    """

    def __init__(self, systems):
        super(HitSystemTracker, self).__init__()
        for system in systems:
            v_hits = system.hits
            for v_hit in v_hits:
                hit = v_hit.hit
                if hit not in self:
                    self[hit] = set()
                self[hit].add(system)


class ClusterSystemTracker(dict):
    """
    track in which system is implied each cluster
    """
    def __init__(self, systems):
        super(ClusterSystemTracker, self).__init__()
        for system in systems:
            clusters = system.clusters
            for clst in clusters:
                if clst not in self:
                    self[clst] = set()
                self[clst].add(system)


class MetaSetOfHits(abc.ABCMeta):

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
        print()
        print("##################### MetaSystem __call__ #################")
        print("### new_system_inst", new_system_inst)
        print("### new_system_inst._supported_status", new_system_inst._supported_status)
        for status in [str(s) for s in new_system_inst._supported_status]:
            # set the private attribute in the Model instance
            setattr(new_system_inst, f"_{status}_occ", {})
            # set the public property in the Model class
            setattr(cls, f"{status}_occ", property(MetaSetOfHits.getter_maker(status)))
        new_system_inst.count()
        return new_system_inst


class AbstractSetOfHits(metaclass=MetaSetOfHits):

    _id = itertools.count(1)


    def __init__(self, model, replicon_name):
        self.id = f"{replicon_name}_{model.name}_{next(self._id)}"
        self.model = model


    @property
    @abc.abstractmethod
    def hits(self):
        pass


    def count(self):
        """
        fill 3 structures one for mandatory, accessory and neutral
        each structure count how many hit for each gene
        :return: None
        """
        for status in [str(s) for s in self._supported_status]:
            setattr(self,
                    f"_{status}_occ",
                    {g.name: [] for g in getattr(self.model, f"{status}_genes")}
                    )

        # all the hits are ValidHit
        for hit in self.hits:
            name = hit.gene_ref.alternate_of().name
            status = str(hit.status)  # transform gene status in lower string
            try:
                getattr(self, f"_{status}_occ")[name].append(hit)
            except AttributeError:
                pass



class System(AbstractSetOfHits):
    """
    Modelize as system. a system is an occurrence of a given model on a replicon.
    """

    _supported_status = _supported_status = (GeneStatus.MANDATORY,
                         GeneStatus.ACCESSORY,
                         GeneStatus.NEUTRAL)

    def __init__(self, model, clusters):
        """

        :param model:  The model which has ben used to build this system
        :type model: :class:`macsypy.model.Model` object
        :param clusters: The list of cluster that form this system
        :type clusters: list of :class:`macsypy.cluster.Cluster` objects
        """
        self._replicon_name = clusters[0].replicon_name
        self.clusters = clusters
        super().__init__(model, self._replicon_name)


    @property
    def replicon_name(self):
        """
        :return: The name of the replicon
        :rtype: str
        """
        return self._replicon_name


    @property
    def wholeness(self):
        """

        :return: a score indicating the genes ratio of the model which have at least one hit
                ('neutral' genes do not count)
        :rtype: float
        """
        # model completude
        # the neutral hit do not participate to the model completude
        score = sum([1 for hits in chain(self._mandatory_occ.values(), self._accessory_occ.values()) if hits]) / \
                   (len(self._mandatory_occ) + len(self._accessory_occ))
        return score


    @property
    def score(self):
        """
        :return: a score take in account
            * if a hit match for the gene or it is an exchangeable gene
            * if a hit is duplicated and already present in the system or the cluster
            * if a hit match for mandatory/accessory gene of the model
        :rtype: float
        """
        score = sum([clst.score for clst in self.clusters])
        for gene in self.model.mandatory_genes + self.model.accessory_genes:
            clst_having_hit = sum([1 for clst in self.clusters if clst.fulfilled_function(gene)])
            if clst_having_hit:
                clst_penalty = (clst_having_hit - 1) * 1.5
                score -= clst_penalty
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
        genes = {g.name: g for g in self.model.genes}
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
        hits = [h for cluster in self.clusters for h in cluster.hits]
        hits.sort(key=attrgetter('position'))
        return hits


    @property
    def loci(self):
        """
        :return: The number of loci of this system (loners are not considered)
        :rtype: int > 0
        """
        # we do not take loners in account
        loci = sum([1 for c in self.clusters if len(c) > 1])
        return loci


    @property
    def multi_loci(self):
        """
        :return: True if the systems is multi_loci. False otherwise
        :rtype: bool
        """
        return self.loci > 1

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


    def is_compatible(self, other):
        """
        :param other: the other systems to test compatibility
        :type other: :class:`macsypy.system.System` object
        :return: True if other system is compatible with this one. False otherwise.
                 Two systems are compatible if they do not share :class:`macsypy.hit.Hit`
                 except hit corresponding to a multi_system gene in the model.

                 .. note::
                    This method is used to compute the best combination of systems.
        """
        other_hits = {vh.hit for vh in other.hits if not vh.multi_system}
        my_hits = {vh.hit for vh in self.hits if not vh.multi_system}
        return not (my_hits & other_hits)



class RejectedClusters(AbstractSetOfHits):
    """
    Handle a set of clusters which has been rejected during the :func:`macsypy.system.match`  step
    This clusters (can be one) does not fill the requirements or contains forbidden genes.
    """
    _supported_status = (GeneStatus.MANDATORY,
                         GeneStatus.ACCESSORY,
                         GeneStatus.NEUTRAL,
                         GeneStatus.FORBIDDEN)

    def __init__(self, model, clusters, reason):
        """
        :param model:
        :type model: :class:`macsypy.model.Model` object
        :param clusters: list of clusters. These Clusters should be created with
                         :class:`macsypy.cluster.Cluster` of :class:`macsypy.hit.ValidHit` objects
        :type clusters: list of :class:`macsypy.cluster.Cluster` objects
        :param str reason: the reason why these clusters have been rejected
        """
        if isinstance(clusters, Cluster):
            self.clusters = [clusters]
        else:
            self.clusters = clusters
        self._replicon_name = clusters[0].replicon_name
        self.reason = reason
        super().__init__(model, self._replicon_name)


    def __str__(self):
        """

        :return: a string representation of this RejectedCluster
        """
        s = ''
        for c in self.clusters:
            s += str(c)
            s += '\n'
        s += f'These clusters has been rejected because:\n{self.reason}'
        return s


class PotentialSystem(AbstractSetOfHits):
    """"
    Handle component that fill the quorum requirements with no idea about
    genetic organization (gene cluster)
    so we cannot take in account forbidden genes

    .. note:
        do not forget that this class inherits from MetaSetOfHits
        so the accessory to mandatory, accessory, neutral, forbidden is dynamically injected
        by the meta class base on  _supported_status
    """

    _supported_status = (GeneStatus.MANDATORY,
                         GeneStatus.ACCESSORY,
                         GeneStatus.NEUTRAL,
                         GeneStatus.FORBIDDEN)


    def __init__(self, model, hits):
        """

        :param model:  The model which has ben used to build this system
        :type model: :class:`macsypy.model.Model` object
        :param hits: The list of hit that form this potentil system
        :type hits: list of :class:`macsypy.hit.ValidHit` objects
        """
        self._replicon_name = hits[0].replicon_name
        self._hits = hits
        super().__init__(model, self._replicon_name)

    @property
    def hits(self):
        return self._hits



class NotPotentialSystem(AbstractSetOfHits):

    _supported_status = (GeneStatus.MANDATORY,
                         GeneStatus.ACCESSORY,
                         GeneStatus.NEUTRAL,
                         GeneStatus.FORBIDDEN)


    def __init__(self, model, hits, reason):
        """

        :param model:  The model which has ben used to build this system
        :type model: :class:`macsypy.model.Model` object
        :param clusters: The list of cluster that form this system
        :type clusters: list of :class:`macsypy.cluster.Cluster` objects
        """
        self._replicon_name = hits[0].replicon_name
        self._hits = hits
        self._reason = reason
        super().__init__(model, self._replicon_name)


   def __str__(self):
        """

        :return: a string representation of this RejectedCluster
        """
        s = ''
        for c in self.clusters:
            s += str(c)
            s += '\n'
        s += f'These clusters has been rejected because:\n{self.reason}'
        return s

