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

import itertools
import logging

import macsypy.gene

from .error import MacsypyError
from .gene import GeneStatus
from .hit import ModelHit, Loner, get_best_hit_4_func

_log = logging.getLogger(__name__)


def build_clusters(hits, rep_info, model, hit_weights):
    """
    From a list of filtered hits, and replicon information (topology, length),
    build all lists of hits that satisfied the constraints:

        * max_gene_inter_space
        * loner
        * multi_system

    If Yes create a cluster
    A cluster contains at least two hits separated by less or equal than max_gene_inter_space
    Except for loner genes which are allowed to be alone in a cluster

    :param hits: list of filtered hits
    :type hits: list of :class:`macsypy.hit.CoreHit` objects
    :param rep_info: the replicon to analyse
    :type rep_info: :class:`macsypy.Indexes.RepliconInfo` object
    :param model: the model to study
    :type model: :class:`macsypy.model.Model` object
    :return: list of clusters
    :rtype: List of :class:`Cluster` objects
    """
    def collocates(h1, h2):
        """
        compute the distance (in number of gene between) between 2 hits

        :param :class:`macsypy.hit.ModelHit` h1: the first hit to compute inter hit distance
        :param :class:`macsypy.hit.ModelHit` h2: the second hit to compute inter hit distance
        :return: True if the 2 hits spaced by lesser or equal genes than inter_gene_max_space.
                 Managed circularity.
        """
        # compute the number of genes between h1 and h2
        dist = h2.get_position() - h1.get_position() - 1
        g1 = h1.gene_ref
        g2 = h2.gene_ref
        inter_gene_max_space = max(g1.inter_gene_max_space, g2.inter_gene_max_space)
        if 0 <= dist <= inter_gene_max_space:
            return True
        elif dist <= 0 and rep_info.topology == 'circular':
            # h1 and h2 overlap the ori
            dist = rep_info.max - h1.get_position() + h2.get_position() - rep_info.min
            return dist <= inter_gene_max_space
        return False

    ##################
    # build clusters #
    ##################
    clusters = []
    cluster_scaffold = []
    # sort hits by increasing position and then descending score
    hits.sort(key=lambda h: (h.position, - h.score))
    # remove duplicates hits (several hits for the same sequence),
    # keep the first one, this with the best score
    # position == sequence rank in replicon
    hits = [next(group) for pos, group in itertools.groupby(hits, lambda h: h.position)]
    if hits:
        hit = hits[0]
        gene = gene = model.get_gene(hit.gene.name)
        cluster_scaffold.append(ModelHit(hit, gene_ref=gene, gene_status=gene.status))
        previous_hit = cluster_scaffold[0]

        for hit in hits[1:]:
            # hit is a CoreHit so hit.gene is a CoreGene
            gene = model.get_gene(hit.gene.name)
            m_hit = ModelHit(hit, gene_ref=gene, gene_status=gene.status)
            if collocates(previous_hit, m_hit):
                cluster_scaffold.append(m_hit)
            else:
                if len(cluster_scaffold) > 1 :
                    # close the current scaffold if it contains at least 2 hits
                    cluster = Cluster(cluster_scaffold, model, hit_weights)
                    clusters.append(cluster)
                elif model.min_genes_required == 1:
                    # close the current scaffold if it contains 1 hit
                    # but it's allowed by the model
                    cluster = Cluster(cluster_scaffold, model, hit_weights)
                    clusters.append(cluster)
                elif model.get_gene(cluster_scaffold[0].gene.name).loner:
                    # close the current scaffold it contains 1 hit
                    # to handle circularity if it's the last cluster
                    cluster = Cluster(cluster_scaffold, model, hit_weights)
                    clusters.append(cluster)
                    # the hit transformation in loner is perform at the end when
                    # circularity and merging is done

                # open new scaffold
                cluster_scaffold = [m_hit]
            previous_hit = m_hit

        # close the last current cluster
        len_scaffold = len(cluster_scaffold)
        if len_scaffold > 1:
            new_cluster = Cluster(cluster_scaffold, model, hit_weights)
            clusters.append(new_cluster)
        elif len_scaffold == 1:
            # handle circularity
            # if there are clusters
            # may be the hit collocate with the first hit of the first cluster
            if clusters and collocates(cluster_scaffold[0], clusters[0].hits[0]):
                new_cluster = Cluster(cluster_scaffold, model, hit_weights)
                clusters[0].merge(new_cluster, before=True)
            elif cluster_scaffold[0].gene_ref.loner:
                # the hit does not collocate but it's a loner
                # handle clusters containing only one loner
                new_cluster = Cluster(cluster_scaffold, model, hit_weights)
                clusters.append(new_cluster)
            elif model.min_genes_required == 1:
                # the hit does not collocate but the model required only one gene
                # handle clusters containing only one gene
                new_cluster = Cluster(cluster_scaffold, model, hit_weights)
                clusters.append(new_cluster)

        # handle circularity
        if len(clusters) > 1:
            if collocates(clusters[-1].hits[-1], clusters[0].hits[0]):
                clusters[0].merge(clusters[-1], before=True)
                clusters = clusters[:-1]

        ###################
        # get True Loners #
        ###################
        # true_loner is a hit which encode for a gene tagged as loner
        # and which does NOT clusterize with some other hits of interest
        true_clusters = []
        true_loners = {}

        for clstr in clusters:
            if len(clstr) > 1:
                true_clusters.append(clstr)
            elif clstr.loner:
                loner = clstr[0]
                func_name = loner.gene_ref.alternate_of().name
                if func_name in true_loners:
                    true_loners[func_name].append(loner)
                else:
                    true_loners[func_name] = [loner]
            else:
                # it's a cluster of 1 hit
                # but it's NOT a loner
                # min_genes_required == 1
                true_clusters.append(clstr)

        for func_name in true_loners:
            # transform ModelHit in Loner
            loners =  true_loners[func_name]
            true_loners[func_name] = []
            for i in range(len(loners)):
                counterpart = loners[:]
                hit = counterpart.pop(i)
                true_loners[func_name].append(Loner(hit, counterpart=counterpart))
            # replace List of Loners by the best Loner
            best_loner = get_best_hit_4_func(func_name, true_loners[func_name], key='score')
            true_loners[func_name] = best_loner

        true_loners = {func_name: Cluster([loner], model, hit_weights) for func_name, loner in true_loners.items()}

        ##########################
        # get multi_systems hits #
        ##########################
        multi_system_hits = {}
        for clstr in true_clusters:
            for hit in clstr.hits:
                if not hit.gene_ref.multi_system:
                    continue
                func_name = hit.gene_ref.alternate_of().name
                if func_name in  multi_system_hits:
                    multi_system_hits[func_name].append(hit)
                else:
                    multi_system_hits[func_name] = [hit]
        for func_name in multi_system_hits:
            best_loner = get_best_hit_4_func(func_name, multi_system_hits[func_name], key='score')
            # replace list of multisystem hit by the best one
            multi_system_hits[func_name] = best_loner

        multi_system_hits = {func_name: Cluster([ms], model, hit_weights) for func_name, ms in multi_system_hits.items()}
    else: #there is not hits
        true_clusters = []
        true_loners = {}
        multi_system_hits = {}
    return true_clusters, true_loners, multi_system_hits



def get_loners(hits, model, hit_weights):
    """
    Create a list of Clusters each cluster is build with one hit matching a loner

    :param hits: The list of hits to filter
    :param model: the model which will used to build the clusters
    :type model: :class:`macsypy.model.Model` object
    :return: The list of cluster which each element is build at least with one loner
    :rtype: [Cluster, ...]
    """
    gene_loners = {g.name for g in model.genes() if g.loner}
    loners = [hit for hit in hits if hit.gene.name in gene_loners]
    loners = [Cluster([hit], model, hit_weights) for hit in loners]
    return loners


def filter_loners(cluster, loners):
    """
    filter loners to remove those which are already in the cluster

    :param cluster: The cluster
    :type cluster: :class:`macsypy.cluster.Cluster` object
    :param loners: the clusters constituted by one loner to filter
    :type loners: list of cluster [Cluster, ...]
    :return: list of loners which are not already in the cluster
    :rtype: [Clsuter, ...]
    """
    cluster_hit_name = {hit.gene.name for hit in cluster.hits}
    filtered_loners = []
    for loner in loners:
        if loner.hits[0].gene.name not in cluster_hit_name:
            filtered_loners.append(loner)
    return filtered_loners


class Cluster:
    """
    Handle hits relative to a model which collocates
    """


    def __init__(self, hits, model, hit_wheights):
        """

        :param hits: the hits constituting this cluster
        :type hits: [ :class:`macsypy.hit.CoreHit` | :class:`macsypy.hit.ModelHit`, ... ]
        :param model: the model associated to this cluster
        :type model: :class:`macsypy.model.Model`
        """
        self.hits = hits
        self.model = model
        self._check_replicon_consistency()
        self._score = None
        self._genes_roles = None
        self._hit_weights = hit_wheights


    def __len__(self):
        return len(self.hits)

    def __getitem__(self, item):
        return self.hits[item]


    @property
    def loner(self):
        """
        :return: True if this cluster is made of only one hit representing a loner gene
                 False otherwise:
                 - contains several hits
                 - contains one hit but gene is not tag as loner (max_gene_required = 1)
        """
        # need this method in build_cluster before to transform ModelHit in Loner
        # so cannot rely on Loner type
        return len(self) == 1 and self.hits[0].gene_ref.loner



    def _check_replicon_consistency(self):
        """
        :raise: MacsypyError if all hits of a cluster are NOT related to the same replicon
        """
        rep_name = self.hits[0].replicon_name
        if not all([h.replicon_name == rep_name for h in self.hits]):
            msg = "Cannot build a cluster from hits coming from different replicons"
            _log.error(msg)
            raise MacsypyError(msg)


    def __contains__(self, v_hit):
        """
        :param v_hit: The hit to test
        :type v_hit: :class:`macsypy.hit.ModelHit` object
        :return: True if the hit is in the cluster hits, False otherwise
        """
        return v_hit in self.hits


    def fulfilled_function(self, gene):
        """

        :param gene: The gene which must be tested.
        :type gene: :class:`macsypy.gene.ModelGene` object
        :return: True if the cluster contains one hit which fulfill the function corresponding to the gene
                 (the gene hitself or an exchageable)
        """
        if self._genes_roles is None:
            self._genes_roles = {h.gene_ref.alternate_of().name for h in self.hits}
        if isinstance(gene, macsypy.gene.ModelGene):
            function = gene.name
        else:
            # gene is a string
            function = gene
        return function in self._genes_roles


    def merge(self, cluster, before=False):
        """
        merge the cluster in this one. (do it in place)

        :param cluster:
        :type cluster: :class:`macsypy.cluster.Cluster` object
        :param bool before: If False the hits of the cluster will be add at the end of this one,
                            Otherwise the cluster hits will be inserted before the hits of this one.
        :return: None
        :raise MacsypyError: if the two clusters have not the same model
        """
        if cluster.model != self.model:
            raise MacsypyError("Try to merge Clusters from different model")
        else:
            if before:
                self.hits = cluster.hits + self.hits
            else:
                self.hits.extend(cluster.hits)

    @property
    def replicon_name(self):
        return self.hits[0].replicon_name


    @property
    def score(self):
        if self._score is not None:
            return self._score
        else:
            seen_hits = {}
            score = 0
            _log.debug(f"===================== compute score for cluster =====================")
            for v_hit in self.hits:
                _log.debug(f"-------------- test hit {v_hit.gene.name} --------------")

                # attribute a score for this hit
                # according to status of the gene_ref in the model: mandatory/accessory
                if v_hit.status == GeneStatus.MANDATORY:
                    hit_score = self._hit_weights.mandatory
                elif v_hit.status == GeneStatus.ACCESSORY:
                    hit_score = self._hit_weights.accessory
                elif v_hit.status == GeneStatus.NEUTRAL:
                    hit_score = self._hit_weights.neutral
                else:
                    raise MacsypyError(f"a Cluster contains hit {v_hit.gene.name} {v_hit.position} which is neither mandatory nor accessory: {v_hit.status}")
                _log.debug(f"{v_hit.id} is {v_hit.status} hit score = {hit_score}")

                # weighted the hit score according to the hit match the gene or
                # is an exchangeable
                if v_hit.gene_ref.is_exchangeable:
                    hit_score *= self._hit_weights.exchangeable
                    _log.debug(f"{v_hit.id} is exchangeable hit score = {hit_score}")
                else:
                    hit_score *= self._hit_weights.itself

                if self.loner and v_hit.multi_system:
                    hit_score *= self._hit_weights.loner_multi_system
                    _log.debug(f"{v_hit.id} is loner and multi_system hit score = {hit_score}")

                # funct is the name of the gene if it code for itself
                # or the name of the reference gene if it's an exchangeable
                funct = v_hit.gene_ref.alternate_of().name
                if funct in seen_hits:
                    # count only one occurrence of each function per cluster
                    # the score use is the max of hit score for this function
                    if hit_score > seen_hits[funct]:
                        seen_hits[funct] = hit_score
                        _log.debug(f"{v_hit.id} code for {funct} update hit_score to {hit_score}")
                    else:
                        _log.debug(f"{v_hit.id} code for {funct} which is already take in count in cluster")
                else:
                    _log.debug(f"{v_hit.id} {v_hit.gene_ref.name} is not already in cluster")
                    seen_hits[funct] = hit_score

            hits_scores = seen_hits.values()
            score = sum(hits_scores)
            _log.debug(f"cluster score = sum({list(hits_scores)}) = {score}")
        _log.debug(f"===============================================================")
        self._score = score
        return score


    def __str__(self):
        """

        :return: a string representation of this cluster
        """
        s = f"""Cluster:
- model = {self.model.name}
- replicon = {self.replicon_name}
- hits = {', '.join([f"({h.id}, {h.gene.name}, {h.position})" for h in self.hits])}"""
        return s
