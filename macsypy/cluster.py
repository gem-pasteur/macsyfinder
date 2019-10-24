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

import itertools
import logging

from .error import MacsypyError
from .hit import hit_weight
from .gene import GeneStatus

_log = logging.getLogger(__name__)


def build_clusters(hits, rep_info, model):
    """
    From a list of filtered hits, and replicon information (topology, length),
    build all lists of hits that satisfied the constraints:

        * max_gene_inter_space
        * loner
        * multi_system

    If Yes create a cluster
    A cluster contains at least to hits separated by less or equal than max_gene_inter_space
    Except for loner genes which are allowed to be alone in a cluster

    :param hits: list of filtered hits
    :type hits: list of :class:`macsypy.report.Hit` objects
    :param rep_info: the replicon to analyse
    :type rep_info: :class:`macsypy.Indexes.RepliconInfo` object
    :param model:
    :type model: :class:`macsypy.model.Model` object
    :return: list of clusters
    :rtype: List of :class:`Cluster` objects
    """
    def collocates(h1, h2, model):
        # compute the number of genes between h1 and h2
        dist = h2.get_position() - h1.get_position() - 1
        g1 = model.get_gene(h1.gene.name)
        g2 = model.get_gene(h2.gene.name)
        inter_gene_max_space = max(g1.inter_gene_max_space, g2.inter_gene_max_space)
        if 0 <= dist <= inter_gene_max_space:
            return True
        elif dist <= 0 and rep_info.topology == 'circular':
            # h1 and h2 overlap the ori
            dist = rep_info.max - h1.get_position() + h2.get_position() - rep_info.min
            return dist <= inter_gene_max_space
        return False

    clusters = []
    cluster_scaffold = []
    # sort hits by increasing position and then descending score
    hits.sort(key=lambda h: (h.position, - h.score))
    # remove duplicates hits (several hits for the same sequence),
    # keep the first one, this with the best score
    # position == sequence rank in replicon
    hits = [next(group) for pos, group in itertools.groupby(hits, lambda h: h.position)]
    if hits:
        cluster_scaffold.append(hits[0])
        previous_hit = hits[0]
        for hit in hits[1:]:
            if collocates(previous_hit, hit, model):
                cluster_scaffold.append(hit)
            else:
                # by definition a loner gene can be alone in a cluster
                is_a_loner = model.get_gene(cluster_scaffold[0].gene.name).loner
                if len(cluster_scaffold) > 1 or is_a_loner:
                    cluster = Cluster(cluster_scaffold, model)
                    clusters.append(cluster)
                cluster_scaffold = [hit]
            previous_hit = hit

        # close the current cluster
        if len(cluster_scaffold) > 1:
            new_cluster = Cluster(cluster_scaffold, model)
            if clusters and collocates(new_cluster.hits[-1], clusters[0].hits[0], model):
                # handle circular replicon
                clusters[0].merge(new_cluster, before=True)
            else:
                clusters.append(new_cluster)
        elif clusters:
            if collocates(previous_hit, clusters[0].hits[0], model):
                clusters[0].merge(Cluster([previous_hit], model), before=True)
            elif model.get_gene(previous_hit.gene.name).loner:
                # this hit is far from other clusters
                # but by definition a loner gene can be alone in a cluster
                cluster = Cluster(cluster_scaffold, model)
                clusters.append(cluster)
    return clusters


def get_loners(hits, model):
    """
    Create a list of Clusters each cluster is build with one hit matching a loner

    :param hits: The list of hits to filter
    :param model: the model which will used to build the clusters
    :type model: :class:`macsypy.model.Model` object
    :return: The list of cluster which each element is build with one loner
    :rtype: [Cluster, ...]
    """
    gene_loners = {g.name for g in model.genes if g.loner}
    loners = [hit for hit in hits if hit.gene.name in gene_loners]
    loners = [Cluster([hit], model) for hit in loners]
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


    def __init__(self, hits, model):
        """

        :param hits: the valid hits constituting this cluster
        :type hits: [ :class:`macsypy.hit.ValidHit` , ... ]
        :param model: the model associated to this cluster
        :type model: :class:`macsypy.model.Model`
        """
        self.hits = hits
        self.model = model
        self._check_replicon_consistency()
        self._score = None
        self._genes_ref = None


    def __len__(self):
        return len(self.hits)


    def _check_replicon_consistency(self):
        rep_name = self.hits[0].replicon_name
        if not all([h.replicon_name == rep_name for h in self.hits]):
            msg = "Cannot build a cluster from hits coming from different replicons"
            _log.error(msg)
            raise MacsypyError(msg)


    def __contains__(self, v_hit):
        """
        :param v_hit: The hit to test
        :type v_hit: :class:`macsypy.hit.ValidHit` object
        :return: True if the hit is in the cluster hits, False otherwise
        """
        return v_hit in self.hits


    def fulfilled_function(self, gene):
        """

        :param gene: The gene which must be tested.
        :type gene: :class:`macsypy.gene.Gene` object
        :return: True if the cluster contains one hit which fulfill the function corresponding to the gene
        """
        if self._genes_ref is None:
            self._genes_ref = {h.gene_ref.name for h in self.hits}
        return gene.name in self._genes_ref


    def merge(self, cluster, before=False):
        """
        merge the cluster in this one. (do it in place)

        :param cluster:
        :type cluster: :class:`macsypy.cluster.Cluster` object
        :param bool before: If False the hits of the cluster will be add at the end of this one,
                            Otherwise the cluster hits will be inserted before the hits of this one.
        :return: None
        :raise MasypyError: if the two clusters have not the same model
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
            seen_hits = set()
            score = 0
            for v_hit in self.hits:
                if v_hit.gene_ref in seen_hits:
                    # count only one occurrence of each hit per cluster
                    continue

                # attribute a score for this hit
                # according to status of the gene_ref in the model: mandatory/accessory
                if v_hit.status == GeneStatus.MANDATORY:
                    hit_score = hit_weight.mandatory
                elif v_hit.status == GeneStatus.ACCESSORY:
                    hit_score = hit_weight.accessory
                elif v_hit.status == GeneStatus.NEUTRAL:
                    hit_score = hit_weight.neutral
                else:
                    raise MacsypyError("a Cluster contains hit which is neither mandatory nor accessory")
                # weighted the hit score according to the hit match the gene or
                # is an analog/homolog
                if v_hit.gene_ref.is_homolog:
                    hit_score *= hit_weight.homolog
                elif v_hit.gene_ref.is_analog:
                    hit_score *= hit_weight.analog
                else:
                    hit_score *= hit_weight.hitself
                score += hit_score
                seen_hits.add(v_hit.gene_ref)
        self._score = score
        return score


    def __str__(self):
        s = """Cluster:
    - model: {}
    - hits: {}""".format(self.model.name, ', '.join(["({}, {}, {})".format(h.id,
                                                                           h.gene.name,
                                                                           h.position) for h in self.hits]))
        return s


class RejectedClusters:
    """
    Handle a set of clusters which has been rejected during the :func:`macsypy.system.match`  step
    This clusters (can be one) does not fill the requirements or contains forbidden genes.
    """

    def __init__(self, model, clusters, reason):
        """
        :param model:
        :type model: :class:`macsypy.model.Model` object
        :param clusters:
        :type clusters: list of :class:`macsypy.cluster.Cluster` objects
        :param str reason: the reason why these clusters have been rejected
        """
        self.model = model
        if isinstance(clusters, Cluster):
            self.clusters = [clusters]
        else:
            self.clusters = clusters
        self.reason = reason

    def __str__(self):
        s = ''
        for c in self.clusters:
            s += str(c)
            s += '\n'
        s += f'These clusters has been rejected because:\n{self.reason}'
        return s
