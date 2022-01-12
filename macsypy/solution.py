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
import networkx as nx

from macsypy.system import RejectedClusters


def find_best_solutions(systems):
    """
    Among the systems choose the combination of systems which does not share :class:`macsypy.hit.CoreHit`
    and maximize the sum of systems scores

    :param systems: the systems to analyse
    :type systems: list of :class:`macsypy.system.System` object
    :return: the list of list of systems which represent one best solution and the it's score
    :rtype: tuple of 2 elements the best solution and it's score
            ([[:class:`macsypy.system.System`, ...], [:class:`macsypy.system.System`, ...]], float score)
            The inner list represent a best solution
    """
    def sort_cliques(clique):
        """
        sort cliques

         - first by the sum of hits of systems composing the solution, most hits in first
         - second by the number of systems, most system in first
         - third by the average of wholeness of the systems
         - and finally by hits position. This criteria is to produce predictable results
           between two runs and to be testable (functional_test gembase)

        :param clique: the solutions to sort
        :type clique: List of :class:`macsypy.system.System` objects
        :return: the clique ordered
        """
        l = []
        for solution in clique:
            hits_pos = {hit.position for syst in solution for hit in syst.hits}
            hits_pos = sorted(list(hits_pos))
            l.append((sorted(solution, key=lambda sys: sys.id), hits_pos))

        sorted_cliques = sorted(l, key=lambda item: (sum([len(sys.hits) for sys in item[0]]),
                                                     len(item[0]),
                                                     item[1],
                                                     sum([sys.wholeness for sys in item[0]]) / len(item[0]),
                                                     '_'.join([sys.id for sys in item[0]])
                                                     ),
                                reverse=True)
        sorted_cliques = [item[0] for item in sorted_cliques]
        return sorted_cliques

    G = nx.Graph()
    # add nodes (vertices)
    G.add_nodes_from(systems)
    # let's create an edges between compatible nodes
    for sys_i, sys_j in itertools.combinations(systems, 2):
        if sys_i.is_compatible(sys_j):
            G.add_edge(sys_i, sys_j)

    cliques = nx.algorithms.clique.find_cliques(G)
    max_score = None
    max_cliques = []
    for c in cliques:
        current_score = sum([s.score for s in c])
        if max_score is None or (current_score > max_score):
            max_score = current_score
            max_cliques = [c]
        elif current_score == max_score:
            max_cliques.append(c)
    # sort the solutions (cliques)
    solutions = sort_cliques(max_cliques)
    return solutions, max_score


def combine_clusters(clusters, true_loners, multi_loci=False):
    """
    generate the combinations of clusters, with loners and multi systems

    :param clusters: the clusters to combines
    :type clusters: list of :class:`macsypy.cluster.Cluster` object
    :param true_loners: the multi-systems hits
    :type true_loners: dict the name of the function code by hit gene_ref.alternate_of as key
                              and 1 :class:`macsypy.cluster.Cluster` with the best a
                              :class:`macsypy.hit.Loner` or
                              :class:`macsypy.hit.LonerMultiSsystem` hit  as value
    :param bool multi_loci: True if the model is multi_loci false otherwise
    :return: all available combination of clusters
    :rtype: List of combination. a combination is a tuple of :class:`macsypy.cluster.Cluster` objects
    """
    if not clusters:
        cluster_combinations = []
    elif multi_loci:
        cluster_combinations = [itertools.combinations(clusters, i) for i in range(1, len(clusters) + 1)]
        cluster_combinations = list(itertools.chain(*cluster_combinations))
    else:
        cluster_combinations = [(clst,) for clst in clusters]

    # add loners
    loners_combinations = true_loners.items()
    loners_combinations = [itertools.combinations(loners_combinations, i) for i in
                           range(1, len(loners_combinations) + 1)]
    loners_combinations = itertools.chain(*loners_combinations)
    combination_w_loners = []
    for loner_comb in loners_combinations:
        loner_functions = [item[0] for item in loner_comb]
        loners = [item[1] for item in loner_comb]
        if cluster_combinations:
            for one_combination in cluster_combinations:
                to_add = True
                for clstr in one_combination:
                    if clstr.fulfilled_function(*loner_functions):
                        to_add = False
                        break
                if to_add:
                    combination_w_loners.append(tuple(list(one_combination) + loners))
        # we always add the loner
        # in case definition may contain only loners
        # or min_gene_required = 1 with one Loner
        combination_w_loners.append(tuple(loners))

    cluster_combinations += combination_w_loners

    return cluster_combinations


def combine_multisystems(rejected_clusters, multi_systems):
    """

    :param rejected_clusters:
    :param multi_systems: sequence of :class:`macsypy.cluster.Cluster`
                          each cluster must be composed of only one :class:`macsypy.hit.MultiSystem` object
    :return: list of cluster combination with teh multisystem
    :rtype: [(:class:`macsypy.cluster.Cluster` cluster1, cluster2, ...),
             (:class:`macsypy.cluster.Cluster` cluster3, cluster4, ...)]
    """
    # instead as clusters with loners
    # we have not to take in account rejected_cluster in new_comb
    # as they already have been challenged and rejected
    if isinstance(rejected_clusters, RejectedClusters):
        rejected_clusters = [rejected_clusters]
    new_comb = []
    ms_combinations = [itertools.combinations(multi_systems, i) for i in range(1, len(multi_systems) + 1)]
    ms_combinations = list(itertools.chain(*ms_combinations))  # tuple of MultiSystem
    for rej_clust in rejected_clusters:
        for one_ms_combination in ms_combinations:
            combination_hits = {h for clst in one_ms_combination for h in clst.hits}
            functions = [h.gene_ref.alternate_of().name for h in combination_hits]
            if not rej_clust.fulfilled_function(*functions):
                new_comb.append(tuple(rej_clust.clusters + list(one_ms_combination)))
    return new_comb

