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
import networkx as nx


def find_best_solutions(systems):
    """
    Among the systems choose the combination of systems which does not share :class:`macsypy.hit.Hit`
    and maximize the sum of systems scores

    :param systems: the systems to analyse
    :type systems: list of :class:`macsypy.system.System` object
    :return: the list of list of systems which represent one best solution and the it's score
    :rtype: tuple of 2 elements
            ([[:class:`macsypy.system.System`, ...], [:class:`macsypy.system.System`, ...]], float score)
            The inner list represent a best solution
    """
    G = nx.Graph()
    # add nodes (vertices)
    G.add_nodes_from(systems)
    # let's create an edges between compatible nodes
    for sys_i, sys_j in itertools.combinations(systems, 2):
        if sys_i.is_compatible(sys_j):
            G.add_edge(sys_i, sys_j)

    cliques = nx.algorithms.clique.find_cliques(G)
    max_score = 0
    max_cliques = []
    for c in cliques:
        current_score = sum([s.score for s in c])
        if current_score > max_score:
            max_score = current_score
            max_cliques = [c]
        elif current_score == max_score:
            max_cliques.append(c)
    return max_cliques, max_score