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


class Solution:
    """
    manage set of systems to build a solution to the systems search
    the best solution should have the highest score and the components (hits)
    should have no overlaps (hit belonging to several system)
    """

    def __init__(self, systems):
        """
        :param systems: the systems which composed the solution
        :type systems: list of :class:`macsypy.system.System` object
        """
        self._score = 0
        self.hits = set()
        self.systems = systems
        for s in systems:
            self.hits.update({vh.hit for vh in s.hits})
            self._score += s.score


    def __iadd__(self, system):
        return Solution(self.systems + [system])


    def __str__(self):
        out = f"Score of the solution = {self.score}\n"
        out += '\n'.join([f"Sys_ID={s.id} Score={s.score}" for s in self.systems])
        return out


    def is_compatible(self, system):
        """
        :param system: the system to compare to the solution
        :type system: :class:`macsypy.system.System` instance
        :return: True if the system does not share components (hits) with the systems that are part of the solution.
                 False otherwise
        """
        s2 = {vh.hit for vh in system.hits}
        return not (self.hits & s2)


    @property
    def score(self):
        """
        The score of the solution
        :rtype: float
        """
        return self._score


def compute_max_bound(systems):
    """
    Computes a grand estimation of the maximal score to be found using the set of systems as ground for a solution

    :param systems:
    :type systems: list of :class:`macsypy.system.System`
    :return: the max score of the systems
    :rtype: float
    """
    return sum([s.score for s in systems])


def find_best_solution(sorted_systems, best_sol, cur_sol):
    """

    :param sorted_systems: the systems to analyse to find the best combination of systems.
                           These systems must be sorted by their score in descending order.
    :type sorted_systems: list of :class:`macsypy.system.System` objets
    :param best_sol: The current best solution
    :type best_sol: :class:`macsypy.solution.Solution` object
    :param cur_sol: The current solution
    :type cur_sol: :class:`macsypy.solution.Solution` object
    :return: The best solution among the sorted_systems
    :rtype: :class:`macsypy.solution.Solution` object
    """
    if not sorted_systems:
        return best_sol

    max_bound = compute_max_bound(sorted_systems)
    if max_bound + cur_sol.score < best_sol.score:
        return best_sol

    system_to_test = sorted_systems.pop(0)
    if cur_sol.is_compatible(system_to_test):
        cur_sol += system_to_test
        if cur_sol.score > best_sol.score:
            best_sol = cur_sol
    return find_best_solution(sorted_systems, best_sol, cur_sol)
