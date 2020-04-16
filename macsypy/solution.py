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

import logging

_log = logging.getLogger(__name__)


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
        self._hits = set()
        self._systems = systems
        for s in systems:
            self._hits.update({vh.hit for vh in s.hits})
            self._score += s.score

    def __eq__(self, other):
        # 2 solutions are equal if they are composed with the same systems
        return {s.id for s in self.systems} == {s.id for s in other.systems}

    def __iadd__(self, system):
        return Solution(self._systems + [system])


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

    @property
    def hits(self):
        """
        :return: The hits composing the solution
                 This is all hits from all systems used to build this solution
        :rtype: set of :class:`macsypy.hit.Hit` object`
        """
        return self._hits

    @property
    def systems(self):
        """
        :return: The systems composing the solution
        :rtype: set of :class:`macsypy.systems.Systems` object`
        """
        return self._systems


def compute_max_bound(systems):
    """
    Computes a grand estimation of the maximal score to be found using the set of systems as ground for a solution

    :param systems:
    :type systems: list of :class:`macsypy.system.System`
    :return: the max score of the systems
    :rtype: float
    """
    return sum([s.score for s in systems])


def find_best_solution(sorted_systems, best_sol, cur_sol, branch=0):
    """

    :param sorted_systems: the systems to analyse to find the best combination of systems.
                           These systems must be sorted by their score in descending order.
    :type sorted_systems: list of :class:`macsypy.system.System` objets
    :param best_sol: The current best solution
    :type best_sol: :class:`macsypy.solution.Solution` object
    :param cur_sol: The current solution
    :type cur_sol: :class:`macsypy.solution.Solution` object
    :param branch: identifier of the solution branch in the solution space
                   This parameter should not be set it is used only for debugging
                   and should be set to 0 when the function is not a recursive call
    :return: The best solution among the sorted_systems
    :rtype: :class:`macsypy.solution.Solution` object
    """
    _log.debug("######################### find_best_solution ###################################")
    _log.debug(f"### {branch} ## sorted_systems {[(s.id, s.score) for s in sorted_systems]}")
    _log.debug(f"### {branch} ## best_sol {[(s.id, s.score) for s in best_sol.systems]} score {best_sol.score}")
    _log.debug(f"### {branch} ## cur_sol {[(s.id, s.score) for s in cur_sol.systems]} score {cur_sol.score}")
    if not sorted_systems:
        _log.debug(f"### {branch} ## NO more systems to explore RETURN best_sol ##\n{best_sol}")
        return best_sol

    # It's IMPORTANT to do a copy of sorted_systems
    # otherwise (for instance using pop to get the first element)
    # the side effect will empty the sorted_systems even outside the function :-(
    system_to_test = sorted_systems[0]
    remaining_syst = sorted_systems[1:]

    max_bound = compute_max_bound(sorted_systems)
    if max_bound + cur_sol.score < best_sol.score:
        _log.debug(f"### {branch} ## max_bound + cur_sol.score < best_sol.score "
                   f"{max_bound} + {cur_sol.score} < {best_sol.score} Stop exploring this branch")
        _log.debug(f"### {branch} ## RETURN best_sol \n{[s.id for s in best_sol.systems]}")
        return best_sol

    if cur_sol.is_compatible(system_to_test):
        _log.debug(f"### {branch} ## cur_sol is compatible with {system_to_test.id}")
        cur_sol += system_to_test
        if cur_sol.score > best_sol.score:
            best_sol = cur_sol
        _log.debug(f"### {branch} ## best sol so far ##\n{best_sol}")
        return find_best_solution(remaining_syst, best_sol, cur_sol, branch=branch)
    else:
        _log.debug(f"### {branch} ## cur_sol is NOT compatible with {system_to_test.id}")
        _log.debug(f"### {branch} ## lets explore new branch of solutions")
        is_the_best = find_best_solution(sorted_systems, best_sol, Solution([]), branch=branch + 1)
        if is_the_best.score > best_sol.score:
            _log.debug(f"### {branch} ## The solution from the new branch become the best solution")
            best_sol = is_the_best
        _log.debug(f"### {branch} ## continue to explore the branch {branch}")
        return find_best_solution(remaining_syst, best_sol, cur_sol, branch=branch)