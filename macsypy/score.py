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
from macsypy.error import MacsypyError


class ComposedScore:
    """
    Modelize a system score based on several dimensions

        - the system score itself
        - the number of overlapping genes (number of genes which are represented several times)
        - the length of overlapping (number of occurrence of overlapping genes)

          *\- H1\- \- H2\- \- H3\- \- H1 \- \- H2' \- \- H1 \-* (H2' code for a exchangeable of H2)

        - overlapping genes = H1 + H2 = 2
        - overlapping length = 3 H1 + 2  H2 = 5
    """

    def __init__(self, system, hit_tracker):
        """
        :param system: the system to compute the composed score
        :type system: :class:`macsypy.system.System` object
        :param hit_tracker: the hit tracker build for this run
        :type hit_tracker: :class:`macsypy.hit.HitTracker` object
        """
        self._system = system
        self._hit_tracker = hit_tracker
        self._sys_score = self._system.score
        self._overlapping_genes = None
        self._overlapping_length = None
        self._overlap()

    def _overlap(self):
        """
        Compute overlapping genes and overlapping length
        """
        used_in_systems = {}
        for vh in self._system.hits:
            used_in_systems[vh] = {s.id for s in self._hit_tracker[vh.hit] if s.model.fqn != self._system.model.fqn}
        self._overlapping_genes = sum([1 for vh in used_in_systems if used_in_systems[vh]])
        self._overlapping_length = sum([1 for used_in in used_in_systems.values() for h in used_in])

    @property
    def system(self):
        """
        :return: The system related to this score
        :rtype: :class:`macsypy.system.System`
        """
        return self._system

    @property
    def sys_score(self):
        """
        :return: the system score itself
        :rtype: float
        """
        return self._sys_score


    @property
    def overlapping_genes(self):
        """
        :return: The number of overlapping genes (number of genes which are represented several times
                 in **different models**)
        :rtype: int
        """
        return self._overlapping_genes

    @property
    def overlapping_length(self):
        """
        :return: The length of overlapping (number of occurrence of overlapping genes)
        :rtype: int
        """
        return self._overlapping_length


class BestSystemSelector:
    """
    Sort systems and filter them to proposed only the best ones
    """

    def __init__(self, systems, hit_tracker):
        """

        :param systems: The systems to sort, filter. The systems must be occurrence of the **same model**
        :type systems: list of :class:`macsypy.system.System` object
        :param hit_tracker: the hit tracker build for this run
        :type hit_tracker: :class:`macsypy.hit.HitTracker` object
        """
        models = {sys.model.fqn for sys in systems}
        if len(models) != 1:
            raise MacsypyError(f"Cannot build Score with system from different models: {','.join(models)}")
        self. systems = systems
        self.hit_tracker = hit_tracker


    def best_system(self):
        """
        :return: the best system
        :rtype: list of :class:`macsypy.system.System` objects

        best system use an heuristic to select the best systems.

            1. select the systems with the highest score
            2. if several systems have same score, choose among them
               the system with the fewest different overlapping genes (:attr:`macsypy.score.ComposedScore.overlapping_genes`)
            3. if there is still several systems, pick the systems
               with the fewest number of genes overlapping (:attr:`macsypy.score.ComposedScore.overlapping_length`)
        """
        systems = sorted(self.systems, key=lambda s: - s.score)
        grouped = itertools.groupby(systems, key=lambda s: s.score)
        score, best_systems = next(grouped)
        best_systems = list(best_systems)
        if len(best_systems) > 1:
            best_score = [ComposedScore(sys, self.hit_tracker) for sys in best_systems]
            criterion = lambda cs: cs.overlapping_genes
            best_score = itertools.groupby(sorted(best_score, key=criterion), key=criterion)
            _, best_score = next(best_score)
            best_score = list(best_score)
            if len(best_score) > 1:
                criterion = lambda cs: cs.overlapping_length
                best_score = itertools.groupby(sorted(best_score, key=criterion), key=criterion)
                _, best_score = next(best_score)
                return [cs.system for cs in best_score]
            else:
                return [cs.system for cs in best_score]
        else:
            return best_systems

