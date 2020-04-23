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
import itertools
from macsypy.utils import no_stdout
from macsypy.error import MacsypyError


with no_stdout():
    # mip display version message @the import :-(
    import mip

_log = logging.getLogger(__name__)


def find_best_solution(systems):
    """
    Among the systems choose the combination of systems which does not share :class:`macsypy.hit.Hit`
    and maximize the sum of systems scores

    :param systems: the systems to analyse
    :type systems: list of :class:`macsypy.system.System` object
    :return: the list of systems which represent the best solution and the it's score
    :rtype: tuple of 2 elements ([:class:`macsypy.system.System`, ...], float score)
    """
    scores = [s.score for s in systems]
    ranks = range(len(systems))
    # create a new Model
    m = mip.Model("")
    if _log.getEffectiveLevel() > 10:
        m.verbose = 0
    # let's create variable which can take 0/1 for each system
    # and add them to the model
    x = [m.add_var(var_type=mip.BINARY) for i in ranks]

    # add constraints
    for i, j in itertools.combinations(ranks, 2):
        if not systems[i].is_compatible(systems[j]):
            m += x[i] + x[j] <= 1, f'{systems[i].id} and {systems[j].id} are not compatible'

    # We want to optimize the scores
    m.objective = mip.maximize(mip.xsum(x[i] * scores[i] for i in ranks))

    status = m.optimize()

    if status == mip.OptimizationStatus.OPTIMAL:
        _log.debug('optimal solution cost {} found'.format(m.objective_value))
        best_ranks = [i for i, v in enumerate(m.vars) if abs(v.x) > 1e-6]
        best_sol = [systems[i] for i in best_ranks]
    else:
        raise MacsypyError("No optimal solution")
    return best_sol, m.objective_value