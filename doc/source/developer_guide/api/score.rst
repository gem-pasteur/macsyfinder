.. MacSyFinder - Detection of macromolecular systems in protein datasets
   using systems modelling and similarity search.
   Authors: Sophie Abby, Bertrand Néron
   Copyright © 2014-2020  Institut Pasteur (Paris), and CNRS.
   See the COPYRIGHT file for details
   MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
   See the COPYING file for details.
    
.. _score:

******
Score
******

It represents a way to score System objects this score is composed of several components,
 - the system score itself
 - the number of different genes composing the system which are used in other systems
 - the number of hit implied in different systems

and select the best systems based on this composed score


ComposedScore
=============
.. autoclass:: macsypy.score.ComposedScore
   :members:
   :private-members:
   :special-members:

BestSystemSelector
==================

The heuristic to pick the best systems based on the :class:`ComposedScore` object

.. autoclass:: macsypy.score.BestSystemSelector
   :members:
   :private-members:
   :special-members: