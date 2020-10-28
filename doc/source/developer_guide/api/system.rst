.. MacSyFinder - Detection of macromolecular systems in protein datasets
   using systems modelling and similarity search.
   Authors: Sophie Abby, Bertrand Néron
   Copyright © 2014-2020  Institut Pasteur (Paris), and CNRS.
   See the COPYRIGHT file for details
   MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
   See the COPYING file for details.
    
.. _system:

******
system
******

This module classes and functions which a given set of hits and a model
compute if this set satisfy the model or not

The object which check the compliance of hits to a model is MatchMaker
which have 2 sub-classes for ordered and unordered replicons

MatchMaker.match method link hit to a model (:class:`macsypy.hit.ValidHit`)
and then check if these valid hit satisfy the quorum constraints defined
in the model. According this it instanciate a :class:`macsypy.system.System`
or :class:`macsypy.system.RejectedClusters` for ordered replicons
or :class:`macsypy.system.LikelySystem` or :class:`macsypy.system.UnlikelySystem`
for unordered replicons

below the inheritance diagram:

.. inheritance-diagram::
      macsypy.system.System
      macsypy.system.RejectedClusters
      macsypy.system.LikelySystem
      macsypy.system.UnlikelySystem
   :parts: 1


.. warning::
   The abstract class :class:`macsypy.system.AbstractSetOfHits` is controlled by the metaclass
   :class:`macsypy.system.MetaSetOfHits` which inject on the fly several private attributes and
   public properties (see more in :class:`macsypy.system.MetaSetOfHits` documentation)

System
=======
.. automodule:: macsypy.system
   :members:
   :private-members:
   :special-members:

