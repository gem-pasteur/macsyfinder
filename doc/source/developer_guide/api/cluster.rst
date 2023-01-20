.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2023  Institut Pasteur (Paris), and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _cluster:


*******
cluster
*******

A cluster is an ordered set of hits related to a model which satisfy the model distance constraints.

.. _cluster_api:

cluster API reference
=====================

cluster
=======
.. autoclass:: macsypy.cluster.Cluster
   :members:
   :private-members:
   :special-members:


build_clusters
==============
.. autofunction:: macsypy.cluster.build_clusters