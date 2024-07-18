.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2023  Institut Pasteur (Paris), and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _hit:


***
hit
***

This module implements class relative to hit and some functions to do some computation on hit objects.

=========================================== =============================================================================
:class:`macsypy.hit.CoreHit`                Modelize a hmm hit on the replicon. There is only one Corehit for a CoreGene.
:class:`macsypy.hit.ModelHit`               Modelize a hit and its relation to the Model.
:class:`macsypy.hit.AbstractCounterpartHit` Parent class of Loner, MultiSystem. It's inherits from ModelHit.
:class:`macsypy.hit.Loner`                  Modelize "true" Loner.
:class:`macsypy.hit.MultiSystem`            Modelize hit which can be used in several Systems (same model)
:class:`macsypy.hit.LonerMultiSystem`       Modelize a hit representing a gene Loner and MultiSystem at same time.
:class:`macsypy.hit.HitWeight`              The weights apply to the hit to compute score
:func:`macsypy.hit.get_best_hit_4_func`     Return the best hit for a given function
:func:`macsypy.hit.sort_model_hits`         Sort hits
:func:`macsypy.hit.compute_best_MSHit`      Choose among svereal multisystem hits the best one
:func:`macsypy.hit.get_best_hits`           If several profile hit the same gene return the best hit
=========================================== =============================================================================

A Hit is created when `hmmsearch` find similarities between a profile and protein of the input dataset

Below the ingheritance diagram of Hits

.. inheritance-diagram::
      macsypy.hit.CoreHit
      macsypy.hit.ModelHit
      macsypy.hit.AbstractCounterpartHit
      macsypy.hit.Loner
      macsypy.hit.MultiSystem
      macsypy.hit.LonerMultiSystem
   :parts: 1


And a diagram showing the interaction between CoreGene, ModelGene, Model, Hit, Loner, ... interactions

.. figure:: ../../_static/gene_obj_interaction.*


    The diagram above represents the models, genes and hit generated from the definitions below.

    .. code-block::

        <model name="A" inter_gene_max_space="2">
            <gene name="abc" presence="mandatory"/>
            <gene name="def" presence="accessory"/>
        </model>

        <model name="B" inter_gene_max_space="5">
            <gene name="def" presence="mandatory"/>
                <exchangeables>
                    <gene name="abc"/>
                </exchangeables>
            <gene name="ghj" presence="accessory"
        </model>





.. _hit_api:

hit API reference
=================

CoreHit
=======
.. autoclass:: macsypy.hit.CoreHit
   :members:
   :private-members:
   :special-members:

ModelHit
========
.. autoclass:: macsypy.hit.ModelHit
   :members:
   :private-members:
   :special-members:

AbstractCounterpartHit
======================
.. autoclass:: macsypy.hit.AbstractCounterpartHit
   :members:
   :private-members:
   :special-members:

Loner
=====
.. autoclass:: macsypy.hit.Loner
   :members:
   :private-members:
   :special-members:

MultiSystem
===========
.. autoclass:: macsypy.hit.MultiSystem
   :members:
   :private-members:
   :special-members:

LonerMultiSystem
================
.. autoclass:: macsypy.hit.LonerMultiSystem
   :members:
   :private-members:
   :special-members:

HitWeight
=========
.. autoclass:: macsypy.hit.HitWeight
   :members:
   :private-members:
   :special-members:

get_best_hit_4_func
===================
.. autofunction:: macsypy.hit.get_best_hit_4_func

sort_model_hits
===============
.. autofunction:: macsypy.hit.sort_model_hits

compute_best_MSHit
==================
.. autofunction:: macsypy.hit.compute_best_MSHit

get_best_hits
=============
.. autofunction:: macsypy.hit.get_best_hits
