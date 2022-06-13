.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2022 Institut Pasteur (Paris), CNRS.
    See the COPYRIGHT file for details
    MacSyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

    Macsyfinder documentation master file, created by sphinx-quickstart

.. image:: _static/logo_macsyfinder.png
   :align: center
   :alt: MacSyFinder logo



Welcome to MacSyFinder's documentation!
=======================================


.. note::

  A **new version of MacSyFinder (v2)** is available, see :ref:`here for an overview of the novelties<new_v2>`. The search engine was changed, and some bugs/unwanted behaviors corrected.
  MacSyFinder's models for v2 are very similar, yet not compatible with those from v1. See here for details on :ref:`how to carry your models to v2<models_v1_v2>`.

  The search engine of v2 being much different from that of v1, we **strongly suggest** to test
  whether the results are relevant by simply "translating" the models from v1 to v2, or if the models need to be adapted to correctly function with v2.


MacSyFinder is a program to **model and detect macromolecular systems, genetic pathways**... in protein datasets.
In prokaryotes, these systems have often evolutionarily conserved properties:

- they are made of **conserved components**,
- they are encoded in **compact loci** (conserved genetic architecture).

The user models these systems with MacSyFinder to reflect these conserved features, and to allow their efficient detection.

Criteria for systems detection include **component content (quorum)**, and **genomic co-localization**.
Each component corresponds to a hidden Markov model (HMM) protein profile to
perform sequence similarity searches with the program Hmmer.

In order to model macromolecular systems, the user:

    - builds or gather from databanks **HMM protein profiles** for components of interest,
    - defines **decision rules** for each system in a dedicated XML grammar (see :ref:`model_definition`).


    .. figure:: _static/figure_main-2.*
        :height: 800px
        :align: center


.. note::

  If you use MacSyFinder, please cite:

  `Abby SS, Néron B, Ménager H, Touchon M, Rocha EPC (2014).
  MacSyFinder: A Program to Mine Genomes for Molecular Systems with an Application to CRISPR-Cas Systems.
  PLoS ONE 9(10): e110726. doi:10.1371/journal.pone.0110726 <http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0110726>`_

==========
User Guide
==========
.. toctree::
   :maxdepth: 2

   user_guide/index


==============
Modeller Guide
==============

.. toctree::
   :maxdepth: 2

   modeler_guide/index


===============
Developer Guide
===============

.. toctree::
   :maxdepth: 2

   developer_guide/index



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
