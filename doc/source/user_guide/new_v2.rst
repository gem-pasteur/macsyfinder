.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2023 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _new_v2:

*****************************
What's new in MacSyFinder v2?
*****************************

For Version 2, MacSyFinder was carried under `Python 3 <https://www.python.org/download/releases/3.0/>`_.


==============================
New features and search engine
==============================

MacSyFinder v2 is a major release. The **search engine** was changed for a more intuitive and comprehensive exploration of putative systems.

The search is now more thorough and avoid undesirable side-effects of the previous search engine. Being more thorough, it now also
includes a **scoring scheme** to build candidate systems from sets of detected components (clusters), and can offer several optimal "solutions" (sets of
detected systems) based on a combinatorial exploration of detected clusters.
See :ref:`here for more details <functioning>`.

.. warning::

  The search engine being different, one might want to check that models carried from v1 to v2 have the expected behaviour.


Several **new features** were added, including:

- a **new type of gene component** "neutral" was added in order to provide more possibilities for systems' modelling in macsy-models. :ref:`See here <components>` for more details.
- a **new component feature** was introduced: "multi-model", that corresponds to components that are allowed to participate in occurrences of systems from different models. :ref:`See here <multi-model-label>` for more.
- more flexibility was introduced in the **search for systems' components using HMMER**. It is now possible to use the `cut_ga` threshold when provided in the HMM profiles used for components' similarity search. This enables to have a search tailored for each HMM profile, and thus component. :ref:`See here <hmmer-options>` for more details.
- a **new file structure** was created to better organize MacSyFinder's packages (i.e. that include systems' models and corresponding HMMER profiles). :ref:`See here <package_structure>` for details.
- a **tool** to easily install and distribute MacSyFinder's packages was created. :ref:`See here <macsydata>` for more details on *macsydata*.
- the **format for MacSyFinder's models** has slightly changed, in order to offer more possibilities, and more readibility. To see **how to carry models from v1 to v2**, :ref:`visit here <models_v1_v2>`.


Also, the search modes corresponding to "unordered" and "unordered_replicon" were merged into the **"unordered"** search mode - as they basically correspond to the same behaviour.

.. note::

 In v2, output files were also re-defined. See :ref:`here for more details <outputs>`.



Dependencies
============

MacSyFinder v2 no longer requires the *formatdb* or *makeblastdb* tools from NCBI.
However, new dependencies are used, but as they are Python libraries, it should be transparent for the user, and not require manual installations. See :ref:`here for details<dependencies>`.





Models are more formalized
==========================

The models data are more formalized, with a well defined structure.
For instance the definitions and profiles must be packed together in what we call a `macsy-model` package
If you intend to model new systems please refer to the :ref:`modeler_guide`.



Models installation
===================

We now provide a new tool to manage the models. See :ref:`macsydata`.


Models configuration
====================

The modeler can provide some spcific configuration values released along the model package. See :ref:`model_configuration`.


Modeller helper tool
====================

To help modellers create new models we provide a new helper tool `macsyprofile`, which analyses HMMER raw output files from
results of a previous MacSyFinder run, to provide information on all hits even if filtered out. See :ref:`macsyprofile`.

:ref:`macsydata` provide also some options to help the modeller as

* **macsydata init** to init a new model package.
* **macsydata check** to check the integrity of a model package, before to use/publish it.

