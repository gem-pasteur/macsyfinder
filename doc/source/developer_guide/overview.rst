.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2025  Institut Pasteur (Paris),and CNRS.
    See the COPYRIGHT file for details
    MacSyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _overview:

MacSyFinder implementation overview
===================================

MacSyFinder is implemented with an object-oriented architecture.
Below a short glossary to fix the vocabulary used in MacSyFinder

.. glossary::
    :sorted:

    Model family
        A set of models, on the same topic.
        It is composed of several definitions which can be sorted in hierachical structure
        and profiles. A profile is a hmm profile file.

    Model
        Is a formal description of a macromolecular system.
        Is composed of a definition and a list of profiles.
        at each gene of the Model must correspond a profile

    ModelDefinition

        Is a definition of model, it's serialize as a xml file

    Cluster

        Is a "contiguous" set of hits.
        two hits are considered contiguous if the number of genes between the 2 genes matching the 2 hits
        in the replicon is lesser than inter-genes-max-space.

    System

        It's an occurrence of a specific Model on a replicon.
        Basically, it's a cluster or set of clusters which satisfy the Model quorum.

    Solution

        It's a systems combination for one replicon.
        The best solution for a replicon, is the combination of all systems found in this replicon which
        maximize the score.


MacSyFinder project structure
-----------------------------

A brief overview of the files and directory constituting the MacSyFinder project

.. glossary::

    doc
        The project is documented using sphinx.
        All sources files needed to generate this documentation is in the directory *doc*

    macsyfinder
        This the MacSyFinder python library.
        Inside macsyfinder there is a subdirectory *scripts* which are the entry points for
        `macsyfinder`, `msf_data`, `msf_profile`, ...

    tests
        The code is tested using `unittests`.
        In *tests* the directory *data* contains all data needed to perform the tests.

    utils
        Contains a binary `setsid` needed macsyfinder to parallelize some steps.
        Usually `setsid` is provides by the system, but some macOS version does not provide it.

    CITATION.yml
        A file indicating how to cite macsyfinder in yaml format.

    CONTRIBUTORS
        A file containing the list of code contributors.

    CONTRIBUTING
        A guide on how to contribute to the project.

    COPYRIGHT
        The macsyfinder copyrights.

    COPYING
        The licencing.
        MacSyFinder is released under GPLv3.

    README.md
        Brief information about the project.

    setup.py
        The installation recipe specific to embeded doc in distrib.

    pyproject.toml
        a configuration file for packaging-related tools (as well as other tools).


MacSyFinder architecture overview
---------------------------------

Starting with `2.1.5` version MacSyFinder is build on top of `MacSyLib` for the internal read
`MacSyLib documentation <https://macsylib.readthedocs.io/en/latest/>`_
