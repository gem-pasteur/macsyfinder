.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2023 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _model_package:

**************
Models Package
**************


MacSyFinder relies on the definition of models of macromolecular systems as a **set of models' components**
to be searched by similarity search, and a **set of rules** regarding their genomic organization and
their requirement level to make a complete system (mandatory, accessory components, number of components required).

See the section :ref:`model-definition-grammar-label` for more details on MacSyFinder's modelling scheme and the section
on :ref:`Functioning <functioning>` for the principles of the MacSyFinder's search engine.


A **MacSyFinder model** (macsy-model for short) is the association of several elements:

    * a **definition** which describes the system to detect with a specific **XML grammar** that is :ref:`described here<model-definition-grammar-label>`.

    * a set of :ref:`HMM profiles <provide-hmm_label>`  (one per component/gene in the model) to enable the similarity search of the systems' components with the HMMER program.

The models are grouped by *family* possibly gathering *sub-families* (multiple levels allowed), for instance *Secretion*, *Cas-proteins*...
A set of models from a same family (coherent set) of systems to detect is called hereafter a **macsy-model package** ``NEW in V2``.



.. _package_structure:


Structure of a macsy-model package
==================================

A macsy-model package follows the following structure: ::

    family_name
        |_______ metadata.yml
        |_______ LICENSE
        |_______ README.md
        |_______ model_conf.xml
        |_______ definitions
        |            |________ model_1.xml
        |            |________ model_2.xml
        |            :
        |
        |_______ profiles
                     |________ geneA.hmm
                     |________ geneB.hmm
                     |________ geneC.hmm.gz


If the package contains sub-families: ::

    family_name
        |_______ metadata.yml
        |_______ LICENSE
        |_______ README.md
        |_______ model_conf.xml
        |_______ definitions
        |            |________ subfamilyA
        |            |            |________ model_1.xml
        |            |            |________ model_2.xml
        |            |
        |            |________ subfamilyB
        |            |            |________ model_3.xml
        |            |            |________ model_4.xml
        |            |
        |            :
        |
        |_______ profiles
                     |________ geneA.hmm
                     |________ geneB.hmm
                     |________ geneC.hmm.gz


For examples of macsy-model packages, please visit https://github.com/macsy-models

You can create a template for your package by using `macsydata init`.
It will create for you:

* the data package directory with the right structure.
* a template of `metadata.yaml` .
* a template of `README.md` file.
* a generic `model_conf.xml` file.
* a LICENSE file if `--license` option is set.
* a COPYRIGHT file if `--holders` option is set.
* a directory `definitions` with an example of model definition (model_example.xml to remove before publishing).
* a directory `profiles` where to put the hmm profiles corresponding to the models genes.

.. note::

    MSF can also read *.gz* compressed files; it will uncompress them on the fly.
    The compressed files must end with the *.gz* extension.
    For the `hmmsearch` step You need to have `gunzip` installed on your system for this to work.


README.md
---------

A description of the package: what kind of systems the package models,
how to use it etc... in `markdown <https://guides.github.com/features/mastering-markdown/>`_ format.
The Readme is displayed to the user on the macsy-models repository on Github.
It is also displayed when the user runs `macsydata help`.


LICENSE
-------

The license is used to protect your work when sharing it.
If you don't know which license to choose, have a look at `CreativeCommons <https://creativecommons.org/share-your-work/>`_
*This file is optional, but highly recommended.*


Metadata file
-------------

The `metadata.yml` file contains some meta information about the package itself.

It is in `YAML <https://en.wikipedia.org/wiki/YAML>`_ format and must have the following structure:

.. code-block:: yaml

    ---
    maintainer:
      name: The name of the person who maintains/to contact for further information. (required)
      email: The email of the maintainer (required)
    short_desc: A one line description of the package (can e.g. be used for *macsydata* searches) (required)
    vers: The package version (DEPRECATED)
    cite: The publication(s) to cite by the user when the package is used (optional, used by `macsydata cite`)
    doc: Where to find extended documentation (optional)
    license: The license under the package is released (optional but highly recommended)
    copyright: The copyright of the package (optional)

For example:

.. code-block:: yaml

    ---
    maintainer:
       name: first name last name
       email: login@my_domain.com
    short_desc: Models for 15 types of secretion systems or bacterial appendages (T1SS, T2SS, T3SS, T4P, pT4SSt, pT4SSi, T5aSS, T5bSS, T5bSS, T6SSi, T6SSii, T6SSiii, Flagellum, Tad, T9SS).
    cite:
       - |
         Abby Sophie S., Cury Jean, Guglielmini Julien, Néron Bertrand, Touchon Marie, Rocha Eduardo P. C. (2016).
         Identification of protein secretion systems in bacterial genomes.
         In Scientific Reports, 6, pp. 23080.
         http://dx.doi.org/10.1038/srep23080
    doc: https://github.com/macsy-models/TXSS
    license: CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/)
    copyright: 2014-2022, Institut Pasteur, CNRS


.. note::

    * *-* specify an item of yaml list
    * *|* is used to specify a single item but over multiple lines.


.. error::
    This `metadata.yml` file is **mandatory**. Without this file your archive/repository will not be considered as a *macsy-model package*.


.. warning::

    The field *vers* (the package version) is deprecated. *macsydata install* rely only on the git tag.


.. _model_configuration:

Model configuration
-------------------

The modeler has the possibility to specify some options that are specific to their package,
different than the MacSyFinder defaults in the `model_conf.xml` file. ``NEW in v2``

These options can be grouped in two families: the scoring weights and filtering options.

Scoring weights:

    * mandatory (*float* default = 1.0)
    * accessory (*float* default = 0.5)
    * exchangeable (*float* default = 0.8)
    * loner_multi_systems (*float* default =  0.7)
    * redundancy_penalty (*float* default = 1.5)

Filtering options:

    * e_value_search (*float* default = 0.1)
    * i_evalue_sel (*float* default = 0.001)
    * profile_coverage (*float* default = 0.5)
    * cut_ga (*bool* default = True)

All these options are optional and can be omitted in the configuration file, **the file itself is optional**.
The precedence rules between the different levels of configuration are:


.. code-block:: text

 system < home < model < project < --cfg-file | --previous-run < command line options

* **system**: the `macsyfinder.conf` file either in /etc/macsyfinder/ or in ${VIRTUAL_ENV}/etc/macsyfinder/
  in case of a *virtualenv* this configuration affects only the MacSyFinder version installed in this virtualenv
* **home**:  the `~/.macsyfinder/macsyfinder.conf` file
* **model**: the `model_conf.xml` file at the root of the model package
* **project**: the `macsyfinder.conf` file found in the directory where the `macsyfinder` command was run
* **cfgfile**: any configuration file specified by the user on the command line (conflicts with the `--previous-run` option)
* **previous-run**: the `macsyfinder.conf` file found in the results directory of the previous run (conflicts with the `--cfg-file` option)
* **command line**: any option specified directly in the command line

The model_conf.xml configuration file is in xml format and must have the following structure:

.. code-block:: xml

    <model_config>
        <weights>
            <mandatory>1</mandatory>
            <accessory>0.5</accessory>
            <exchangeable>0.8</exchangeable>
            <redundancy_penalty>1.5</redundancy_penalty>
            <out_of_cluster>0.7</out_of_cluster>
        </weights>
        <filtering>
            <e_value_search>0.1</e_value_search>
            <i_evalue_sel>0.01</i_evalue_sel>
            <coverage_profile>0.5</coverage_profile>
            <cut_ga>True</cut_ga>
        </filtering>
    </model_config>


:ref:`Details about the scoring method can be obtained here <combinatorial-exploration>`.
