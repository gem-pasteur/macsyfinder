.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2020  Institut Pasteur, Paris.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _model_definition:

********************************
Macromolecular models definition
********************************

Principles
==========

MacSyFinder relies on the definition of models of macromolecular systems with an **XML grammar**
that is described :ref:`below<model-definition-grammar-label>`.

MacSyFinder relies on models of macromolecular systems.

A model is the association of several elements:
    * a definition which describe the systems to capture.
      the definitions use a xml grammar. for more information about definition see :ref:`below<model-definition-grammar-label>`.
    * a set of profile for each gene belonging the model

The models are grouped by family and eventually subfamily, for instance *secretion* or *cas protein*, ...
A models family are a model package.

.. _package_structure:
A package model follow the structure below ::

    family_name
        |_______ metadata.yml
        |_______ LICENCE
        |_______ README.md
        |_______ definitions
        |            |________ model_1.xml
        |            |________ model_2.xml
        |            :
        |
        |_______ profiles
                     |________ geneA.hmm
                     |________ geneB.hmm


If the package contains sub-family ::

    family_name
        |_______ metadata.yml
        |_______ LICENCE
        |_______ README.md
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


for concrete examples of macsy-models package visit https://github.com/macsy-models


How to install new models
=========================

MacSyFinder does not provide models. You must install Models before to use MacSyFinder.
There is a small utility tool shipped with `MacSyFinder` to deal with macsy-models: ``macsydata``.


macsydata <subcommand> [options]

The main sub-command are

* ``macsydata available`` to get the models family available
* ``macsydata search`` to search a model given its name or a pattern or in description
* ``macsydata install`` to install a package (the installed version can be set see --help)
* ``macsydata cite`` how to cite the model
* ``macsydata --help`` will give you the extended list of available subcommand

and macsydata subcommand --help to have specific help about the subcommand


Where the models are located
============================

MacSyFinder look at several locations to find macsy-models.

system wide
-----------

By default macsydata installed models in shared location (set by --install-data option)
By default in `/usr/share/macsyfinder/` or `/usr/local/share/macsyfinder` depending on your distribution.
If you use a virtualenv the shared resources are located in the `<virtualenv>/share/macsyfinder` directory.


user wide
---------

If you have not rights to installed in system wide you can installed models in the macsyfinder cache
located in your home `$HOME/.macsyfinder/data/`.
macsydata install packages in this location when you use the `--user` option.
The packages installed in user land is added to the system wide packages.
If two packages have the same name the package in user land superseded the system wide package.

project wide
------------

If you cannot install macsy_models package in system or user land locations you can specify a
specific location with the ``--models-dir`` command line option. The path must point a directory
that contains models packages.

 .. _model_package:

Write my own model-package
==========================

The whole package structure is described `above <package_structure>`_

metadata file
-------------

This file contains some meta information about the package itself.
This file is in `YAML <https://en.wikipedia.org/wiki/YAML>`_ format.
This file must have the following structure: ::

    ---
    maintainer:
      name: The name of the person who maintain/to contact for further informations. (required)
      email: The email of the maintainer (required)
    short_desc: A one line description of the package (can be used with macsydata search) (required)
    vers: the package version (required)
    cite: The list of publication to mentioned when the user have to cite the package (optional used by `macsydata cite`)
    doc: where to find extended documentation (optional)
    licence: The licence under the package is released (optional but highly recommended)
    copyright: The copyright of the package (optional)

for example ::

    ---
    maintainer:
      name: first name last name
      email: login@my_domain.com
    short_desc: Models for 15 types of secretion systems or bacterial appendages (T1SS, T2SS, T3SS, T4P, pT4SSt, pT4SSi, T5aSS, T5bSS, T5bSS, T6SSi, T6SSii, T6SSiii, Flagellum, Tad, T9SS).
    vers: 0.0a1
    cite:
      - |
        Abby Sophie S., Cury Jean, Guglielmini Julien, Néron Bertrand, Touchon Marie, Rocha Eduardo P. C. (2016).
        Identification of protein secretion systems in bacterial genomes.
        In Scientific Reports, 6, pp. 23080.
        http://dx.doi.org/10.1038/srep23080
    doc: https://github.com/macsy-models/TXSS
    licence: CC BY-NC-SA 4.0 (https://creativecommons.org/licenses/by-nc-sa/4.0/)
    copyright: 2014-2020, Institut Pasteur, CNRS

*This file is **mandatory**. Without this file your archive/repository will not considered as a package.*

.. note::

    * *-* specify an item of yaml list
    * *|* is used to specify a single item but over multi lines.


README.md
---------

A description of the package, what kind of systems the package models.
How to use it, ... in `markdown <https://guides.github.com/features/mastering-markdown/>`_ format.

LICENCE
-------

The licence use to protect and share your work.
If you don't know which licence to choose.
Have a look on `CreativeCommons <https://creativecommons.org/share-your-work/>`_
*This file is optional, but highly recommended.*

Write my definitions
--------------------


(*e.g.*, 'T1SS.xml' for T1SS, the Type 1 Secretion System) by a set of **components**
(*i.e.* proteins, or protein-coding genes given the context) with different attributes and that are used
for **content description**.
Features regarding **co-localization** parameters for system detection are also defined in this system-specific file.

Three distinct types of components can be used to model a given system content,
and which corresponds to Gene objects, and the corresponding HMM protein profiles.

* **Mandatory** components represent essential components to be found to infer the System presence.
* **Accessory** components correspond to components that can be found in some systems occurrence,
  or fastly evolving components that are hard to detect with a single profile.
* **neutral** components used to build the clusters but not take in account to build the system.
* **Forbidden** components are components which presence is eliminatory for the System assessment. 


    .. image:: ../_static/Figure1_figure_system_no_mb-new3_2col.*
     :height: 500px
     :align: left


.. _model-definition-grammar-label:

The XML hierarchy
"""""""""""""""""

* The element root is "model".

  * It has a mandatory attribute: "inter_gene_max_space", an integer representing the maximal number of components
    without a match between two components with a match for a component profile.
  * the version of the xml grammar (the actual version is "2.0")
  * The element "model" may have attributes:
  
     * **min_mandatory_genes_required**: an integer representing the minimal number of mandatory genes required
       to infer the system presence.
     * **min_genes_required**: an integer representing the minimal number of mandatory or accessory genes
       (whose corresponding proteins match a profile of the model) required to infer the system presence.
     * **max_nb_genes**: an integer representing the maximal number of mandatory or accessory genes in the system.
     * **multi_loci**: a boolean set to True ("1", "true" or "True") to allow the definition of "scattered" systems
       (systems encoded by different loci). If not specified, *default value is false*.
     
  * The model contains one or more element "gene".
  
* The element "gene" has several mandatory attributes: 

   * **name**: which must match to a profile in the profile directory.
   * **presence**: which can take three values "mandatory", "accessory", "neutral", "forbidden".

 The element "gene" may have other attributes: 

   * **loner**: which is a boolean. If a gene is loner that means this gene can be isolated on the genome ( *default false* ).
   * **exchangeable**: which is a boolean. If a gene is exchangeable (value set to "1", "true" or "True") that
     means this gene or one of its homologs or analogs can be interchanged for the assessment of the presence
     of the macromolecular system ( *default false* ).
   * **multi_system**: which is a boolean. If a gene is "multi_system" (value set to "1", "true" or "True"),
     it means that it can be used to fill by multiple systems occurrences. ( *default false* ).
   * **inter_gene_max_space**: an integer that defines gene-wise value of system's "inter_gene_max_space" parameter (see above).

 The element "gene" may have one "exchangeables" child element:

* The element "exchangeables" can contains one or more elements "gene".

Example of a model definition in XML: ::
  
  <model inter_gene_max_space="5" ver="2.0">
    <gene name="gspD" presence="mandatory">
       <exchangeables>
           <gene name="sctC"/>
       </exchangeables>
    </gene>
    <gene name="sctN_FLG" presence="mandatory" loner="1"/>
       <exchangeables>
           <gene name="gspE"/>
           <gene name="pilT"/>
       </exchangeables>
    <gene name="sctV_FLG" presence="mandatory"/>
    <gene name="flp" presence="accessory"/>
  </model>

.. warning::
  
    * a gene is identified by its name.
    * this name is case sensitive.
    * this name must be unique inside a family of models.
    * a hmm profile with the same name must be exists in the `profiles` directory


Provide hmm profiles
--------------------

For each gene mentioned in each model you have to provide a hmm profile
that capture this gene. The hmm profile must be created from specific alignment with `hmmbuild`
from `HMMER package <http://hmmer.org/>`_.
This profile *MUST* have the same name as the name of the gene mentioned in the definition.
The names are Case sensitive. All the profile must be placed in the `profiles` directory.


Share your models
=================

If you want to share your models you can create a :ref:`macsy-model package <model_package>` in your github repository
check the validity of your package with the ``macsydata check`` command
create a tag and the submit a pull request to https://github.com/macsy-models organization.
So once your PR will be accepted the model package will be automatically available through the macsydata tool.
If you don't want to submit a PR you can provide the tag release tarball (tar.gz) as is to your collaborators.
The archive is also usable with the `macsydata` tool.