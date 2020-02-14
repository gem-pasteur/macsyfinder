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
There is a small utility tool shiped with macsyfinder to deal with macsy-models: `macsydata`


macsydata <subcommand> [options]

The main sub-command are

* ``macsydata available`` to get the models family available
* ``macsydata search`` to search a model given its name or a pattern or in description
* ``macsydata install`` to install a package (we can fix version see --help)
* ``macsydata cite`` how to cite the model
* ``macsydata --help`` will give you the extende list of available subcommand

and macsydata subcommand --help to have specific help about the subcommand


Where the models are located
============================

system wide
-----------

user wide
---------

project wide
------------


Write my own model-package
==========================

metadata file
-------------

README.md
---------

LICENCE
-------

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


provide hmm profiles
--------------------


share your models
=================

