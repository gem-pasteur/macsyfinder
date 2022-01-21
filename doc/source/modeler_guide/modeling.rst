.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2022 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _modeling:

*********************
Macromolecular models
*********************


MacSyFinder relies on the definition of models of macromolecular systems as a **set of models' components** 
to be searched by similarity search, and a **set of rules** regarding their genomic organization and 
their requirement level to make a complete system (mandatory, accessory components, number of components required). 

See :ref:`below<model-definition-grammar-label>` for more details on MacSyFinder's modelling scheme and the section 
on :ref:`Functioning <functioning>` for the principles of the MacSyFinder's search engine.


A **MacSyFinder model** (macsy-model for short) is the association of several elements:

    * a **definition** which describes the system to detect with a specific **XML grammar** that is described :ref:`below<model-definition-grammar-label>`.
    
    * a set of :ref:`HMM profiles <provide-hmm_label>`  (one per component/gene in the model) to enable the similarity search of the systems' components with the HMMER program.

The models are grouped by *family* possibly gathering *sub-families* (multiple levels allowed), for instance *Secretion*, *Cas-proteins*...
A set of models from a same family (coherent set) of systems to detect is called hereafter a **macsy-model package** ``NEW in V2``.



.. _writing-models:

Principles, and how to write macsy-models definitions
=====================================================

Macsy-models are written as XML files, and should be named with the name of the system to detect as a prefix,
and the XML file extension as a suffix. For example, 'T1SS.xml' for T1SS (Type I Secretion System).

A macsy-model defines a macromolecular System as:

* A set of **components** (*i.e.* proteins, or protein-coding genes given the context) with different attributes that are used for system's **content description**.
* Features regarding the **genomic architecture** of the systems' components for system detection.
* Rules for **quorum** specifying how many components are required to infer the presence of a complete system.


.. _components:

Macsy-model Components
----------------------

Four distinct **types of components** can be used to model the System's content.
Components correspond to Gene objects in MacSyFinder's implementation, and point to corresponding HMM protein profiles.

* **mandatory** components represent components that are essential to be found to infer the system's presence.
* **accessory** components correspond to components that can be found in some systems' occurrence
  (or quickly evolving components that are hard to detect with a single HMM profile and thus can be missed along similarity search).
* **neutral** components are used to build/extend clusters of proximal genes/components on the replicon analysed, but are not part of the quorum (i.e., not taken into account to assess the system's presence). ``NEW in V2``
* **forbidden** components are components which presence is eliminatory for the system's presence assessment.


.. _model-definition-genomic-orga:

Specifying a genomic organization
---------------------------------

Beyond its list of Components, a MacSyFinder's model of a System is defined by the genomic organization of its components.
This genomic organization can be defined in several ways:

* the general System's architecture, whether it is `single-locus` or `multi-loci` (encoded at one or several loci)
* the co-localization criteria defined either at the System level or at the Gene (component) level:

    * the `inter-gene-max-space` parameter (system- or gene- wise)
    * the `loner` parameter (gene- wise)


See :ref:`below<model-definition-grammar-label>` for more details on how to specify these parameters in a macsy-model.


.. _model-definition-grammar-label:

The XML hierarchy
-----------------

A System's model is defined using a specific XML grammar that is hereby described.
It consists in a hierarchic view of a Model that has specific features described through parameters, and is made of a set of Genes that have specific features themselves.
All these elements and corresponding parameters will parametrize the search of Systems matching the search by MacSyFinder, in terms of Gene content and genomic architecture criteria.


.. image:: ../_static/MSF_modelling.*
    :height: 1000px
    :align: left


* The element root of a System's model is "model".

  * It has a mandatory attribute: "inter_gene_max_space", an integer representing the maximal number of components
    without a match between two components with a match for a component profile in order to consider them contiguous (part of a same *Cluster*).
  * The version of the XML grammar (the actual version is "2.0")
  * The element "model" may have attributes:

     * **min_mandatory_genes_required**: an *integer* representing the minimal number of mandatory genes required
       to infer the system's presence.
     * **min_genes_required**: an *integer* representing the minimal number of mandatory or accessory genes
       (whose corresponding proteins match a profile of the model) required to infer the system's presence.
     * **multi_loci**: a *boolean* set to True ("1", "true" or "True") to allow the definition of "scattered" systems
       (i.e., systems encoded at different genomic loci or by different gene *clusters*). If not specified, *default value is false*.
     * **max_nb_genes** define how many genes is necessary to consider a system as full. By default it is the sum of mandatory and accessory genes.
       But sometimes in special cases, there is 2 profiles, so 2 *msf* genes in model for one real gene. So in system only one gene can be detected
       and the whaleness is false.

  * The model contains one or more element(s) "gene" that correspond(s) to the genetic components of the macromolecular system.

* The element "gene" has several mandatory attributes:

   * **name**: a *string* representing the name of the component/gene which must match that of a profile enclosed in the profile directory of the macsy-model package (see :ref:`below <provide-hmm_label>`).
   * **presence**: a *string* representing the status of the gene's presence in the system. It can take four values among "mandatory", "accessory", "neutral", "forbidden" (see above).

 The element "gene" may have other attributes:

   * **loner**: a *boolean*. A *loner* gene can be isolated on the genome and does not have to be part of a cluster of genes to be considered for system's assessment ( *default false* ).
   * **multi_system**: a *boolean*. If a gene has the feature "multi_system" (value set to "1", "true" or "True"),
     it means that it can be used to fill multiple systems' occurrences - and thus be considered part of several systems. ( *default false* ).
   * **inter_gene_max_space**: an *integer* that defines gene-wise value of system's "inter_gene_max_space" parameter (see above). It supersedes the system-wise parameter to give the gene a specific co-localization parameter.

.. _exchangeables_label:

 The element "gene" may have one "exchangeables" child element:

   * The element "exchangeables" can contain one or more elements "gene".

   For a Gene to have "exchangeables" Genes listed, means that this Gene can be replaced *in the quorum* by the listed child Genes.



.. note::

  If not specified by the user, several features will have their values assigned **by default**:

  * the **genomic architecture** of the System being searched will consist in a **single locus**. If a System may be made of Genes from multiple loci, consider setting the `multi_loci` parameter to `True`.
  * the **quorum parameters** `min_mandatory_genes_required` and `min_genes_required` will be set to the number of mandatory Genes listed - the `accessory` Genes being deemed not required to infer a complete System.




Example of a macsy-model definition in XML:

.. code-block:: xml

  <model inter_gene_max_space="5" vers="2.0">
    <gene name="gspD" presence="mandatory">
       <exchangeables>
           <gene name="sctC"/>
       </exchangeables>
    </gene>
    <gene name="sctN_FLG" presence="mandatory" loner="1">
       <exchangeables>
           <gene name="gspE"/>
           <gene name="pilT"/>
       </exchangeables>
    </gene>
    <gene name="sctV_FLG" presence="mandatory"/>
    <gene name="flp" presence="accessory"/>
  </model>



In this example, the described System consists of three mandatory and one accessory components:

  * Two components, the Gene "GspD" and the Gene "sctN_FLG" can respectively be replaced by sctC, and gspE and pilT genes in the quorum.
  * To be considered as part of such System, the components should be co-localized in loci (Clusters of Genes), which in this case would amount to being located from each other at a distance of 5-Genes maximum, except for the Gene "sctN_FLG" that is allowed to be located "alone" in the genome being investigated, by a `loner` parameter being set to True. As the `multi_loci` parameter is not set, by default the System should be made of a single locus (Cluster of co-localized Genes - except for the ones listed as `loners`).
  * To be considered a complete System, the quorum of Genes should be reached. In this case, the `min_genes_required` and `min_mandatory_genes_required` are not specified and therefore assigned to their default values: `min_mandatory_genes_required` is set to the number of mandatory Genes listed as well as the `min_genes_required` parameter (see above).


.. warning::

    * a gene is identified by its name.
    * this name is case sensitive.
    * this name must be unique inside a family of models.
    * a HMM profile with a gene-based name must exist in the `profiles` directory of the macsy-model package (see :ref:`below <provide-hmm_label>`).


.. _provide-hmm_label:

Providing HMM profiles
----------------------

For each gene mentioned in each model you have to provide **a HMM profile**
to enable the similarity search of this gene. The HMM profile must have been created by the user from a curated multiple sequence alignment with the `hmmbuild` program
from the `HMMER package <http://hmmer.org/>`_, or can have been obtained from HMM profiles' databases such as `TIGRFAM <https://dx.doi.org/10.1093%2Fnar%2Fgkg128>`_ or `PFAM <https://pfam.xfam.org/>`_ .

This profile *MUST* have the same name as the name of the gene mentioned in the definition.
For instance, a component named "GeneA" in the macsy-model would correspond by default to a HMM profile "GeneA.hmm" enclosed in the macsy-model package.
The names are **case-sensitive**. All HMM profiles must be placed in the `profiles` directory of the macsy-model package.


.. note::
	For a detailed tutorial on how to define your macsy-model's features, parameters and HMM profiles, you can have a look at our cookbook in `this book chapter <https://link.springer.com/protocol/10.1007/978-1-4939-7033-9_1>`_ .




