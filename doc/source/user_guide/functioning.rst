.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2020 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _functioning:


MacSyFinder functioning
=======================

********************
Functioning overview
********************

    .. image:: ../_static/MSF_functioning.svg
     :height: 500px
     :align: left

.. A. MacSyFinder is ran from the command-line using a variety of input files and options.
   See :ref:`input-dataset-label` for more details.

.. B. Depending on the input dataset type ("ordered" or "unordered"),
   the hits detected are processed using their contiguity or not.
   More details are provided in the :ref:`section below<system_assessment>`


MacSyFinder is run from the command-line using a variety of input files and options.
See :ref:`input-dataset-label` for more details. Below follows a description of its overall functioning. 


************************************
A. Searching for systems' components
************************************

Initially, MacSyFinder **searches for the components** of the system(s) to detect by sequence similarity search.

From the list of system(s) to detect, a **non-redundant list of components to search** is built.
For each system, the list can include:

    - mandatory components
    - accessory components
    - neutral components
    - forbidden components
    - exchangeable components that can be functionally replaced by other components (usually by analogs or homologs). These other components are thus also added to the list of components to search.

HMMER is run on the corresponding set of components' HMM profiles, and the hits are filtered according to the criteria defined
by the user or by default (see :ref:`Hmmer options <hmmer-options>` and :ref:`HMMReport`).
This step, and the extraction of significant hits can be performed in parallel (`-w` command-line option).
See the :ref:`command-line-label`, and the :ref:`search_genes API <search_genes>` for more details.

.. _system_assessment:

******************************
B. Assessing systems' presence
******************************

The following steps depend on whether the input dataset is *ordered* (complete or nearly complete genome(s)),
or *unordered*  (metagenomes, or unassembled genome) (see :ref:`input-dataset-label`).

In the case of **ordered datasets**, the hits are filtered to keep only hits related to the system's model we are looking for.
These hits are used to build *clusters of co-localized genes* as defined in the macsy-model files.
These clusters are then screened to check for the model specifications such as the minimal quorum of
"Mandatory" or "Accessory" genes, or the absence of "Forbidden" components.

When the **gene order is unknown** the power of the analysis is more **limited**.
In this case, the presence of systems can only be suggested on the basis of
the quorum of genes - and not based on genomic context information. 

.. _note:
    The `neutral` components are used to build clusters of co-localized genes.
    They do not play any role in components' quorum assessment.


For *ordered* datasets:
-----------------------

1. The search starts first with the hits filtering to keep all hits related to a model (mandatory, accessory, neutral,
   forbidden, exchangeable)

2.  We looking for the formation of clusters of contiguous hits
    **(co-localization criterion)** for each replicon.
    Two hits are said contiguous if their genomic location is separated by less than D protein-encoding gene D
    being the maximum of the parameter "inter_gene_max_space"
    from the two genes with hits (system-specific, of gene-specific parameter).
    The `loner` components may form a cluster on their own.

3. Clusters are used to fill system according the quorum (:func:`macsypy.system.match`).
   In this step, the **quorum** criteria for the system assessment are checked according to the model's definition:
   min_genes_required, min_mandatory_genes_required.
   In the case of single loci (default) each clusters + loners are evaluated for quorum separately.
   In the case of multi loci (``multi_loci=True``) each clusters and all clusters combination are evaluated for the quorum.
   The clusters which fill the quorum are reported the the `all_systems.txt` and `all_systems.tsv` files see :ref:`outputs`.
   The clusters which not fulfill the quorum are reported in the `rejected_clusters.txt` file.
   
Then the :ref:`next step of the search <combinatorial-exploration>` is performed to compute global solutions for the replicon(s) analysed (sets of compatible systems). 

For *unordered* datasets: 
-------------------------

1. The Hits are filtered by model.
2. They are used to check if they fill the quorum (in other words the clustering step is skipped).

.. note::
    The "unordered" mode of detection is less powerful, as a single occurrence of a given model is filled for
    an entire dataset with hits that origin is unknown. Please consider the assessment of systems with caution in this mode.

For unordered datasets, the search so ends, and MacSyFinder generates the final :ref:`output files <outputs>`. 


.. _combinatorial-exploration:

**************************************************************
C. Computing possible solutions, defining the best one(s)
**************************************************************

This step only applies to the most powerful search mode, i.e., on **ordered datasets**. ``NEW in V2``

The **new search engine** implemented since version 2 of MacSyFinder better explores the space of possible solutions regarding the presence of systems in replicons analysed. 
It creates clusters of hits for systems' components separately for each system searched, and therefore might find candidate occurrences of systems that overlap. 
Moreover, if a system is possibly encoded at several locations on the replicon analysed (option multi_loci set to "True" in the model), this calls for a combinatorial analysis of the different clusters to assemble them into coherent systems regarding the macsy-models.
We therefore introduced a **scoring scheme for candidate systems**, to easily separate combinations of clusters that are readily more similar to a system's model than others.  

The assumptions behind this scoring scheme are the following:
	* we set a score for the different types of genes/components:
		- +X is added when a mandatory gene is present 
		- +X is added when a accessory gene is present 
		
	* when combinations of clusters are explored in order to fulfill macsy-models' requirements ("multi_loci" mode, several clusters can make up a complete System), we want to **favor concise sets of clusters** to fulfill a System's model. We thus **penalize the adjunction of a cluster** to a candidate System when this cluster does not bring any new components to the System's quorum, or when it brings **redundant components**. Thus:
		- -X is added when a redundant mandatory gene is added when adjuncting the cluster to a candidate System
		- -X is added when a redundant accessory gene is added when adjuncting the cluster to a candidate System

	* overall, only candidate sets of clusters that fulfill a macsy-model and that are thus designated candidate Systems, obtain a **System's score**

This search for candidate Systems results in a number of possible Solutions representing combinations of putative sets of Systems in the analysed dataset. 
We define a Solution as being a set of compatible Systems, since we do not allow to have overalps between inferred Systems, i.e. components part of several Systems.
All possible Solutions are combinatorially explored and consist in all possible sets of compatible Systems. 

A scoring system also enables to separate between sets of Solution. It is basically the **sum of the Systems' scores**.  
The overall procedure of exploring the space of possible Solutions while finding the optimal one, i.e. that with the maximal score, is performed at once using a graph solution to this problem, implemented in the networkx package. 
This allows to provide the user with one, or multiple Solutions that have the best score possible among all combinations of compatible Systems. 

******************************
C. Selecting the best solution
******************************

1. At the end of the previous step MacSyFinder has computed all potential systems present in the replicon.
   But all this systems does not exists at the same time.
   Sometimes for a given model MacSyFinder found several potential systems. But these systems share lot of components.
   So only ones of these systems is really present. So to choose the most system probable system, we compute for each one
   a score, based on the model wholeness, the number of loci, if each function is coded by a gene or it's exchangeable, ...


2. So Model also share some components for instance let's consider a Model B with a gene G4 as in Model A
   and the cluster C5 contains the same hit H4 than in cluster C2.
   So we had to choose to attribute the hit H4 to systems A or B.
   To do that we consider all systems combination, the systems which share components are considering incompatible.
   For instance SA_1 and SB_1 share the Hit H4 (respectively in cluster C2 and C5) so these two systems cannot exists together
   Then we choose the largest systems combination which maximize the score.
   So in our example the system SA_2 and SB_3
   The results of this step are reported in `best_systems.tsv` file see :ref:`outputs`.

