.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014  Institut Pasteur, Paris.                           
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _functioning:


MacSyFinder functioning
=======================

*******************************
Searching for "gene" components
*******************************

Initially, MacSyFinder **searches for the components** of a system using Hmmer with the corresponding protein profiles. This step, and the extraction of significant hits can be performed in parallel (`-w` command-line option). See :ref:`the search_genes API <search_genes>`. 

From the list of systems to detect, a non-redundant list of components (and corresponding HMM protein profiles) to search with Hmmer is built. For each system, the search is performed for:

    - mandatory components
    - accessory components
    - forbidden components
    - homologs of these three types of components in the case they are *"exchangeable"*  


Hmmer is run on the corresponding set of HMM profiles, and the hits are filtered according to criteria defined by the user (see :ref:`Hmmer options <hmmer-options>` and :ref:`HMMReport`). 


**************************
Assessing systems presence
**************************

The following steps depend on whether the input dataset is *ordered* (complete or nearly complete genome(s)), or *unordered*  (metagenomes, or unassembled genome) (see :ref:`input-dataset-label`). 
In the case of **ordered datasets**, the hits of the previous analysis are used to build *clusters of co-localized genes* as defined in the XML files. These clusters are then scanned to check for the model specifications like minimal quorum of "Mandatory" or "Accessory" genes or the absence of "Forbidden" components. 
When the gene order is unknown the power of the analysis is more limited. In this case, and depending on the type of dataset, the presence of systems can be suggested only on the basis of the quorum of genes. The results are outputted in a tabular and graphical form (see :ref:`outputs`). 


For *ordered* datasets: 
-----------------------

1. The search starts first with the formation of clusters of contiguous hits **(co-localization criterion)** for each replicon. Two hits are said contiguous if their genomic location is separated by less than D protein-encoding gene, D being the maximum of the parameter "inter_gene_max_space" from the two genes with hits (system-specific, of gene-specific parameter). 

2. Clusters are then scanned, and those containing only genes from a single system are kept for further analyses (step 4.), wether those with multiple systems represented are analysed with a disambiguation step (step 3.).

3. The disambiguation step consists in splitting clusters that contains genes from different systems into sub-clusters that contain genes from a single system. Valid sub-clusters are then analysed like other clusters (step 2.). In the complex cases where genes from a same system are scattered into the cluster, then corresponding sub-clusters will not be further analysed for system detection.

4. Valid clusters are used to fill system occurrences (:class:`macsypy.search_systems.SystemOccurence`). In this step, the **quorum** criteria for the system assessment are checked according to the system's definition. In the case a single locus fills a complete system occurrence, it is stored and reported in the output files (**"single-locus"** occurrence). Otherwise, if the cluster corresponds to a valid but incomplete system, it is stored for inclusion in a scattered system occurrence (**"multi-loci"** occurrence).

5. When all clusters, "loner" genes and "multi_system" genes were scanned for inclusion in system occurrences, a decision is made for every system occurrence regarding the **quorum rules** defined for the corresponding system. 

.. note::
   When the "multi_loci" option is turned on, a single "multi-loci" system is assessed per replicon, even if it could correspond to multiple scattered systems. Thus, the "single-locus" systems correspond to a more powerful mode of detection.

.. warning::
    Cases where systems are consecutive will be treated, and separate systems will be detected, but complex cases of detection, *i.e.* when systems' components are intermingled will not be considered.


For *unordered* datasets: 
-------------------------

1. The Hits are grouped by system. 
2. They are used to fill a single system occurrence (:class:`macsypy.search_systems.SystemOccurence`) per system type.

.. note::
    The "unordered" mode of detection is less powerful, as a single occurrence of a given system is filled for an entire dataset with hits that origin is unknown. Please consider "systems assessments" with caution in this mode. 



