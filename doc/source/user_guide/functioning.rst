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


************************************
A. Searching for systems' components
************************************
MacSyFinder is run from the command-line using a variety of input files and options.
See :ref:`input-dataset-label` for more details.

Initially, MacSyFinder **searches for the components** of a system by sequence similarity search. 

From the list of systems to detect, a non-redundant list of components to search is built.
For each system, the list includes:

    - mandatory components
    - accessory components
    - neutral components
    - forbidden components
    - exchangeables these components can be functionally replaced by other components (usually by analogs or homologs)

Hmmer is run on the corresponding set of HMM profiles, and the hits are filtered according to criteria defined
by the user (see :ref:`Hmmer options <hmmer-options>` and :ref:`HMMReport`).
This step, and the extraction of significant hits can be performed in parallel (`-w` command-line option).
See the :ref:`command-line-label`, and the :ref:`search_genes API <search_genes>` for more details.

.. _system_assessment:

*****************************
B. Assessing systems presence
*****************************

The following steps depend on whether the input dataset is *ordered* (complete or nearly complete genome(s)),
or *unordered*  (metagenomes, or unassembled genome) (see :ref:`input-dataset-label`).
In the case of **ordered datasets**, the hits are filtered to keep only hits related to the model we looking for.
The these hits are used to build *clusters of co-localized genes* as defined in the Model files.
These clusters are then scanned to check for the model specifications like minimal quorum of
"Mandatory" or "Accessory" genes or the absence of "Forbidden" components.
When the gene order is unknown the power of the analysis is more limited.
In this case, and depending on the type of dataset, the presence of systems can be suggested only on the basis of
the quorum of genes. The results are outputted in a tabular and graphical form (see :ref:`outputs`).

.. _note:
    The `neutral` components are used to build clusters of co-localized genes.
    But does not play any role in quorum assessment.


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
   The clusters which fill the quorum are reported the the `systems.txt` and `systems.tsv` files see :ref:`outputs`.
   The clusters which not fulfill the quorum are reported in the `rejected_clusters.txt` file.

For *unordered* datasets: 
-------------------------

1. The Hits are filtered by model.
2. They are used to check if they fill the quorum (in other words the clustering step is skipped).

.. note::
    The "unordered" mode of detection is less powerful, as a single occurrence of a given model is filled for
    an entire dataset with hits that origin is unknown. Please consider "systems assessments" with caution in this mode.



