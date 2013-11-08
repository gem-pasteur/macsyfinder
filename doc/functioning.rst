.. _functioning:

*******************
TXSScan functioning
*******************

Initially, TXSScan **searches for the components** of a system using Hmmer with the corresponding protein profiles. The following steps depend on whether the input dataset is ordered (complete or nearly complete genomes), unordered replicon (partial genomes) or unordered dataset (metagenomes) (see :ref:`database`). 
In case of **ordered datasets**, the hits of the previous analysis are used to build *clusters of co-localized genes* as defined in the XML files. These clusters are then scanned to check for the model specifications like minimal quorum of "Mandatory" or "Allowed" genes or the absence of "Forbidden" components. 
When the gene order is unknown the power of the analysis is more limited. In this case, and depending on the type of dataset, the presence of systems can be suggested only on the basis of the quorum of genes. The results are outputted in a tabular and graphical form (see :ref:`outputs`). 

.. toctree::
   :maxdepth: 2

   search_genes
   search_systems

