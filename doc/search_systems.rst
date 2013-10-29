.. _search_systems:

******************
Search for systems
******************
The search for systems starts with the extraction of the list of components to be detected with HMM protein profiles. Then Hmmer is run on this set of profiles, and hits are filtered. 

For *ordered* datasets: 

1. The search starts first with the formation of clusters of contiguous hits **(co-localization criterion)**. Two hits are said contiguous if their genomic location is separated by less than D protein-encoding gene, D being the maximum of the parameter "inter_gene_max_space" from the two genes with hits (system-specific, of gene-specific parameter). 
2. Clusters are then scanned, and those containing only genes from a single system are kept for further analyses (step 4.), wether those with multiple systems represented are analysed with a disambiguation step (step 3.).
3. The disambiguation step consists in ....
4. Valid clusters are used to fill system occurrences.... 
5. When all clusters or loner genes were scanned for inclusion in system occurrences, a decision is made for every system occurrence regarding the **quorum rules** defined for the corresponding system. 

 
search_systems API reference
============================

.. automodule:: txsscanlib.search_systems
   :members: 
   :private-members:
   :special-members:
  
