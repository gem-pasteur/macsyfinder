.. _search_systems:

**************************
Assessing systems presence
**************************
The search for systems starts with the extraction of the list of components to be detected with HMM protein profiles. Then Hmmer is run on this set of profiles, and hits are filtered. 

For *ordered* datasets: 

1. The search starts first with the formation of clusters of contiguous hits **(co-localization criterion)**. Two hits are said contiguous if their genomic location is separated by less than D protein-encoding gene, D being the maximum of the parameter "inter_gene_max_space" from the two genes with hits (system-specific, of gene-specific parameter). 
2. Clusters are then scanned, and those containing only genes from a single system are kept for further analyses (step 4.), wether those with multiple systems represented are analysed with a disambiguation step (step 3.).
3. The disambiguation step consists in splitting clusters that contains genes from different systems into sub-clusters that contain genes from a single system. Valid sub-clusters are then analysed like other clusters (step 2.). In the complex cases where genes from a same system are scattered into the cluster, then corresponding sub-clusters will not be further analysed for system detection.
4. Valid clusters are used to fill system occurrences. In this step, the **quorum** criteria for the system assessment are checked according to the system's definition. In the case a system occurrence corresponds to a complete system, it is stored and reported in the output files. Otherwise, if the cluster corresponds to a valid but incomplete system, it is stored for inclusion in a scattered system occurrence.
5. When all clusters, "loner" genes and "multi_system" genes were scanned for inclusion in system occurrences, a decision is made for every system occurrence regarding the **quorum rules** defined for the corresponding system. 

 
search_systems API reference
============================

.. automodule:: txsscanlib.search_systems
   :members: 
   :private-members:
   :special-members:
  
