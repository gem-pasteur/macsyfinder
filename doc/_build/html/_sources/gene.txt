.. _gene:

****
gene
****

the Gene class represents the genes defined by a System
each gene must have a profile (the profiles nust be stored in ... ). 
a gene may have homologs. an homologs is a Gene with 2 additional properties:

* a reference gene, which is a reference toward the gene which is the homolog
* and if the profile is aligned on the gene reference or not.

 
 
Gene API reference
==================

.. automodule:: txsscanlib.gene
   :members: Gene
   :private-members:
   :special-members:
  
Homolog API reference  
=====================

.. automodule:: txsscanlib.gene
   :members: Homolog
   :private-members:
   :special-members: