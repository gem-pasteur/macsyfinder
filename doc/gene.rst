.. _gene:

****
gene
****

the Gene class represents the genes defined by a System
each gene must have a profile (the profiles must be stored in ... ). 
a gene may have homologs.

 homolog have to properties :

* a reference gene, which is a reference toward the gene which is the homolog
* and if the profile is aligned on the gene reference or not.

and is composed of a gene. All attribute/methods which are not directly implemented in 
Homolog are redirected to this gene.

To optimize computation and to avoid concurrency problems when we search several systems.
Each gene must be instanciated only one time,and stored in gene_bank.
gene_bank is a  :class:`txsscanlib.gene.GeneBank` object. 
the gene_bank and system bank are filled by a system_parser ( :class:`txsscanlib.system_parser.SystemParser` )

example to get a gene object: ::
  
    from txsscanlib.system import system_bank
    from txsscanlib.config import Config
    from txsscanlib.gene import gene_bank
      
      
    #get a system  
    t2ss =  system_bank["T2SS"]
    
    #get of a gene
    pilO = gene_bank["pilO"]
 
GeneBank API reference
=========================
 .. automodule:: txsscanlib.gene
   :members: GeneBank
   :private-members:
   :special-members:

.. note::

   Don't instanciate your own GeneFactory use the gene_bank at the top level of the module. ::
     
     from txsscanlib.gene import gene_bank

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
   
   