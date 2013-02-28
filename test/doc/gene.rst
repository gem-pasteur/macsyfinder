.. _gene:

****
gene
****

the Gene class represents the genes defined by a System
each gene must have a profile (the profiles nust be stored in ... ). 
a gene may have homologs.

 homolog have to properties :

* a reference gene, which is a reference toward the gene which is the homolog
* and if the profile is aligned on the gene reference or not.

and is composed of a gene. All attribute/methods which are not directly implemented in 
Homolog are redireted to this gene.

To optimize computation and to avoid concurrency problems when we search several systems.
Each gene must be instanciate only one time. This is ensure by using the gene_factory.
gene_factory is a  :class:`txsscanlib.gene.GeneFactory` object. 

example to get a gene object: ::
  
    from txsscanlib.system import system_factory
    from txsscanlib.config import Config
    from txsscanlib.gene import gene_factory
      
    config = Config()
      
    #instanciation of a system  
    t2ss =  system_factory.get_system("T2SS", config)
    
    #instanciation of a gene
    pilO = gene_factory.get_gene("pilO", t2ss, config)
 
GeneFactory API reference
=========================
 .. automodule:: txsscanlib.gene
   :members: GeneFactory
   :private-members:
   :special-members:

.. note::

   Don't instanciate your own GeneFactory use the gene_factory at the top level of the module. ::
     
     from txsscanlib.gene import gene_factory

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
   
   