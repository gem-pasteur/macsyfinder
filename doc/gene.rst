.. _gene:

********
Gene API
********

The Gene object represents genes encoding the protein components of a System. See :ref:`gene-implementation` for an overview of the implementation. 


.. warning::
    To optimize computation and to avoid concurrency problems when we search several systems, each gene must be instanciated only once, and stored in gene_bank.
    gene_bank is a :class:`txsscanlib.gene.GeneBank` object. 
    The gene_bank and system bank are filled by a system_parser (:class:`txsscanlib.system_parser.SystemParser`)

Example to get a gene object: ::
  
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



.. note::

    All attributes/methods which are not directly implemented in Homolog are redirected to that of the encapsulated Gene.

.. _homolog-API:
  
Homolog API reference  
=====================

.. automodule:: txsscanlib.gene
   :members: Homolog
   :private-members:
   :special-members:
   
   
