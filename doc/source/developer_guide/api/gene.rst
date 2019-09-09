.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014  Institut Pasteur, Paris.                           
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _gene:

********
Gene API
********

The Gene object represents genes encoding the protein components of a System. See :ref:`gene-implementation` for an overview of the implementation. 


.. warning::
    To optimize computation and to avoid concurrency problems when we search several systems,
    each gene must be instanciated only once, and stored in gene_bank.
    gene_bank is a :class:`macsypy.gene.GeneBank` object. 
    The gene_bank and system_bank (:class:`macsypy.model.ModelBank` object)
    are filled by a system_parser (:class:`macsypy.defintion_parser.DefinitionParser`)

Example to get a gene object: ::
  
    from macsypy.system import system_bank
    from macsypy.config import Config
    from macsypy.gene import gene_bank
      
      
    #get a system  
    t2ss =  system_bank["T2SS"]
    
    #get of a gene
    pilO = gene_bank["pilO"]

 
GeneBank API reference
=========================
 .. automodule:: macsypy.gene
   :members: GeneBank
   :private-members:
   :special-members:

.. note::

   Don't instanciate your own GeneFactory use the gene_bank at the top level of the module. ::
     
     from macsypy.gene import gene_bank

Gene API reference
==================

.. automodule:: macsypy.gene
   :members: Gene
   :private-members:
   :special-members:



.. note::

    All attributes/methods which are not directly implemented in Homolog are redirected to that of the encapsulated Gene.

.. _homolog-API:
  
Homolog API reference  
=====================

.. automodule:: macsypy.gene
   :members: Homolog
   :private-members:
   :special-members:


.. _analog-API:   

Analog API reference  
====================

.. automodule:: macsypy.gene
   :members: Analog
   :private-members:
   :special-members:
   
 
