.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2020  Institut Pasteur (Paris), and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _gene:

****
gene
****

The Gene object represents genes encoding the protein components of a System.
See :ref:`gene-implementation` for an overview of the implementation.


.. warning::
    To optimize computation and to avoid concurrency problems when we search several models,
    each gene must be instantiated only once, and stored in gene_bank.
    gene_bank is a :class:`macsypy.gene.GeneBank` object. 
    The gene_bank and model_bank (:class:`macsypy.model.ModelBank` object)
    are instantiated in :func:`macsypy.scripts.macsyfinder.main` function
    and filled by a definition_parser (:class:`macsypy.defintion_parser.DefinitionParser`)

Example to get a gene object: ::
  
    #get a model
    t2ss =  modelb_ank["T2SS"]
    
    #get of a gene
    pilO = gene_bank["pilO"]


Homolog/Analog are Composition with Gene.
Then a gene in some model is seen as a Gene, in some other models as an Homolg or Analog.
But there only one instance of this gene.::

    sctn = Gene(name="sctN", model_A)
    sctn_flg = Gene(name="sctN_FLG", model_A)
    sctn_homolog = Homolg(sctn, sctn_flg)
    sctn_flg.add_homolg(sctn_homolog)

which means that sctn_flg has an homolog which is sctn
sctn appear as a gene `sctn` and as homolog `sctn_homolg`


GeneBank
========

 .. automodule:: macsypy.gene
   :members: GeneBank
   :private-members:
   :special-members:

.. _gene_api:

Gene
====

.. automodule:: macsypy.gene
   :members: Gene
   :private-members:
   :special-members:

.. _homolog_api:

.. note::

    All attributes/methods which are not directly implemented in Homolog are redirected to that of the encapsulated Gene.

Homolog
=======

.. automodule:: macsypy.gene
   :members: Homolog
   :private-members:
   :special-members:


.. _analog_api:

.. note::

    All attributes/methods which are not directly implemented in Analog are redirected to that of the encapsulated Gene.

Analog
======

.. automodule:: macsypy.gene
   :members: Analog
   :private-members:
   :special-members:
   
 
