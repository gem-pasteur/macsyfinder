.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2020  Institut Pasteur (Paris), and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _gene_module:

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


Exchangeable is a Composition with Gene.
Then a gene in some model is seen as a Gene, in some other models as an Exchangeable.
But there only one instance of this gene.::

    sctn = Gene(name="sctN", model_A)
    sctn_flg = Gene(name="sctN_FLG", model_A)
    sctn_ex = Exchangeable(sctn, sctn_flg)

which means that sctn_flg can replaced by sctn
sctn appear as a gene `sctn` and as exchangeable `sctn_ex`


GeneBank
========

 .. autoclass:: macsypy.gene.GeneBank
   :members:
   :private-members:
   :special-members:


Gene
====

There is two classes to modelize a gene: :class:`macsypy.gene.CoreGene` and :class:`macsypy.gene.ModelGene`.
The CoreGene are created using the :class:`macsypy.gene.GeneBank` factory and there is only one instance
of a CoreGene with a given name.
Whereas several ModelGene with the same name can appear in different model and can have differents properties,
`loner` in one model and not in an other, have different `inter_gene_max_space` ...
The ModelGene is attached to the model and is composed of a CoreGene.

.. note::
    The :class:`macsypy.hit.Hit` object are link to a CoreGene, whereas the :class:`macsypy.hit.ValidHit` `ref_gene`
    attribute reference a :class:`macsypy.gene.ModelGene`


.. _core_gene_api:

CoreGene
========

.. autoclass:: macsypy.gene.CoreGene
   :members:
   :private-members:
   :special-members:

.. _model_gene_api:

ModelGene
=========

.. autoclass:: macsypy.gene.ModelGene
   :members:
   :private-members:
   :special-members:

.. _exchangeable_api:

.. note::

    All attributes/methods which are not directly implemented in Exchangeable are redirected
    to that of the encapsulated ModelGene.

Exchangeable
============

.. autoclass:: macsypy.gene.Exchangeable
   :members:
   :private-members:
   :special-members:


.. _gene_status_api:

GeneStatus
==========

.. autoclass:: macsypy.gene.GeneStatus
   :members:
   :private-members:
   :special-members:
