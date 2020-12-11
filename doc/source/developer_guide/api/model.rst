.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2020  Institut Pasteur (Paris), and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _model:


*****
model
*****

The model is a formal representation of system.
The model is describe in terms of components.
There are 3 component classes:

    * genes which are mandatory
    * genes which are accessory
    * genes which are forbiden

each genes can have homolgs or analogs.
Analogs and Homolgs refer to genes described in an other model.

The models describe also distance constraints between genes:

    * inter_gene_max_space
    * loner
    * multi_loci

and quorum constraints

    * min_mandatory_genes_required
    * min_genes_required


.. _model_bank_api:

ModelBank
=========
 .. autoclass:: macsypy.model.ModelBank
   :members:
   :private-members:
   :special-members:


.. _model_api:

Model
=====

.. autoclass:: macsypy.model.Model
   :members:
   :private-members:
   :special-members:
