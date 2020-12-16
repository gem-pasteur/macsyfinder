.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2020  Institut Pasteur (Paris), and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.

.. _hit:



***
hit
***

A Hit is created when `hmmsearch` find similarities between a profile and protein of the input dataset


.. figure:: ../../_static/gene_obj_interaction.*

    A diagram showing the interaction between CoreGene, ModelGene, Model, HIt, ValidHit interactions
    The diagram above represents the models, genes and hit generated from the definitions below.

    .. code-block::

        <model name="A" inter_gene_max_space="2">
            <gene name="abc" presence="mandatory"/>
            <gene name="def" presence="accessory"/>
        </model>

        <model name="B" inter_gene_max_space="5">
            <gene name="def" presence="mandatory"/>
                <exchangeables>
                    <gene name="abc"/>
                </exchangeables>
            <gene name="ghj" presence="accessory"
        </model>



.. _hit_api:

hit
===
.. automodule:: macsypy.hit
   :members:
   :private-members:
   :special-members:
