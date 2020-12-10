.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.
    Authors: Sophie Abby, Bertrand Néron
    Copyright © 2014-2020 Institut Pasteur (Paris), and CNRS.
    See the COPYRIGHT file for details
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3).
    See the COPYING file for details.
    
.. _definition_parser:

*****************
definition_parser
*****************

The model definition parser object "DefinitionParser" instantiates Models and Genes objects from
XML model definitions (see :ref:`model_definition`).
The parsing consists in three phases.

Phase 1. 

   * For each model to parse
   
     * create the Model
     * add this Model to the model_bank
     * findall genes defined in this model what are the level in the model definition.
     * create the CoreGene (a Gene which is not bind to a model).
       For each gene name there is only one instance of CoreGene
     * add these CoreGene in the gene_bank
    
Phase 2.

   * For each model to search
   
     * For each Gene defined in this System:
     
         * link the gene to the model. Create a ModelGene by encapsulating CoreGene from the gene_bank
           It can exists at each run several ModelGene for one CoreGene
         * If a gene has exhangeables create them (an Exchangeable inherits from ModeleGene)
           and add them to the current ModelGene


For instance::

    Syst_1
    <system inter_gene_max_space="10">
        <gene name=”A” mandatory=”1” loner="1">
            <exchangeables>
                <gene name=”B”>
            </exchangeables>
        </gene>
    <system>
    
    Syst_2
    <system inter_gene_max_space="15">
        <gene name=”B” mandatory=”1”>
            <exchangeables>
                <gene name=”C”>
            </exchangeables>
        </gene>
    <system>
    
    Syst_3
    <system inter_gene_max_space="20">
        <gene name=”c” mandatory=”1” />
    <system>


With the example above:

* the CoreGene A, B, C will be created
* the ModelGene (Syst_1, A)  (Syst_1, B), (Syst_2, B), (Syst_2, C), (Syst_3, C)
* The ModeleGene (Syst_1, A), (Syst_2, B) and (Syst_3, C) are directly link to their respective Models
* and where (Syst_1, B) (Syst_2, C) are exchangeables and link respectively to (Syst_1, A) and (Syst_2, B)
* the ModelGene has attributes defined in the model where they appear
  (Syst_1, B) inter_gene_max_space="10"
  (Syst_2, B) inter_gene_max_space="15"

.. note::
    The only "full" Systems (*i.e.,* with all corresponding Genes created) are those to detect.

.. _defintion_parser_api:

DefinitionParser
================
.. automodule:: macsypy.definition_parser
   :members:
   :private-members:
   :special-members:


