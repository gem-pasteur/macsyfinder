.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014  Institut Pasteur, Paris.                           
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _system_parser:

*********************************
The Parser of Systems definitions
*********************************

The system parser object "SystemParser" instanciates Systems and Genes objects from XML system definitions (see :ref:`system_definition`).
The parsing consists in three phases.

Phase 1. 

   * each Gene is parsed from the System it is defined
   * From the list of System to detect, the list of Systems to parse is established

Phase 2.

   * For each system to parse 
   
     * create the System
     * add this System to the system_bank
     * create the Genes defined in this System with their attributes but not their Homologs
     * add these Genes in the gene_bank
    
Phase 3. 

   * For each System to search
   
     * For each Gene defined in this System:
     
         * create the Homologs by encapsulating Genes from the gene_bank
         * add the Gene to the System


For instance::

    Syst_1
    <system inter_gene_max_space="10">
        <gene name=”A” mandatory=”1” loner="1">
            <homologs>
                <gene name=”B” sys_ref=”Syst_2”>
            </homologs>
        </gene>
    <system>
    
    Syst_2
    <system inter_gene_max_space="15">
        <gene name=”B” mandatory=”1”>
            <homologs>
                <gene name=”B” sys_ref=”Syst_1”
                <gene name=”C” sys_ref=”Syst_3”>
            </homologs>
        </gene>
    <system>
    
    Syst_3
    <system inter_gene_max_space="20">
        <gene name=”c” mandatory=”1” />
    <system>


With the example above:

* the Syst_1 has a gene_A 
* the gene_A has homolog gene_B
* the gene_B has a reference to Syst_2
* gene_B attributes from the Syst_2 are used to build the Gene
* the Syst_2 has attributes as defined in the corresponding XML file (inter_gene_max_space ,...)

Contrariwise: 

* the gene_B has no Homologs
* the Syst_2 has no Genes

.. note::
    The only "full" Systems (*i.e.,* with all corresponding Genes created) are those to detect.


SystemParser API reference
==========================
.. automodule:: macsypy.system_parser
   :members:
   :private-members:
   :special-members:


