.. _system_parser:

*************
system_parser
*************

the system parser create system and genes from the xml system definition.
the parsing in 3 phases

for instance::

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


1. phase: 

   * each gene must be parse from the system where it is defined
   * so from the list of system to search establish the list of system to parse

2. phase:

   * for each system to parse 
   
     * create the system
     * add this system in the system_bank
     * create of the genes defined in this system with its attribute but not its homologs
     * add these genes in the gene_bank
    
3. phase: 

   * for each system to search
   
     * for each gene defined in this system:
     
         * create the homologs getting the gene from the gene_bank
         * add the gene in the system

at the end with the example above:

* the Syst_1 have a gene_A 
* the gene_A have homolog gene_B
* the gene_B have a reference to Syst_2
* the gene_B attributes from the Syst_2 are use to build the gene
* the Syst_2 have attributes as define in xml (inter_gene_max_space ,...)

countrariwise: 

* the gene_B have no homologs
* the Syst_2 have not gene

only the systems in the list to the systems to search are full.

SystemParser API reference
==========================
.. automodule:: txsscanlib.system_parser
   :members:
   :private-members:
   :special-members:


