.. _system:

*****************
The system object
*****************

It represents a macromolecular system to detect. 
It is defined by a definition file in XML stored in 
The XML syntax is

for instance::

    <system inter_gene_max_space="20"> 
        <gene name="sctJ_FLG" presence="mandatory">
           <homologs>
               <gene name="sctJ" system_ref="T2SS"/>
               <gene name="pilO" system_ref="T2SS" aligned="true"/>
           </homologs>
        </gene>
        <gene name="sctN_FLG" presence="mandatory"/>
        <gene name="sctQ_FLG" presence="mandatory"/>
        <gene name="sctR_FLG" presence="mandatory"/>
        <gene name="sctS_FLG" presence="mandatory"/>
        <gene name="sctT_FLG" presence="mandatory"/>
        <gene name="sctU_FLG" presence="mandatory"/>
        <gene name="sctV_FLG" presence="mandatory"/>
        <gene name="flgB" presence="allowed"/>
        <gene name="flgC" presence="allowed"/>
        <gene name="fliE" presence="allowed"/>
        <gene name="sctC" presence="forbidden"/>
    </system>

 
An object system_parser is used to build a system object from its XML definition file.

A system is named after the file name of its XML definition.
A system has an attribute inter_gene_max_space which is an integer,
and three kind of components are listed in function of their presence in the system:

* The genes that must be present in the genome to defined this system ("mandatory").
* The genes that can be present, but do not have to be found in every cases ("allowed").
* The genes that must not be present in the system ("forbiden").

.. note:: 
    
    a complete description of the secretion system grammar is available here :ref:`system-definition-grammar-label`

SystemBank API reference
===========================
 .. automodule:: txsscanlib.system
   :members: SystemBank
   :private-members:
   :special-members:

 
.. note::

   Don't instanciate your own SystemBank use the system_bank at the top level of the module. ::
     
     from txsscanlib.system import system_bank
 
System API reference
====================

.. automodule:: txsscanlib.system
   :members: System
   :private-members:
   :special-members:
