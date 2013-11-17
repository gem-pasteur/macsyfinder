.. _system:

*****************
The System object
*****************

It represents a macromolecular system to detect. 
It is defined by a definition file in XML stored in a dedicated location that can be specified *via* the configuration file, or the command-line (`-d` parameter). See :ref:`system-definition-grammar-label` for more details on the XML grammar. 
The XML syntax is
for instance::

    <system inter_gene_max_space="15"> 
        <gene name="sctC" presence="mandatory">
           <homologs>
               <gene name="gspD" system_ref="T2SS"/>
               <gene name="pilQ" system_ref="T4P"/>
               <gene name="rcpA" system_ref="Tad"/>
           </homologs>
        </gene>
        <gene name="sctJ" presence="mandatory"/>
        <gene name="sctN" presence="mandatory"/>
        <gene name="sctQ" presence="mandatory"/>
        <gene name="sctR" presence="mandatory"/>
        <gene name="sctS" presence="mandatory"/>
        <gene name="sctT" presence="mandatory"/>
        <gene name="sctU" presence="mandatory"/>
        <gene name="sctV" presence="mandatory"/>
        <gene name="flgB" presence="forbidden"/>
        <gene name="flgC" presence="forbidden"/>
        <gene name="fliE" presence="forbidden"/>
    </system>

 
An object system_parser is used to build a system object from its XML definition file.

A system is named after the file name of its XML definition.
A system has an attribute inter_gene_max_space which is an integer,
and three kind of components are listed in function of their presence in the system:

* The genes that must be present in the genome to defined this system ("mandatory").
* The genes that can be present, but do not have to be found in every case ("allowed").
* The genes that must not be present in the system ("forbiden").

.. note:: 
    
    a complete description of the secretion system modelling is available in the section :ref:`system_definition`

SystemBank API reference
========================
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
