.. _system:

****************
secretion system
****************

Represents a secretion system. 
It is defined by a definition file in xml stored in 
The xml syntax is

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


a system_parser is used to build a system object from it's xml definition

a system as a name based on the file name of the xml definition
a system as an attribute inter_gene_max_space which is an integer
and 3 kind of genes lists in function of their presence:

* The genes which must be present in the genome to defined this system.
* The genes which can be presents but not in all cases.
* The genes which must not be present.

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
