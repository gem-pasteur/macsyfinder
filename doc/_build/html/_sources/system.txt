.. _system:

****************
secretion system
****************

represent a secretion system. 
it's defined by a definition file in xml stored in 
the xml syntax is

for instance::

    <system> 
        <gene name="sctJ_FLG" presence="mandatory">
           <homologs>
               <gene name="sctJ" />
               <gene name="pilO" aligned="true"/>
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
and 3 kind of genes lists in function of their presence:

* The genes which must be present in the genome to defined this system.
* The genes which can be presents but not in all cases.
* The genes which must not be present.
 
 
System API reference
====================

.. automodule:: txsscanlib.system
   :members: System
   :private-members:
   :special-members:
