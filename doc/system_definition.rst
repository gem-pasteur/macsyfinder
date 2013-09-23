.. _system_definition:


.. _system-definition-grammar-label:

******************
Systems definition
******************

TXSScan allows to fully customize the definition of macromolecular systems via a **XML grammar** that is described below. 
A system is defined in a dedicated file named after the system (*e.g.*, 'T1SS.xml' for T1SS, the Type 1 Secretion System) by a set of **components** (*i.e.* proteins, or protein-coding genes given the context) with different attributes and that are used for **content description**. Some components are specific to the system, and some are possibly from other systems. In the latter case, the full description of the gene with its attributes must be defined in the XML file of the original system. 
Features regarding **co-localization** parameters for system detection are also defined in this system-specific file.

Three distinct types of components can be used to model a given system content, and which corresponds to Gene objects, and the corresponding HMM protein profile. 

* **Mandatory** components represent essential components to be found to infer the System presence.
* **Allowed** components correspond to components that can be found in some systems occurrence, or fastly evolving components that are hard to detect with a single profile. 
* **Forbidden** components are components which presence is eliminatory for the System assessment. 

Here is the XML hierarchy:

* The element root is "system". 

  * The element "system" must have an attribute "inter_gene_max_space" which is an integer representing the maximal number of components without a match between two components with a match for a component profile.
  * The element "system" may have attributes:
  
     * "min_mandatory_genes_required" which is an integer representing the minimal number of mandatory genes required to infer the system presence.
     * "min_genes_required" which is an integer representing the minimal number of mandatory or allowed genes (whose corresponding proteins match a profile of the system) required to infer the system presence.
     
  * The system contains one or more element "gene".
  
* The element "gene". 

   * must have an attribute "name" which must match to a profile in the profile directory.
   * must have an attribute "presence" which can take 3 values "mandatory", "allowed", "forbidden".
   * may have an attribute "system_ref" which is a reference to the secretion system where the gene 
      comes from (this attribute is used for forbidden gene and homologs gene). 
      If system_ref is not specified that mean the gene is from the current system.
   * may have attribute "loner" which is a boolean. If a gene is loner that means this gene can be isolated on the genome ( *default false* ).
   * may have attribute "exchangeable" which is a boolean. If a gene is exchangeable that means this gene or one of the homologs can be found without
      impacts on the secretion system.
   * "aligned" which is a boolean (this attribute is used only for homologs).
   * "inter_gene_max_space" which is an integer. 
   * may have one element "homologs".
   
* The element "homologs" contains one or more element "gene".

Example of a system definition in XML: ::
  
  <system> 
    <gene name="gspD" presence="mandatory" exchangeable="1">
       <homologs>
           <gene name="sctJ" system_ref="T3SS"/>
       </homologs>
    </gene>
    <gene name="sctN_FLG" presence="mandatory" loner="1"/>
    <gene name="sctV_FLG" presence="mandatory"/>
    <gene name="flp" presence="allowed"/>
    <gene name="sctC" presence="forbidden" system_ref="T3SS"/>
  </system>

.. note::
  
* a gene is identified by its name.
* a gene can be defined only one time (in one system).
* the other occurences of this gene must be references ("system_ref" attribute).
* if a gene specify the attribute "system_ref" it means that this gene is define in this system.
* if a gene not specify the attribute "system_ref" it means that this gene is define in the current system.
    
