.. _system_definition:


.. _system-definition-grammar-label:

*****************
System definition
*****************

The secretion systems are defined in xml format. Each file describe **one** system.
The name of the file corresponds to the name of the secretion system.

* The element root is "system". 

  * The element "system" must have an attribute "inter_gene_max_space" which is an integer representing
  * The element "system" may have attributes:
  
     * "min_mandatory_genes_required" which is an integer representing
     * "min_genes_required" which is an integer representing
     
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

Example of secretion system defintion: ::
  
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
    