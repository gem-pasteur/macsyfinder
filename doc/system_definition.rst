.. _system_definition:

*****************
system definition
*****************

The secretion system are defined in xml format. Each file describe **one** system.
The name of the file correspond to the name of the secretion system.

* The element root is "system". 
* The system contains one or more element "gene".
* The element "gene". 
   * must have an attribute "name" which must match to a profile in profile directory.
   * must have an attribute "presence" chich can take 3 values "mandatory", "allowed", "forbidden".
   * may have an attribute "system_ref" which is a reference to the secretion system where the gene 
     come from (this attribute is mainly used for forbidden gene and homologs gene). 
     If system_ref is not specified that mean the gene is from the current system.
   * may have attribute "loner" which is a boolean. If a gene is loner that means this gene can be isolated on the genome ( *default false* ).
   * may have attibute "exhangeable" which is a boolean. If a gene is exchangeable that means this gene or one of the homologs can be found without
     impacts on the secretion system.
   * may have one element "homologs".
* The element homologs contains one or more element "gene".
 

example of secretion system defintion: ::
  
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