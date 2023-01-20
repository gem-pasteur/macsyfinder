.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2023 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _models_v1_v2:

=============================
Carrying models from v1 to v2 
=============================

Models from v1 are not compatible straight away with v2. 
For those who had designed MacSyFinder's models for Version 1 and would like to carry them for Version 2, here are the changes to consider:  

- the keyword "system" was changed:
  `<system>` ::arrow:: `<model>`
- the keyword `<system_ref>` was removed.
  For a given systems' package, each gene has to be defined only once in a macsy-model.
  There is no need anymore to reference which model it is from, when used as a component in another system's model.
- now the version of the macsy-models' type has to be documented as a feature of the "model" keyword, like this: `vers = "2.0"`
- the following keywords have been replaced (but see :ref:`below<ex2>` for more details):

  * homologs => exchangeables
  * analogs => exchangeables

.. note::
 
 "exchangeable" is not a feature anymore, but is replaced by the keyword "exchangeables". 


.. note::
 
 These changes in the grammar used to specify model is also accompanied by a change on how to organize folders with models and profiles.
 In particular, the new file architecture enables an :ref:`easier shipping<publish_package>` of the developed macsy-models. See :ref:`here<model_package>` for more details.
  

Here follow some examples of updates from v1 to v2.


1. A very simple model.
-----------------------

`T1SS.xml` under **v1**::

	<system inter_gene_max_space="5" min_mandatory_genes_required="3" min_genes_required="3">
	    <gene name="T1SS_abc" presence="mandatory"/>
	    <gene name="T1SS_mfp" presence="mandatory"/>
	    <gene name="T1SS_omf" presence="mandatory" loner="1" multi_system="1"/>
	</system>


`T1SS.xml` under **v2**::

	<model inter_gene_max_space="5" min_mandatory_genes_required="3" min_genes_required="3" vers = "2.0">
	    <gene name="T1SS_abc" presence="mandatory"/>
	    <gene name="T1SS_mfp" presence="mandatory"/>
	    <gene name="T1SS_omf" presence="mandatory" loner="1" multi_system="1"/>
	</model>


.. note::

	In a nutshell, the minimal changes from v1 to v2 for a simple macsy-model listing components are the following:
	
	  - <system> => <model>
	  - `vers = "2.0"`

.. _ex2:

2. A model with homologs.
-------------------------

`Tad.xml` under **v1**::

	<system inter_gene_max_space="5" min_mandatory_genes_required="4" min_genes_required="6" multi_loci="0">
	    <gene name="Tad_rcpA" presence="mandatory">    	
	    	<homologs>
	        	<gene name="T2SS_gspD" system_ref="T2SS"/>
	    	        <gene name="T4P_pilQ" system_ref="T4P"/>
		        <gene name="T3SS_sctC" system_ref="T3SS"/>
		</homologs> 
	    </gene>
	    <gene name="Tad_tadA" presence="mandatory"/>
	    <gene name="Tad_tadB" presence="mandatory"/>
	    <gene name="Tad_tadC" presence="mandatory"/>
	    <gene name="Tad_tadV" presence="mandatory"/>
	    <gene name="Tad_tadZ" presence="mandatory"/>
	    <gene name="Tad_flp" presence="accessory"/>
	    <gene name="Tad_tadE" presence="accessory"/>
	    <gene name="Tad_tadF" presence="accessory"/>
	</system>


`Tad.xml` under **v2**::

	<model inter_gene_max_space="5" min_mandatory_genes_required="4" min_genes_required="6" multi_loci="0" vers="2.0">

	    <gene name="Tad_rcpA" presence="mandatory"/>    	
	    <gene name="Tad_tadA" presence="mandatory"/>
	    <gene name="Tad_tadB" presence="mandatory"/>
	    <gene name="Tad_tadC" presence="mandatory"/>
	    <gene name="Tad_tadV" presence="mandatory"/>
	    <gene name="Tad_tadZ" presence="mandatory"/>
	    <gene name="Tad_flp" presence="accessory"/>
	    <gene name="Tad_tadE" presence="accessory"/>
	    <gene name="Tad_tadF" presence="accessory"/>

	</model>

.. note::

	The `homologs` and `analogs` keyword having disappeared, it is not necessary anymore to list homologous components (e.g., those that may match several HMM profiles during the sequence similarity search), unless they are `exchangeables`. 
	
	

3. A model with exchangeable homologs.
--------------------------------------

`T3SS.xml` under **v1**::

	<system inter_gene_max_space="10" min_mandatory_genes_required="7" min_genes_required="7" multi_loci="1">
	    <gene name="T3SS_sctC" presence="mandatory" exchangeable="1">        
	        <homologs>
	    	    <gene name="T2SS_gspD" system_ref="T2SS"/>
	    	    <gene name="T4P_pilQ" system_ref="T4P"/>
	    	    <gene name="Tad_rcpA" system_ref="Tad"/>
	        </homologs>
	    </gene>
	    <gene name="T3SS_sctJ" presence="mandatory">       
	        <homologs>
	    	    <gene name="Flg_sctJ_FLG" system_ref="Flagellum"/>
	        </homologs>
	    </gene>
	    <gene name="T3SS_sctN" presence="mandatory">       
	    	<homologs>
	    	    <gene name="Flg_sctN_FLG" system_ref="Flagellum"/>
	    	</homologs>
	    </gene>
	    <gene name="T3SS_sctQ" presence="mandatory">  
	        <homologs>
	    	    <gene name="Flg_sctQ_FLG" system_ref="Flagellum"/>
	    	</homologs>
	    </gene>
	    <gene name="T3SS_sctR" presence="mandatory">    
	    	<homologs>
	            <gene name="Flg_sctR_FLG" system_ref="Flagellum"/>
	        </homologs>
	    </gene>
	    <gene name="T3SS_sctS" presence="mandatory">    
    		<homologs>
	            <gene name="Flg_sctS_FLG" system_ref="Flagellum"/>
	        </homologs>
	    </gene>
	    <gene name="T3SS_sctT" presence="mandatory">    
	        <homologs>
	    	    <gene name="Flg_sctT_FLG" system_ref="Flagellum"/>
	        </homologs>
	    </gene>
	    <gene name="T3SS_sctU" presence="mandatory">    
	        <homologs>
	            <gene name="Flg_sctU_FLG" system_ref="Flagellum"/>
	    	</homologs>
	    </gene>
	    <gene name="T3SS_sctV" presence="mandatory">    
	    	<homologs>
	            <gene name="Flg_sctV_FLG" system_ref="Flagellum"/>
	        </homologs>
	    </gene>
	    <gene name="Flg_fliE" presence="forbidden" system_ref="Flagellum"/>
	    <gene name="Flg_flgB" presence="forbidden" system_ref="Flagellum"/>
	    <gene name="Flg_flgC" presence="forbidden" system_ref="Flagellum"/>
	</system>



`T3SS.xml` under **v2**::

	<model inter_gene_max_space="10" min_mandatory_genes_required="7" min_genes_required="7" multi_loci="1" vers="2.0">
	    <gene name="T3SS_sctC" presence="mandatory">        
	        <exchangeables>
	    	    <gene name="T2SS_gspD"/>
	    	    <gene name="T4P_pilQ"/>
	    	    <gene name="Tad_rcpA"/>
	        </exchangeables>
	    </gene>
	    <gene name="T3SS_sctJ" presence="mandatory"/>       
	    <gene name="T3SS_sctN" presence="mandatory"/>       
	    <gene name="T3SS_sctQ" presence="mandatory"/>  
	    <gene name="T3SS_sctR" presence="mandatory"/>    
	    <gene name="T3SS_sctS" presence="mandatory"/>    
	    <gene name="T3SS_sctT" presence="mandatory"/>    
	    <gene name="T3SS_sctU" presence="mandatory"/>    
	    <gene name="T3SS_sctV" presence="mandatory"/>    
	    <gene name="Flg_fliE" presence="forbidden"/>
	    <gene name="Flg_flgB" presence="forbidden"/>
	    <gene name="Flg_flgC" presence="forbidden"/>
	</model>

.. note::

	- As only the secretin component 'T3SS_sctC' was exchangeable in its role within T3SS with its homologs T2SS_gspD, T4P_pilQ and Tad_rcpA, these three components are now set as `exchangeables` (they can functionally *replace* the component 'T3SS_sctC'), and all other `homologs` do not need to be listed anymore.  
	- The keyword `system_ref` is not needed anymore. Therefore, the **v2** definition of T3SS is way more compact than that for **v1**.
