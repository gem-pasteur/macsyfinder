.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2021  Institut Pasteur (Paris),and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _gallery_models:

*******************************************
Gallery of examples of MacSyFinder's models 
*******************************************

.. contents:: Table of contents of gallery
	:local: 
        :depth: 1 


Here follows a "gallery" of MacSyFinder models we have developed over the years, attempting to describe the reasoning behind the modeling process. 

These examples are extracted from published work, see the following references (they include more examples):

- Abby et al. 2016, *Scientific Reports* for the description of T1SS and T3SS models
- Denise et al. 2019, *Plos Biology* for the description of T2SS and type IV-filament super-family models
- Abby et al. 2014, *Plos ONE* and Bernheim et al. 2018, for the description of the Cas systems models



.. _T1SS:

Getting started with a (not-so-)simple example: modelling the T1SS
==================================================================


1. Identifying genetic components
---------------------------------

The type I secretion system consists in three conserved components: 

- an ABC transporter (ABC)
- a membrane-fusion protein (MFP)
- an outer membrane protein (OMF)

For their detection, we therefore need to provide HMM profiles for each component, for example: "abc.hmm", "mfp.hmm" and "omf.hmm". 
These can be specifically designed, or taken from HMM profiles databanks such as PFAM or TIGRFAM. 


2. Determining the role of the components
-----------------------------------------

From litterature, the three components listed above *must* be present to have a viable T1SS. Therefore, these are all deemed *mandatory* in the model of the T1SS. 


3. Describing their genetic architecture
----------------------------------------

According to the litteraure, the genes encoding the three components listed above are generally found lying next to each other in genomes. Therefore, these are considered as "single-locus" system. In addition, there is the particular case of the OMF component. It can either be found:
- next to the two other components, as explained just below
- in some other cases, it can be involved in other cellular machineries functioning, and thus be encoded some place else that at the main T1SS' locus (in this case, made of ABC+MFP). 

Therefore, we can attribute the `loner` feature to the OMF component. 

In addition to the latter exception described, it means that this OMF component can also be involved in the functioning of not a single, but several machineries at the same time. In practice, this would mean that two full sets of T1SS components can be inferred with a single OMF component found in the genome. This corresponds to the `multi-system` feature. 



.. image:: ../_static/T1SS_example.*
    :height: 4000px
    :align: center



4. Writing down the model
-------------------------

Now that all elements of the model are listed, the model for the T1SS can be written using the dedicated MacSyFinder XML grammar:


.. code-block:: xml

  <model inter_gene_max_space="5" min_mandatory_genes_required="3" min_genes_required="3" vers="2.0">
      <gene name="T1SS_abc" presence="mandatory"/>
      <gene name="T1SS_mfp" presence="mandatory"/>
      <gene name="T1SS_omf" presence="mandatory" loner="1" multi_system="1"/>
  </model>




.. _T3SS:

The case of T3SS and the bacterial flagellum, or how to distinguish homologous cellular machineries
===================================================================================================


A toy example on how to model similar yet distinct machineries :ref:`here<model-definition-grammar-label>`. 






.. image:: ../_static/T3SS_example.*
    :height: 4000px
    :align: center




Model of the T3SS: 

.. code-block:: xml

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


Model of the Flagellum:

.. code-block:: xml

  <model inter_gene_max_space="20" min_mandatory_genes_required="9" min_genes_required="10" multi_loci="1" vers="2.0">
    <gene name="Flg_sctJ_FLG" presence="mandatory"/>
    <gene name="Flg_sctN_FLG" presence="mandatory"/>
    <gene name="Flg_sctQ_FLG" presence="mandatory"/>
    <gene name="Flg_sctR_FLG" presence="mandatory"/>
    <gene name="Flg_sctS_FLG" presence="mandatory"/>
    <gene name="Flg_sctT_FLG" presence="mandatory"/>
    <gene name="Flg_sctU_FLG" presence="mandatory"/>
    <gene name="Flg_sctV_FLG" presence="mandatory"/>
    <gene name="Flg_flgB" presence="mandatory"/>
    <gene name="Flg_flgC" presence="mandatory"/>
    <gene name="Flg_fliE" presence="mandatory"/>
    <gene name="T3SS_sctC" presence="forbidden"/>
 </model>
