.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2021  Institut Pasteur (Paris),and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _gallery_models:

********************************
Examples of MacSyFinder's models 
********************************

.. contents:: Table of contents of the examples' gallery
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


4. Writing down the model
-------------------------

Now that all elements of the model are listed, the model for the T1SS can be written using the dedicated MacSyFinder XML grammar:



.. _T3SS:

The case of T3SS and the bacterial flagellum, or how to distinguish homologous cellular machineries
===================================================================================================
