.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2020 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _new_v2:

*****************************
What's new in MacSyFinder V2? 
*****************************

For Version 2, MacSyFinder was carried under **Python 3** https://www.python.org/download/releases/3.0/

==============================
New features and search engine 
==============================

The **search engine** was changed for a more intuitive and comprehensive exploration of putative systems. 

The search is now more thorough and avoid undesirable side-effects of the previous search engine. Being more thorough, it now also 
includes a scoring scheme to build candidate systems from sets of detected components (clusters), and can offer several optimal "solutions" (sets of 
detected systems) based on a combinatorial exploration of detected clusters. 
See :ref:`here for more details <functioning>`.

Several **new features** were added, including:

- a new type of gene component "neutral" was added in order to provide more possibilities for systems' modelling in macsy-models. :ref:`See here <hmmer-options>` for more details.
- more flexibility was introduced in the search for systems' components using HMMER. It is now possible to use the `cut_ga` threshold when provided in the HMM profiles used for components' similarity search. This enables to have a search tailored for each HMM profile, and thus component. :ref:`See here <hmmer-options>` for more details.
- a new file structure was created to better organize MacSyFinder's packages (i.e. that include systems' models and corresponding HMME profiles). :ref:`See here <package_structure>` for details.
- a tool to easily install and distribute MacSyFinder's packages was created. :ref:`See here <macsydata>` for more details on *macsydata*.
- the format for MacSyFinder has slightly changed, in order to offer more possibilities, and more readibility. To see how to carry models from V1 to V2, :ref:`see below <models_v1_v2>`. 


.. note::
 
 In v2, output files were also re-defined. See :ref:`here for more details <outputs>`.


============
Dependencies
============

MacSyFinder v2 no longer requires the *formatdb* or *makeblastdb* tools from NCBI. 
However, new dependencies are used, but as they are Python libraries, it should be transparent for the user, and not require manual installations. See :ref:`here for details<dependencies>`.



.. _models_v1_v2:

=============================
Carrying models from V1 to V2 
=============================

Models from V1 are not compatible straight away with V2.
For those who had designed MacSyFinder's models for Version 1 and would like to carry them for Version 2, here are the changes to consider:


* the keyword "system" was changed:
<system> => <model>

* the keyword `<system_ref>` was removed. 
For a given systems' package, each gene has to be defined only once in a macsy-model. There is no need anymore to reference which model it is from, when used as a component in another system's model. 

* now the version of the macsy-models' type have to be documented as a feature of the "model" keyword, like this: `BLBLBLA="V2"` 


* the following keywords have been replaced:
homologs => exchangeables
analogs => exchangeables

.. note::
 
 "exchangeable" is not a feature anymore, but is replaced by the keyword "exchangeables". 
