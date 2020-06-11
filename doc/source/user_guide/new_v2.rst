.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2020 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _new_v2:


What's new in MacSyFinder V2? 
=============================

For Version 2, MacSyFinder was carried under **Python 3** https://www.python.org/download/releases/3.0/

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
