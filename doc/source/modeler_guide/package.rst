.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2020 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _package:

**************
Models Package
**************


MacSyFinder relies on the definition of models of macromolecular systems as a **set of models' components** 
to be searched by similarity search, and a **set of rules** regarding their genomic organization and 
their requirement level to make a complete system (mandatory, accessory components, number of components required). 

See :ref:`below<model-definition-grammar-label>` for more details on MacSyFinder's modelling scheme and the section 
on :ref:`Functioning <functioning>` for the principles of the MacSyFinder's search engine.


A **MacSyFinder model** (macsy-model for short) is the association of several elements:

    * a **definition** which describes the system to detect with a specific **XML grammar** that is described :ref:`below<model-definition-grammar-label>`.
    
    * a set of :ref:`HMM profiles <provide-hmm_label>`  (one per component/gene in the model) to enable the similarity search of the systems' components with the HMMER program.

The models are grouped by *family* possibly gathering *sub-families* (multiple levels allowed), for instance *Secretion*, *Cas-proteins*...
A set of models from a same family (coherent set) of systems to detect is called hereafter a **macsy-model package** ``NEW in V2``.



.. _package_structure:


Structure of a macsy-model package
==================================

A macsy-model package follows the following structure ::

    family_name
        |_______ metadata.yml
        |_______ LICENCE
        |_______ README.md
        |_______ definitions
        |            |________ model_1.xml
        |            |________ model_2.xml
        |            :
        |
        |_______ profiles
                     |________ geneA.hmm
                     |________ geneB.hmm


If the package contains sub-families ::

    family_name
        |_______ metadata.yml
        |_______ LICENCE
        |_______ README.md
        |_______ definitions
        |            |________ subfamilyA
        |            |            |________ model_1.xml
        |            |            |________ model_2.xml
        |            |
        |            |________ subfamilyB
        |            |            |________ model_3.xml
        |            |            |________ model_4.xml
        |            |
        |            :
        |
        |_______ profiles
                     |________ geneA.hmm
                     |________ geneB.hmm


For examples of macsy-model packages, please visit https://github.com/macsy-models




README.md
---------

A description of the package: what kind of systems the package models, how to use it etc... in `markdown <https://guides.github.com/features/mastering-markdown/>`_ format.

LICENCE
-------

The licence use to protect and share your work.
If you don't know which licence to choose, have a look at `CreativeCommons <https://creativecommons.org/share-your-work/>`_
*This file is optional, but highly recommended.*

