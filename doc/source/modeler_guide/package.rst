.. MacSyFinder - Detection of macromolecular systems in protein datasets
    using systems modelling and similarity search.            
    Authors: Sophie Abby, Bertrand Néron                                 
    Copyright © 2014-2021 Institut Pasteur (Paris) and CNRS.
    See the COPYRIGHT file for details                                    
    MacsyFinder is distributed under the terms of the GNU General Public License (GPLv3). 
    See the COPYING file for details.  
    
.. _model_package:

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
        |_______ model_conf.xml
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
        |_______ model_conf.xml
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

A description of the package: what kind of systems the package models,
how to use it etc... in `markdown <https://guides.github.com/features/mastering-markdown/>`_ format.
The Readme is display to the user on the macsy-models repository on github.
It is also display whe the user run `macsydata help`.


LICENCE
-------

The licence use to protect and share your work.
If you don't know which licence to choose, have a look at `CreativeCommons <https://creativecommons.org/share-your-work/>`_
*This file is optional, but highly recommended.*


metadata file
-------------

This file contains some meta information about the package itself.


model configuration
-------------------

The modeler have the possibility to specify some options specific for its package
different than the masyfinder defaults (new in v2).

This options can be grouped in two families: the scoring weights and filtering options.

scoring weights:

    * mandatory (*float* default = 1.0)
    * accessory (*float* default = 0.5)
    * exchangeable (*float* default = 0.8)
    * loner_multi_systems (*float* default =  0.7)
    * redundancy_penalty (*float* default = 1.5)

filtering options:

    * e_value_search (*float* default = 0.1)
    * i_evalue_sel (*float* default = 0.001)
    * profile_coverage (*float* default = 0.5)
    * cut_ga (*bool* default = True)

All this options are optional and can be omitted in the configuration file, the file itself is optional.
The precedence rules between the different level of configuration are:

 system < home < model < project < --cfg-file | --previous-run < command line options


 * **system**: file in /etc/macsyfinder/macsyfinder.conf on in virtalenv/etc/macsyfinder/macsyfinder.conf
   in case of virtualenv this configuration affect only the macsyfinder installed in this virtualenv
 * **home**:  ~/.macsyfinder/macsyfinder.conf
 * **model**: file model_conf.xml at the root of model package
 * **project**: a file macsyfinder.conf in the directory where is run the macsyfinder command
 * **cfgfile**: any configuration file specify by the user on the command line (conflict with --previous-run opt)
 * **previous-run**: the macsyfinder.comf find in the results directory of the previous run (conflict with --cfg-file opt)
 * **command line**: any option specify directly on the command line

The model_conf.xml configuration file is in xml format and must have the following structure

.. code-block:: yaml

    <model_config>
        <weights>
            <mandatory>1</mandatory>
            <accessory>0.5</accessory>
            <exchangeable>0.8</exchangeable>
            <redundancy_penalty>1.5</redundancy_penalty>
            <loner_multi_system>0.7</loner_multi_system>
        </weights>
        <filtering>
            <e_value_search>0.1</e_value_search>
            <i_evalue_sel>0.01</i_evalue_sel>
            <coverage_profile>0.5</coverage_profile>
            <cut_ga>True</cut_ga>
        </filtering>
    </model_config>


:ref:`Details about scoring method <combinatorial-exploration>`