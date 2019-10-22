#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2020  Institut Pasteur (Paris) and CNRS.           #
# See the COPYRIGHT file for details                                    #
#                                                                       #
# This file is part of MacSyFinder package.                             #
#                                                                       #
# MacSyFinder is free software: you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or     #
# (at your option) any later version.                                   #
#                                                                       #
# MacSyFinder is distributed in the hope that it will be useful,        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
# GNU General Public License for more details .                         #
#                                                                       #
# You should have received a copy of the GNU General Public License     #
# along with MacSyFinder (COPYING).                                     #
# If not, see <https://www.gnu.org/licenses/>.                          #
#########################################################################

import os.path
import xml.etree.ElementTree as Et
import logging
_log = logging.getLogger(__name__)

from .model import Model
from .gene import CoreGene, ModelGene
from .gene import Homolog, Analog
from .registries import split_def_name, join_def_path
from .error import MacsypyError, ModelInconsistencyError


class DefinitionParser:
    """
    Build a Model instance from the corresponding model definition described in the XML file.
    """

    def __init__(self, cfg, model_bank, gene_bank, model_registry, profile_factory):
        """
        :param cfg: the configuration object of this run
        :type cfg: :class:`macsypy.config.Config` object
        :param model_bank: the model factory
        :type model_bank: :class:`macsypy.model.ModelBank` object
        :param gene_bank: the gene factory
        :type gene_bank: :class:`macsypy.gene.GeneBank` object
        """
        self.cfg = cfg
        self.model_bank = model_bank
        self.gene_bank = gene_bank
        self.model_registry = model_registry
        self.profile_factory = profile_factory


    def parse(self, models_2_detect):
        """
        Parse models definition in XML format to build the corresponding Model objects,
        and add them to the model factory after checking its consistency.
        To get the model ask it to model_bank

        :param models_2_detect: a list with the fully qualified names of the models to parse
                               (eg ['TXSS/T2SS', 'CRISPR-Cas/typing/CAS-TypeII', ...])
        :type models_2_detect: list of string
        """
        models_2_check = []
        _log.info("Models Parsing")
        for def_loc in models_2_detect:
            path = def_loc.path
            if path is None:
                raise MacsypyError(f"{path}: No such model definitions")
            model_node = self._get_model_node(def_loc)
            model = self._create_model(def_loc, model_node)
            self.model_bank.add_model(model)
            model_location = self.model_registry[def_loc.family_name]
            self._fill_gene_bank(model_node, model_location, def_loc)

            self._parse_genes(model, model_node)
            models_2_check.append(model)
        self.check_consistency(models_2_check)


    def _get_model_node(self, def_loc):
        """
        :param def_2_parse: the set of definitions fqn to parse ([TXSS/T6SS TXSS/T3SS, ...])
        :type def_2_parse: set of strings {string, ...}
        :param parsed_models: the fully qualified name of the models already build.
        :type parsed_models: set of strings {string, ...}
        :return: a set of definitions' fully qualified names to parse.
                 Scan the whole chain of 'model_ref' in a recursive way.
        :rtype: {string, ...}
        :raises: MacsypyError when Model definition is not found.
        """
        path = def_loc.path
        try:
            tree = Et.parse(path)
            model_node = tree.getroot()

            vers = model_node.get('vers')
            if vers != '2.0':
                msg = f"The model defintion {os.path.basename(path)} is not versioned. " \
                      f"Please update your model."
                raise ModelInconsistencyError(msg)

            if model_node.tag == 'system':
                msg = f"The model defintion {os.path.basename(path)} is obsolete. Please update your model."
                raise ModelInconsistencyError(msg)

            # get all genes which are define in an other model
            # and add these models to the list of models to parse
            sys_ref = model_node.findall(".//gene[@system_ref]")
            if sys_ref:
                msg = f"The model defintion {os.path.basename(path)} is obsolete. Please update your model."
                raise ModelInconsistencyError(msg)
        except Exception as err:
            msg = f"unable to parse model definition '{def_loc.fqn}' : {err}"
            _log.critical(msg)
            raise MacsypyError(msg) from None
        return model_node


    def _create_model(self, def_loc, model_node):
        """
        :param def_loc: the definition location to parse.
        :type def_fqn: :class:`macsypy.registries.DefinitionLocation` object
        :param model_node: the node corresponding to the model.
        :type model_node: :class"`Et.ElementTree` object.

        :return: the model corresponding to the definition location.
        :rtype: :class:`macsypy.model.Model` object.
        """

        inter_gene_max_space = model_node.get('inter_gene_max_space')
        if inter_gene_max_space is None:
            msg = f"Invalid model definition ({def_loc.path}): inter_gene_max_space must be defined"
            _log.critical(msg)
            raise SyntaxError(msg)
        try:
            inter_gene_max_space = int(inter_gene_max_space)
        except ValueError:
            msg = f"Invalid model definition ({def_loc.path}): " \
                  f"inter_gene_max_space must be an integer: {inter_gene_max_space}"
            _log.critical(msg)
            raise SyntaxError(msg)
        min_mandatory_genes_required = model_node.get('min_mandatory_genes_required')
        if min_mandatory_genes_required is not None:
            try:
                min_mandatory_genes_required = int(min_mandatory_genes_required)
            except ValueError:
                msg = f"Invalid model definition ({def_loc.path}): " \
                      f"min_mandatory_genes_required must be an integer: {min_mandatory_genes_required}"
                _log.critical(msg)
                raise SyntaxError(msg)

        min_genes_required = model_node.get('min_genes_required')
        if min_genes_required is not None:
            try:
                min_genes_required = int(min_genes_required)
            except ValueError:
                msg = f"Invalid model definition ({def_loc.path}):\
 min_genes_required must be an integer: {min_genes_required}"
                _log.critical(msg)
                raise SyntaxError(msg)

        max_nb_genes = model_node.get('max_nb_genes')
        if max_nb_genes is not None:
            try:
                max_nb_genes = int(max_nb_genes)
            except ValueError:
                msg = f"Invalid model definition ({def_loc.path}): max_nb_genes must be an integer: {max_nb_genes}"
                _log.critical(msg)
                raise SyntaxError(msg)

        multi_loci = model_node.get('multi_loci')
        if multi_loci is not None:
            multi_loci = multi_loci.lower() in ("1", "true")
        else:
            multi_loci = False

        # overload value get from xml
        # by these read from configuration (file or command line)
        cfg_inter_gene_max_space = self.cfg.inter_gene_max_space(def_loc.fqn)
        if cfg_inter_gene_max_space is not None:
            inter_gene_max_space = cfg_inter_gene_max_space

        cfg_min_mandatory_genes_required = self.cfg.min_mandatory_genes_required(def_loc.fqn)
        if cfg_min_mandatory_genes_required is not None:
            min_mandatory_genes_required = cfg_min_mandatory_genes_required

        cfg_min_genes_required = self.cfg.min_genes_required(def_loc.fqn)
        if cfg_min_genes_required is not None:
            min_genes_required = cfg_min_genes_required

        cfg_max_nb_genes = self.cfg.max_nb_genes(def_loc.fqn)
        if cfg_max_nb_genes:
            max_nb_genes = cfg_max_nb_genes

        cfg_multi_loci = self.cfg.multi_loci(def_loc.fqn)
        if cfg_multi_loci:
            multi_loci = cfg_multi_loci

        model = Model(def_loc.fqn,
                      inter_gene_max_space,
                      min_mandatory_genes_required=min_mandatory_genes_required,
                      min_genes_required=min_genes_required,
                      max_nb_genes=max_nb_genes,
                      multi_loci=multi_loci)
        return model


    def _fill_gene_bank(self, model_node, model_location, def_loc):
        gene_nodes = model_node.findall(".//gene")
        for gene_node in gene_nodes:
            gene_name = gene_node.get("name")
            if not gene_name:
                msg = f"Invalid model definition '{def_loc.fqn}': gene without name"
                _log.error(msg)
                raise SyntaxError(msg)
            # self.gene_bank.add_new_gene if we add twice same (model_location, gene_name)
            # the second time is NOOP
            self.gene_bank.add_new_gene(model_location, gene_name, self.profile_factory)


    def _parse_genes(self, model, model_node):
        """
        Create genes belonging to the models. Be careful, the returned genes have not their homologs/analogs set yet.
        all genes belonging to an other model (model_ref) are ignored

        :param model: the Model currently parsing
        :type model: :class:`macsypy.model.Model` object
        :param def_node: the element gene
        :type def_node: :class"`Et.ElementTree` object
        :return: a list of the genes belonging to the model.
        :rtype: [:class:`macsypy.gene.Gene`, ...]
        """
        gene_nodes = model_node.findall("./gene")
        for gene_node in gene_nodes:
            name = gene_node.get("name")
            attrs = {}
            for attr in ('loner', 'exchangeable', 'multi_system'):
                val = gene_node.get(attr)
                if val in ("1", "true", "True"):
                    val = True
                elif val in (None, "0", "false", "False"):
                    val = False
                attrs[attr] = val
            inter_gene_max_space = gene_node.get("inter_gene_max_space")
            try:
                inter_gene_max_space = int(inter_gene_max_space)
            except ValueError:
                msg = f"Invalid model definition '{model.name}': " \
                      f"inter_gene_max_space must be an integer: {inter_gene_max_space}"
                _log.critical(msg)
                raise SyntaxError(msg)
            except TypeError:
                pass
            else:
                attrs['inter_gene_max_space'] = inter_gene_max_space
            new_gene = ModelGene(self.gene_bank[(model.family_name, name)], model, **attrs)

            for homolog_node in gene_node.findall("homologs/gene"):
                new_gene.add_homolog(self._parse_homolog(homolog_node, new_gene, model))
            for analog_node in gene_node.findall("analogs/gene"):
                new_gene.add_analog(self._parse_analog(analog_node, new_gene, model))

            presence = gene_node.get("presence")
            if not presence:
                msg = f"Invalid model definition '{model.fqn}': gene '{name}' without presence"
                _log.error(msg)
                raise SyntaxError(msg)
            if presence in model.gene_category:
                getattr(model, f'add_{presence}_gene')(new_gene)
            else:
                msg = f"Invalid model '{model.name}' definition: presence value must be either: " \
                      f"""{', '.join(["'{}'".format(c) for c in model.gene_category])} not {presence}"""
                _log.error(msg)
                raise SyntaxError(msg)


    def _parse_homolog(self, node, gene_ref, curr_model):
        """
        Parse a xml element gene and build the corresponding object

        :param node: a "node" corresponding to the gene element in the XML hierarchy
        :type node: :class:`xml.etree.ElementTree.Element` object.
        :param gene_ref: the gene which this gene is homolog to
        :type gene_ref: class:`macsypy.gene.Gene` object
        :param curr_model: the model being parsed .
        :type curr_model: :class:`macsypy.model.Model` object
        :return: the gene object corresponding to the node
        :rtype: :class:`macsypy.gene.Homolog` object
        :raise SyntaxError: if the model definition does not follow the grammar.
        """
        name = node.get("name")
        model_name = split_def_name(curr_model.fqn)[0]
        key = (model_name, name)
        # It cannot fail
        # all genes in the xml are created and insert in GeneBank before this step
        gene = self.gene_bank[key]
        homolog = Homolog(gene, gene_ref)
        for homolog_node in node.findall("homologs/gene"):
            h2 = self._parse_homolog(homolog_node, gene, curr_model)
            homolog.add_homolog(h2)
        return homolog


    def _parse_analog(self, node, gene_ref, curr_model):
        """
        Parse a xml element gene and build the corresponding object

        :param node: a "node" corresponding to the gene element in the XML hierarchy
        :type node: :class:`xml.etree.ElementTree.Element` object.
        :param gene_ref: the gene which this gene is homolog to
        :type gene_ref: class:`macsypy.gene.Gene` object.
        :type curr_model: :class:`macsypy.model.Model` object
        :return: the gene object corresponding to the node
        :return: the gene object corresponding to the node
        :rtype: :class:`macsypy.gene.Analog` object 
        """
        name = node.get("name")
        if not name:
            msg = f"Invalid model definition '{curr_model.name}': gene without name"
            _log.error(msg)
            raise SyntaxError(msg)
        try:
            model_name = split_def_name(curr_model.fqn)[0]
            key = (model_name, name)
            gene = self.gene_bank[key]
        except KeyError:
            msg = f"Invalid model definition '{curr_model.name}': The gene '{name}' described as analog of " \
                  f"'{gene_ref.name}' in model '{curr_model.name}' is not in the 'GeneBank' gene factory"
            _log.critical(msg)
            raise ModelInconsistencyError(msg)
        model_ref = node.get("model_ref") or node.get("system_ref")
        if model_ref is not None and model_ref != gene.model.name:
            msg = f"Inconsistency in models definitions: the gene '{name}' described as analog of\
 '{gene_ref.name}' with model_ref '{model_ref}' has an other model in bank ({gene.model.name})"
            _log.critical(msg)
            raise ModelInconsistencyError(msg)
        analog = Analog(gene, gene_ref)
        for analog_node in node.findall("analogs/gene"):
            h2 = self._parse_analog(analog_node, gene, curr_model)
            analog.add_analog(h2)
        return analog


    def check_consistency(self, models):
        """
        Check the consistency of the co-localization features between the different values given as an input: 
        between XML definitions, configuration file, and command-line options.
        
        :param models: the list of models to check
        :type models: list of `class:macsypy.model.Model` object
        :raise: :class:`macsypy.error.ModelInconsistencyError` if one test fails

        (see `feature <https://projets.pasteur.fr/issues/1850>`_)
          
        In the different possible situations, different requirements need to be fulfilled
        ("mandatory_genes" and "accessory_genes" consist of lists of genes defined as such in the model definition):
          
          - **If:** min_mandatory_genes_required = None  ; min_genes_required = None
          - **Then:** min_mandatory_genes_required = min_genes_required = len(mandatory_genes)
          
          *always True by Models design*

          - **If:** min_mandatory_genes_required = value  ; min_genes_required = None
          - **Then:** min_mandatory_genes_required <= len(mandatory_genes) 
          - AND min_genes_required = min_mandatory_genes_required
          
          *always True by design*

          - **If:** min_mandatory_genes_required =  None ; min_genes_required = Value
          - **Then:** min_mandatory_genes_required = len(mandatory_genes)
          - AND min_genes_required >= min_mandatory_genes_required 
          - AND min_genes_required <= len(mandatory_genes+accessory_genes)
          
          *to be checked*

          - **If:** min_mandatory_genes_required =  Value ; min_genes_required = Value
          - **Then:** min_genes_required <= len(accessory_genes+mandatory_genes)
          - AND min_genes_required >= min_mandatory_genes_required
          - AND min_mandatory_genes_required <= len(mandatory_genes)
           
          *to be checked*   
                 
        """
        for model in models:
            len_accessory_genes = len(model.accessory_genes)
            len_mandatory_genes = len(model.mandatory_genes)
            if not (model.min_genes_required <= (len_accessory_genes + len_mandatory_genes)):
                msg = f"model '{model.name}' is not consistent: min_genes_required {model.min_genes_required:d} " \
                      f"must be lesser or equal than the number of \"accessory\" and \"mandatory\" components " \
                      f"in the model: {len_accessory_genes + len_mandatory_genes:d}"
                _log.critical(msg)
                raise ModelInconsistencyError(msg)

            if not (model.min_mandatory_genes_required <= len_mandatory_genes):
                msg = f"model '{model.name}' is not consistent: 'min_mandatory_genes_required':" \
                      f" {model.min_mandatory_genes_required:d} must be lesser or equal than the number of " \
                      f"'mandatory' components in the model: {len_mandatory_genes:d}"
                _log.critical(msg)
                raise ModelInconsistencyError(msg)
            # the following test
            # model.min_mandatory_genes_required <= model.min_genes_required
            # is done during the model.__init__


