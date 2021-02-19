#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2021  Institut Pasteur (Paris) and CNRS.           #
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
from .gene import ModelGene
from .gene import Exchangeable
from .registries import split_def_name
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
        :param model_registry: The registry with all model location
        :type model_registry: :class:`macsypy.registry.ModelRegistry` object
        :param profile_factory: The profile factory
        :type profile_factory: :class:`macsypy.profil.ProfilFactory` object
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

        :param models_2_detect: a list of model definition to parse.
        :type models_2_detect: list of :class:`macsypy.registry.DefinitionLocation`
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
        :param def_loc: a definition location to parse.
        :type def_loc: return the node corresponding to the 'model' tag

        """
        path = def_loc.path
        try:
            tree = Et.parse(path)
            model_node = tree.getroot()
            self._check_syntax(model_node, path)
        except Exception as err:
            msg = f"unable to parse model definition '{def_loc.fqn}' : {err}"
            _log.critical(msg)
            raise MacsypyError(msg) from None
        return model_node


    def _check_syntax(self, model_node, path):
        """
        Check if the definition does not contains logical error which is allow by syntax
        and absence of explicit grammar.

        :param model_node: the node correponding to the model
        :type model_node:  :class:`Et.Element` object
        :param str path: the path of the definition.
        :return: None
        :raises ModelInconsistencyError: if an error is encountered in the document.
        """
        vers = model_node.get('vers')
        msg = None
        if vers is None:
            msg = f"The model definition {os.path.basename(path)} is not versioned. " \
                  f"Please update your model."
        elif vers != '2.0':
            msg = f"The model definition {os.path.basename(path)} has not the right version. " \
                  f"version supported is '2.0'. Please update your model."
            raise ModelInconsistencyError(msg)
        elif model_node.tag == 'system':
            msg = f"The model definition {os.path.basename(path)} is obsolete. Please update your model."

        if msg:
            raise ModelInconsistencyError(msg)

        # get all genes which are define in an other model
        # and add these models to the list of models to parse
        sys_ref = model_node.findall(".//gene[@system_ref]")
        if sys_ref:
            msg = f"The model definition {os.path.basename(path)} is obsolete. Please update your model."
            raise ModelInconsistencyError(msg)
        # Since ElementTree from stdlib provides only limited xpath support,
        # you can use | xpath OR operator only if you are using lxml
        if model_node.findall(".//homologs") + model_node.findall(".//analogs"):
            msg = f"The model definition {os.path.basename(path)} is obsolete. Please update your model."
            raise ModelInconsistencyError(msg)


    def _create_model(self, def_loc, model_node):
        """
        :param def_loc: the definition location to parse.
        :type def_loc: :class:`macsypy.registries.DefinitionLocation` object
        :param model_node: the node corresponding to the model.
        :type model_node: :class:`Et.ElementTree` object.

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

        cfg_max_nb_genes =  self.cfg.max_nb_genes(def_loc.fqn)
        if cfg_max_nb_genes is not None:
            max_nb_genes = cfg_max_nb_genes
        else:
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
        """
        find all gene node and add them to the gene_bank

        :param model_node: :param model_node: the node corresponding to the model.
        :type model_node: :class:`Et.ElementTree` object.
        :param model_location:
        :type model_location: class:`macsypy.registries.ModelLocation` object.
        :param def_loc: a definition location to parse.
        :type def_loc: the node corresponding to the 'model' tag
        :return: None
        """
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
        Create genes belonging to the models.
        Each gene is directly added to the model in it's right category ('mandatory, accessory, ...)

        :param model: the Model currently parsing
        :type model: :class:`macsypy.model.Model` object
        :param model_node: the element 'model'
        :type model_node: :class"`Et.ElementTree` object
        """
        gene_nodes = model_node.findall("./gene")
        for gene_node in gene_nodes:
            name = gene_node.get("name")
            attrs = {}
            for attr in ('loner', 'multi_system'):
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

            for exchangeable_node in gene_node.findall("exchangeables/gene"):
                new_gene.add_exchangeable(self._parse_exchangeable(exchangeable_node, new_gene, model))

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


    def _parse_exchangeable(self, node, gene_ref, curr_model):
        """
        Parse a xml element gene child of exchangeable and build the corresponding object

        :param node: a "node" corresponding to the gene element in the XML hierarchy
        :type node: :class:`xml.etree.ElementTree.Element` object.
        :param gene_ref: the gene which this gene is homolog to
        :type gene_ref: class:`macsypy.gene.ModelGene` object
        :param curr_model: the model being parsed .
        :type curr_model: :class:`macsypy.model.Model` object
        :return: the gene object corresponding to the node
        :rtype: :class:`macsypy.gene.Exchangeable` object
        """
        name = node.get("name")
        model_name = split_def_name(curr_model.fqn)[0]
        key = (model_name, name)
        # It cannot fail
        # all genes in the xml are created and insert in GeneBank before this step
        c_gene = self.gene_bank[key]
        ex = Exchangeable(c_gene, gene_ref)
        return ex


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
