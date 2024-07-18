#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2024  Institut Pasteur (Paris) and CNRS.           #
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

"""
Module use to parse XML model defintion and create a python Model and Genes, ...
"""

import os.path
import xml.etree.ElementTree as Et
import logging

from .model import Model, ModelBank
from .gene import ModelGene, GeneBank
from .gene import Exchangeable
from .error import MacsypyError, ModelInconsistencyError

from typing import Any
from .config import Config, NoneConfig
from .registries import ModelRegistry, DefinitionLocation, ModelLocation
from .profile import ProfileFactory

_log = logging.getLogger(__name__)


class DefinitionParser:
    """
    Build a Model instance from the corresponding model definition described in the XML file.
    """

    def __init__(self,
                 cfg: Config | NoneConfig,
                 model_bank: ModelBank,
                 gene_bank: GeneBank,
                 model_registry: ModelRegistry,
                 profile_factory: ProfileFactory) -> None:
        """
        :param cfg: the configuration object of this run
        :param model_bank: the model factory
        :param gene_bank: the gene factory
        :param model_registry: The registry with all model location
        :param profile_factory: The profile factory
        """
        self.cfg = cfg
        self.model_bank = model_bank
        self.gene_bank = gene_bank
        self.model_registry = model_registry
        self.profile_factory = profile_factory


    def parse(self, models_2_detect: list[DefinitionLocation]) -> None:
        """
        Parse models definition in XML format to build the corresponding Model objects,
        and add them to the model factory after checking its consistency.
        To get the model ask it to model_bank

        :param models_2_detect: a list of model definition to parse.
        """
        models_2_check = []
        _log.info("Models Parsing")
        for def_loc in models_2_detect:
            path = def_loc.path
            if path is None:
                raise MacsypyError(f"{path}: No such model definitions")
            model_location = self.model_registry[def_loc.root_name(def_loc.fqn)]
            model_node = self._get_model_node(def_loc)
            model = self._create_model(def_loc, model_node)
            self.model_bank.add_model(model)
            self._fill_gene_bank(model_node, model_location, def_loc)

            self._parse_genes(model, model_node)
            models_2_check.append(model)
        self.check_consistency(models_2_check)


    def _get_model_node(self, def_loc: DefinitionLocation) -> Et.ElementTree:
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


    def _check_syntax(self, model_node: Et.ElementTree, path: str) -> None:
        """
        Check if the definition does not contain logical error which is allowed by syntax
        and absence of explicit grammar.

        :param model_node: the node corresponding to the model
        :param path: the path of the definition.
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

        model_allowed_attributes = {'inter_gene_max_space',
                                    'min_mandatory_genes_required',
                                    'min_genes_required',
                                    'max_nb_genes',
                                    'multi_loci',
                                    'vers'}
        model_all_attributes = set(model_node.attrib.keys())
        model_unallowed_attribute = model_all_attributes - model_allowed_attributes
        if model_unallowed_attribute:
            msg = f"The model definition {os.path.basename(path)} has an unknow attribute " \
                  f"'{', '.join(model_unallowed_attribute)}'. Please fix the definition."
            raise ModelInconsistencyError(msg)

        gene_allowed_attributes = {'name', 'presence', 'loner', 'multi_system', 'multi_model', 'inter_gene_max_space'}
        gene_all_attributes = set()
        for gene in model_node.iter('gene'):
            gene_all_attributes |= set(gene.attrib.keys())
        gene_unallowed_attribute = gene_all_attributes - gene_allowed_attributes
        if gene_unallowed_attribute:
            msg = f"The model definition {os.path.basename(path)} has an unknown attribute " \
                  f"'{', '.join(gene_unallowed_attribute)}' for a gene." \
                  f" Please fix the definition."
            raise ModelInconsistencyError(msg)


    def _create_model(self, def_loc: DefinitionLocation, model_node: Et.ElementTree) -> Model:
        """
        :param def_loc: the definition location to parse.
        :param model_node: the node corresponding to the model.
        :return: the model corresponding to the definition location.
        """

        inter_gene_max_space = model_node.get('inter_gene_max_space')
        if inter_gene_max_space is None:
            msg = f"Invalid model definition ({def_loc.path}): inter_gene_max_space must be defined"
            _log.critical(msg)
            raise SyntaxError(msg)
        try:
            inter_gene_max_space = int(inter_gene_max_space)
        except ValueError as err:
            msg = f"Invalid model definition ({def_loc.path}): " \
                  f"inter_gene_max_space must be an integer: {inter_gene_max_space}"
            _log.critical(msg)
            raise SyntaxError(msg) from err
        min_mandatory_genes_required = model_node.get('min_mandatory_genes_required')
        if min_mandatory_genes_required is not None:
            try:
                min_mandatory_genes_required = int(min_mandatory_genes_required)
            except ValueError as err:
                msg = f"Invalid model definition ({def_loc.path}): " \
                      f"min_mandatory_genes_required must be an integer: {min_mandatory_genes_required}"
                _log.critical(msg)
                raise SyntaxError(msg) from err

        min_genes_required = model_node.get('min_genes_required')
        if min_genes_required is not None:
            try:
                min_genes_required = int(min_genes_required)
            except ValueError as err:
                msg = f"Invalid model definition ({def_loc.path}):\
 min_genes_required must be an integer: {min_genes_required}"
                _log.critical(msg)
                raise SyntaxError(msg) from err

        cfg_max_nb_genes = self.cfg.max_nb_genes(def_loc.fqn)
        if cfg_max_nb_genes is not None:
            max_nb_genes = cfg_max_nb_genes
        else:
            max_nb_genes = model_node.get('max_nb_genes')
            if max_nb_genes is not None:
                try:
                    max_nb_genes = int(max_nb_genes)
                except ValueError as err:
                    msg = f"Invalid model definition ({def_loc.path}): max_nb_genes must be an integer: {max_nb_genes}"
                    _log.critical(msg)
                    raise SyntaxError(msg) from err
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


    def _fill_gene_bank(self,
                        model_node: Et.ElementTree,
                        model_location: ModelLocation,
                        def_loc: DefinitionLocation) -> None:
        """
        find all gene node and add them to the gene_bank

        :param model_node: the node corresponding to the model.
        :param model_location:
        :param def_loc: a definition location corresponding to the 'model' to parse.
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


    def _parse_gene_attrs(self, gene_node: Et.ElementTree) -> dict[str: Any]:
        attrs = {}
        for attr in ('loner', 'multi_system', 'multi_model'):
            val = gene_node.get(attr)
            if val in ("1", "true", "True"):
                val = True
            elif val in ("0", "false", "False"):
                val = False
            attrs[attr] = val
        inter_gene_max_space = gene_node.get("inter_gene_max_space")
        try:
            inter_gene_max_space = int(inter_gene_max_space)
        except ValueError:
            msg = f"inter_gene_max_space must be an integer: {inter_gene_max_space}"
            raise SyntaxError(msg)
        except TypeError:
            # no attribute inter_gene_max_space
            pass
        else:
            attrs['inter_gene_max_space'] = inter_gene_max_space
        return attrs


    def _parse_genes(self, model: Model, model_node: Et.ElementTree) -> None:
        """
        Create genes belonging to the models.
        Each gene is directly added to the model in its right category ('mandatory, accessory, ...)

        :param model: the Model currently parsing
        :param model_node: the element 'model'
        """

        gene_nodes = model_node.findall("./gene")
        for gene_node in gene_nodes:
            name = gene_node.get("name")
            try:
                attrs = self._parse_gene_attrs(gene_node)
            except SyntaxError as err:
                msg = f"Invalid model definition '{model.fqn}': {err}"
                _log.critical(msg)
                raise SyntaxError(msg)
            new_gene = ModelGene(self.gene_bank[(model.family_name, name)], model, **attrs)

            for exchangeable_node in gene_node.findall("exchangeables/gene"):
                ex = self._parse_exchangeable(exchangeable_node, new_gene, model)
                new_gene.add_exchangeable(ex)

            presence = gene_node.get("presence")
            if not presence:
                msg = f"Invalid model definition '{model.fqn}': gene '{name}' without presence"
                _log.error(msg)
                raise SyntaxError(msg)
            if presence in model.gene_category:
                getattr(model, f'add_{presence}_gene')(new_gene)
            else:
                msg = f"Invalid model '{model.fqn}' definition: presence value must be either: " \
                      f"""{', '.join(["'{}'".format(c) for c in model.gene_category])} not {presence}"""
                _log.error(msg)
                raise SyntaxError(msg)


    def _parse_exchangeable(self, gene_node: Et.ElementTree, gene_ref: ModelGene, curr_model: Model) -> Exchangeable:
        """
        Parse a xml element gene child of exchangeable and build the corresponding object

        :param gene_node: a "node" corresponding to the gene element in the XML hierarchy
        :param gene_ref: the gene which this gene is homolog to
        :param curr_model: the model being parsed .
        :return: the gene object corresponding to the node
        """
        name = gene_node.get("name")
        family_name = curr_model.family_name
        try:
            attrs = self._parse_gene_attrs(gene_node)
        except SyntaxError as err:
            msg = f"Invalid model definition '{curr_model.fqn}': {err}"
            _log.critical(msg)
            raise SyntaxError(msg)

        key = (family_name, name)
        # It cannot fail
        # all genes in the xml are created and insert in GeneBank before this step
        c_gene = self.gene_bank[key]
        ex = Exchangeable(c_gene, gene_ref, **attrs)
        return ex


    def check_consistency(self, models: list[Model]) -> None:
        """
        Check the consistency of the co-localization features between the different values given as an input:
        between XML definitions, configuration file, and command-line options.

        :param models: the list of models to check
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

            if len_mandatory_genes == 0 and len_accessory_genes == 1:
                msg = f"model '{model.name}' is not consistent: there is only one gene in your model. " \
                      f"So its status should be 'mandatory'."
                _log.critical(msg)
                raise ModelInconsistencyError(msg)
