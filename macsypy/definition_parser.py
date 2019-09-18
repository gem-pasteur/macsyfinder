# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Neron                                         #
# Copyright (c) 2014-2019  Institut Pasteur (Paris) and CNRS.                  #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################

import os.path
import xml.etree.ElementTree as Et
import logging
_log = logging.getLogger(__name__)

from .model import Model
from .gene import Gene
from .gene import Homolog, Analog
from .registries import split_def_name, join_def_path
from .error import MacsypyError, ModelInconsistencyError


class DefinitionParser:
    """
    Build a Model instance from the corresponding model definition described in the XML file.
    """

    def __init__(self, cfg, model_bank, gene_bank, profile_factory, model_registry):
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
        self.profile_factory = profile_factory
        self.model_registry = model_registry


    def definition_to_parse(self, def_2_parse, parsed_models):
        """
        :param def_2_parse: the set of definitions fqn to parse ([TXSS/T6SS TXSS/T3SS, ...])
        :type def_2_parse: set of strings {string, ...}
        :param parsed_models: the fully qualified name of the models already build.
        :type parsed_models: set of strings {string, ...}
        :return: a set of definitions' fully qualified names to parse.
                 Scan the whole chain of 'system_ref' in a recursive way.
        :rtype: {string, ...}
        :raises: MacsypyError when Model definition is not found.
        """
        diff_def = parsed_models ^ def_2_parse
        if not diff_def:
            return def_2_parse
        else:            
            for def_fqn in diff_def:
                parsed_models.add(def_fqn)
                def_path = split_def_name(def_fqn)
                model_name = def_path[0]
                try:
                    model_location = self.model_registry[model_name]
                    definition_location = model_location.get_definition(def_fqn)
                except KeyError:
                    raise MacsypyError("{}: No such Models in {}".format(model_name, self.cfg.models_dir()))
                except ValueError:
                    raise MacsypyError("{}: No such definition".format(def_fqn))

                path = definition_location.path
                try:
                    tree = Et.parse(path)
                    root = tree.getroot()
                    if root.tag == 'system':
                        _log.warning(f"'system' is deprecated as xml root. "
                                     f"Migrate {os.path.basename(path)} with macsydef_1to2 script.")
                    # get all genes which are define in an other model
                    # and add these models to the list of models to parse
                    model_ref = root.findall(".//gene[@model_ref]")
                    sys_ref = root.findall(".//gene[@system_ref]")
                    if sys_ref:
                        _log.warning(f"'system_ref' is deprecated. "
                                     f"Migrate {os.path.basename(path)} with macsydef_1to2 script.")
                        model_ref += sys_ref
                    for gene_node in model_ref:
                        def_ref = gene_node.get("model_ref")
                        if not def_ref:
                            def_ref = gene_node.get("system_ref")
                        def_ref_fqn = def_path[:-1]
                        def_ref_fqn.append(def_ref)
                        def_ref_fqn = join_def_path(*def_ref_fqn)
                        def_2_parse.add(def_ref_fqn)
                except Exception as err:
                    msg = "unable to parse model definition \"{}\" : {}".format(def_fqn, err)
                    _log.critical(msg)
                    raise MacsypyError(msg)
            return self.definition_to_parse(def_2_parse, parsed_models)


    def _create_model(self, def_fqn, def_node):
        """
        :param def_fqn: the fully qualified name of the definition.\
          This name must match the path to a definition file.
        :type def_fqn: string
        :param def_node: the node corresponding to the model.
        :type def_node: :class"`Et.ElementTree` object.

        :return: the model corresponding to the definition fully qualified name.
        :rtype: :class:`macsypy.model.Model` object.
        """
        model_name = split_def_name(def_fqn)[0]
        model_location = self.model_registry[model_name]
        definition_location = model_location.get_definition(def_fqn)
        path = definition_location.path

        inter_gene_max_space = def_node.get('inter_gene_max_space')
        if inter_gene_max_space is None:
            msg = "Invalid model definition ({}): inter_gene_max_space must be defined".format(path)
            _log.critical(msg)
            raise SyntaxError(msg)
        try:
            inter_gene_max_space = int(inter_gene_max_space)
        except ValueError:
            msg = "Invalid model definition ({}): " \
                  "inter_gene_max_space must be an integer: {}".format(path, inter_gene_max_space)
            _log.critical(msg)
            raise SyntaxError(msg)
        min_mandatory_genes_required = def_node.get('min_mandatory_genes_required')
        if min_mandatory_genes_required is not None:
            try:
                min_mandatory_genes_required = int(min_mandatory_genes_required)
            except ValueError:
                msg = "Invalid model definition ({}): " \
                      "min_mandatory_genes_required must be an integer: {}".format(path, min_mandatory_genes_required)
                _log.critical(msg)
                raise SyntaxError(msg)

        min_genes_required = def_node.get('min_genes_required')
        if min_genes_required is not None:
            try:
                min_genes_required = int(min_genes_required)
            except ValueError:
                msg = "Invalid model definition ({}):\
 min_genes_required must be an integer: {}".format(path, min_genes_required)
                _log.critical(msg)
                raise SyntaxError(msg)

        max_nb_genes = def_node.get('max_nb_genes')
        if max_nb_genes is not None:
            try:
                max_nb_genes = int(max_nb_genes)
            except ValueError:
                msg = "Invalid model definition ({}): max_nb_genes must be an integer: {}".format(path, max_nb_genes)
                _log.critical(msg)
                raise SyntaxError(msg)

        multi_loci = def_node.get('multi_loci')
        if multi_loci is not None:
            multi_loci = multi_loci.lower() in ("1", "true")
        else:
            multi_loci = False

        # overload value get from xml
        # by these read from configuration (file or command line)
        cfg_inter_gene_max_space = self.cfg.inter_gene_max_space(def_fqn)
        if cfg_inter_gene_max_space is not None:
            inter_gene_max_space = cfg_inter_gene_max_space

        cfg_min_mandatory_genes_required = self.cfg.min_mandatory_genes_required(def_fqn)
        if cfg_min_mandatory_genes_required is not None:
            min_mandatory_genes_required = cfg_min_mandatory_genes_required

        cfg_min_genes_required = self.cfg.min_genes_required(def_fqn)
        if cfg_min_genes_required is not None:
            min_genes_required = cfg_min_genes_required

        cfg_max_nb_genes = self.cfg.max_nb_genes(def_fqn)
        if cfg_max_nb_genes:
            max_nb_genes = cfg_max_nb_genes

        cfg_multi_loci = self.cfg.multi_loci(def_fqn)
        if cfg_multi_loci:
            multi_loci = cfg_multi_loci

        model = Model(def_fqn,
                      inter_gene_max_space,
                      min_mandatory_genes_required,
                      min_genes_required,
                      max_nb_genes,
                      multi_loci)
        return model


    def _create_genes(self, model, def_node):
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
        genes = []
        created_genes = set()
        gene_nodes = def_node.findall(".//gene")
        gene_nodes = [gene_node for gene_node in gene_nodes if gene_node.get("system_ref") is None]
        for node in gene_nodes:
            name = node.get("name")
            if name in created_genes:
                continue
            if not name:
                msg = "Invalid model definition '{}': gene without a name".format(model.name)
                _log.error(msg)
                raise SyntaxError(msg)

            attrs = {}
            for attr in ('loner', 'exchangeable', 'multi_system'):
                val = node.get(attr)
                if val in ("1", "true", "True"):
                    val = True
                elif val in (None, "0", "false", "False"):
                    val = False
                attrs[attr] = val
            inter_gene_max_space = node.get("inter_gene_max_space")
            try:
                inter_gene_max_space = int(inter_gene_max_space)
            except ValueError:
                msg = "Invalid model definition '{}': " \
                      "inter_gene_max_space must be an integer: {}".format(model.name, inter_gene_max_space)
                _log.critical(msg)
                raise SyntaxError(msg)
            except TypeError:
                pass
            else:
                attrs['inter_gene_max_space'] = inter_gene_max_space
            model_name = split_def_name(model.fqn)[0]
            model_location = self.model_registry[model_name]
            new_gene = Gene(self.profile_factory, name, model, model_location, **attrs)
            genes.append(new_gene)
            created_genes.add(new_gene.name)
        return genes


    def _fill(self, model, def_node):
        """
        Fill the model with genes found in this model definition. Add homologs to the genes if necessary.
        
        :param model: the model to fill
        :type model: :class:`macsypy.model.Model` object
        :param def_node: the "node" in the XML hierarchy corresponding to the model
        :type def_node: :class"`Et.ElementTree` object
        """
        genes_nodes = def_node.findall("gene")
        for gene_node in genes_nodes:
            presence = gene_node.get("presence")
            if not presence:
                msg = "Invalid model definition '{}': gene without presence".format(model.name)
                _log.error(msg)
                raise SyntaxError(msg)
            gene_name = gene_node.get('name')
            model_name = split_def_name(model.fqn)[0]
            key = (model_name, gene_name)
            gene = self.gene_bank[key]
            for homolog_node in gene_node.findall("homologs/gene"):
                gene.add_homolog(self._parse_homolog(homolog_node, gene, model))
            for analog_node in gene_node.findall("analogs/gene"):
                gene.add_analog(self._parse_analog(analog_node, gene, model))
            if presence == 'mandatory':
                model.add_mandatory_gene(gene)
            elif presence == 'accessory':
                model.add_accessory_gene(gene)
            elif presence == 'forbidden':
                model.add_forbidden_gene(gene)
            else:
                msg = "Invalid model '{}' definition: presence value must be either\
 [mandatory, accessory, forbidden] not {}".format(model.name, presence)
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
        if not name:
            msg = "Invalid model definition: gene without name"
            _log.error(msg)
            raise SyntaxError(msg)
        try:
            model_name = split_def_name(curr_model.fqn)[0]
            key = (model_name, name)
            gene = self.gene_bank[key]
        except KeyError:
            msg = "Invalid model definition '{}': The gene '{}' described as homolog of\
 '{}' in model '{}' is not in the 'GeneBank' gene factory".format(curr_model.name, name,
                                                                  gene_ref.name, curr_model.name)
            _log.critical(msg)
            raise ModelInconsistencyError(msg)

        model_ref = node.get("system_ref")
        if model_ref is not None and model_ref != gene.model.name:
            msg = "Inconsistency in models definitions: the gene '{}' described as homolog of '{}'\
 with system_ref '{}' has an other model in bank ({})".format(name, gene_ref.name, model_ref, gene.model.name)
            _log.critical(msg)
            raise ModelInconsistencyError(msg)
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
            msg = "Invalid model definition '{}': gene without name".format(curr_model.name)
            _log.error(msg)
            raise SyntaxError(msg)
        try:
            model_name = split_def_name(curr_model.fqn)[0]
            key = (model_name, name)
            gene = self.gene_bank[key]
        except KeyError:
            msg = "Invalid model definition '{}': The gene '{}' described as analog of '{}' in model '{}'\
 is not in the 'GeneBank' gene factory".format(curr_model.name, name, gene_ref.name, curr_model.name)
            _log.critical(msg)
            raise ModelInconsistencyError(msg)
        model_ref = node.get("system_ref")
        if model_ref is not None and model_ref != gene.model.name:
            msg = f"Inconsistency in models definitions: the gene '{name}' described as analog of\
 '{gene_ref.name}' with system_ref '{model_ref}' has an other model in bank ({gene.model.name})"
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
                msg = "model '{}' is not consistent: min_genes_required {:d} must be\
 lesser or equal than the number of \"accessory\" and \"mandatory\" components\
 in the model: {:d}".format(model.name, model.min_genes_required, len_accessory_genes + len_mandatory_genes)
                _log.critical(msg)
                raise ModelInconsistencyError(msg)

            if not (model.min_mandatory_genes_required <= len_mandatory_genes):
                msg = "model '{}' is not consistent: 'min_mandatory_genes_required': {:d} " \
                      "must be lesser or equal than the number of 'mandatory' components " \
                      "in the model: {:d}".format(model.name, model.min_mandatory_genes_required, len_mandatory_genes)
                _log.critical(msg)
                raise ModelInconsistencyError(msg)
            # the following test
            # model.min_mandatory_genes_required <= model.min_genes_required
            # is done during the model.__init__


    def parse(self, models_2_detect):
        """
        Parse models definition in XML format to build the corresponding Model objects,
        and add them to the model factory after checking its consistency.
        To get the model ask it to model_bank

        :param models_2_detect: a list with the fully qualified names of the models to parse
                               (eg ['TXSS/T2SS', 'CRISPR-Cas/typing/CAS-TypeII', ...])
        :type models_2_detect: list of string
        """
        # one opening/closing file / definition
        parsed_defs = set()
        models_2_detect = {s for s in models_2_detect}
        # one opening /closing file /definition
        defs_2_parse = self.definition_to_parse(models_2_detect, parsed_defs)
        msg = "\nModel(s) to parse (recursive inclusion of 'system_ref'):"
        
        for s in defs_2_parse:
            msg += "\n\t-{}".format(s)
        _log.info(msg)
        
        for def_fqn in defs_2_parse:
            model_name = split_def_name(def_fqn)[0]
            model_location = self.model_registry[model_name]
            definition_location = model_location.get_definition(def_fqn)
            path = definition_location.path
            if path is None:
                raise MacsypyError("{}: No such model definitions".format(path))
            tree = Et.parse(path)
            model_node = tree.getroot()
            model = self._create_model(def_fqn, model_node)  # one opening&closing file /definition
            self.model_bank.add_model(model)
            genes = self._create_genes(model, model_node)
            for g in genes:
                try:
                    self.gene_bank.add_gene(g)
                except KeyError:
                    msg = f"gene '{g.name}' define in '{def_fqn}' model is already defined in an another model"
                    _log.error(msg)
                    raise MacsypyError(msg)

        # Now, all model definition related (e.g. via system_ref) to the one to detect are filled appropriately.
        for def_fqn in defs_2_parse:
            model = self.model_bank[def_fqn]
            model_name = split_def_name(def_fqn)[0]
            model_location = self.model_registry[model_name]
            definition = model_location.get_definition(def_fqn)
            path = definition.path

            tree = Et.parse(path)
            def_node = tree.getroot()
            self._fill(model, def_node)
        model_2_check = {self.model_bank[s] for s in models_2_detect}
        self.check_consistency(model_2_check)


