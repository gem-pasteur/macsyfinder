# -*- coding: utf-8 -*-

#################################
# Created on Jul 12, 2013
#
# @author: bneron
# @contact: user_email
# @organization: Institut Pasteur
# @license: license
###################################
import os
import xml.etree.ElementTree as ET
import logging
_log = logging.getLogger('txsscan.' + __name__)

from system import System
from gene import Gene
from gene import Homolog
from registries import ProfilesRegistry, DefinitionsRegistry
from txsscan_error import TxsscanError, SystemInconsistencyError


class SystemParser(object):
    """
    Build a System instance from the corresponding System definition described in the XML file (named after the system's name) found at the dedicated location ("-d" command-line option). 
    """

    def __init__(self, cfg, system_bank, gene_bank):
        """
        Constructor

        :param cfg: the configuration object of this run
        :type cfg: :class:`txsscanlib.config.Config` object
        :param cfg: the system factory
        :type cfg: :class:`txsscanlib.system.SystemBank` object
        :param cfg: the gene factory
        :type cfg: :class:`txsscanlib.gene.GeneBank` object
        """
        self.cfg = cfg
        self.system_bank = system_bank
        self.gene_bank = gene_bank
        self.profiles_registry = ProfilesRegistry(cfg)
        self.definitions_registry = DefinitionsRegistry(cfg)
        
        
    def system_to_parse(self, sys_2_detect):
        """
        :param sys_2_detect: the list of systems to detect
        :type sys_2_detect: [string, ...]
        :return: the list of systems' names to parse
        :rtype: [string, ...]
        """
        systems_2_parse = {}
        for system_name in sys_2_detect:
            systems_2_parse[system_name] = None
            path = os.path.join(self.cfg.def_dir, system_name + ".xml")
            if not os.path.exists(path):
                raise TxsscanError("%s: No such system definitions" % path)
            try:
                tree = ET.parse(path)
                root = tree.getroot()
                sys_ref = root.findall(".//gene[@system_ref]")
                for gene_node in sys_ref:
                    systems_2_parse[gene_node.get("system_ref")] = None
            except Exception, err:
                msg = "unable to parse system definition \"{0}\" : {1}".format(system_name, err)
                _log.critical(msg)
                raise TxsscanError(msg)
            
        return systems_2_parse.keys()

    def _create_system(self, system_name, system_node):
        """
        :param system_name: the name of the system to create. This name must match a XML file in the definition directory ("-d" option in the command-line)
        :type param: string
        :return: the system corresponding to the name
        :rtype: :class:`txsscanlib.system.System` object 
        """
        path = os.path.join(self.cfg.def_dir, system_name + ".xml")
        if not os.path.exists(path):
            raise Exception("%s: No such system definitions" % path)
        tree = ET.parse(path)
        root = tree.getroot()
        inter_gene_max_space = root.get('inter_gene_max_space')
        if inter_gene_max_space is None:
            msg = "Invalid system definition (%s): inter_gene_max_space must be defined" % path
            _log.critical(msg)
            raise SyntaxError(msg)
        try:
            inter_gene_max_space = int(inter_gene_max_space)
        except ValueError:
            msg = "Invalid system definition: inter_gene_max_space must be an integer: %s" % inter_gene_max_space
            _log.critical(msg)
            raise SyntaxError(msg)
        min_mandatory_genes_required = root.get('min_mandatory_genes_required')
        if min_mandatory_genes_required is not None:
            try:
                min_mandatory_genes_required = int(min_mandatory_genes_required)
            except ValueError:
                msg = "Invalid system definition: min_mandatory_genes_required must be an integer: %s" % min_mandatory_genes_required
                _log.critical(msg)
                raise SyntaxError(msg)

        min_genes_required = root.get('min_genes_required')
        if min_genes_required is not None:
            try:
                min_genes_required = int(min_genes_required)
            except ValueError:
                msg = "Invalid system definition: min_genes_required must be an integer: %s" % min_genes_required
                _log.critical(msg)
                raise SyntaxError(msg)

        max_nb_genes = root.get('max_nb_genes')
        if max_nb_genes is not None:
            try:
                max_nb_genes = int(max_nb_genes)
            except ValueError:
                msg = "Invalid system definition: max_nb_genes must be an integer: %s" % max_nb_genes
                _log.critical(msg)
                raise SyntaxError(msg)

        multi_loci = root.get('multi_loci')
        if multi_loci is not None:
            multi_loci = multi_loci.lower() in ("1", "true")
        else:
            multi_loci = False
        system = System(self.cfg,
                        system_name,
                        inter_gene_max_space,
                        min_mandatory_genes_required,
                        min_genes_required,
                        max_nb_genes,
                        multi_loci)
        return system

    def _create_genes(self, system, system_node):
        """
        Create genes belonging to the systems. Be careful, the returned genes have not their homologs set yet.
        
        :param gene_name: 
        :type gene_name: string
        :param system:
        :type system: :class:`txsscanlib.system.System` object
        :return: a list of the genes belonging to the system.
        :rtype: [:class:`txsscanlib.gene.Gene`, ...]
        """
        genes = []
        #gene_nodes = system_node.findall(".//gene[not(@system_ref)]")
        gene_nodes = system_node.findall(".//gene")
        gene_nodes = [gene_node for gene_node in gene_nodes if gene_node.get("system_ref") is None ]
        for node in gene_nodes:
            name = node.get("name")
            if not name:
                msg = "Invalid system definition: gene without a name"
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
                msg = "Invalid system definition: inter_gene_max_space must be an integer: %s" % inter_gene_max_space
                _log.critical(msg)
                raise SyntaxError(msg)
            except TypeError:#if None
                pass
            else:
                attrs['inter_gene_max_space'] = inter_gene_max_space
            genes.append(Gene(self.cfg, name, system, self.profiles_registry, **attrs))
        return genes

    def _fill(self, system, system_node):
        """
        Fill the system with genes found in this system definition. Add homolgs to the genes if necessary.
        
        :param system: the system to fill
        :type system: :class:`txsscanlib.system.System` object
        :param system_node: the "node" in the XML hierarchy corresponding to the system
        :type system_node: :class"`ET.ElementTree` object
        """
        genes_nodes = system_node.findall("gene")
        for gene_node in genes_nodes:
            presence = gene_node.get("presence")
            if not presence:
                msg = "Invalid system definition: gene without presence"
                _log.error(msg)
                raise SyntaxError(msg)
            gene_name = gene_node.get('name')
            gene = self.gene_bank[gene_name]
            for homolog_node in gene_node.findall("homologs/gene"):
                gene.add_homolog(self._parse_homolog(homolog_node, gene, system))
            if presence == 'mandatory':
                system.add_mandatory_gene(gene)
            elif presence == 'allowed':
                system.add_allowed_gene(gene)
            elif presence == 'forbidden':
                system.add_forbidden_gene(gene)
            else:
                msg = "Invalid system definition: presence value must be either [mandatory, allowed, forbidden] not %s" % presence
                _log.error(msg)
                raise SyntaxError(msg)


    def _parse_homolog(self, node, gene_ref, curr_system):
        """
        Parse a xml element gene and build the corresponding object

        :param node: a "node" corresponding to the gene element in the XML hierarchy
        :type node: :class:`xml.etree.ElementTree.Element` object 
        :return: the gene object corresponding to the node
        :rtype: :class:`txsscanlib.gene.Homolog` object 
        """
        name = node.get("name")
        if not name:
            msg = "Invalid system definition: gene without name"
            _log.error(msg)
            raise SyntaxError(msg)
        aligned = node.get("aligned")
        if aligned in ("1", "true", "True"):
            aligned = True
        elif aligned in (None, "0", "false", "False"):
            aligned = False
        else:
            msg = 'Invalid system definition: invalid value for an attribute of gene %s: %s allowed values are "1", "true", "True", "0", "false", "False"'% (aligned, name)
            _log.error(msg)
            raise SyntaxError(msg)
        try:
            gene = self.gene_bank[name]
        except KeyError:
            msg = "The gene %s described as homolog of %s in system %s is not in the \"GeneBank\" gene factory" % (name, gene_ref.name, curr_system.name)
            _log.critical(msg)
            raise SystemInconsistencyError(msg)
        system_ref = node.get("system_ref")       
        #if system_ref != gene.system.name:     
        if system_ref != None and system_ref != gene.system.name:
            msg = "Inconsistency in systems definitions: the gene %s described as homolog of %s with system_ref %s has an other system in bank (%s)" % (name, gene_ref.name, system_ref, gene.system.name)
            _log.critical(msg)
            raise SystemInconsistencyError(msg)
        homolog = Homolog(gene, gene_ref, aligned)
        for homolog_node in node.findall("homologs/gene"):
            h2 = self._parse_homolog(homolog_node , gene, curr_system)
            homolog.add_homolog(h2)
        return homolog


    def check_consitency(self, systems, cfg):
        """
        Check the consistency of the co-localization features between the different values given as an input: 
        between XML definitions, configuration file, and command-line options.
        
        :param systems: the list of systems to check
        :type systems: list of `class:txssnalib.system.System` object
        :param cfg: the configuration object
        :type cfg: a :class:`txsscanlib.config.Config` object
        :raise: :class:`txsscanlib.txsscan_error.SystemInconsistencyError` if one test fails

        (see `feature <https://projets.pasteur.fr/issues/1850>`_)
          
        In the different possible situations, different requirements need to be fulfilled ("mandatory_genes" and "allowed_genes" consist of lists of genes defined as such in the system definition):
          
          - **If:** min_mandatory_genes_required = None  ; min_genes_required = None
          - **Then:** min_mandatory_genes_required = min_genes_required = len(mandatory_genes)
          
          *always True by Systems design*

          - **If:** min_mandatory_genes_required = value  ; min_genes_required = None
          - **Then:** min_mandatory_genes_required <= len(mandatory_genes) 
          - AND min_genes_required = min_mandatory_genes_required
          
          *always True by design*

          - **If:** min_mandatory_genes_required =  None ; min_genes_required = Value
          - **Then:** min_mandatory_genes_required = len(mandatory_genes)
          - AND min_genes_required >= min_mandatory_genes_required 
          - AND min_genes_required <= len(mandatory_genes+allowed_genes)
          
          *to be checked*

          - **If:** min_mandatory_genes_required =  Value ; min_genes_required = Value
          - **Then:** min_genes_required <= len(allowed_genes+mandatory_genes)
          - AND min_genes_required >= min_mandatory_genes_required
          - AND min_mandatory_genes_required <= len(mandatory_genes)
           
          *to be checked*   
                 
        """
        for system in systems:
            len_allowed_genes = len(system.allowed_genes)
            len_mandatory_genes = len(system.mandatory_genes)
            if not (system.min_genes_required <= (len_allowed_genes + len_mandatory_genes)) :
                msg = "system %s is not consistent: min_genes_required %d must be lesser or equal than the number of \"allowed\" and \"mandatory\" components in the system: %d" %(system.name, 
                                                                                                                                      system.min_genes_required, 
                                                                                                                                      len_allowed_genes + len_mandatory_genes)
                _log.critical(msg)
                raise SystemInconsistencyError(msg)

            if not (system.min_mandatory_genes_required <= len_mandatory_genes):
                msg = "system %s is not consistent: min_mandatory_genes_required %d must be lesser or equal than the number of \"mandatory\" components in the system: %d" %(system.name, 
                                                                                                                      system.min_mandatory_genes_required, 
                                                                                                                      len_mandatory_genes)
                _log.critical(msg)
                raise SystemInconsistencyError(msg)
            if not(system.min_mandatory_genes_required <= system.min_genes_required):
                msg = "system %s is not consistent: min_mandatory_genes_required %d must be lesser or equal than min_genes_required %d" %(system.name, 
                                                                                                                      system.min_mandatory_genes_required, 
                                                                                                                      system.min_genes_required)
                _log.critical(msg)
                raise SystemInconsistencyError(msg)


    def parse(self, systems_2_detect):
        """
        Parse systems definition in XML format to build the corresponding system objects, and add them to the system factory.

        :param systems_2_detect: a list with the names of the systems to parse
        :type systems_2_detect: list of string
        """
        systems_2_parse = self.system_to_parse(systems_2_detect)  # une ouverture fermeture de fichier /systeme
        for system_name in systems_2_parse:
            path = self.definitions_registry.get(system_name)
            if path is None:
                raise TxsscanError("%s: No such system definitions" % path)
            tree = ET.parse(path)
            system_node = tree.getroot()
            sys = self._create_system(system_name, system_node)  # une ouverture par fichier
            self.system_bank.add_system(sys)
            genes = self._create_genes(sys, system_node)
            for g in genes:
                self.gene_bank.add_gene(g)
        for system_name in systems_2_detect: #une ouverture par fichier
            system = self.system_bank[system_name]
            path = os.path.join(self.cfg.def_dir, system_name + ".xml")
            tree = ET.parse(path)
            system_node = tree.getroot()
            self._fill(system, system_node)
        system_2_check = [self.system_bank[s] for s in systems_2_detect]
        self.check_consitency(system_2_check, self.cfg)

