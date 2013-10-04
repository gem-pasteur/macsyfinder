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
from txsscan_error import TxsscanError, SystemInconsistencyError


class SystemParser(object):
    """
    build a System instance from System definition write in XML and build a 
    """

    def __init__(self, cfg, system_bank, gene_bank):
        """
        Constructor

        :param cfg: the configuration of this run
        :type cfg: :class:`txsscanlib.config.Config` object
        """
        self.cfg = cfg
        self.system_bank = system_bank
        self.gene_bank = gene_bank

    def system_to_parse(self, sys_2_detect):
        """
        :param sys_2_detect: the list of systems to detect
        :type sys_2_detect: [string, ...]
        :return: the list of systems name to parse
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

    def create_system(self, system_name, system_node):
        """
        :param system_name: the name of system. This name must match a file in Definitions directory
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
        system = System(self.cfg, system_name, inter_gene_max_space, min_mandatory_genes_required, min_genes_required)
        return system

    def create_genes(self, system, system_node):
        """
        create genes belonging to the systems. be carefull the return genes have not their homologs.
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
                msg = "Invalid system definition: gene without name"
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
            genes.append(Gene(self.cfg, name, system, **attrs))
        return genes

    def fill(self, system, system_node):
        """
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
        parse a xml element gene and build the corresponding object

        :param node: a node corresponding to gene element
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
            msg = 'Invalid system definition: invalid value for attribute type for gene %s: %s allowed values are "1","true", "True","0" , "false" , "False" '% (aligned, name)
            _log.error(msg)
            raise SyntaxError(msg)
        try:
            gene = self.gene_bank[name]
        except KeyError:
            msg = "the gene %s describe as homolog of %s in system %s in is not in the bank" % (name, gene_ref.name, gene.system.name)
            _log.critical(msg)
            raise SystemInconsistencyError(msg)
        system_ref = node.get("system_ref")       
        if system_ref is None:
            system_ref = curr_system.name
        if system_ref != gene.system.name:
            msg = "inconsitency in systems definitions: the gene %s describe as homolog of %s with system_ref %s has an other system in bank(%s)" % (name, gene_ref.name, system_ref, gene.system.name)
            _log.critical(msg)
            raise SystemInconsistencyError(msg)
        homolog = Homolog(gene, gene_ref, aligned)
        for homolog_node in node.findall("homologs/gene"):
            h2 = self._parse_homolog(homolog_node , gene, curr_system)
            homolog.add_homolog(h2)
            #h2.add_homolog(self._parse_homolog(homolog_node , gene) )
        return homolog


    def check_consitency(self, systems, cfg):
        """
        :param systems: the list of systems to check
        :type systems: list of `class:txssnalib.system.System` object
        :param cfg: the configuration
        :type cfg: a `class:txsscan.config.Config` object
        :raise: SystemInconsistencyError if one test fail

        see `feature <https://projets.pasteur.fr/issues/1850>`_

          - min_mandatory_genes_required = None  ; min_genes_required = None
          - min_mandatory_genes_required = min_genes_required = len(mandatory_genes)
          always True by Systems design

          - min_mandatory_genes_required = value  ; min_genes_required = None
          - min_mandatory_genes_required = len(mandatory_genes) 
          - AND len(allowed_genes + mandatory_genes) >= min_genes_required >= len(mandatory_genes)
           always True By design

          - min_mandatory_genes_required =  None ; min_genes_required = Value
          - min_genes_required = min_mandatory_genes_required 
          - AND min_mandatory_genes_required <= len(mandatory_genes)
          - min_mandatory_genes_required <= len(mandatory_genes)

          - min_mandatory_genes_required =  Value ; min_genes_required = Value
          - len(allowed_genes+mandatory_genes) >= min_genes_required 
          - AND min_mandatory_genes_required <= len(mandatory_genes) 
          - AND min_genes_required >= min_mandatory_genes_required
        """
        for system in systems:
            # feature https://projets.pasteur.fr/issues/1850

            # min_mandatory_genes_required = None  ; min_genes_required = None
            # min_mandatory_genes_required = min_genes_required = len(mandatory_genes)
            # always True by Systems design

            # min_mandatory_genes_required = value  ; min_genes_required = None
            # min_mandatory_genes_required = len(mandatory_genes) 
            # AND len(allowed_genes+mandatory_genes) >= min_genes_required >= len(mandatory_genes)
            # allways True 

            # min_mandatory_genes_required =  None ; min_genes_required = Value
            # min_genes_required = min_mandatory_genes_required 
            # AND min_mandatory_genes_required <= len(mandatory_genes)
            # min_mandatory_genes_required <= len(mandatory_genes)

            # min_mandatory_genes_required =  Value ; min_genes_required = Value
            # len(allowed_genes+mandatory_genes) >= min_genes_required 
            # AND min_mandatory_genes_required <= len(mandatory_genes) 
            # AND min_genes_required >= min_mandatory_genes_required

            len_allowed_genes = len(system.allowed_genes)
            len_mandatory_genes = len(system.mandatory_genes)
            if system.min_genes_required > (len_allowed_genes + len_mandatory_genes) :
                msg = "systems %s is not consistent min_genes_required %d must be lesser or equal than allowed_genes + mandatory_genes %d" %(system.name, 
                                                                                                                                      system.min_genes_required, 
                                                                                                                                      len_allowed_genes + len_mandatory_genes)
                _log.critical(msg)
                raise SystemInconsistencyError(msg)

            if system.min_mandatory_genes_required > len_mandatory_genes:
                msg = "systems %s is not consistent min_mandatory_genes_required %d must be lesser or equal than mandatory_genes %d" %(system.name, 
                                                                                                                      system.min_mandatory_genes_required, 
                                                                                                                      len_mandatory_genes)
                _log.critical(msg)
                raise SystemInconsistencyError(msg)
            if system.min_mandatory_genes_required > system.min_genes_required:
                msg = "systems %s is not consistent min_mandatory_genes_required %d must be lesser or equal than min_genes_required %d" %(system.name, 
                                                                                                                      system.min_mandatory_genes_required, 
                                                                                                                      system.min_genes_required)
                _log.critical(msg)
                raise SystemInconsistencyError(msg)


    def parse(self, systems_2_detect):
        """
        parse  systems definition in xml format to build the corresponding objects

        :param systems_name: the name of the system to parse
        :type systems_name: string
        :return: the system corresponding to the name 
        :rtype: :class:`txsscanlib.secretion.System` object 
        """
        systems_2_parse = self.system_to_parse(systems_2_detect)  # une ouverture fermeture de fichier /systeme
        for system_name in systems_2_parse:
            path = os.path.join(self.cfg.def_dir, system_name + ".xml")
            if not os.path.exists(path):
                raise Exception("%s: No such system definitions" % path)
            tree = ET.parse(path)
            system_node = tree.getroot()
            sys = self.create_system(system_name, system_node)  # une ouverture par fichier
            self.system_bank.add_system(sys)
            genes = self.create_genes(sys, system_node)
            for g in genes:
                self.gene_bank.add_gene(g)
        for system_name in systems_2_detect: #une ouverture par fichier
            system = self.system_bank[system_name]
            path = os.path.join(self.cfg.def_dir, system_name + ".xml")
            tree = ET.parse(path)
            system_node = tree.getroot()
            self.fill(system, system_node)
        self.check_consitency(self.system_bank, self.cfg)

