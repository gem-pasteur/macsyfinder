# -*- coding: utf-8 -*-

###############################
# Created on Feb 14, 2013
#
# @author: sabby
# @contact: sabby@pasteur.fr
# @organization: Institut Pasteur
# @license: GPLv3
################################

import threading
import logging
from report import Hit # required? 
import os.path
from collections import namedtuple, Counter, OrderedDict
import itertools, operator
import json
from operator import attrgetter # To be used with "sorted"

from txsscan_error import TxsscanError, SystemDetectionError
from database import RepliconDB
from system import system_bank

_log = logging.getLogger('txsscan.' + __name__)

class ClustersHandler(object):
    """
    Deals with sets of clusters found in a dataset. Conceived to store only clusters from a same replicon.
    """

    def __init__(self):
        """
        :param cfg: The configuration object built from default and user parameters.
        :type cfg: :class:`txsscanlib.config.Config` 
        """
        self.clusters = []
        self.replicon_name = ""

    def add(self, cluster):
        if not self.replicon_name:
            self.replicon_name = cluster.replicon_name
            cluster.save()
            self.clusters.append(cluster)
        elif self.replicon_name == cluster.replicon_name:
            cluster.save()
            self.clusters.append(cluster)
        else:
            msg = "Attempting to add a cluster in a ClustersHandler dedicated to another replicon !"
            for c in self.clusters:
                msg+=str(c)
            msg+="To add: %s"%str(cluster)
            _log.critical(msg)
            #raise Exception(msg)
            raise SystemDetectionError(msg)

    def __str__(self):
        to_print=""
        for cluster in self.clusters:
            to_print+=str(cluster)

        return to_print

    def circularize(self, rep_info):
        """
        This function takes into account the circularity of the replicon by merging clusters when appropriate (typically at replicon's ends). 
        It has to be called only if the replicon_topology is set to \"circular\".
        """
        # We assume this function is called when appropriate (i.e. for circular replicons)
        if len(self.clusters) > 1:
            clust_first = self.clusters[0]
            clust_last = self.clusters[len(self.clusters)-1]

            pos_min = rep_info.min
            pos_max = rep_info.max
            dist_clust = clust_first.begin - pos_min + pos_max - clust_last.end

            #if (dist_clust <= max(clust_first.hits[0].get_syst_inter_gene_max_space(), clust_last.hits[len(clust_first.hits)-1].get_syst_inter_gene_max_space())):
            if (dist_clust <= max(clust_first.hits[0].get_syst_inter_gene_max_space(), clust_last.hits[len(clust_last.hits)-1].get_syst_inter_gene_max_space())):
                # Need to circularize ! 
                print " A cluster needs to be \"circularized\" ! "
                self.clusters.pop(0)
                for h in clust_first.hits:
                    clust_last.add(h)
                clust_last.save(True) # Force to re-save the updated cluster
                print clust_last


class Cluster(object):
    """
    Stores a set of contiguous hits. The Cluster object can have different states regarding 
    its content in different genes' systems: 
    
      - ineligible: not a cluster to analyze
      - clear: a single system is represented in the cluster
      - ambiguous: several systems are represented in the cluster => might need a disambiguation
      
    """
    
    #def __init__(self):
    def __init__(self, systems_to_detect):
        self.hits = []
        self.systems_to_detect = systems_to_detect # NEW
        self.systems = {}
        self.replicon_name = ""
        self.begin = 0
        self.end = 0
        self._state = ""
        self._putative_system = ""
        self._compatible_systems = [] # NEW!

    def __len__(self):
        return len(self.hits)

    def __str__(self):
        pos = []
        seq_ids = []
        gene_names = [] 
        for h in self.hits:
            pos.append(h.position)
            seq_ids.append(h.id)
            gene_names.append(h.gene.name)

        if self.state == "clear":
            return "--- Cluster %s ---\n%s\n%s\n%s"%(self.putative_system, str(seq_ids), str(gene_names), str(pos))
        else:
            return "--- Cluster %s ? ---\n%s\n%s\n%s"%(self.putative_system, str(seq_ids), str(gene_names), str(pos))

    @property
    def state(self):
        """
        :return: the state of the cluster of hits
        :rtype: string
        """
        return self._state

    @property
    def putative_system(self):
        """
        :return: the name of the putative system represented by the cluster
        :rtype: string
        """
        return self._putative_system
    
    @property
    def compatible_systems(self):
        """
        :return: the list of the names of compatible systems represented by the cluster
        :rtype: string
        """
        return self._compatible_systems

    def add(self, hit):
        """
        Hits are always added at the end of the cluster (appended to the list of hits). Thus, 'begin' and 'end' positions of the Cluster are always the position of the 1st and of the last hit respectively.
        """
        # need to update cluster bounds
        if len(self.hits) == 0:
            self.begin = hit.get_position()
            self.end = self.begin
            self.replicon_name = hit.replicon_name
            self.hits.append(hit)
        else:
            if(self.replicon_name == hit.replicon_name):
                # To be updated !! make this work also with "circularized" clusters
                """
                if hit.get_position() < self.begin:
                    self.begin = hit.get_position()
                elif hit.get_position() > self.end:
                    self.end = hit.get_position()
                else:
                    if hit.get_position() > self.begin and hit.get_position() < self.end:
                        _log.debug("Weird cluster inclusion hit : %s"%hit)
                """
                # The end position only is updated, as Hits are always appended.
                self.end = hit.get_position()
                self.hits.append(hit)
            else:
                msg = "Attempting to gather in a cluster hits from different replicons ! "
                _log.critical(msg)
                #raise Exception(msg)
                raise SystemDetectionError(msg)


    def save(self, force=False):
        """
        Check the status of the cluster regarding systems which have hits in it. 
        Update systems represented, and assign a putative system (self._putative_system), which is the system with most hits in the cluster. 
        The systems represented are stored in a dictionary in the self.systems variable. 
        The execution of this function can be forced, even if it has already run for the cluster with the option force=True.
        """
        if not self.putative_system or force:
            # First compute the "Majoritary" system
            systems={} # Counter of occcurrences of systems in the cluster (keys are systems' names)
            genes=[]
            #systems_object={} # Store systems object. # To be replaced by "system_bank"? Yep ! done
            
            # Version with hits "reference" systems
            for h in self.hits:
                syst=h.system.name
                if not systems.has_key(syst):
                    systems[syst]=1
                    #systems_object[syst]=h.system
                else:
                    systems[syst]+=1
                if genes.count(h.gene.name) == 0:
                    genes.append(h.gene.name)
            
            self.systems = systems
            
            # Useless ?! Or change? 
            max_syst=0
            tmp_syst_name=""
            for x,y in systems.iteritems():
                if y>=max_syst:
                    tmp_syst_name = x
                    max_syst = y

            self._putative_system = tmp_syst_name # Remove cause useless?
                    
            # NEW Version with hits all "compatible" systems
            systems_compat={} # Counter of occcurrences of COMPATIBLE (-extended- list of) systems in the cluster. Keys are systems' names
            for h in self.hits:
                #syst_list=h.gene.get_compatible_systems(self.systems_to_detect) # Need the list of systems (obj!) to be detected... in the cfg?
                # Now exclude forbidden genes from those that define the list of compatible systems
                #syst_list=h.gene.get_compatible_systems(self.systems_to_detect, False) # Need the list of systems (obj!) to be detected... in the cfg?
                syst_list=h.gene.get_compatible_systems(self.systems_to_detect, True) # Need the list of systems (obj!) to be detected... in the cfg?
                for syst in syst_list:
                    syst_name=syst.name
                    if not systems_compat.has_key(syst_name):
                        systems_compat[syst_name]=1
                        #systems_object[syst]=h.system
                    else:
                        systems_compat[syst_name]+=1
                    if genes.count(h.gene.name) == 0:
                        genes.append(h.gene.name)
            
            #print systems_compat
            # We sort the list of compatible systems per decreasing nb of systems occurrences
            systems_compat=OrderedDict(sorted(systems_compat.items(), reverse=True, key=lambda t: t[1]))
            #print systems_compat
            
            if len(genes) == 1 and self.hits[0].gene.loner == False:                            
                self._state = "ineligible"
            else:
                # Check for foreign "allowed" genes... Might increase nb of systems predicted in the cluster, even if they are tolerated in the cluster. 
                # Also deal with foreign "exchangeable" genes for the same reasons... NB !! Maybe just not add the system to the list if exchangeable?  
                if len(systems.keys()) == 1:
                    self._state = "clear"
                    syst = systems.keys()[0]
                    #print syst
                    self._putative_system = syst
                    #self._compatible_systems.append(system_bank[syst])
                    # Store only compatible systems that are searched for !!
                    self._compatible_systems.append(syst)
                    
                else:
                    # Check for foreign "allowed" genes regarding the majoritary system... They might increase nb of systems predicted in the cluster artificially, even if they are tolerated in the cluster. For that need to scan again all hits and ask wether they are allowed foreign genes. 
                    #def try_system(hits, putative_system, systems):
                    def try_system(hits, putative_system, counter_systems_in_clust):   
                        """
                        Test if the putative_system is compatible with the systems of hits (counter_systems_in_clust) 
                        
                        :param hits: a list of hits
                        :type hits: a list of :class:`txsscanlib.report.Hit`
                        :param putative_system: the name of a putative system to consider
                        :type putative_system: string
                        :param counter_systems_in_clust: a dictionary with systems occurrences when exploring the hits
                        :type counter_systems_in_clust: a dictionary with systems' names as keys, and counts as values
                        :return: the "state" of the set of hits regarding the putative system
                        :rtype: string
                        
                        """
                        foreign_allowed = 0
                        auth = 0 # counts nb of hits that are authorized in the putative system
                        #print "Unclear state with multiple systems to deal with..."
                        #print systems
                        for h in hits:
                            #print h
                            #print putative_system
                            #if h.system.name != putative_system and h.gene.is_authorized(systems_object[putative_system]):
                            #if h.system.name != putative_system and h.gene.is_authorized(system_bank[putative_system.name]):
                            # Too complex !
                            #if h.system.name != putative_system and h.gene.is_authorized(system_bank[putative_system]): 
                            #    foreign_allowed+=1
                            #if h.gene.is_authorized(system_bank[putative_system]):
                            # Exclude the consideration of "forbidden" genes !
                            #if h.gene.is_authorized(system_bank[putative_system], False):
                            if h.gene.is_authorized(system_bank[putative_system], True):
                                auth+=1
                    	#if foreign_allowed == sum(counter_systems_in_clust.values())-systems[putative_system]:
                    	#if foreign_allowed == sum(counter_systems_in_clust.values())-counter_systems_in_clust[putative_system]: # Nope ! Too much complicated !
                    	if auth == len(hits):
                            # Case where all foreign genes are allowed in the majoritary system => considered as a clear case, does not need disambiguation.
                            state = "clear"
                    	else:
                            state = "ambiguous"
                        return state

                    #for putative_system in systems.keys():
                    # Sort systems to consider by decreasing counts.
                    cluster_compatible_systems = []
                    for putative_system in systems_compat.keys():
                        # Add that it has to be done first from the most rep systems by decreasing order of systems.values.
                        #state=try_system(self.hits, putative_system, systems)
                        state=try_system(self.hits, putative_system, systems_compat)
                        #state=try_system(self.hits, putative_system.name, systems_compat)
                        if state == "clear":
                            #print "BUENO SYSTEMO %s"%putative_system
                            #self._state="clear" 
                            #self._putative_system=putative_system 
                            #break
                            cluster_compatible_systems.append(putative_system)
                        #else:
                        #    self._state="ambiguous"
                        #    # Aoutch in this case no putative_system?!
                        #    #print "YAPABON...%s"%putative_system
                    if len(cluster_compatible_systems) >= 1:
                        self._state="clear"
                        self._putative_system = cluster_compatible_systems[0]
                        self._compatible_systems = cluster_compatible_systems
                    else:
                        self._state="ambiguous"
                            

class SystemNameGenerator(object):
    """
    Creates and stores the names of detected systems. Ensures the uniqueness of the names.  
    """
    name_bank={}

    def getSystemName(self, replicon, system):
        basename = self._computeBasename(replicon, system)
        if basename in self.name_bank:
            self.name_bank[basename]+=1
        else:
            self.name_bank[basename] = 1

        system_name =  basename+str(self.name_bank[basename])   
        return system_name

    def _computeBasename(self, replicon, system):
        return replicon+"_"+system+"_"

system_name_generator = SystemNameGenerator()



class SystemOccurence(object):
    """
    This class is instantiated for a specific system that has been asked for detection. It can be filled step by step with hits. 
    A decision can then be made according to the parameters defined *e.g.* quorum of genes. 

    The SystemOccurence object has a "state" parameter, with the possible following values: 
      - "empty" if the SystemOccurence has not yet been filled with genes of the decision rule of the system
      - "no_decision" if the filling process has started but the decision rule has not yet been applied to this occurence
      - "single_locus" 
      - "multi_loci" 
      - "uncomplete"
        
    """
    def __init__(self, system):
        """
        :param system: the system to \"fill\" with hits.
        :type system: :class:`txsscanlib.system.System` 
        """
        self.system = system
        self.system_name = system.name

        # Variables to be updated during the system detection 
        self.valid_hits = [] # validSystemHit are stored with the "fill_with" function, and ready for extraction in case of a positive detection

        self.loci_positions = [] # list of tuples

        self._state = "empty"
        self.nb_cluster = 0
        self._nb_syst_genes = 0
        self.unique_name = ""

        # System definition
        # Make those attributes non modifiable?
        self.mandatory_genes = {}
        self.exmandatory_genes = {} # List of 'exchanged' mandatory genes

        # New ! Add of a list of "multi_system" genes, fed only from mandatory and allowed genes from the actual system (and not 'exchanged')
        self.multi_syst_genes = {}
        
        for g in system.mandatory_genes:
            self.mandatory_genes[g.name] = 0
            if g.exchangeable:
                homologs=g.get_homologs()
                for h in homologs:
                    self.exmandatory_genes[h.name] = g.name
            if g.multi_system:
                self.multi_syst_genes[g.name] = 0

        self.allowed_genes = {}
        self.exallowed_genes = {} # List of 'exchanged' allowed genes
        for g in system.allowed_genes:
            self.allowed_genes[g.name] = 0
            if g.exchangeable:
                homologs=g.get_homologs()
                for h in homologs:
                    self.exallowed_genes[h.name] = g.name
            if g.multi_system:
                self.multi_syst_genes[g.name] = 0

        self.forbidden_genes = {}
        for g in system.forbidden_genes:
            self.forbidden_genes[g.name] = 0 

    def __str__(self):
        out=""
        if self.mandatory_genes: 
            out+="Mandatory genes: \n"
            for k, g in self.mandatory_genes.iteritems():
                out+="%s\t%d\n"%(k, g) 
        if self.allowed_genes:
            out+="Allowed genes: \n"  
            for k, g in self.allowed_genes.iteritems():
                out+="%s\t%d\n"%(k, g)  
        if self.forbidden_genes:
            out+="Forbidden genes: \n"  
            for k, g in self.forbidden_genes.iteritems():
                out+="%s\t%d\n"%(k, g)
        # NEW  
        if self.multi_syst_genes:
            out+="Multi_syst genes:\n"
            for k, g in self.multi_syst_genes.iteritems():
                out+="%s\t%d\n"%(k, g)
          
        return out

    #def get_gene_counter_output(self):
    def get_gene_counter_output(self, forbid_exclude=False):
        """
        Returns a dictionary ready for printing in system summary, with genes (mandatory, allowed and forbidden if specified) occurences in the system occurrence        
        """
        out=""
        out+=str(self.mandatory_genes)
        out+="\t%s"%str(self.allowed_genes)
        if not forbid_exclude:
            out+="\t%s"%str(self.forbidden_genes)
        else:
            out+="\t{}"
        return out

    @property
    def state(self):
        """
        :return: the state of a system occurence
        :rtype: string
        """
        return self._state

    def get_system_unique_name(self, replicon_name):
        """
        Attributes unique name to the system occurrence with the class :class:`txsscanlib.search_systems.SystemNameGenerator`.
        Generates the name if not already set. 
        
        :return: the unique name of the :class:`txsscanlib.search_systems.SystemOccurence`
        :rtype: string
        """
        if not self.unique_name:
            self.unique_name = system_name_generator.getSystemName(replicon_name, self.system_name)
        return self.unique_name

    def get_system_name_unordered(self, suffix="_putative"):
        """
        Attributes a name to the system occurrence for an "unordered" dataset => generating a generic name based on the system name and the suffix given. 
        
        :return: a name for a system in an "unordered" dataset to the :class:`txsscanlib.search_systems.SystemOccurence`
        :rtype: string
        """
        return self.system_name+suffix


    def compute_system_length(self, rep_info):
        """
        Returns the length of the system, all loci gathered, in terms of protein number (even those non matching any system gene)
        """
        length=0
        # To be updated to deal with "circular" clusters
        for(begin, end) in self.loci_positions:
            if begin<=end:
                length+=(end-begin+1)
            elif rep_info.topology == "circular":
                locus_length=end-begin+rep_info.max-rep_info.min+2
                length+=locus_length
            else:
                msg="Inconsistency in locus positions in the case of a linear replicon. The begin position of a locus cannot be higher than the end position. \n"
                msg+="Problem with locus found with positions begin: %d end: %d"%(begin, end)
                _log.critical(msg)
                #raise Exception(msg)
                raise SystemDetectionError(msg)
        return length

    @property
    def nb_syst_genes(self):
        """
        This value is set after a decision was made on the system in :func:`txsscanlib.search_systems.SystemOccurence:decision_rule`
        
        :return: the number of mandatory and allowed genes with at least one occurence (number of different allowed genes)
        :rtype: integer
        """
        return self._nb_syst_genes

    def compute_nb_syst_genes(self):
        return self.count_genes(self.mandatory_genes)+self.count_genes(self.allowed_genes)

    def compute_nb_syst_genes_tot(self):
        return self.count_genes_tot(self.mandatory_genes)+self.count_genes_tot(self.allowed_genes)

    def count_genes(self, gene_dict):
        """
        Counts the nb of genes with at least one occurrence in a dictionary with a counter of genes. 
        """
        total = 0
        for v in gene_dict.values():
            if v>0:
                total+=1
        return total

    def count_genes_tot(self, gene_dict):
        """
        Counts the nb of matches in a dictionary with a counter of genes, independently of the nb of genes matched.
        """
        total = 0
        for v in gene_dict.values():
            total+=v
        return total

    def compute_missing_genes_list(self, gene_dict):
        """
        :returns: the list of genes with no occurence in the gene counter. 
        :rtype: list
        """
        missing=[]
        for k,v in gene_dict.iteritems():
            if v==0:
                missing.append(k)
        return missing


    def count_missing_genes(self, gene_dict):
        """
        Counts the number of genes with no occurence in the gene counter.

        :rtype: integer
        """
        return len(self.compute_missing_genes_list(gene_dict))


    def is_complete(self):
        if self.state == "single_locus" or self.state == "multi_loci":
            return True
        else:
            return False 

    def get_summary_header(self):
        """
        Returns a string with the description of the summary returned by self.get_summary()

        :rtype: string
        """
        return "#Replicon_name\tSystem_Id\tReference_system\tSystem_status\tNb_loci\tNb_Ref_mandatory\tNb_Ref_allowed\tNb_Ref_Genes_detected_NR\tNb_Genes_with_match\tSystem_length\tNb_Mandatory_NR\tNb_Allowed_NR\tNb_missing_mandatory\tNb_missing_allowed\tList_missing_mandatory\tList_missing_allowed\tLoci_positions\tOccur_Mandatory\tOccur_Allowed\tOccur_Forbidden"


    def get_summary(self, replicon_name, rep_info):
        """
        Gives a summary of the system occurrence in terms of gene content and localization.

        :return: a tabulated summary of the :class:`txsscanlib.search_systems.SystemOccurence`
        :rtype: string
        """

        report_str = replicon_name+"\t"+self.get_system_unique_name(replicon_name)
        report_str+="\t%s"%self.system_name
        report_str+="\t%s"%self.state
        report_str+="\t%d"%self.nb_cluster # Nb of loci included to fill the system occurrence
        report_str+="\t%d"%len(self.mandatory_genes) # Nb mandatory_genes in the definition of the system
        report_str+="\t%d"%len(self.allowed_genes) # Nb allowed_genes in the definition of the system
        report_str+="\t%d"%self.nb_syst_genes # Nb syst genes NR
        report_str+="\t%d"%self.compute_nb_syst_genes_tot() # Nb syst genes matched
        #report_str+="\t%d"%self.compute_system_length() # The total length of the locus in protein number, delimited by hits for profiles of the system.
        report_str+="\t%d"%self.compute_system_length(rep_info) # The total length of the locus in protein number, delimited by hits for profiles of the system.

        report_str+="\t%d"%self.count_genes(self.mandatory_genes) # Nb mandatory_genes matched at least once
        report_str+="\t%d"%self.count_genes(self.allowed_genes) # Nb allowed_genes matched at least once

        missing_mandatory = self.compute_missing_genes_list(self.mandatory_genes)        
        missing_allowed = self.compute_missing_genes_list(self.allowed_genes)

        report_str+="\t%d"%len(missing_mandatory) # Nb mandatory_genes with no occurrence in the system
        report_str+="\t%d"%len(missing_allowed) # Nb allowed_genes with no occurrence in the system
        report_str+="\t%s"%str(missing_mandatory) # List of mandatory genes with no occurrence in the system
        report_str+="\t%s"%str(missing_allowed) # List of allowed genes with no occurrence in the system

        report_str+="\t%s"%self.loci_positions # The positions of the loci (begin, end) as delimited by hits for profiles of the system.
        report_str+="\t%s"%self.get_gene_counter_output() # A dico per type of gene 'Mandatory, Allowed, Forbidden' with gene occurrences in the system

        return report_str

    def get_summary_unordered(self, replicon_name):
        """
        Gives a summary of the system occurrence in terms of gene content only (specific of "unordered" datasets).

        :return: a tabulated summary of the :class:`txsscanlib.search_systems.SystemOccurence`
        :rtype: string
        """

        #report_str = replicon_name+"\t"+self.get_system_unique_name(replicon_name)
        # No replicon name for unordered... get it from the config object in future developments... 
        report_str = replicon_name+"\t"+self.get_system_name_unordered()
        report_str+="\t%s"%self.system_name
        report_str+="\t%s"%self.state
        
        #report_str+="\t%d"%self.nb_cluster # Nb of loci included to fill the system occurrence
        report_str+="\tNone"# No loci in unordered
        report_str+="\t%d"%len(self.mandatory_genes) # Nb mandatory_genes in the definition of the system
        report_str+="\t%d"%len(self.allowed_genes) # Nb allowed_genes in the definition of the system
        report_str+="\t%d"%self.nb_syst_genes # Nb syst genes NR
        report_str+="\t%d"%self.compute_nb_syst_genes_tot() # Nb syst genes matched
        
        #report_str+="\t%d"%self.compute_system_length(rep_info) # The total length of the locus in protein number, delimited by hits for profiles of the system.
        report_str+="\tNone" # No loci in unordered

        report_str+="\t%d"%self.count_genes(self.mandatory_genes) # Nb mandatory_genes matched at least once
        report_str+="\t%d"%self.count_genes(self.allowed_genes) # Nb allowed_genes matched at least once

        missing_mandatory = self.compute_missing_genes_list(self.mandatory_genes)        
        missing_allowed = self.compute_missing_genes_list(self.allowed_genes)

        report_str+="\t%d"%len(missing_mandatory) # Nb mandatory_genes with no occurrence in the system
        report_str+="\t%d"%len(missing_allowed) # Nb allowed_genes with no occurrence in the system
        report_str+="\t%s"%str(missing_mandatory) # List of mandatory genes with no occurrence in the system
        report_str+="\t%s"%str(missing_allowed) # List of allowed genes with no occurrence in the system

        #report_str+="\t%s"%self.loci_positions # The positions of the loci (begin, end) as delimited by hits for profiles of the system.
        report_str+="\tNone" # No loci in unordered
        report_str+="\t%s"%self.get_gene_counter_output(True) # A dico per type of gene 'Mandatory, Allowed, Forbidden' with gene occurrences in the system

        return report_str

    def fill_with_cluster(self, cluster):
        """
        Adds hits from a cluster to a system occurence, and check which are their status according to the system definition.
        Set the system occurence state to "no_decision" after calling of this function.

        :param cluster: the set of contiguous genes to treat for :class:`txsscanlib.search_systems.SystemOccurence` inclusion. 
        :type cluster: :class:`txsscanlib.search_systems.Cluster`
        """
        included = True
        self._state = "no_decision"
        for hit in cluster.hits:
            # Need to check first that this cluster is eligible for system inclusion
            # Stores hits for system extraction (positions, sequences) when complete. 

            if hit.gene.is_mandatory(self.system):
                self.mandatory_genes[hit.gene.name]+=1
                valid_hit=validSystemHit(hit, self.system_name, "mandatory")
                self.valid_hits.append(valid_hit)
                # NEW
                if hit.gene.multi_system:
                    self.multi_syst_genes[hit.gene.name]+=1                    
            elif hit.gene.is_allowed(self.system):
                self.allowed_genes[hit.gene.name]+=1
                valid_hit=validSystemHit(hit, self.system_name, "allowed")
                self.valid_hits.append(valid_hit)
                # NEW
                if hit.gene.multi_system:
                    self.multi_syst_genes[hit.gene.name]+=1  
            elif hit.gene.is_forbidden(self.system):
                self.forbidden_genes[hit.gene.name]+=1
                included=False
            else:
                if hit.gene.name in self.exmandatory_genes.keys():
                    self.mandatory_genes[self.exmandatory_genes[hit.gene.name]]+=1
                    valid_hit=validSystemHit(hit, self.system_name, "mandatory")
                    self.valid_hits.append(valid_hit)
                elif hit.gene.name in self.exallowed_genes.keys():
                    self.allowed_genes[self.exallowed_genes[hit.gene.name]]+=1
                    valid_hit=validSystemHit(hit, self.system_name, "allowed")
                    self.valid_hits.append(valid_hit)
                else:
                    msg="Foreign gene %s in cluster %s"%(hit.gene.name, self.system_name)
                    print msg
                    #_log.info(msg)

        if included:
            # Update the number of loci included in the system
            self.nb_cluster += 1            
            # Update the positions of the system
            self.loci_positions.append((cluster.begin, cluster.end))

    def fill_with_hits(self, hits):
        """
        Adds hits to a system occurence, and check which are their status according to the system definition.
        Set the system occurence state to "no_decision" after calling of this function. 
        
        .. note::
            Forbidden genes will only be included if they do belong to the current system (and not to another specified with "system_ref" in the current system's definition). 

        :param hits: a list of Hits to treat for :class:`txsscanlib.search_systems.SystemOccurence` inclusion. 
        :type list of: :class:`txsscanlib.report.Hit`
        """
        included = True
        self._state = "no_decision"
        for hit in hits:
            # Need to check first that this cluster is eligible for system inclusion
            # Stores hits for system extraction (positions, sequences) when complete. 

            if hit.gene.is_mandatory(self.system):
                self.mandatory_genes[hit.gene.name]+=1
                valid_hit=validSystemHit(hit, self.system_name, "mandatory")
                self.valid_hits.append(valid_hit)
            elif hit.gene.is_allowed(self.system):
                self.allowed_genes[hit.gene.name]+=1
                valid_hit=validSystemHit(hit, self.system_name, "allowed")
                self.valid_hits.append(valid_hit)
            elif hit.gene.is_forbidden(self.system):
                self.forbidden_genes[hit.gene.name]+=1
                included=False
            else:
                if hit.gene.name in self.exmandatory_genes.keys():
                    self.mandatory_genes[self.exmandatory_genes[hit.gene.name]]+=1
                    valid_hit=validSystemHit(hit, self.system_name, "mandatory")
                    self.valid_hits.append(valid_hit)
                elif hit.gene.name in self.exallowed_genes.keys():
                    self.allowed_genes[self.exallowed_genes[hit.gene.name]]+=1
                    valid_hit=validSystemHit(hit, self.system_name, "allowed")
                    self.valid_hits.append(valid_hit)
                else:
                    msg="Foreign gene %s in cluster %s"%(hit.gene.name, self.system_name)
                    print msg
                    #_log.info(msg)
                    
    def fill_with_multi_systems_genes(self, multi_systems_hits):
        """
        This function fills the SystemOccurrence with genes putatively coming from other systems (feature "multi_system").
        Those genes are used only if the occurrence of the corresponding gene was not yet filled with a gene from a cluster of the system. 
        
        :param multi_systems_hits: a list of hits of genes that are "multi_system" which correspond to mandatory or allowed genes from the current system for which to fill a SystemOccurrence 
        :type list of: :class:`txsscanlib.report.Hit`
        
        """
        # For each "multi_system" gene missing:
        for g in self.multi_syst_genes:
            if self.multi_syst_genes[g] == 0:
                #multi_systems_hits should be a dico gene.name-wise?
                # We check wether this missing "multi_system" gene was found elsewhere:
                #if g in multi_gene_names in [multi_gene.name for multi_gene in [hit.gene for hit in multi_systems_hits]]:
                if g in [multi_gene.name for multi_gene in [hit.gene for hit in multi_systems_hits]]:
                    # If so, then the SystemOccurrence is filled with this:
                    if g in self.allowed_genes.keys():
                        self.allowed_genes[g]+=1
                        # Add a valid_hit with a special status? e.g "allowed_multi_system"?
                        #self.allowed_genes[hit.gene.name]+=1
                        #valid_hit=validSystemHit(hit, self.system_name, "allowed")
                        #self.valid_hits.append(valid_hit)
                        
                    elif g in self.mandatory_genes.keys():
                        self.mandatory_genes[g]+=1
                        # Add a valid_hit with a special status? e.g "mandatory_multi_system"?
                        #self.mandatory_genes[hit.gene.name]+=1
                        #valid_hit=validSystemHit(hit, self.system_name, "mandatory")
                        #self.valid_hits.append(valid_hit)
                        
                    print "Gene %s supplied from a multi_system gene"%g
        #all_hits = [hit for subl in [report.hits for report in all_reports ] for hit in subl]
    
    
    def decision_rule(self):
        """
        This function applies the decision rules for system assessment in terms of quorum:
            - the absence of forbidden genes is checked
            - the minimal number of mandatory genes is checked (\"min_mandatory_genes_required\")
            - the minimal number of genes in the system is checked (\"min_genes_required\")
            
        When a decision is made, the status (self.status) of the 
        :class:`txsscanlib.search_systems.SystemOccurence` is set either to:
            - "\single_locus\" when a complete system in the form of a single cluster was found
            - "\multi_loci\" when a complete system in the form of several clusters was found
            - "\uncomplete\" when no system was assessed (quorum not reached)
            - "\empty\" when no gene for this system was found
            - "\exclude\" when no system was assessed (at least one forbidden gene was found)            
        """
        nb_forbid = self.count_genes(self.forbidden_genes)
        nb_mandat = self.count_genes(self.mandatory_genes)
        nb_allowed = self.count_genes(self.allowed_genes)
        self._nb_syst_genes = self.compute_nb_syst_genes()

        msg = "====> Decision rule for putative system %s\n"%self.system_name
        msg += str(self)
        msg += "\nnb_forbid : %d\nnb_mandat : %d\nnb_allowed : %d"%(nb_forbid, nb_mandat, nb_allowed)

        if ( nb_forbid == 0 ):
            if (nb_mandat >= self.system.min_mandatory_genes_required) and (self.nb_syst_genes >= self.system.min_genes_required) and (self.nb_syst_genes  >= 1):
                if self.nb_cluster == 1: 
                    self._state = "single_locus"
                else:
                    self._state = "multi_loci"

                msg += "\nYeah complete system \"%s\"."%self.state
                msg += "\n******************************************\n"
                #print msg
                #_log.info(msg)

            elif self.nb_syst_genes > 0:
                msg += "\nuncomplete system."
                msg += "\n******************************************\n"
                #print msg
                #_log.info(msg)
                self._state = "uncomplete"

            else:
                msg += "\nempty system."
                msg += "\n******************************************\n"
                #print msg
                #_log.info(msg)
                self._state = "empty"
        else:
            msg += "\nexclude."
            msg += "\n******************************************\n"
            #print msg
            #_log.info(msg)
            self._state = "exclude"

class validSystemHit(object):
    """
    Encapsulates a :class:`txsscanlib.report.Hit`
    This class stores a Hit that has been attributed to a detected system. Thus, it also stores:  
    
    - the system, 
    - the status of the gene in this system,
    
    It also aims at storing information for results extraction:
    
    - system extraction (e.g. genomic positions)
    - sequence extraction
        
    """
    def __init__(self, hit, detected_system, gene_status):
        self._hit = hit
        self.predicted_system = detected_system
        self.reference_system = hit.system.name
        self.gene_status = gene_status

    def __getattr__(self, attr_name):
        return getattr(self._hit, attr_name)
    
    def __str__(self):
        return "%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%d\t%d\n" % (self.id,
                                                     self.replicon_name,
                                                     self.position,
                                                     self.seq_length,
                                                     self.gene.name,
                                                     self.reference_system,
                                                     self.predicted_system,
                                                     self.gene_status,
                                                     self.i_eval,
                                                     self.score,
                                                     self.profile_coverage,
                                                     self.sequence_coverage,
                                                     self.begin_match,
                                                     self.end_match)

    def output_system(self, system_name, system_status):
        return "%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%d\t%d\n" % (self.id,
                                                     self.replicon_name,
                                                     self.position,
                                                     self.seq_length,
                                                     self.gene.name,
                                                     self.reference_system,
                                                     self.predicted_system,
                                                     system_name,
                                                     system_status,
                                                     self.gene_status,
                                                     self.i_eval,
                                                     self.score,
                                                     self.profile_coverage,
                                                     self.sequence_coverage,
                                                     self.begin_match,
                                                     self.end_match)


    def output_system_header(self):
        return "#Hit_Id\tReplicon_name\tPosition\tSequence_length\tGene\tReference_system\tPredicted_system\tSystem_Id\tSystem_status\tGene_status\ti-evalue\tScore\tProfile_coverage\tSequence_coverage\tBegin_match\tEnd_match\n"


class systemDetectionReport(object):
    """
    Stores the detected systems to report for each replicon: 
        - by system name, 
        - by state of the systems (single vs multi loci)
    
    """
    
    def __init__(self, replicon_name, systems_occurences_list, systems):
        self._systems_occurences_list = systems_occurences_list
        #self._system_textlist = []
        self.replicon_name = replicon_name
            
    def counter_output(self):
        """
        Builds a counter of systems per replicon, with different "states" separated (single-locus vs multi-loci systems)
        """
        system_textlist=[]
        for so in self._systems_occurences_list:
            system_textlist.append(so.system_name+"_"+so.state)
               
        return Counter(system_textlist)
 
    def tabulated_output_header(self, system_occurence_states, system_names):
        """
        Returns a string containing the header of the tabulated output
        """
        # Can be done intra-class 
        header = "#Replicon"
        for syst_name in system_names:
            for state in system_occurence_states:
                header += "\t"+syst_name+"_"+state
    
        header+="\n"        
        
        return header

    def tabulated_output(self, system_occurence_states, system_names, reportfilename, print_header = False):
        """
        Write a tabulated output with number of detected systems for each replicon. 
        """
        system_counter=self.counter_output()
        print system_counter    
        report_str = self.replicon_name
        for s in system_names:
            for o in system_occurence_states:
                index=s+"_"+str(o)
                if system_counter.has_key(index):
                    report_str+="\t"
                    report_str+=str(system_counter[index])
                else:
                    report_str+="\t0"
        report_str+="\n"      
            
        with open(reportfilename, 'a') as _file:
            if print_header:
                _file.write(self.tabulated_output_header(system_occurence_states, system_names))
            _file.write(report_str)
    
           
    def report_output(self, reportfilename, print_header = False):
        """
        Writes a report of sequences forming the detected systems, with information in their status in the system, 
        their localization on replicons, and statistics on the Hits. 
        """
        report_str=""
        for so in self._systems_occurences_list:
            so_unique_name = so.get_system_unique_name(self.replicon_name)
            for hit in so.valid_hits:
                if print_header:
                    report_str+=hit.output_system_header()
                    print_header = False
                report_str+=hit.output_system(so_unique_name, so.state)
        
        with open(reportfilename, 'a') as _file:
            _file.write(report_str)    

    def summary_output(self, reportfilename, rep_info, print_header = False):
        """
        Writes a report with the summary of systems detected in replicons. For each system, a summary is done including: 
                   
            - the number of mandatory/allowed genes in the reference system (as defined in XML files)
            - the number of mandatory/allowed genes detected
            - the number and list of missing genes
            - the number of loci encoding the system
            
        """
        
        report_str = ""
        for so in self._systems_occurences_list:
            if print_header:
                report_str+="%s\n"%so.get_summary_header()
                print_header=False

            report_str+="%s\n"%so.get_summary(self.replicon_name, rep_info)

        with open(reportfilename, 'a') as _file:
            _file.write(report_str)

    def json_output(self, path, rep_db):
        """
        Generates the report in json format

        :param path: the path to a file where to write the report in json format
        :type path: string
        :param rep_db: the replicons database
        :type rep_db: a class:`txsscanlib.database.RepliconDB` object
        """
        with open(path, 'w') as _file:
            all_systems_occurences = []
            for so in self._systems_occurences_list:
                system = {}
                system['name'] = so.unique_name
                system['replicon'] = {}
                system['replicon']['name'] = so.valid_hits[0].replicon_name
                rep_info = rep_db[system['replicon']['name']]
                system['replicon']['length'] = rep_info.max - rep_info.min
                system['replicon']['topology'] = rep_info.topology
                system['genes'] = []
                for valid_hit in so.valid_hits:
                    gene = {}
                    gene['id'] = valid_hit.id
                    gene['position'] = valid_hit.position
                    gene['sequence_length'] = valid_hit.seq_length
                    gene['system'] = valid_hit.reference_system
                    gene['match'] = valid_hit.gene.name
                    gene['gene_status'] = valid_hit.gene_status
                    gene['i_eval'] = valid_hit.i_eval
                    gene['score'] = valid_hit.score
                    gene['profile_coverage'] = valid_hit.profile_coverage
                    gene['sequence_coverage'] = valid_hit.sequence_coverage
                    gene['begin_match'] = valid_hit.begin_match
                    gene['end_match'] = valid_hit.end_match
                    system['genes'].append(gene)
                system['summary'] = {}
                system['summary']['mandatory'] = so.mandatory_genes
                system['summary']['exmandatory_genes'] = so.exmandatory_genes
                system['summary']['allowed'] = so.allowed_genes
                system['summary']['exallowed_genes'] = so.exallowed_genes
                system['summary']['forbiden'] = so.forbidden_genes
                system['summary']['state'] = so._state
                all_systems_occurences.append(system)
            json.dump(all_systems_occurences, _file, indent = 2)

class systemDetectionReportUnordered(object):
    """
    Stores a report for putative detected systems gathering all hits from a search in an unordered dataset: 
        - by system.
    
    Mandatory and allowed genes only are reported in the "json" and "report" output, but all hits matching a system component are reported in the "summary".  
    
    """
    
    def __init__(self, systems_occurences_list, systems):
        self._systems_occurences_list = systems_occurences_list
        #self.replicon_name = replicon_name            

             
    def report_output(self, reportfilename, print_header = False):
        """
        Writes a report of sequences forming the detected systems, with information in their status in the system, 
        their localization on replicons, and statistics on the Hits. 
        """
        report_str=""
        for so in self._systems_occurences_list:
            #so_unique_name = so.get_system_unique_name(self.replicon_name)
            so_unique_name = so.get_system_name_unordered()
            #so_unique_name = so.system_name+"_putative"
            for hit in so.valid_hits:
                if print_header:
                    report_str+=hit.output_system_header()
                    print_header = False
                report_str+=hit.output_system(so_unique_name, so.state)
        
        with open(reportfilename, 'a') as _file:
            _file.write(report_str)    

    def summary_output(self, reportfilename, print_header = False):
        """
        Writes a report with the summary for putative systems in an unordered dataset. For each system, a summary is done including: 
                   
            - the number of mandatory/allowed genes in the reference system (as defined in XML files)
            - the number of mandatory/allowed genes detected
            
        """
        
        report_str = ""
        for so in self._systems_occurences_list:
            if print_header:
                report_str+="%s\n"%so.get_summary_header()
                print_header=False

            #report_str+="%s\n"%so.get_summary(self.replicon_name, rep_info)
            # Get a fake "replicon_name" from the config object in future devt.
            report_str+="%s\n"%so.get_summary_unordered("Unordered")

        with open(reportfilename, 'a') as _file:
            _file.write(report_str)



def disambiguate_cluster_old(cluster):
    """
    This disambiguation step is used on clusters with hits for multiple systems (when cluster.state is set to "ambiguous"). 
    It returns a "cleansed" list of clusters, ready to use for system occurence detection (and that are "clear" cases). It: 

    - splits the cluster in two if it seems that two systems are nearby
    - removes single hits that are not forbidden for the "main" system and that are at one end of the current cluster in this case, check that they are not "loners", cause "loners" can be stored.

    """
    res_clusters = []               
    syst_dico = cluster.systems
    print "Disambiguation step"
    print syst_dico
     
    cur_syst = cluster.hits[0].system.name
    cur_nb_syst_genes_tot = syst_dico[cur_syst]
    cur_nb_syst_genes = 1
             
    #cur_cluster=Cluster() 
    cur_cluster=Cluster(cluster.systems_to_detect) # New
    cur_cluster.add(cluster.hits[0])       
    for h in cluster.hits[1:]:
        syst = h.system.name
        if syst == cur_syst:
            cur_nb_syst_genes+=1
            cur_cluster.add(h)
        else:
            # Deal with "allowed foreign genes", and system attribution when the current gene can be in both aside systems! 
            #if h.gene.is_authorized(cur_syst):
            #    if not h.gene.is_authorized(syst):
            #        cur_cluster.add(h)
            # ==> Part done elsewhere ! now systems of foreign genes are not counted.

            # Case 1: the current gene can not be found in the last system
            if cur_nb_syst_genes == cur_nb_syst_genes_tot:
                cur_cluster.save()
                # Check cluster status before storing it or not:
                if cur_cluster.state == "clear":
                    res_clusters.append(cur_cluster)
           
            cur_syst = syst
            cur_nb_syst_genes = 1
            cur_nb_syst_genes_tot = syst_dico[cur_syst]
            #cur_cluster = Cluster()
            cur_cluster = Cluster(cluster.systems_to_detect) # NEW
            cur_cluster.add(h)
            
    if cur_nb_syst_genes == cur_nb_syst_genes_tot:
        #print cur_nb_syst_genes
        #print cur_nb_syst_genes_tot
        #print cur_syst
        cur_cluster.save()
        # Check cluster status before storing it or not:
        if cur_cluster.state == "clear":
            res_clusters.append(cur_cluster) 
    
    for r in res_clusters:
        print r
                                   
    return res_clusters


def get_compatible_systems(systems_liste1, systems_liste2):
    """
    Returns the intersection of the two input lists.s
    """
    inter=[]
    #print "\nliste 1"
    #print [ system.name for system in systems_liste1]
    #print "liste 2"
    #print [ system.name for system in systems_liste2]
    for el in systems_liste1:
        if el in systems_liste2:
            if not el in inter:
                inter.append(el)
    #print "=> inter"
    #print [ system.name for system in inter]
    return inter


def disambiguate_cluster(cluster):
    """
    This disambiguation step is used on clusters with hits for multiple systems (when cluster.state is set to "ambiguous"). 
    It returns a "cleansed" list of clusters, ready to use for system occurence detection (and that are "clear" cases). It: 

    - splits the cluster in two if it seems that two systems are nearby
    - removes single hits that are not forbidden for the "main" system and that are at one end of the current cluster in this case, check that they are not "loners", cause "loners" can be stored.

    """
    res_clusters = []               
    syst_dico = cluster.systems
    print "Disambiguation step"
    print syst_dico
     
    cur_syst = cluster.hits[0].system.name # to rm
    cur_nb_syst_genes_tot = syst_dico[cur_syst]
    cur_nb_syst_genes = 1
             
    #cur_cluster=Cluster() 
    cur_cluster=Cluster(cluster.systems_to_detect) # New
    cur_cluster.add(cluster.hits[0])
    # Now more complex, deals with compatible systems also for disambiguation.
    cur_compatible=cluster.hits[0].gene.get_compatible_systems(cluster.systems_to_detect)
    
    for h in cluster.hits[1:]:
        print h.gene.name
        compatible_systems=h.gene.get_compatible_systems(cluster.systems_to_detect)
        compat_list=get_compatible_systems(cur_compatible, compatible_systems)

        syst = h.system.name # to rm
        #if syst == cur_syst: # to rm
        if compat_list:
            cur_nb_syst_genes+=1
            cur_cluster.add(h)
            cur_compatible=compat_list
        else:
            # Deal with "allowed foreign genes", and system attribution when the current gene can be in both aside systems! 
            #if h.gene.is_authorized(cur_syst):
            #    if not h.gene.is_authorized(syst):
            #        cur_cluster.add(h)
            # ==> Part done elsewhere ! now systems of foreign genes are not counted.

            # Case 1: the current gene can not be found in the last system
            if cur_nb_syst_genes == cur_nb_syst_genes_tot:
                cur_cluster.save()
                # Check cluster status before storing it or not:
                if cur_cluster.state == "clear":
                    res_clusters.append(cur_cluster)
           
            cur_syst = syst
            cur_compatible=compatible_systems
            #cur_compatible=compat_list
            
            cur_nb_syst_genes = 1
            cur_nb_syst_genes_tot = syst_dico[cur_syst]
            #cur_cluster = Cluster()
            cur_cluster = Cluster(cluster.systems_to_detect) # NEW
            cur_cluster.add(h)
            
    if cur_nb_syst_genes == cur_nb_syst_genes_tot:
        #print cur_nb_syst_genes
        #print cur_nb_syst_genes_tot
        #print cur_syst
        cur_cluster.save()
        # Check cluster status before storing it or not:
        if cur_cluster.state == "clear":
            res_clusters.append(cur_cluster) 
    
    for r in res_clusters:
        print r
                                   
    return res_clusters

    
def analyze_clusters_replicon(clusters, systems, multi_systems_genes):
    """
    Analyzes sets of contiguous hits (clusters) stored in a ClustersHandler for system detection:
        
    - split clusters if needed
    - delete them if they are not relevant
    - add eventual genes from other systems "multi_system" genes
    - check the QUORUM for each system to detect, *i.e.* mandatory + allowed - forbidden
          
    Only for \"ordered\" datasets representing a whole replicon. 
    Reports systems occurence. 
    
    :param clusters: the set of clusters to analyze
    :type clusters: :class:`txsscanlib.search_systems.ClustersHandler` 
    :param systems: the set of systems to detect
    :type systems: a list of :class:`txsscanlib.system.System`
    :param multi_systems_genes: a dictionary with genes that could belong to multiple systems (keys are system names)
    :return: a set of systems occurence filled with hits found in clusters
    :rtype: a list of :class:`txsscanlib.search_systems.SystemOccurence` 
    
    """
    
    # Global Hits collectors, for uncomplete cluster Hits
    systems_occurences_scattered = {}
    systems_occurences_list = []
    
    syst_dict = {}
    for system in systems:
        syst_dict[system.name] = system
        systems_occurences_scattered[system.name] = SystemOccurence(system)
    
    for clust in clusters.clusters:
        print "\n%s"%str(clust)
        #if clust.state == "clear": 
        systems_to_consider =  get_compatible_systems([system_bank[s] for s in clust.compatible_systems], clust.systems_to_detect)
        if clust.state == "clear" and len(systems_to_consider) > 0:           
            # Local Hits collector
            # Check the putative system belongs to the list of systems to detect !! If it does not, do not go further with this cluster of genes.
            # New! different compatible systems are tested: then update cluster._putative_system w the good one?
            #print clust.compatible_systems
            # Arbitratily, if none of the set of compatible_systems pass the decision rule step, then the 1st system will store a scattered version of this...
            first=True
            #store_scattered=True
            store_scattered=False
            store_clust=None
            store_so=None
            exclude=False
            #for putative_system in clust.compatible_systems:
            for putative_system in [s.name for s in systems_to_consider]:
                print "Considering %s - compatible system"%putative_system
                if putative_system in syst_dict.keys():
                    so = SystemOccurence(syst_dict[putative_system])
                    so.fill_with_cluster(clust)
                    # NEW!
                    if putative_system in multi_systems_genes.keys():
                        so.fill_with_multi_systems_genes(multi_systems_genes[putative_system])
                    so.decision_rule()
                    so_state = so.state
                    print so_state
                    if so_state != "exclude":
                        if so_state != "single_locus":            
                            # Store it to pool genes found with genes from other clusters.
                            # Do not do it if the so has a forbidden gene !!!
                            # NEW !! Now do not do this at the 1st try ! only if not complete is stored in the loop ! 
                            if first:
                                store_clust=clust
                                store_so=so
                                store_scattered=True
                            #print "...\nStored for later treatment of scattered systems.\n"
                            #systems_occurences_scattered[putative_system].fill_with_cluster(clust)
                        else:
                            print so 
                            #print so_state
                            print "...\nComplete %s %s system stored.\n"%(putative_system, so_state)
                            systems_occurences_list.append(so)
                            store_scattered=False
                            break
                    else:
                        exclude=True
                        #store_scattered=True
            #if store_scattered and not exclude:
            if store_scattered:
                print store_so
                print "...\nPutative %s locus stored for later treatment of scattered systems.\n"%clust.compatible_systems[0]
                systems_occurences_scattered[clust.compatible_systems[0]].fill_with_cluster(store_clust)
                
            """ OLD CODE
            if clust.putative_system in syst_dict.keys():
                so = SystemOccurence(syst_dict[clust.putative_system])
                so.fill_with_cluster(clust)
                # NEW!
                if clust.putative_system in multi_systems_genes.keys():
                    so.fill_with_multi_systems_genes(multi_systems_genes[clust.putative_system])
                so.decision_rule()
                so_state = so.state
                if so_state != "exclude":
                    if so_state != "single_locus":            
                        # Store it to pool genes found with genes from other clusters.
                        # Do not do it if the so has a forbidden gene !!!
                        print "...\nStored for later treatment of scattered systems.\n"
                        systems_occurences_scattered[clust.putative_system].fill_with_cluster(clust)
                    else: 
                        print "...\nComplete system stored.\n"
                        systems_occurences_list.append(so)
            """
                    
        elif clust.state == "ambiguous":
            # Implement a way to "clean" the clusters. For instance :
            # - split the cluster in two if it seems that two systems are nearby
            # - remove single hits that are not forbidden for the "main" system and that are at one end of the current cluster
            # in this case, check that they are not "loners", cause "loners" can be stored.
            disamb_clusters = disambiguate_cluster(clust)
            # Add those new clusters to the set of clusters to fill system_occurrences?
            for c in disamb_clusters:
                clusters.add(c)
            if disamb_clusters:
                print "=> disambiguated cluster(s) stored for later treatment"
            else:
                print "=> none kept"
                
        else:
            print "------- next -------"
    print "\n\n***************************************************\n******* Report scattered/uncomplete systems *******\n***************************************************\n"
    for system in systems:
        #print systems_occurences[system]
        so = systems_occurences_scattered[system.name]
        so.decision_rule()
        so_state=so.state
        if so.is_complete():
            systems_occurences_list.append(so)
    print "******************************************\n"
    
    # Stores results in this list? Or code a new object : systemDetectionReport ? 
    return systems_occurences_list



#def build_clusters(hits, rep_info):
def build_clusters(hits, systems_to_detect, rep_info):
    """
    Gets sets of contiguous hits according to the minimal inter_gene_max_space between two genes. Only for \"ordered\" datasets.
     
    :param hits: a list of Hmmer hits to analyze 
    :type hits: a list of :class:`txsscanlib.report.Hit`
    :param systems_to_detect: the list of systems to detect 
    :type systems_to_detect: a list of :class:`txsscanlib.system.System`
    :param cfg: the configuration object built from default and user parameters.
    :type cfg: :class:`txsscanlib.config.Config`
    :return: a set of clusters and a dictionary with \"multi_system\" genes stored in a system-wise way for further utilization.
    :rtype: :class:`txsscanlib.search_systems.ClustersHandler`
    """    
    
    _log.debug("Starting cluster detection with build_clusters... ")
    
    # Deals with different dataset types using Pipeline ?? 
    clusters = ClustersHandler()
    prev = hits[0]
    #cur_cluster = Cluster()
    cur_cluster = Cluster(systems_to_detect)
    positions = []
    loner_state=False
    
    # New: storage of multi_system genes:
    multi_system_genes_system_wise={}
    
    tmp=""
    for cur in hits[1:]:
        
        _log.debug("Hit %s"%str(cur))
        prev_max_dist = prev.get_syst_inter_gene_max_space()
        cur_max_dist = cur.get_syst_inter_gene_max_space()
        inter_gene = cur.get_position() - prev.get_position() - 1
        
        tmp="\n****\n"
        tmp+="prev_max_dist : %d\n"%(prev_max_dist)
        tmp+="cur_max_dist : %d\n"%(cur_max_dist)
        tmp+="Intergene space : %d\n"%(inter_gene)
        tmp+="Cur : %s"%cur
        tmp+="Prev : %s"%prev
        tmp+="Len cluster: %d\n"%len(cur_cluster)
        #print tmp
        
        # First condition removes duplicates (hits for the same sequence)
        # the two others takes into account either system1 parameter or system2 parameter

        #smaller_dist = min(prev_max_dist, cur_max_dist)    
        if(inter_gene <= prev_max_dist or inter_gene <= cur_max_dist ):
        #if(inter_gene <= smaller_dist ):
            #print "zero"
            if positions.count(prev.position) == 0:
                #print "un - ADD prev in cur_cluster"
                cur_cluster.add(prev)
                positions.append(prev.position)
                # New : Storage of multi_system genes:
                if prev.gene.multi_system:
                    if not prev.system.name in multi_system_genes_system_wise.keys():
                        multi_system_genes_system_wise[prev.system.name]=[]
                    multi_system_genes_system_wise[prev.system.name].append(prev)
                
            if positions.count(cur.position) == 0:
                #print "deux - ADD cur in cur_cluster"
                cur_cluster.add(cur)
                positions.append(cur.position)
                # New : Storage of multi_system genes:
                if cur.gene.multi_system:
                    if not cur.system.name in multi_system_genes_system_wise.keys():
                        multi_system_genes_system_wise[cur.system.name]=[]
                    multi_system_genes_system_wise[cur.system.name].append(cur)
                
            if prev.gene.loner:
                #print "trois - loner_state"
                #print "--- PREVLONER %s %s"%(prev.id, prev.gene.name)
                loner_state = True
        else:
            # Storage of the previous cluster
            if len(cur_cluster)>1:
                #print cur_cluster
                #print "quatre - ADD cur_cluster"
                clusters.add(cur_cluster) # Add an in-depth copy of the object? 
                #print(cur_cluster)
                 
                #cur_cluster = Cluster()
                cur_cluster = Cluster(systems_to_detect)
                loner_state = False
                
            elif len(cur_cluster) == 1 and loner_state == True: # WTF?
                #print cur_cluster
                #print "cinq - ADD cur_cluster"
                #print "PREVLONER %s %s"%(prev.id, prev.gene.name)
                clusters.add(cur_cluster) # Add an in-depth copy of the object? 
                #print(cur_cluster)
                     
                #cur_cluster = Cluster()     
                cur_cluster = Cluster(systems_to_detect)
                loner_state = False
            
            if prev.gene.loner:
                #print "six - check"
                #print "PREVLONER ?? %s %s"%(prev.id, prev.gene.name)
                
                if positions.count(prev.position) == 0:
                    #print "six - ADD prev in cur_cluster, ADD cur_cluster"
                    cur_cluster.add(prev) 
                    clusters.add(cur_cluster) 
                    #print(cur_cluster)
                
                    # New : Storage of multi_system genes:
                    if prev.gene.multi_system:
                        if not prev.system.name in multi_system_genes_system_wise.keys():
                            multi_system_genes_system_wise[prev.system.name]=[]
                        multi_system_genes_system_wise[prev.system.name].append(prev)
                    
                    positions.append(prev.position)
                    loner_state = False  
                    #cur_cluster = Cluster() 
                    cur_cluster = Cluster(systems_to_detect) 
                
            #cur_cluster = Cluster() 
            cur_cluster = Cluster(systems_to_detect)     
            
        prev=cur
    
    if len(cur_cluster)>1 or (len(cur_cluster)==1 and prev.gene.loner):
        clusters.add(cur_cluster)
    
    if rep_info.topology == "circular":
        #print "Circ TEST"
        clusters.circularize(rep_info)
        
    return (clusters,multi_system_genes_system_wise)

def get_best_hits(hits, tosort=False, criterion="score"):
    """
    Returns from a putatively redundant list of hits a list of best matching hits.
    Analyzes quorum and co-localization if required for system detection. 
    By default, hits are already sorted by position, and the hit with the best score is kept. Possible criteria are:
        
    - maximal score (criterion=\"score\")
    - minimal i-evalue (criterion=\"i_eval\")
    - maximal percentage of the profile covered by the alignment with the query sequence (criterion=\"profile_coverage\")
    
    """
    if tosort:
        hits = sorted(hits, key=attrgetter('position'))
    best_hits=[]

    prev_hit=hits[0]
    prev_pos=prev_hit.get_position()
    #print hits[0]
    for h in hits[1:]:
        pos=h.get_position()
        if pos !=prev_pos:
            best_hits.append(prev_hit)
            #print "******* no comp ****"
            #print prev_hit
            #print "******* ****** ****"
            prev_hit=h
            prev_pos=pos
        else:
            #print "******* COMP ****"
            #print h 
            #print prev_hit
            if criterion == "score":
                if prev_hit.score<h.score:
                    prev_hit=h
            elif criterion == "i_eval":
                if getattr(prev_hit, 'i_eval') > getattr(h, 'i_eval'):
                    prev_hit=h
            elif criterion == "profile_coverage":
                if getattr(prev_hit, 'profile_coverage') < getattr(h, 'profile_coverage'):
                    prev_hit=h
            else:
                raise TxsscanError("The criterion for Hits comparison % does not exist or is not available. \nIt must be either \"score\", \"i_eval\" or \"profile_coverage\"."%criterion)
            
            #print "BEST"
            #print prev_hit
            #print "******* ****** ****"
    
    best_hits.append(prev_hit)
    
    #settupkes=[(h, pos) for h in hits for pos in [getattr(z, 'position') for z in hits]]
        
    return best_hits


 
def search_systems(hits, systems, cfg):
    """
    Runs search of systems from a set of hits. Criteria for system assessment will depend on the kind of input dataset provided: 
    
      - analyze **quorum and co-localization** for "ordered_replicon" and "gembase" datasets.
      - analyze **quorum only** (and in a limited way) for "unordered_replicon" and "unordered" datasets.
    
    :param hits: the list of hits for input systems components
    :type hits: list of :class:`txsscanlib.report.Hit`
    :param systems: the list of systems asked for detection
    :type systems: list of :class:`txsscanlib.system.System`
    :param cfg: the configuration object
    :type cfg: :class:`txsscanlib.config.Config`
    
    """
    
    tabfilename = os.path.join(cfg.working_dir, 'txsscan.tab')
    reportfilename = os.path.join(cfg.working_dir, 'txsscan.report')
    summaryfilename = os.path.join(cfg.working_dir, 'txsscan.summary')
    json_filename = os.path.join(cfg.working_dir, 'txsscan.json')
    
    # For the headers of the output files: no report so far ! print them in the loop at the 1st round ! 
    system_occurences_states = ['single_locus', 'multi_loci']
    system_names = []
    for s in systems:
        syst_name = s.name
        system_names.append(syst_name)
    
    # Specify to build_clusters the rep_info (min, max positions), and replicon_type... 
    # Test with psae_circular_test.prt: pos_min = 1 , pos_max = 5569
    #RepInfo= namedtuple('RepInfo', ['topology', 'min', 'max'])
    #rep_info=RepInfo("circular", 1, 5569)
    
    header_print = True
    if cfg.db_type == 'gembase':
        # Construction of the replicon database storing info on replicons: 
        rep_db = RepliconDB(cfg)
        
        # Use of the groupby() function from itertools : allows to group Hits by replicon_name, 
        # and then apply the same build_clusters functions to replicons from "gembase" and "ordered_replicon" types of databases.
        for k, g in itertools.groupby(hits, operator.attrgetter('replicon_name')):
            sub_hits=list(g)
            rep_info=rep_db[k]
            print rep_info
            
            # The following applies to any "replicon"
            #print "\n************\nBuilding clusters for %s \n************\n"%k
            #(clusters, multi_syst_genes)=build_clusters(sub_hits, rep_info)         
            (clusters, multi_syst_genes)=build_clusters(sub_hits, systems, rep_info)          
            print "\n************************************\n Analyzing clusters for %s \n************************************\n"%k
            # Make analyze_clusters_replicon return an object systemOccurenceReport?
            # Note: at this stage, ther is no control of which systems are looked for... But systemsOccurrence do not have to be created for systems not searched. 
            # 
            #systems_occurences_list = analyze_clusters_replicon(clusters, systems)
            systems_occurences_list = analyze_clusters_replicon(clusters, systems, multi_syst_genes)  
            
            print "******************************************"
            #print "Reporting systems for %s : \n"%k
            print " Building reports for %s: \n"%k
            report = systemDetectionReport(k, systems_occurences_list, systems)
                
            # TO DO: Add replicons with no hits in tabulated_output!!! But where?! No trace of these replicons as replicons are taken from hits. 
            report.tabulated_output(system_occurences_states, system_names, tabfilename, header_print)
            report.report_output(reportfilename, header_print)
            report.summary_output(summaryfilename, rep_info, header_print)
            #report.json_output(json_filename, rep_db)
            print "******************************************"
            
            header_print = False
            
    elif cfg.db_type == 'ordered_replicon':
        # Basically the same as for 'gembase' (except the loop on replicons)
        rep_db = RepliconDB(cfg)
        rep_info = rep_db[RepliconDB.ordered_replicon_name]
        
        #(clusters, multi_syst_genes)=build_clusters(hits, rep_info) 
        (clusters, multi_syst_genes)=build_clusters(hits, systems, rep_info) 
        #for syst in multi_syst_genes:
        #    for g in multi_syst_genes[syst]:
        #        print g
        print "\n************************************\n Analyzing clusters \n************************************\n"
        #systems_occurences_list = analyze_clusters_replicon(clusters, systems)              
        systems_occurences_list = analyze_clusters_replicon(clusters, systems, multi_syst_genes)                    
        print "******************************************"
        #print "Reporting detected systems : \n"
        print " Building reports of detected systems\n "
        report = systemDetectionReport(RepliconDB.ordered_replicon_name, systems_occurences_list, systems)            
        report.tabulated_output(system_occurences_states, system_names, tabfilename, header_print)
        report.report_output(reportfilename, header_print)
        report.summary_output(summaryfilename, rep_info, header_print)
        #report.json_output(json_filename, rep_db)
        print "******************************************"
    
    elif cfg.db_type == 'unordered_replicon' or cfg.db_type == 'unordered':
        
        # implement a new function "analyze_cluster" => Fills a systemOccurence per system
        systems_occurences_list=[]
        # Hits with best score are first selected. 
        hits=get_best_hits(hits, True)
        # Then system-wise treatment:
        hits=sorted(hits, key=attrgetter('system'))
        for k, g in itertools.groupby(hits, operator.attrgetter('system')):
            if k in systems:
                sub_hits=list(g)
                so=SystemOccurence(k)
                #resy=so.fill_with_hits(sub_hits) # does not return anything
                so.fill_with_hits(sub_hits)
                print "******************************************"
                print k.name
                print "******************************************"
                print so
                systems_occurences_list.append(so)
        print "******************************************"
        print " Building reports of detected systems "
        report = systemDetectionReportUnordered(systems_occurences_list, systems)
        report.report_output(reportfilename, header_print)
        report.summary_output(summaryfilename, header_print)
        #report.json_output(json_filename)
        print "******************************************"
       
    else:
        raise ValueError("Invalid database type. ")





