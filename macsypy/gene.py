# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Néron                                         #
# Copyright © 2014  Institut Pasteur (Paris) and CNRS.                         #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################


import os
import logging
_log = logging.getLogger(__name__)

from subprocess import Popen
from threading import Lock
from .report import GembaseHMMReport, GeneralHMMReport, OrderedHMMReport
from .error import MacsypyError
from . import registries

class GeneBank(object):
    """
    Store all Gene objects. Ensure that genes are instanciated only once.
    """

    def __init__(self):
        self._genes_bank = {}

    def __getitem__(self, key):
        """
        :param key: The key to retrieve a gene.
        The key is composed of the name of models family and the gene name.
        for instance CRISPR-Cas/cas9_TypeIIB or TXSS/T6SS_tssH
        :type name: tuple (string, string)
        :return: return the Gene corresponding to the key.
        :rtype: :class:`macsypy.gene.Gene` object
        :raise KeyError: if the key does not exist in GeneBank.
        """
        try:
            return self._genes_bank[key]
        except KeyError:
            raise KeyError("No such gene {} in this bank".format(key))


    def __contains__(self, gene):
        """
        Implement the membership test operator

        :param gene: the gene to test
        :type gene: :class:`macsypy.gene.Gene` object
        :return: True if the gene is in, False otherwise
        :rtype: boolean
        """
        return gene in list(self._genes_bank.values())


    def __iter__(self):
        """
        Return an iterator object on the genes contained in the bank
        """
        return iter(self._genes_bank.values())


    def add_gene(self, gene):
        """
        Add a gene in the bank

        :param gene: the gene to add
        :type gene: :class:`macsypy.gene.Gene` object
        :raise: KeyError if a gene with the same name is already registered
        """
        model_name = registries.split_def_name(gene.model.fqn)[0]
        key = (model_name, gene.name)
        if key in self._genes_bank:
            raise KeyError("a gene named '{0}/{1}' is already registered".format(model_name, gene.name))
        else:
            self._genes_bank[key] = gene


gene_bank = GeneBank()


class Gene(object):
    """
    Handle Gene of a (secretion) System

    """

    def __init__(self, cfg, name,
                 model,
                 model_location,
                 loner=False,
                 exchangeable=False,
                 multi_system=False,
                 inter_gene_max_space=None):
        """
        handle gene

        :param cfg: the configuration object. 
        :type cfg: :class:`macsypy.config.Config` object
        :param name: the name of the Gene.
        :type name: string.
        :param model: the model that owns this Gene
        :type model: :class:`macsypy.model.Model` object.
        :param model_loc: where all the paths profiles and definitions are register for a kind of model.
        :type model_loc: :class:`macsypy.registries.ModelLocation` object.
        :param loner: True if the Gene can be isolated on the genome (with no contiguous genes), False otherwise.
        :type loner: boolean.
        :param exchangeable: True if this Gene can be replaced with one of its homologs or analogs
          without any effects on the model assessment, False otherwise.
        :type exchangeable: boolean.
        :param multi_system: True if this Gene can belong to different occurrences of this System. 
        :type multi_system: boolean.
        :param inter_gene_max_space: the maximum space between this Gene and another gene of the System.
        :type inter_gene_max_space: integer
        """
        self.name = name 
        self.profile = profile_factory.get_profile(self, cfg, model_location)
        """:ivar profile: The HMM protein Profile corresponding to this gene :class:`macsypy.gene.Profile` object"""

        self.homologs = []
        self.analogs = []
        self._model = model
        self._loner = loner
        self._exchangeable = exchangeable
        self._multi_system = multi_system
        self._inter_gene_max_space = inter_gene_max_space


    def __str__(self):
        """
        Print the name of the gene and of its homologs/analogs.
        """
        s = "name : {0}".format(self.name)
        s += "\ninter_gene_max_space: {:d}".format(self.inter_gene_max_space)
        if self.loner:
            s += "\nloner"
        if self.multi_system:
            s += "\nmulti_system"
        if self.exchangeable:
            s += "\nexchangeable"
        if self.homologs:
            s += "\n    homologs: "
            for h in self.homologs:
                s += h.name + ", "
            s = s[:-2]
        if self.analogs:
            s += "\n    analogs: "
            for a in self.analogs:
                s += a.name + ", "
            s = s[:-2]
        return s


    @property
    def model(self):
        """
        :return: the Model that owns this Gene
        :rtype: :class:`macsypy.model.Model` object
        """
        return self._model


    @property
    def loner(self):
        """
        :return: True if the gene can be isolated on the genome, False otherwise
        :rtype: boolean
        """
        return self._loner


    @property
    def exchangeable(self):
        """
        :return: True if this gene can be replaced with one of its homologs or analogs whithout any effects on the model, False otherwise.
        :rtype: boolean.
        """
        return self._exchangeable


    @property
    def multi_system(self):
        """
        :return: True if this Gene can belong to different occurrences of **this System** (and can be used for multiple System assessments), False otherwise.
        :rtype: boolean.
        """
        return self._multi_system


    @property
    def inter_gene_max_space(self):
        """
        :return: The maximum distance allowed between this gene and another gene for them to be considered co-localized. 
                 If the value is not set at the Gene level, return the value set at the System level.
        :rtype: integer.
        """
        if self._inter_gene_max_space is not None:
            return self._inter_gene_max_space
        else:
            return self._model.inter_gene_max_space


    def add_homolog(self, homolog):
        """
        Add a homolog gene to the Gene

        :param homolog: homolog to add
        :type homolog:  :class:`macsypy.gene.Homolog` object 
        """
        self.homologs.append(homolog)


    def get_homologs(self):
        """
        :return: the Gene homologs
        :type: list of :class:`macsypy.gene.Homolog` object
        """
        return self.homologs


    def add_analog(self, analog):
        """
        Add an analogous gene to the Gene

        :param analog: analog to add
        :type analog:  :class:`macsypy.gene.Analog` object 
        """
        self.analogs.append(analog)


    def get_analogs(self):
        """
        :return: the Gene analogs
        :type: list of :class:`macsypy.gene.Analog` object
        """
        return self.analogs


    def __eq__(self, gene):
        """
        :return: True if the gene names (gene.name) are the same, False otherwise.
        :param gene: the query of the test
        :type gene: :class:`macsypy.gene.Gene` object.
        :rtype: boolean.
        """
        return self.name == gene.name


    def __hash__(self):
        # needed to be hashable in Py3 when __eq__ is defined
        # see https://stackoverflow.com/questions/1608842/types-that-define-eq-are-unhashable  
        
        return id(self)


    def is_homolog(self, gene):
        """
        :return: True if the two genes are homologs, False otherwise.
        :param gene: the query of the test
        :type gene: :class:`macsypy.gene.Gene` object.
        :rtype: boolean.
        """

        if self == gene:
            return True
        else:
            for h in self.homologs:
                if gene == h.gene:
                    return True
        return False


    def is_analog(self, gene):
        """
        :return: True if the two genes are analogs, False otherwise.
        :param gene: the query of the test
        :type gene: :class:`macsypy.gene.Gene` object.
        :rtype: boolean.
        """

        if self == gene:
            return True
        else:
            for h in self.analogs:
                if gene == h.gene:
                    return True
        return False


    def is_mandatory(self, model):
        """
        :return: True if the gene is within the *mandatory* genes of the model, False otherwise.
        :param model: the query of the test
        :type model: :class:`macsypy.model.Model` object.
        :rtype: boolean.
        """
        if self in model.mandatory_genes:
            return True
        else:
            return False


    def is_accessory(self, model):
        """
        :return: True if the gene is within the *accessory* genes of the model, False otherwise.
        :param model: the query of the test
        :type model: :class:`macsypy.model.Model` object.
        :rtype: boolean.
        """
        if self in model.accessory_genes:
            return True
        else:
            return False


    def is_forbidden(self, model):
        """
        :return: True if the gene is within the *forbidden* genes of the model, False otherwise.
        :param model: the query of the test
        :type model: :class:`macsypy.model.Model` object.
        :rtype: boolean.
        """
        if self in model.forbidden_genes:
            return True
        else:
            return False


    def is_authorized(self, model, include_forbidden=True):
        """
        :return: True if this gene is found in the Model, False otherwise.
        :param model: the query of the test
        :type model: :class:`macsypy.model.Model` object.
        :param include_forbidden: tells if forbidden genes should be considered as "authorized" or not
        :type include_forbidden: boolean
        :rtype: boolean.
        """
        genes = model.mandatory_genes + model.accessory_genes
        if include_forbidden:
            genes = genes + model.forbidden_genes
        for g in genes:
            if self == g:
                return True
            if g.exchangeable and (g.is_homolog(self) or g.is_analog(self)):
                return True
        return False


    def get_compatible_models(self, model_list, include_forbidden=True):
        """
        Test every model in model_list for compatibility with the gene using the is_authorized function.

        :param model_list: a list of model names to test
        :type model_list: list of strings
        :param include_forbidden: tells if forbidden genes should be considered as defining a compatible models or not
        :type include_forbidden: boolean
        :return: the list of compatible models
        :rtype: list of :class:`macsypy.model.Model` objects, or void list if none compatible
        """
        compatibles = [model for model in model_list if self.is_authorized(model, include_forbidden=include_forbidden)]
        return compatibles


class Homolog(object):
    """
    Handle homologs, encapsulate a Gene
    """

    def __init__(self, gene, gene_ref, aligned=False):
        """
        :param gene: the gene
        :type gene: :class:`macsypy.gene.Gene` object.
        :param gene_ref: the gene to which the current is homolog.
        :type gene_ref: :class:`macsypy.gene.Gene` object.
        :param aligned: if True, the profile of this gene overlaps totally the sequence of the reference gene profile.
                        Otherwise, only partial overlapping between the profiles.
        :type aligned: boolean
        """
        self.gene = gene 
        """:ivar gene: gene """

        self.ref = gene_ref 
        self.aligned = aligned

    def __getattr__(self, name):
        return getattr(self.gene, name)


    def is_aligned(self):
        """
        :return: True if this gene homolog is aligned to its homolog, False otherwise.
        :rtype: boolean
        """
        return self.aligned


    @property
    def gene_ref(self):
        """
        :return: the gene to which this one is homolog to (reference gene)
        :rtype: :class:`macsypy.gene.Gene` object
        """
        return self.ref


class Analog(object):
    """
    Handle analogs, encapsulate a Gene
    """

    def __init__(self, gene, gene_ref):
        """
        :param gene: the gene
        :type gene: :class:`macsypy.gene.Gene` object.
        :param gene_ref: the gene to which the current is analog.
        :type gene_ref: :class:`macsypy.gene.Gene` object.
        """
        self.gene = gene 
        """:ivar gene: gene """

        self.ref = gene_ref 


    def __getattr__(self, name):
        return getattr(self.gene, name)


    @property
    def gene_ref(self):
        """
        :return: the gene to which this one is analog to (reference gene)
        :rtype: :class:`macsypy.gene.Gene` object
        """
        return self.ref


class ProfileFactory(object):
    """
    Build and store all Profile objects. Profiles must not be instanciated directly.
    The profile_factory must be used. The profile_factory ensures there is only one instance
    of profile for a given name.
    To get a profile, use the method get_profile. If the profile is already cached, this instance is returned.
    Otherwise a new profile is built, stored in the profile_factory and then returned.

    """
    _profiles = {}

    def get_profile(self, gene, cfg, model_location):
        """
        :param gene: the gene associated to this profile
        :type gene: :class:`macsypy.gene.Gene` or :class:`macsypy.gene.Homolog` or :class:`macsypy.gene.Analog` object
        :param profiles_registry: the registry where are stored the path of the profiles
        :type profiles_registry: the registry of profiles
        :param profiles_registry: :class:`macsypy.registries.ProfilesRegistry` instance.
        :return: the profile corresponding to the name.
                 If the profile already exists, return it. Otherwise build it, store it and return it.
        :rtype: :class:`macsypy.gene.Profile` object
        """
        key = (model_location.name, gene.name)
        if key in self._profiles:
            profile = self._profiles[key]
        else:
            try:
                path = model_location.get_profile(gene.name)
            except KeyError:
                raise MacsypyError("{0}: No such profile".format(gene.name))
            profile = Profile(gene, cfg, path)
            self._profiles[key] = profile
        return profile


profile_factory = ProfileFactory()


class Profile(object):
    """
    Handle a HMM protein profile
    """

    def __init__(self, gene, cfg, path):
        """

        :param gene: the gene corresponding to this profile
        :type gene_name: :class:`macsypy.secretion.Gene` object
        :param cfg: the configuration 
        :type cfg: :class:`macsypy.config.Config` object
        :param path: the path to the hmm profile.
        :type path: string
        """
        self.gene = gene 
        self.path = path
        self.len = self._len()
        self.cfg = cfg 
        self.hmm_raw_output = None
        self._report = None
        self._lock = Lock()


    def __len__(self):
        """
        :return: the length of the HMM protein profile
        :rtype: int
        """
        return self.len


    def _len(self):
        """
        Parse the HMM profile file to get and store the length.
        This private method is called at the Profile init.
        """
        length = None
        with open(self.path) as f:
            for l in f:
                if l.startswith("LENG"):
                    length = int(l.split()[1])
                    break
        return length


    def __str__(self):
        """
        Print the name of the corresponding gene and the path to the HMM profile.
        """
        return "{0} : {1}".format(self.gene.name, self.path)


    def execute(self):
        """
        Launch the Hmmer search (hmmsearch executable) with this profile

        :return: an object storing information on the results of the HMM search (HMMReport)
        :rtype:  :class:`macsypy.report.HMMReport` object
        """
        with self._lock:
            # the results of HMM is cached 
            # so HMMsearch is executed only once per run
            # if this method is called several times,
            # the first call induce the execution of HMMsearch and generate a report
            # the other calls return directly this report
            if self._report is not None:
                return self._report
            hmmer_dir = os.path.join(self.cfg.working_dir(), self.cfg.hmmer_dir())
            if not os.path.exists(hmmer_dir):
                os.mkdir(hmmer_dir)
            output_path = os.path.join(hmmer_dir,  self.gene.name + self.cfg.res_search_suffix())
            err_path = os.path.join(hmmer_dir,
                                    self.gene.name + os.path.splitext(self.cfg.res_search_suffix())[0] + ".err")

            with open(err_path, 'w') as err_file:
                options = {"hmmer_exe": self.cfg.hmmer(),
                           "output_file": output_path,
                           "e_value_res": self.cfg.e_value_search(),
                           "profile": self.path,
                           "sequence_db": self.cfg.sequence_db(),
                           }
                command = "{hmmer_exe} --cpu 0 -o {output_file} -E {e_value_res:f} {profile} {sequence_db} ".format(**options)
                _log.debug("{0} Hmmer command line : {1}".format(self.gene.name, command))
                try:
                    hmmer = Popen(command,
                                  shell=True,
                                  stdout=None,
                                  stdin=None,
                                  stderr=err_file,
                                  close_fds=False,
                                 )
                except Exception as err:
                    msg = "Hmmer execution failed: command = {0} : {1}".format(command, err)
                    _log.critical(msg, exc_info=True)
                    raise err
                hmmer.wait()

            if hmmer.returncode != 0:
                if hmmer.returncode == -15:
                    msg = "The Hmmer execution was aborted: command = {0} : return code = {1:d} check {2}".format(command, hmmer.returncode, err_path)
                    _log.critical(msg)
                    return
                else:
                    msg = "an error occurred during Hmmer execution: command = {0} : return code = {1:d} check {2}".format(command, hmmer.returncode, err_path)
                    _log.debug(msg, exc_info=True)
                    _log.critical(msg)
                    raise RuntimeError(msg)
            self.hmm_raw_output = output_path
            db_type = self.cfg.db_type()
            if db_type == 'gembase':
                report = GembaseHMMReport(self.gene, output_path, self.cfg)
            elif db_type == 'ordered_replicon':
                report = OrderedHMMReport(self.gene, output_path, self.cfg)
            else:
                report = GeneralHMMReport(self.gene, output_path, self.cfg)
            self._report = report
            return report


