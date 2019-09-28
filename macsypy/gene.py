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


import os
import logging
_log = logging.getLogger(__name__)
from enum import Enum

from . import registries


class GeneBank:
    """
    Store all Gene objects. Ensure that genes are instanciated only once.
    """

    def __init__(self):
        self._genes_bank = {}

    def __getitem__(self, key):
        """
        :param key: The key to retrieve a gene.
                    The key is composed of the name of models family and the gene name.
                    for instance CRISPR-Cas/cas9_TypeIIB ('CRISPR-Cas' , 'cas9_TypeIIB') or
                    TXSS/T6SS_tssH ('TXSS', 'T6SS_tssH')
        :type key: tuple (string, string)
        :return: return the Gene corresponding to the key.
        :rtype: :class:`macsypy.gene.Gene` object
        :raise KeyError: if the key does not exist in GeneBank.
        """
        try:
            return self._genes_bank[key]
        except KeyError:
            raise KeyError(f"No such gene {key} in this bank")


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
            raise KeyError(f"a gene named '{model_name}/{gene.name}' is already registered")
        else:
            self._genes_bank[key] = gene


class Gene:
    """
    Handle Gene of a (secretion) System

    """

    def __init__(self, profile_factory, name, model, model_location, loner=False, exchangeable=False,
                 multi_system=False, inter_gene_max_space=None):
        """
        handle gene

        :param name: the name of the Gene.
        :type name: string.
        :param model: the model that owns this Gene
        :type model: :class:`macsypy.model.Model` object.
        :param model_location: where all the paths profiles and definitions are register for a kind of model.
        :type model_location: :class:`macsypy.registries.ModelLocation` object.
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
        self.profile = profile_factory.get_profile(self, model_location)
        """:ivar profile: The HMM protein Profile corresponding to this gene :class:`macsypy.profile.Profile` object"""

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
        s = f"name : {self.name}"
        s += f"\ninter_gene_max_space: {self.inter_gene_max_space:d}"
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
        :return: True if this gene can be replaced with one of its homologs or analogs without any effects on the model,
                 False otherwise.
        :rtype: boolean.
        """
        return self._exchangeable


    @property
    def multi_system(self):
        """
        :return: True if this Gene can belong to different occurrences of **the model**
                (and can be used for multiple System assessments), False otherwise.
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


class Homolog:
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
        self._gene = gene
        self._ref = gene_ref
        self.aligned = aligned

    def __getattr__(self, name):
        return getattr(self._gene, name)


    def is_aligned(self):
        """
        :return: True if this gene homolog is aligned to its homolog, False otherwise.
        :rtype: boolean
        """
        return self.aligned

    @property
    def gene(self):
        """

        :return: the gene which encapsulated in this homolog
        """
        return self._gene

    @property
    def gene_ref(self):
        """
        :return: the gene to which this one is homolog to (reference gene)
        :rtype: :class:`macsypy.gene.Gene` object
        """
        return self._ref


class Analog:
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
        self._gene = gene
        self._ref = gene_ref


    def __getattr__(self, name):
        return getattr(self.gene, name)


    @property
    def gene(self):
        """

        :return: the gene which encapsulated in this Analog
        """
        return self._gene

    @property
    def gene_ref(self):
        """
        :return: the gene to which this one is analog to (reference gene)
        :rtype: :class:`macsypy.gene.Gene` object
        """
        return self._ref


class GeneStatus(Enum):

    MANDATORY = 1
    ACCESSORY = 2
    FORBIDDEN = 3
    NEUTRAL = 4

    def __str__(self):
        return self.name.lower()
