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


import logging
_log = logging.getLogger(__name__)
from enum import Enum

from .error import MacsypyError


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
        :rtype: :class:`macsypy.gene.CoreGene` object
        :raise KeyError: if the key does not exist in GeneBank.
        """
        try:
            return self._genes_bank[key]
        except KeyError:
            raise KeyError(f"No such gene '{key}' in this bank")


    def __len__(self):
        return len(self._genes_bank)


    def __contains__(self, gene):
        """
        Implement the membership test operator

        :param gene: the gene to test
        :type gene: :class:`macsypy.gene.CoreGene` object
        :return: True if the gene is in, False otherwise
        :rtype: boolean
        """
        return gene in set(self._genes_bank.values())


    def __iter__(self):
        """
        Return an iterator object on the genes contained in the bank
        """
        return iter(self._genes_bank.values())


    def add_new_gene(self, model_location, name, profile_factory):
        """
        Create a gene and store it in the bank. If the same gene (same name) is add twice,
        it is created only the first time.

        :param model_location: the location where the model family can be found.
        :type model_location: :class:`macsypy.registry.ModelLocation` object
        :param name: the name of the gene to add
        :type name: str
        :param profile_factory: The Profile factory
        :type profile_factory: :class:`profile.ProfileFactory` object.
        """
        key = (model_location.name, name)
        if key not in self._genes_bank:
            gene = CoreGene(model_location, name, profile_factory)
            self._genes_bank[key] = gene


class CoreGene:
    """
    Modelize gene attach to a profile.
    It can be only one instance with the the same name (familly name, gene name)
    """
    def __init__(self, model_location, name, profile_factory):
        self._name = name
        self._model_family_name = model_location.name
        self._profile = profile_factory.get_profile(self, model_location)

    def __hash__(self):
        return hash((self._name, self._model_family_name))

    @property
    def name(self):
        """
        The name of the gene a hmm profile with the same name must exists.
        """
        return self._name

    @property
    def model_family_name(self):
        """
        The name of the model family for instance 'CRISPRCas' or 'TXSS'
        """
        return self._model_family_name

    @property
    def profile(self):
        """
        The HMM protein Profile corresponding to this gene :class:`macsypy.profile.Profile` object
        """

        return self._profile


class ModelGene:
    """
    Handle Gene describe in a Model
    """

    def __init__(self, gene, model, loner=False, multi_system=False, inter_gene_max_space=None):
        """
        Handle gene described in a Model

        :param gene: a gene link to a profile
        :type gene: a :class:`macsypy.gene.CoreGene` object.
        :param model: the model that owns this Gene
        :type model: :class:`macsypy.model.Model` object.
        :param loner: True if the Gene can be isolated on the genome (with no contiguous genes), False otherwise.
        :type loner: boolean.
        :param multi_system: True if this Gene can belong to different occurrences of this System.
        :type multi_system: boolean.
        :param inter_gene_max_space: the maximum space between this Gene and another gene of the System.
        :type inter_gene_max_space: integer
        """
        if not isinstance(gene, CoreGene):
            raise MacsypyError(f"The ModeleGene gene argument must be a CoreGene not {type(gene)}.")
        self._gene = gene
        self._exchangeables = []
        self._model = model
        self._loner = loner
        self._multi_system = multi_system
        self._inter_gene_max_space = inter_gene_max_space


    def __getattr__(self, item):
        try:
            return getattr(self._gene, item)
        except AttributeError:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{item}'")

    def __str__(self):
        """
        Print the name of the gene and of its exchangeable genes.
        """
        s = f"name : {self.name}"
        s += f"\ninter_gene_max_space: {self.inter_gene_max_space:d}"
        if self.loner:
            s += "\nloner"
        if self.multi_system:
            s += "\nmulti_system"
        if self._exchangeables:
            s += "\n    exchangeables: "
            for h in self.exchangeables:
                s += h.name + ", "
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
    def exchangeables(self):
        """
        :return: the list of genes which can replace this one without any effect on the model
        :rtype: list of :class:`macsypy.gene.ModelGene` objects
        """
        return self._exchangeables[:]


    @property
    def is_exchangeable(self):
        """
        :return: True if this gene is describe in the model as an exchangeable.
                 False if ot is describe as first level gene.
        """
        return False


    def alternate_of(self):
        """
        :return: the gene to which this one is an exchangeable to (reference gene),
                 or itself if it is a first level gene.
        :rtype: :class:`macsypy.gene.ModelGene` object
        """
        return self


    def add_exchangeable(self, exchangeable):
        """
        Add a exchangeable gene to this Gene

        :param exchangeable: the exchangeable to add
        :type exchangeable:  :class:`macsypy.gene.Exchangeable` object
        """
        self._exchangeables.append(exchangeable)


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


    def __hash__(self):
        # needed to be hashable in Py3 when __eq__ is defined
        # see https://stackoverflow.com/questions/1608842/types-that-define-eq-are-unhashable
        return id(self)


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


class Exchangeable(ModelGene):
    """
    Handle Exchangeables. Exchangeable are ModelGene which can replaced functionally an other ModelGene.
    Biologically it can be Homolog or Analog
    """

    def __init__(self, c_gene, gene_ref):
        """
        :param c_gene: the gene
        :type c_gene: :class:`macsypy.gene.CoreGene` object.
        :param gene_ref: the gene to which the current can replace it.
        :type gene_ref: :class:`macsypy.gene.ModelGene` object.
        """
        super().__init__(c_gene, gene_ref.model,
                         loner=gene_ref.loner,
                         multi_system=gene_ref.multi_system,
                         inter_gene_max_space=gene_ref.inter_gene_max_space)
        self._ref = gene_ref


    @property
    def is_exchangeable(self):
        """
        :return: True
        """
        return True


    def alternate_of(self):
        """
        :return: the gene to which this one is an exchangeable to (reference gene)
        :rtype: :class:`macsypy.gene.ModelGene` object
        """
        return self._ref


    def add_exchangeable(self, exchangeable):
        """
        This method should never be called, it's a security to avoid to add exchangeable to an exchangeable.

        :param exchangeable:
        :type exchangeable: :class:`macsypy.gene.Exchangeable`
        :raise MacsypyError:
        """
        raise MacsypyError("Cannot add 'Exchangeable' to an Exchangeable")


class GeneStatus(Enum):
    """
    Handle status of Gene
    GeneStatus can take 4 value:

    * MANDATORY
    * ACCESSORY
    * FORBIDDEN
    * NEUTRAL

    """

    MANDATORY = 1
    ACCESSORY = 2
    FORBIDDEN = 3
    NEUTRAL = 4

    def __str__(self):
        return self.name.lower()
