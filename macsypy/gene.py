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
from __future__ import annotations

from enum import Enum
import logging

from .error import MacsypyError

from typing import Iterator, Any, TYPE_CHECKING
from .registries import ModelLocation
if TYPE_CHECKING:
    from .model import Model
    from .profile import ProfileFactory, Profile

_log = logging.getLogger(__name__)


class GeneBank:
    """
    Store all Gene objects. Ensure that genes are instanciated only once.
    """

    def __init__(self) -> None:
        self._genes_bank = {}

    def __getitem__(self, key: tuple[str, str]) -> CoreGene:
        """
        :param key: The key to retrieve a gene.
                    The key is composed of the name of models family and the gene name.
                    for instance CRISPR-Cas/cas9_TypeIIB ('CRISPR-Cas' , 'cas9_TypeIIB') or
                    TXSS/T6SS_tssH ('TXSS', 'T6SS_tssH')
        :return: return the Gene corresponding to the key.
        :raise KeyError: if the key does not exist in GeneBank.
        """
        try:
            return self._genes_bank[key]
        except KeyError:
            raise KeyError(f"No such gene '{key}' in this bank")


    def __len__(self) -> int:
        return len(self._genes_bank)


    def __contains__(self, gene: CoreGene) -> bool:
        """
        Implement the membership test operator

        :param gene: the gene to test
        :return: True if the gene is in, False otherwise
        :rtype: boolean
        """
        return gene in set(self._genes_bank.values())


    def __iter__(self) -> Iterator:
        """
        Return an iterator object on the genes contained in the bank
        """
        return iter(self._genes_bank.values())


    def genes_fqn(self) -> list[str]:
        """
        :return: the fully qualified name for all genes in the bank
        """
        return [f"{fam}/{gen_nam}" for fam, gen_nam in self._genes_bank]


    def add_new_gene(self, model_location: ModelLocation, name: str, profile_factory: ProfileFactory) -> None:
        """
        Create a gene and store it in the bank. If the same gene (same name) is add twice,
        it is created only the first time.

        :param model_location: the location where the model family can be found.
        :param name: the name of the gene to add
        :param profile_factory: The Profile factory
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
    def __init__(self, model_location: ModelLocation, name: str, profile_factory: ProfileFactory) -> None:
        self._name = name
        self._model_family_name = model_location.name
        self._profile = profile_factory.get_profile(self, model_location)

    def __hash__(self) -> int:
        return hash((self._name, self._model_family_name))

    @property
    def name(self) -> str:
        """
        The name of the gene a hmm profile with the same name must exists.
        """
        return self._name

    @property
    def model_family_name(self) -> str:
        """
        The name of the model family for instance 'CRISPRCas' or 'TXSS'
        """
        return self._model_family_name

    @property
    def profile(self) -> Profile:
        """
        The HMM protein Profile corresponding to this gene
        """

        return self._profile


class ModelGene:
    """
    Handle Gene describe in a Model
    """

    def __init__(self,
                 gene: CoreGene,
                 model: Model,
                 loner: bool = False,
                 multi_system: bool = False,
                 inter_gene_max_space: int = None,
                 multi_model: bool = False):
        """
        Handle gene described in a Model

        :param gene: a gene link to a profile
        :param model: the model that owns this Gene
        :param loner: True if the Gene can be isolated on the genome (with no contiguous genes), False otherwise.
        :param multi_system: True if this Gene can belong to different occurrences of this System.
        :param inter_gene_max_space: the maximum space between this Gene and another gene of the System.
        :param multi_model: True if this Gene is allowing to appear in several system occurence from diferent model.
        """
        if not isinstance(gene, CoreGene):
            raise MacsypyError(f"The ModeleGene gene argument must be a CoreGene not {type(gene)}.")
        self._gene = gene
        self._exchangeables = []
        self._model = model
        self._loner = loner if loner else False
        self._multi_system = multi_system if multi_system else False
        self._multi_model = multi_model if multi_model else False
        self._inter_gene_max_space = inter_gene_max_space
        self._status = None


    def __getattr__(self, item: str) -> Any:
        try:
            return getattr(self._gene, item)
        except AttributeError as err:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{item}'") from err


    def __str__(self) -> str:
        """
        Print the name of the gene and of its exchangeable genes.
        """
        rep = f"name : {self.name}"
        rep += f"\ninter_gene_max_space: {self.inter_gene_max_space}"
        if self.loner:
            rep += "\nloner"
        if self.multi_system:
            rep += "\nmulti_system"
        if self.multi_model:
            rep += "\nmulti_model"
        if self._exchangeables:
            rep += "\n    exchangeables: "
            for m_hit in self.exchangeables:
                rep += m_hit.name + ", "
            rep = rep[:-2]
        return rep


    @property
    def status(self) -> GeneStatus:
        """
        :return: The status of this gene
        """
        return self._status


    def set_status(self, status: GeneStatus) -> None:
        """
        Set the status for this gene

        :param status: the status of this gene
        """
        self._status = status


    @property
    def model(self) -> Model:
        """
        :return: the Model that owns this Gene
        """
        return self._model

    @property
    def core_gene(self) -> CoreGene:
        """
        :return: The CoreGene associated to this ModelGene
        """
        return self._gene


    @property
    def loner(self) -> bool:
        """
        :return: True if the gene can be isolated on the genome, False otherwise
        """
        return self._loner


    @property
    def exchangeables(self) -> list[ModelGene]:
        """
        :return: the list of genes which can replace this one without any effect on the model
        """
        return self._exchangeables[:]


    @property
    def is_exchangeable(self) -> bool:
        """
        :return: True if this gene is described in the model as an exchangeable.
                 False if it is described as first level gene.
        """
        return False


    def alternate_of(self) -> ModelGene:
        """
        :return: the gene to which this one is an exchangeable to (reference gene),
                 or itself if it is a first level gene.
        """
        return self


    def add_exchangeable(self, exchangeable: Exchangeable):
        """
        Add an exchangeable gene to this Gene

        :param exchangeable: the exchangeable to add
        """
        self._exchangeables.append(exchangeable)


    @property
    def multi_system(self) -> bool:
        """
        :return: True if this Gene can belong to different occurrences of **the model**
                (and can be used for multiple System assessments), False otherwise.
        """
        return self._multi_system


    @property
    def multi_model(self) -> bool:
        """
        :return: True if this Gene can belong to different occurrences of systems from different
                 model :class:`macsypy.model.Model`
                (and can be used for multiple System assessments), False otherwise.
        :rtype: boolean.
        """
        return self._multi_model


    @property
    def inter_gene_max_space(self) -> int | None:
        """
        :return: The maximum distance allowed between this gene and another gene for them to be considered co-localized.
                 If the value is not set at the Gene level, return None.
        """
        return self._inter_gene_max_space


    def __hash__(self) -> int:
        # needed to be hashable in Py3 when __eq__ is defined
        # see https://stackoverflow.com/questions/1608842/types-that-define-eq-are-unhashable
        return id(self)


    def is_mandatory(self, model: Model) -> bool:
        """
        :return: True if the gene is within the *mandatory* genes of the model, False otherwise.
        :param model: the query of the test
        """
        return self in model.mandatory_genes


    def is_accessory(self, model: Model) -> bool:
        """
        :return: True if the gene is within the *accessory* genes of the model, False otherwise.
        :param model: the query of the test
        """
        return self in model.accessory_genes


    def is_forbidden(self, model: Model) -> bool:
        """
        :return: True if the gene is within the *forbidden* genes of the model, False otherwise.
        :param model: the query of the test
        """
        return self in model.forbidden_genes


class Exchangeable(ModelGene):
    """
    Handle Exchangeables. Exchangeable are ModelGene which can replaced functionally an other ModelGene.
    Biologically it can be Homolog or Analog
    """

    def __init__(self,
                 c_gene: CoreGene,
                 gene_ref: ModelGene,
                 loner: bool | None = None,
                 multi_system: bool | None = None,
                 multi_model: bool | None = None,
                 inter_gene_max_space: int | None = None) -> None:
        """
        :param c_gene: the gene
        :param gene_ref: the gene to which the current can replace it.
        """
        super().__init__(c_gene, gene_ref.model,
                         loner=loner if loner is not None else gene_ref.loner,
                         multi_system=multi_system if multi_system is not None else gene_ref.multi_system,
                         multi_model=multi_model if multi_model is not None else gene_ref.multi_model,
                         inter_gene_max_space=inter_gene_max_space if inter_gene_max_space is not None
                         else gene_ref.inter_gene_max_space)
        self._ref = gene_ref


    @property
    def is_exchangeable(self) -> bool:
        """
        :return: True
        """
        return True


    def alternate_of(self) -> ModelGene:
        """
        :return: the gene to which this one is an exchangeable to (reference gene)
        """
        return self._ref


    def add_exchangeable(self, exchangeable: Exchangeable) -> None:
        """
        This method should never be called, it's a security to avoid to add exchangeable to an exchangeable.

        :param exchangeable: the exchangeable gene to add
        :raise MacsypyError:
        """
        raise MacsypyError("Cannot add 'Exchangeable' to an Exchangeable")

    @property
    def status(self) -> GeneStatus:
        """
        :return: The status of this gene. if the status is not define for this gene itself,
                 return the status of the reference gene.
        """
        if self._status:
            return self.status
        return self._ref.status


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

    def __str__(self) -> str:
        return self.name.lower()
