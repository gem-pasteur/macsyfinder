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

import abc
from string import Template


class SystemSerializer(metaclass=abc.ABCMeta):
    """
    handle the different way to serialize a system
    """

    @abc.abstractmethod
    def serialize(self, system, hit_system_tracker):
        pass


class TxtSystemSerializer(SystemSerializer):
    """
    Handle System serialization in text
    """


    def serialize(self, system, hit_system_tracker):
        """
        :return: a string representation of system readable by human
        """
        clst = ", ".join(["[" + ", ".join([str((v_h.id, v_h.gene.name, v_h.position)) for v_h in cluster.hits]) + "]"
                          for cluster in system.clusters])

        s = f"""system id = {system.id}
model = {system.model.fqn}
replicon = {system.replicon_name}
clusters = {clst}
occ = {system.occurrence()}
wholeness = {system.wholeness:.3f}
loci nb = {system.loci}
score = {system.score:.3f}
"""
        for title, genes in (("mandatory", system.mandatory_occ),
                             ("accessory", system.accessory_occ),
                             ("neutral", system.neutral_occ)):
            s += f"\n{title} genes:\n"
            for g_name, hits in genes.items():
                s += f"\t- {g_name}: {len(hits)} "
                all_hits_str = []
                for h in hits:
                    used_in_systems = [s.id for s in hit_system_tracker[h.hit]
                                       if s.model.fqn != system.model.fqn]
                    used_in_systems.sort()
                    if used_in_systems:
                        hit_str = f"{h.gene.name} [{', '.join(used_in_systems)}]"
                    else:
                        hit_str = f"{h.gene.name}"
                    all_hits_str.append(hit_str)
                s += f'({", ".join(all_hits_str)})\n'

        return s


class TsvSystemSerializer(SystemSerializer):
    """
    Handle System serialization in tsv format
    """

    header = "replicon\thit_id\tgene_name\thit_pos\tmodel_fqn\tsys_id\tsys_loci\tsys_wholeness\tsys_score\tsys_occ" \
             "\thit_gene_ref\thit_status\thit_seq_len\thit_i_eval\thit_score\thit_profile_cov\thit_seq_cov\t" \
             "hit_begin_match\thit_end_match\tused_in"

    template = Template("$sys_replicon_name\t$vh_id\t$vh_gene_name\t$vh_position\t$sys_model_fqn\t"
                        "$sys_id\t$sys_loci\t$sys_wholeness\t$sys_score\t"
                        "$sys_occurrence\t$vh_gene_role\t$vh_status\t$vh_seq_length\t$vh_i_eval\t"
                        "$vh_score\t$vh_profile_coverage\t$vh_sequence_coverage\t$vh_begin_match"
                        "\t$vh_end_match\t$used_in_systems\n")


    def serialize(self, system, hit_system_tracker):
        r"""

        :return: a serialisation of this system in tabulated separated value format
                 each line represent a hit and have the following structure:
                     replicon\\thit_id\\tgene_name\\thit_pos\\tmodel_fqn\\tsys_id\\tsys_loci\\tsys_wholeness\\tsys_score
                     \\tsys_occ\\thit_gene_ref.alternate_of\\thit_status\\thit_seq_len\\thit_i_eval\\thit_score\\thit_profile_cov
                     \\thit_seq_cov\\tit_begin_match\\thit_end_match

        :rtype: str
        """
        tsv = ''
        for cluster in system.clusters:
            for vh in cluster.hits:
                used_in_systems = [s.id for s in hit_system_tracker[vh.hit] if s.model.fqn != system.model.fqn]
                tsv += self.template.substitute(
                    sys_replicon_name=system.replicon_name,
                    vh_id=vh.id,
                    vh_gene_name=vh.gene.name,
                    vh_position=vh.position,
                    sys_model_fqn=system.model.fqn,
                    sys_id=system.id,
                    sys_loci=system.loci,
                    sys_wholeness=f"{system.wholeness:.3f}",
                    sys_score=f"{system.score:.3f}",
                    sys_occurrence=system.occurrence(),
                    vh_gene_role=vh.gene_ref.alternate_of().name,
                    vh_status=vh.status,
                    vh_seq_length=vh.seq_length,
                    vh_i_eval=vh.i_eval,
                    vh_score=f"{vh.score:.3f}",
                    vh_profile_coverage=f"{vh.profile_coverage:.3f}",
                    vh_sequence_coverage=f"{vh.sequence_coverage:.3f}",
                    vh_begin_match=vh.begin_match,
                    vh_end_match=vh.end_match,
                    used_in_systems=','.join(used_in_systems)
                )

        return tsv


class TsvSolutionSerializer:
    """
    Handle Solution serialization in tsv format
    """

    def __init__(self):
        """Constructor """
        __class__.header = 'sol_id\t' + TsvSystemSerializer.header
        __class__.template = Template(f"$$sol_id\t{TsvSystemSerializer.template.template}")


    def serialize(self, solution, sol_id, hit_system_tracker):
        """
        :param solution: the solution to serialize
        :type solution: list of :class:`macsypy.system.sSystem` object
        :param hit_system_tracker:
        :type hit_system_tracker: :class:`macsypy.system.HitSystemTracker` object
        :return:
        """
        tsv = ''
        sys_ser = TsvSystemSerializer()
        sys_ser.template = self.template

        for system in solution:
            sol_temp = Template(sys_ser.serialize(system, hit_system_tracker))
            tsv += f"{sol_temp.substitute(sol_id=sol_id)}\n"
        return tsv
