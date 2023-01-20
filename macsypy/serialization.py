#########################################################################
# MacSyFinder - Detection of macromolecular systems in protein dataset  #
#               using systems modelling and similarity search.          #
# Authors: Sophie Abby, Bertrand Neron                                  #
# Copyright (c) 2014-2023  Institut Pasteur (Paris) and CNRS.           #
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
This module focus on the way to serialize the different macsyfinder modules
"""

import abc
from string import Template

from macsypy.gene import GeneStatus


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
loci nb = {system.loci_nb}
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

    header = "replicon\thit_id\tgene_name\thit_pos\tmodel_fqn" \
             "\tsys_id\tsys_loci\tlocus_num\tsys_wholeness\tsys_score\tsys_occ" \
             "\thit_gene_ref\thit_status\thit_seq_len\thit_i_eval\thit_score\thit_profile_cov\thit_seq_cov\t" \
             "hit_begin_match\thit_end_match\tcounterpart\tused_in"

    template = Template("$sys_replicon_name\t$mh_id\t$mh_gene_name\t$mh_position\t$sys_model_fqn\t"
                        "$sys_id\t$sys_loci\t$locus_num\t$sys_wholeness\t$sys_score\t"
                        "$sys_occurrence\t$mh_gene_role\t$mh_status\t$mh_seq_length\t$mh_i_eval\t"
                        "$mh_score\t$mh_profile_coverage\t$mh_sequence_coverage\t$mh_begin_match"
                        "\t$mh_end_match\t$mh_counterpart\t$used_in_systems\n")


    def serialize(self, system, hit_system_tracker):
        r"""
        :param :class:`macsypy.system.System` system: The system to serialize.
        :param hit_system_tracker: The hit_system_tracker which allow to know for each hit
               in which system it is implied.
        :type hit_system_tracker: :class:`macsypy.system.HitSystemTracker` object
        :return: a serialisation of this system in tabulated separated value format
                 each line represent a hit and have the following structure:

        .. code-block:: text

            replicon\thit_id\tgene_name\thit_pos\tmodel_fqn\tsys_id\tsys_loci\tlocus_num\tsys_wholeness\tsys_score
            \tsys_occ\thit_gene_ref.alternate_of\thit_status\thit_seq_len\thit_i_eval\thit_score\thit_profile_cov
            \thit_seq_cov\tit_begin_match\thit_end_match\tcounterpart\tused_in_systems

        :rtype: str
        """
        tsv = ''
        loci_num = system.loci_num
        for locus_num, cluster in zip(loci_num, system.clusters):
            for mh in sorted(cluster.hits, key=lambda mh: mh.position):
                used_in_systems = [s.id for s in hit_system_tracker[mh.hit] if s.model.fqn != system.model.fqn]
                used_in_systems.sort()
                tsv += self.template.substitute(
                    sys_replicon_name=system.replicon_name,
                    mh_id=mh.id,
                    mh_gene_name=mh.gene.name,
                    mh_position=mh.position,
                    sys_model_fqn=system.model.fqn,
                    sys_id=system.id,
                    sys_loci=system.loci_nb,
                    locus_num=locus_num,
                    sys_wholeness=f"{system.wholeness:.3f}",
                    sys_score=f"{system.score:.3f}",
                    sys_occurrence=system.occurrence(),
                    mh_gene_role=mh.gene_ref.alternate_of().name,
                    mh_status=mh.status,
                    mh_seq_length=mh.seq_length,
                    mh_i_eval=mh.i_eval,
                    mh_score=f"{mh.score:.3f}",
                    mh_profile_coverage=f"{mh.profile_coverage:.3f}",
                    mh_sequence_coverage=f"{mh.sequence_coverage:.3f}",
                    mh_begin_match=mh.begin_match,
                    mh_end_match=mh.end_match,
                    mh_counterpart=','.join([h.id for h in mh.counterpart]),
                    used_in_systems=','.join(used_in_systems)
                )
        return tsv


class TsvSolutionSerializer:
    """
    Handle Solution (list of Systems) serialization in tsv format
    """

    def __init__(self):
        """Constructor """
        __class__.header = 'sol_id\t' + TsvSystemSerializer.header
        __class__.template = Template(f"$$sol_id\t{TsvSystemSerializer.template.template}")


    def serialize(self, solution, sol_id, hit_system_tracker):
        """
        :param solution: the solution to serialize
        :type solution: list of :class:`macsypy.system.System` object
        :param int sol_id: the solution identifier
        :param hit_system_tracker:
        :type hit_system_tracker: :class:`macsypy.system.HitSystemTracker` object
        :return: a serialisation of this solution (a list of systems) in tabulated separated value format
                 each line represent a hit and have the same structure as system serialization
                 :meth:`macsypy.serialization.TsvSystemSerializer.serialize` but with an extra column
                 sol_id which is a technical id to identified the different solutions.
        """
        tsv = ''
        sys_ser = TsvSystemSerializer()
        sys_ser.template = self.template

        for system in solution:
            sol_temp = Template(sys_ser.serialize(system, hit_system_tracker))
            tsv += f"{sol_temp.substitute(sol_id=sol_id)}\n"
        return tsv


class TxtLikelySystemSerializer(SystemSerializer):
    """
    Handle System serialization in text
    """


    def serialize(self, system, hit_system_tracker):
        """
        :param :class:`macsypy.system.LikelySystem` system: The likely system to serialize.
                                                            Use only for unordered db-type
        :param hit_system_tracker: The hit_system_tracker which allow to know for each hit
               in which system it is implied.
        :type hit_system_tracker: :class:`macsypy.system.HitSystemTracker` object
        :return: a string representation of system readable by human
        """
        hits = ", ".join([str((h.id, h.gene.name, h.position)) for h in system.hits])
        if system.forbidden_hits:
            warning = "WARNING there quorum is reached but there is also some forbidden genes.\n"
        else:
            warning = '\n'

        s = f"""This replicon contains genetic materials needed for system {system.model.fqn}
{warning}
system id = {system.id}
model = {system.model.fqn}
replicon = {system.replicon_name}
hits = [{hits}]
wholeness = {system.wholeness:.3f}
"""
        for title, genes in (("mandatory", system.mandatory_occ),
                             ("accessory", system.accessory_occ),
                             ("neutral", system.neutral_occ),
                             ("forbidden", system.forbidden_occ)):
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

        s += "\nUse ordered replicon to have better prediction.\n"
        return s


class TsvLikelySystemSerializer(SystemSerializer):
    """
    Handle potential System from unordered replicon
    serialization in tsv format
    """

    header = "replicon\thit_id\tgene_name\thit_pos\tmodel_fqn\tsys_id\tsys_wholeness" \
             "\thit_gene_ref\thit_status\thit_seq_len\thit_i_eval\thit_score\thit_profile_cov\thit_seq_cov\t" \
             "hit_begin_match\thit_end_match\tused_in"

    template = Template("$sys_replicon_name\t$mh_id\t$mh_gene_name\t$mh_position\t$sys_model_fqn\t"
                        "$sys_id\t$sys_wholeness\t"
                        "$mh_gene_role\t$mh_status\t$mh_seq_length\t$mh_i_eval\t"
                        "$mh_score\t$mh_profile_coverage\t$mh_sequence_coverage\t$mh_begin_match"
                        "\t$mh_end_match\t$used_in_systems\n")


    def serialize(self, system, hit_system_tracker):
        r"""
        :param :class:`macsypy.system.LikelySystem` system: The likely system to serialize.
                                                            Use only for unordered db-type
        :param hit_system_tracker: The hit_system_tracker which allow to know for each hit
               in which system it is implied.
        :type hit_system_tracker: :class:`macsypy.system.HitSystemTracker` object
        :return: a serialisation of this system in tabulated separated value format
                 each line represent a hit and have the following structure:

        .. code-block:: text

            replicon\thit_id\tgene_name\thit_pos\tmodel_fqn\tsys_id\tsys_wholeness
            \thit_gene_ref.alternate_of\thit_status\thit_seq_len\thit_i_eval\thit_score\thit_profile_cov
            \thit_seq_cov\tit_begin_match\thit_end_match\t$used_in_systems

        :rtype: str
        """
        tsv = ''
        for status in (s.lower() for s in GeneStatus.__members__):
            try:
                hits = getattr(system, f"{status}_hits")
                hits = sorted(hits, key=lambda mh: mh.gene.name)
            except AttributeError:
                continue
            for mh in hits:
                used_in_systems = [s.id for s in hit_system_tracker[mh.hit] if s.model.fqn != system.model.fqn]
                used_in_systems.sort()
                tsv += self.template.substitute(
                    sys_replicon_name=system.replicon_name,
                    mh_id=mh.id,
                    mh_gene_name=mh.gene.name,
                    mh_position=mh.position,
                    sys_model_fqn=system.model.fqn,
                    sys_id=system.id,
                    sys_wholeness=f"{system.wholeness:.3f}",
                    mh_gene_role=mh.gene_ref.alternate_of().name,
                    mh_status=mh.status,
                    mh_seq_length=mh.seq_length,
                    mh_i_eval=mh.i_eval,
                    mh_score=f"{mh.score:.3f}",
                    mh_profile_coverage=f"{mh.profile_coverage:.3f}",
                    mh_sequence_coverage=f"{mh.sequence_coverage:.3f}",
                    mh_begin_match=mh.begin_match,
                    mh_end_match=mh.end_match,
                    used_in_systems=','.join(used_in_systems)
                )

        return tsv


class TxtUnikelySystemSerializer(SystemSerializer):
    """
    Handle System serialization in text
    """


    def serialize(self, system):
        """
        :param system: The unlikely system to serialize. (used only if db-type is "unordered_replicon")
        :type system: :class:`macsypy.system.UnlikelySystem` object
        :return: a string representation of system readable by human
        """
        hits = ", ".join([str((h.id, h.gene.name, h.position)) for h in system.hits])
        reasons = '\n'.join(system.reasons)
        s = f"""This replicon probably not contains a system {system.model.fqn}:
{reasons}

system id = {system.id}
model = {system.model.fqn}
replicon = {system.replicon_name}
hits = [{hits}]
wholeness = {system.wholeness:.3f}
"""
        for title, genes in (("mandatory", system.mandatory_occ),
                             ("accessory", system.accessory_occ),
                             ("neutral", system.neutral_occ),
                             ("forbidden", system.forbidden_occ)):
            s += f"\n{title} genes:\n"
            for g_name, hits in genes.items():
                s += f"\t- {g_name}: {len(hits)} "
                all_hits_str = [f"{h.gene.name}" for h in hits]
                s += f'({", ".join(all_hits_str)})\n'

        s += "\nUse ordered replicon to have better prediction.\n"
        return s


class TsvSpecialHitSerializer:
    """
    Serialize special hits: :class:`macsypy.hit.Loner` and :class:`macsypy.hit.MultiSystem` in tsv format
    """

    def serialize(self, best_hits):
        """
        :param best_hits: the special hits to serialized
        :type best_hits: sequence of :class:`macsypy.hit.Loner` or :class:`macsypy.hit.MultiSystem` objects
        """
        s = ""
        if best_hits:
            header = "replicon\tmodel_fqn\tfunction\tgene_name\t" \
                     "hit_id\thit_pos\thit_status\thit_seq_len\t" \
                     "hit_i_eval\thit_score\thit_profile_cov\t" \
                     "hit_seq_cov\thit_begin_match\thit_end_match\n"
            s += header
            special_hits = set(best_hits)
            for best_hit in best_hits:
                special_hits.update(best_hit.counterpart)
            special_hits = list(special_hits)
            special_hits.sort(key=lambda h: h.position)
            for one_hit in special_hits:
                row = f"{one_hit.replicon_name}\t{one_hit.gene_ref.model.fqn}\t{one_hit.gene_ref.alternate_of().name}\t" \
                      f"{one_hit.gene_ref.name}\t{one_hit.id}\t{one_hit.position:d}\t{one_hit.status}\t" \
                      f"{one_hit.seq_length:d}\t{one_hit.i_eval:.3e}\t{one_hit.score:.3f}\t" \
                      f"{one_hit.profile_coverage:.3f}\t{one_hit.sequence_coverage:.3f}\t" \
                      f"{one_hit.begin_match:d}\t{one_hit.end_match:d}\n"
                s += row
        return s


class TsvRejectedCandidatesSerializer:
    """
    Serialize Rejected Cluster in tsv format
    """

    def serialize(self, candidates):
        """
        :param candidates: list of rejected candidates to serialize
        :type candidates: [ :class:`macsypy.system.RejectedCandidate` object, ...]
        """
        s = ""
        if candidates:
            header = "candidate_id\treplicon\tmodel_fqn\tcluster_id\thit_id\thit_pos\tgene_name\tfunction\treasons\n"
            s += header
            for candidate in candidates:
                reasons = '/'.join(candidate.reasons)
                for cluster in candidate.clusters:
                    for hit in cluster.hits:
                        row = f"{candidate.id}\t{candidate.replicon_name}\t{candidate.model.fqn}\t" \
                              f"{cluster.id}\t{hit.id}\t{hit.position}\t{hit.gene_ref.name}\t{hit.gene_ref.alternate_of().name}\t" \
                              f"{reasons}\n"
                        s += row
                s += '\n'
        return s