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

import os
import argparse
import random

from macsypy.hit import CoreHit, ModelHit, Loner, MultiSystem, LonerMultiSystem,  HitWeight
from macsypy.config import Config, MacsyDefaults
from macsypy.gene import CoreGene, ModelGene, Exchangeable, GeneStatus
from macsypy.profile import ProfileFactory
from macsypy.model import Model
from macsypy.registries import ModelLocation
from macsypy.cluster import Cluster
from macsypy.system import System, RejectedCandidate
from macsypy.solution import find_best_solutions, combine_clusters, combine_multisystems, Solution
from tests import MacsyTest

    
def _build_clusters(cfg, profile_factory):
    
    model_name = 'foo'
    model_location = ModelLocation(path=os.path.join(cfg.models_dir()[0], model_name))

    models = {}
    cg_sctn_flg = CoreGene(model_location, "sctN_FLG", profile_factory)
    cg_sctj_flg = CoreGene(model_location, "sctJ_FLG", profile_factory)
    cg_flgB = CoreGene(model_location, "flgB", profile_factory)
    cg_tadZ = CoreGene(model_location, "tadZ", profile_factory)
    cg_sctn = CoreGene(model_location, "sctN", profile_factory)
    cg_sctj = CoreGene(model_location, "sctJ", profile_factory)
    cg_gspd = CoreGene(model_location, "gspD", profile_factory)
    cg_abc = CoreGene(model_location, "abc", profile_factory)
    cg_sctc = CoreGene(model_location, "sctC", profile_factory)

    ###########
    # Model A #
    ###########
    models['A'] = Model("foo/A", 10)
    mgA_sctn = ModelGene(cg_sctn, models['A'])
    mgA_sctn_hom = Exchangeable(cg_sctn_flg, mgA_sctn)
    mgA_sctn.add_exchangeable(mgA_sctn_hom)
    mgA_sctj = ModelGene(cg_sctj, models['A'])
    mgA_sctj_an = Exchangeable(cg_sctj_flg, mgA_sctj)
    mgA_sctj.add_exchangeable(mgA_sctj_an)
    mgA_gspd = ModelGene(cg_gspd, models['A'])
    mgA_gspd_an = Exchangeable(cg_flgB, mgA_gspd)
    mgA_gspd.add_exchangeable(mgA_gspd_an)
    mgA_abc = ModelGene(cg_abc, models['A'])
    mgA_abc_ho = Exchangeable(cg_tadZ, mgA_abc)
    mgA_abc.add_exchangeable(mgA_abc_ho)

    models['A'].add_mandatory_gene(mgA_sctn)
    models['A'].add_mandatory_gene(mgA_sctj)
    models['A'].add_accessory_gene(mgA_gspd)
    models['A'].add_forbidden_gene(mgA_abc)

    models['A']._min_mandatory_genes_required = 2
    models['A']._min_genes_required = 2

    ###########
    # Model B #
    ###########
    models['B'] = Model("foo/B", 10)
    mgB_sctn_flg = ModelGene(cg_sctn_flg, models['B'])
    mgB_sctj_flg = ModelGene(cg_sctj_flg, models['B'])
    mgB_flgB = ModelGene(cg_flgB, models['B'])
    mgB_tadZ = ModelGene(cg_tadZ, models['B'])

    models['B'].add_mandatory_gene(mgB_sctn_flg)
    models['B'].add_mandatory_gene(mgB_sctj_flg)
    models['B'].add_accessory_gene(mgB_flgB)
    models['B'].add_accessory_gene(mgB_tadZ)

    models['B']._min_mandatory_genes_required = 1
    models['B']._min_genes_required = 2

    ###########
    # Model C #
    ###########
    models['C'] = Model("foo/C", 10)
    mgC_sctn_flg = ModelGene(cg_sctn_flg, models['C'])
    mgC_sctj_flg = ModelGene(cg_sctj_flg, models['C'])
    mgC_flgB = ModelGene(cg_flgB, models['C'])
    mgC_tadZ = ModelGene(cg_tadZ, models['C'])
    mgC_gspd = ModelGene(cg_gspd, models['C'])

    models['C'].add_mandatory_gene(mgC_sctn_flg)
    models['C'].add_mandatory_gene(mgC_sctj_flg)
    models['C'].add_mandatory_gene(mgC_flgB)
    models['C'].add_accessory_gene(mgC_tadZ)
    models['C'].add_accessory_gene(mgC_gspd)

    models['C']._min_mandatory_genes_required = 1
    models['C']._min_genes_required = 2

    ###########
    # Model D #
    ###########
    models['D'] = Model("foo/D", 10)
    mgD_abc = ModelGene(cg_abc, models['D'])
    mgD_sctn = ModelGene(cg_sctn, models['D'])
    models['D'].add_mandatory_gene(mgD_abc)
    models['D'].add_accessory_gene(mgD_sctn)

    models['D']._min_mandatory_genes_required = 1
    models['D']._min_genes_required = 1
    ###########
    # Model E #
    ###########
    models['E'] = Model("foo/E", 10)
    mgE_gspd = ModelGene(cg_gspd, models['E'])
    models['E'].add_accessory_gene(mgE_gspd)

    models['E']._min_mandatory_genes_required = 0
    models['E']._min_genes_required = 1

    ###########
    # Model F #
    ###########
    models['F'] = Model("foo/F", 10)
    mgF_abc = ModelGene(cg_abc, models['F'])
    models['F'].add_mandatory_gene(mgF_abc)

    models['F']._min_mandatory_genes_required = 1
    models['F']._min_genes_required = 1

    #####################
    # Model G idem as C #
    #####################
    models['G'] = Model("foo/G", 10)
    mgG_sctn_flg = ModelGene(cg_sctn_flg, models['G'])
    mgG_sctj_flg = ModelGene(cg_sctj_flg, models['G'])
    mgG_flgB = ModelGene(cg_flgB, models['G'])
    mgG_tadZ = ModelGene(cg_tadZ, models['G'])
    mgG_gspd = ModelGene(cg_gspd, models['G'])
    models['G'].add_mandatory_gene(mgG_sctn_flg)
    models['G'].add_mandatory_gene(mgG_sctj_flg)
    models['G'].add_mandatory_gene(mgG_flgB)
    models['G'].add_accessory_gene(mgG_tadZ)
    models['G'].add_accessory_gene(mgG_gspd)

    #####################
    # Model H idem as D #
    #####################
    models['H'] = Model("foo/H", 10)
    mgH_abc = ModelGene(cg_abc, models['H'])
    mgH_sctn = ModelGene(cg_sctn, models['H'])
    models['H'].add_mandatory_gene(mgH_abc)
    models['H'].add_accessory_gene(mgH_sctn)

    models['H']._min_mandatory_genes_required = 1
    models['H']._min_genes_required = 1

    ###########
    # Model I #
    ###########
    models['I'] = Model("foo/I", 10)
    mgI_abc = ModelGene(cg_abc, models['I'])
    mgI_flgB = ModelGene(cg_flgB, models['I'])
    mgI_tadZ = ModelGene(cg_tadZ, models['I'])
    models['I'].add_mandatory_gene(mgI_abc)
    models['I'].add_mandatory_gene(mgI_flgB)
    models['I'].add_accessory_gene(mgI_tadZ)

    models['I']._min_mandatory_genes_required = 1
    models['I']._min_genes_required = 1

    ###########
    # model J #
    ###########
    models['J'] = Model("foo/J", 10)
    mgJ_abc = ModelGene(cg_abc, models['J'])
    mgJ_gspd = ModelGene(cg_gspd, models['J'])
    mgJ_tadZ = ModelGene(cg_tadZ, models['J'])
    mgJ_sctc = ModelGene(cg_sctc, models['J'])
    models['J'].add_mandatory_gene(mgJ_abc)
    models['J'].add_mandatory_gene(mgJ_gspd)
    models['J'].add_accessory_gene(mgJ_tadZ)
    models['J'].add_accessory_gene(mgJ_sctc)

    models['J']._min_mandatory_genes_required = 1
    models['J']._min_genes_required = 1

    ###########
    # model K #
    ###########
    models['K'] = Model("foo/K", 10)
    mgK_flgB = ModelGene(cg_flgB, models['K'])
    mgK_sctn_flg = ModelGene(cg_sctn_flg, models['K'])
    mgK_sctj_flg = ModelGene(cg_sctj_flg, models['K'])
    mgK_sctn = ModelGene(cg_sctn, models['K'])
    models['K'].add_mandatory_gene(mgK_flgB)
    models['K'].add_mandatory_gene(mgK_sctn_flg)
    models['K'].add_accessory_gene(mgK_sctj_flg)
    models['K'].add_accessory_gene(mgK_sctn)

    models['K']._min_mandatory_genes_required = 1
    models['K']._min_genes_required = 1

    ###########
    # model L #
    ###########
    models['L'] = Model("foo/L", 10)
    mgL_flgB = ModelGene(cg_flgB, models['L'])
    mgL_sctn_flg = ModelGene(cg_sctn_flg, models['L'])
    mgL_sctj_flg = ModelGene(cg_sctj_flg, models['L'])
    mgL_sctn = ModelGene(cg_sctn, models['L'], loner=True)
    models['L'].add_mandatory_gene(mgL_flgB)
    models['L'].add_mandatory_gene(mgL_sctn_flg)
    models['L'].add_accessory_gene(mgL_sctj_flg)
    models['L'].add_accessory_gene(mgL_sctn)

    ###########
    # model M #
    ###########
    models['M'] = Model("foo/L", 10)
    mgM_sctj = ModelGene(cg_sctj, models['M'])
    mgM_gspd = ModelGene(cg_gspd, models['M'])
    mgM_sctn = ModelGene(cg_sctn, models['M'], multi_system=True)
    mgM_tadZ = ModelGene(cg_tadZ, models['M'])
    mgM_abc = ModelGene(cg_abc, models['M'])
    models['M'].add_mandatory_gene(mgM_sctj)
    models['M'].add_mandatory_gene(mgM_gspd)
    models['M'].add_accessory_gene(mgM_sctn)
    models['M'].add_accessory_gene(mgM_tadZ)
    models['M'].add_accessory_gene(mgM_abc)

    ###########
    # model N #
    ###########
    models['N'] = Model("foo/N", 10)
    mgN_flgB = ModelGene(cg_flgB, models['N'])
    mgN_sctn_flg = ModelGene(cg_sctn_flg, models['N'])
    mgN_sctj = ModelGene(cg_sctj, models['N'])
    mgN_sctj_flg = ModelGene(cg_sctj_flg, models['N'])
    mgN_sctn = ModelGene(cg_sctn, models['N'], loner=True)
    mgN_tadZ = ModelGene(cg_tadZ, models['N'], loner=True)
    models['N'].add_mandatory_gene(mgN_flgB)
    models['N'].add_mandatory_gene(mgN_sctn_flg)
    models['N'].add_accessory_gene(mgN_sctj)
    models['N'].add_accessory_gene(mgN_sctj_flg)
    models['N'].add_accessory_gene(mgN_sctn)
    models['N'].add_accessory_gene(mgN_tadZ)

    ###########
    # model O #
    ###########
    models['O'] = Model("foo/O", 10)
    mgO_sctj = ModelGene(cg_sctj, models['O'], multi_system=True)
    mgO_sctj_flg = Exchangeable(cg_sctj_flg, mgO_sctj)
    mgO_sctj.add_exchangeable(mgO_sctj_flg)
    mgO_gspd = ModelGene(cg_gspd, models['O'], loner=True, multi_system=True)
    mgO_sctn = ModelGene(cg_sctn, models['O'], multi_system=True)
    mgO_sctn_flg = Exchangeable(cg_sctn_flg, mgO_sctn)
    mgO_sctn.add_exchangeable(mgO_sctn_flg)
    mgO_tadZ = ModelGene(cg_tadZ, models['O'], loner=True)
    mgO_abc = ModelGene(cg_abc, models['O'])

    models['O'].add_mandatory_gene(mgO_sctj)
    models['O'].add_mandatory_gene(mgO_gspd)
    models['O'].add_accessory_gene(mgO_sctn)
    models['O'].add_accessory_gene(mgO_tadZ)
    models['O'].add_neutral_gene(mgO_abc)

    ch_sctj = CoreHit(cg_sctj, "hit_sctj", 803, "replicon_id", 1, 1.0, 1.0, 1.0, 1.0, 10, 20)
    ch_sctn = CoreHit(cg_sctn, "hit_sctn", 803, "replicon_id", 2, 1.0, 1.0, 1.0, 1.0, 10, 20)
    ch_gspd = CoreHit(cg_gspd, "hit_gspd", 803, "replicon_id", 3, 1.0, 1.0, 1.0, 1.0, 10, 20)
    ch_sctn_flg = CoreHit(cg_sctn_flg, "hit_sctn_flg", 803, "replicon_id", 4, 1.0, 1.0, 1.0, 1.0, 10, 20)
    ch_sctj = CoreHit(cg_sctj, "hit_sctj", 803, "replicon_id", 5, 1.0, 1.0, 1.0, 1.0, 10, 20)
    ch_sctj_flg = CoreHit(cg_sctj_flg, "hit_sctj_flg", 803, "replicon_id", 6, 1.0, 1.0, 1.0, 1.0, 10, 20)
    ch_flgB = CoreHit(cg_flgB, "hit_flgB", 803, "replicon_id", 7, 1.0, 1.0, 1.0, 1.0, 10, 20)
    ch_tadZ = CoreHit(cg_tadZ, "hit_tadZ", 803, "replicon_id", 8, 1.0, 1.0, 1.0, 1.0, 10, 20)
    ch_abc = CoreHit(cg_abc, "hit_abc", 803, "replicon_id", 9, 1.0, 1.0, 1.0, 1.0, 10, 20)

    hit_weights = HitWeight(**cfg.hit_weights())

    clusters = {}
    clusters['c1'] = Cluster([ModelHit(ch_sctj, gene_ref=mgA_sctj, gene_status=GeneStatus.MANDATORY),
                              ModelHit(ch_sctn, gene_ref=mgA_sctn, gene_status=GeneStatus.MANDATORY),
                              ModelHit(ch_gspd, gene_ref=mgA_gspd, gene_status=GeneStatus.ACCESSORY)
                              ],
                              models['A'], hit_weights)
    clusters['c2'] = Cluster([ModelHit(ch_sctj, gene_ref=mgA_sctj, gene_status=GeneStatus.MANDATORY),
                              ModelHit(ch_sctn, gene_ref=mgA_sctn, gene_status=GeneStatus.MANDATORY)],
                              models['A'], hit_weights)

    clusters['c3'] = Cluster([ModelHit(ch_sctj_flg, gene_ref=mgB_sctj_flg, gene_status=GeneStatus.MANDATORY),
                              ModelHit(ch_tadZ, gene_ref=mgB_tadZ, gene_status=GeneStatus.ACCESSORY),
                              ModelHit(ch_flgB, gene_ref=mgB_flgB, gene_status=GeneStatus.ACCESSORY)],
                              models['B'], hit_weights)

    clusters['c4'] = Cluster([ModelHit(ch_sctj_flg, gene_ref=mgC_sctj_flg, gene_status=GeneStatus.MANDATORY),
                              ModelHit(ch_tadZ, gene_ref=mgC_tadZ, gene_status=GeneStatus.ACCESSORY),
                              ModelHit(ch_flgB, gene_ref=mgC_flgB, gene_status=GeneStatus.MANDATORY),
                              ModelHit(ch_gspd, gene_ref=mgC_gspd, gene_status=GeneStatus.ACCESSORY)],
                              models['C'], hit_weights)

    clusters['c5'] = Cluster([ModelHit(ch_abc, gene_ref=mgD_abc, gene_status=GeneStatus.MANDATORY),
                              ModelHit(ch_sctn, gene_ref=mgD_sctn, gene_status=GeneStatus.ACCESSORY)],
                              models['D'], hit_weights)

    clusters['c6'] = Cluster([ModelHit(ch_gspd, gene_ref=mgE_gspd, gene_status=GeneStatus.ACCESSORY)],
                              models['E'], hit_weights)

    clusters['c7'] = Cluster([ModelHit(ch_abc, gene_ref=mgF_abc, gene_status=GeneStatus.MANDATORY)],
                              models['F'], hit_weights)

    clusters['c8'] = Cluster([ModelHit(ch_flgB, gene_ref=mgI_flgB, gene_status=GeneStatus.MANDATORY),
                              ModelHit(ch_tadZ, gene_ref=mgI_tadZ, gene_status=GeneStatus.ACCESSORY)],
                              models['I'], hit_weights)

    clusters['c9'] = Cluster([ModelHit(ch_abc, gene_ref=mgJ_abc, gene_status=GeneStatus.MANDATORY),
                              ModelHit(ch_tadZ, gene_ref=mgJ_tadZ, gene_status=GeneStatus.ACCESSORY)],
                              models['J'], hit_weights)

    clusters['c10'] = Cluster([ModelHit(ch_flgB, gene_ref=mgK_flgB, gene_status=GeneStatus.MANDATORY),
                               ModelHit(ch_sctn, gene_ref=mgK_sctn, gene_status=GeneStatus.ACCESSORY)],
                               models['K'], hit_weights)
    clusters['c11'] = Cluster([ModelHit(ch_flgB, gene_ref=mgL_flgB, gene_status=GeneStatus.MANDATORY),
                               ModelHit(ch_sctn_flg, gene_ref=mgL_sctn_flg, gene_status=GeneStatus.MANDATORY)],
                               models['L'], hit_weights)
    clusters['c12'] = Cluster([ModelHit(ch_sctj_flg, gene_ref=mgL_sctj_flg, gene_status=GeneStatus.ACCESSORY),
                               ModelHit(ch_sctn, gene_ref=mgL_sctn, gene_status=GeneStatus.ACCESSORY)],
                              models['L'], hit_weights)
    clusters['c13'] = Cluster([Loner(ch_sctn, gene_ref=mgL_sctn, gene_status=GeneStatus.ACCESSORY)],
                              models['L'], hit_weights)

    clusters['c14'] = Cluster([ModelHit(ch_sctj, mgM_sctj, gene_status=GeneStatus.MANDATORY),
                              MultiSystem(ch_sctn, gene_ref=mgM_sctn, gene_status=GeneStatus.ACCESSORY),
                              ModelHit(ch_gspd, gene_ref=mgM_gspd, gene_status=GeneStatus.ACCESSORY)
                              ],
                              models['M'], hit_weights)
    clusters['c15'] = Cluster([ModelHit(ch_tadZ, gene_ref=mgM_tadZ, gene_status=GeneStatus.ACCESSORY),
                               ModelHit(ch_abc, gene_ref=mgM_abc, gene_status=GeneStatus.ACCESSORY)
                               ],
                              models['M'], hit_weights)
    clusters['c16'] = Cluster([MultiSystem(ch_sctn, gene_ref=mgM_sctn, gene_status=GeneStatus.ACCESSORY)],
                              models['M'], hit_weights)

    clusters['c17'] = Cluster([ModelHit(ch_flgB, mgL_flgB, GeneStatus.MANDATORY),
                               ModelHit(ch_sctn_flg, mgL_sctn_flg, GeneStatus.MANDATORY)],
                               models['N'], hit_weights)
    clusters['c18'] = Cluster([ModelHit(ch_sctj, mgN_sctj, GeneStatus.MANDATORY),
                               ModelHit(ch_sctj_flg, mgL_sctj_flg, GeneStatus.MANDATORY)],
                               models['N'], hit_weights)
    clusters['c19'] = Cluster([Loner(ch_sctn, mgL_sctn, GeneStatus.ACCESSORY)],
                              models['N'], hit_weights)
    clusters['c20'] = Cluster([Loner(ch_tadZ, mgN_tadZ, GeneStatus.ACCESSORY)],
                              models['N'], hit_weights)
    clusters['c21'] = Cluster([ModelHit(ch_sctj, mgO_sctj, GeneStatus.MANDATORY),
                               ModelHit(ch_abc, mgO_abc, GeneStatus.NEUTRAL),
                               ModelHit(ch_tadZ, mgO_tadZ, GeneStatus.ACCESSORY)],
                              models['O'], hit_weights)
    clusters['c22'] = Cluster([ModelHit(ch_sctn_flg, mgO_sctn_flg, GeneStatus.ACCESSORY),
                               ModelHit(ch_gspd, mgO_gspd, GeneStatus.MANDATORY),
                               ModelHit(ch_tadZ, mgO_tadZ, GeneStatus.ACCESSORY)],
                              models['O'], hit_weights)
    clusters['c23'] = Cluster([Loner(ch_gspd, mgO_gspd, gene_status=GeneStatus.MANDATORY)],
                             models['O'], hit_weights)
    clusters['c24'] = Cluster([MultiSystem(ch_gspd, mgO_gspd, gene_status=GeneStatus.MANDATORY)],
                              models['O'], hit_weights)
    clusters['c25'] = Cluster([MultiSystem(ch_sctn, mgO_sctn, gene_status=GeneStatus.ACCESSORY)],
                              models['O'], hit_weights)
    clusters['c26'] = Cluster([MultiSystem(ch_sctj_flg, mgO_sctj_flg, gene_status=GeneStatus.MANDATORY)],
                              models['O'], hit_weights)
    return models, clusters


def _build_systems(models, clusters, cfg):

    systems = {}
    # we need to tweek the replicon_id to have stable ressults
    # whatever the number of tests ran
    # or the tests order
    systems['A'] = System(models['A'], [clusters['c1'], clusters['c2']], cfg.redundancy_penalty())  # 5 hits
    systems['A'].id = "replicon_id_A"
    systems['B'] = System(models['B'], [clusters['c3']], cfg.redundancy_penalty())  # 3 hits
    systems['B'].id = "replicon_id_B"
    systems['C'] = System(models['C'], [clusters['c4']], cfg.redundancy_penalty())  # 4 hits
    systems['C'].id = "replicon_id_C"
    systems['D'] = System(models['D'], [clusters['c5']], cfg.redundancy_penalty())  # 2 hits
    systems['D'].id = "replicon_id_D"
    systems['E'] = System(models['E'], [clusters['c6']], cfg.redundancy_penalty())  # 1 hit
    systems['E'].id = "replicon_id_E"
    systems['F'] = System(models['F'], [clusters['c7']], cfg.redundancy_penalty())  # 1 hit
    systems['F'].id = "replicon_id_F"
    systems['G'] = System(models['G'], [clusters['c4']], cfg.redundancy_penalty())  # 4 hits
    systems['G'].id = "replicon_id_G"
    systems['H'] = System(models['H'], [clusters['c5']], cfg.redundancy_penalty())  # 2 hits
    systems['H'].id = "replicon_id_H"
    systems['I'] = System(models['I'], [clusters['c8']], cfg.redundancy_penalty())  # 2 hits
    systems['I'].id = "replicon_id_I"
    systems['J'] = System(models['J'], [clusters['c9']], cfg.redundancy_penalty())  # 2 hits
    systems['J'].id = "replicon_id_J"
    systems['K'] = System(models['K'], [clusters['c10']], cfg.redundancy_penalty())  # 2 hits
    systems['K'].id = "replicon_id_K"
    return systems


class SolutionTest(MacsyTest):

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        self.cfg = Config(MacsyDefaults(), args)
        # we need to reset the ProfileFactory
        # because it's a like a singleton
        # so other tests are influenced by ProfileFactory and it's configuration
        # for instance search_genes get profile without hmmer_exe
        self.profile_factory = ProfileFactory(self.cfg)
        self.models, self.clusters = _build_clusters(self.cfg, self.profile_factory)
        self.systems = _build_systems(self.models, self.clusters, self.cfg)


    def test_gt(self):
        s1 = Solution([self.systems[k] for k in 'AB'])   # 5 + 3     = 8 hits, 2 syst, 0.875, [2, 2, 3, 5, 5, 6, 7, 8]
        s2 = Solution([self.systems[k] for k in 'HIJ'])  # 2 + 2 + 2 = 6 hits, 3 syst, 0.722, [2, 9, 7, 8, 8, 9]
        s3 = Solution([self.systems[k] for k in 'DC'])   # 4 + 2     = 6 hits, 2 syst  0.900, [2, 9, 3, 6, 7, 8]
        s4 = Solution([self.systems[k] for k in 'CF'])   # 4 + 1     = 5 hits, 2 syst, 0.900, [3, 6, 7, 8, 9]
        s5 = Solution([self.systems[k] for k in 'EC'])   # 1 + 4     = 5 hits, 2 syst, 0.900, [3, 3, 6, 7, 8]
        s6 = Solution([self.systems[k] for k in 'DB'])   # 2 + 3     = 5 hits, 2 syst, 0.875, [2, 9, 6, 7, 8]
        self.assertGreater(s1, s2)  # s1 more hits than s2
        self.assertGreater(s2, s3)  # s2 more systems than s3
        self.assertGreater(s5, s6)  # s5 greater awholeness than s6
        self.assertGreater(s4, s5)  # s4 "bigger positions than s5

    def test_lt(self):
        s1 = Solution([self.systems[k] for k in 'AB'])   # 5 + 3     = 8 hits, 2 syst, 0.875, [2, 2, 3, 5, 5, 6, 7, 8]
        s2 = Solution([self.systems[k] for k in 'HIJ'])  # 2 + 2 + 2 = 6 hits, 3 syst, 0.722, [2, 9, 7, 8, 8, 9]
        s3 = Solution([self.systems[k] for k in 'DC'])   # 4 + 2     = 6 hits, 2 syst  0.900, [2, 9, 3, 6, 7, 8]
        s4 = Solution([self.systems[k] for k in 'CF'])   # 4 + 1     = 5 hits, 2 syst, 0.900, [3, 6, 7, 8, 9]
        s5 = Solution([self.systems[k] for k in 'EC'])   # 1 + 4     = 5 hits, 2 syst, 0.900, [3, 3, 6, 7, 8]
        s6 = Solution([self.systems[k] for k in 'DB'])   # 2 + 3     = 5 hits, 2 syst, 0.875, [2, 9, 6, 7, 8]
        self.assertLess(s2, s1)  # s1 more hits than s2
        self.assertLess(s3, s2)  # s2 more systems than s3
        self.assertLess(s6, s5)  # s5 greater awholeness than s6
        self.assertLess(s5, s4)  # s4 "bigger positions than s5


    def test_sorting(self):
        s1 = Solution([self.systems[k] for k in 'AB'])   # 5 + 3     = 8 hits, 2 syst, 0.875, [2, 2, 3, 5, 5, 6, 7, 8]
        s2 = Solution([self.systems[k] for k in 'HIJ'])  # 2 + 2 + 2 = 6 hits, 3 syst, 0.722, [2, 9, 7, 8, 8, 9]
        s3 = Solution([self.systems[k] for k in 'DC'])   # 4 + 2     = 6 hits, 2 syst  0.900, [2, 9, 3, 6, 7, 8]
        s4 = Solution([self.systems[k] for k in 'CF'])   # 4 + 1     = 5 hits, 2 syst, 0.900, [3, 6, 7, 8, 9]
        s5 = Solution([self.systems[k] for k in 'EC'])   # 1 + 4     = 5 hits, 2 syst, 0.900, [3, 3, 6, 7, 8]
        s6 = Solution([self.systems[k] for k in 'DB'])   # 2 + 3     = 5 hits, 2 syst, 0.875, [2, 9, 6, 7, 8]

        expected_order = [s6, s5, s4, s3, s2, s1]
        shuffled_sol = expected_order[:]
        random.shuffle(shuffled_sol)
        sorted_sol = sorted(shuffled_sol)
        self.assertEqual(expected_order, sorted_sol)

    def test_systems(self):
        s1 = Solution([self.systems[k] for k in 'CD'])
        self.assertEqual([self.systems['D'], self.systems['C']],
                         s1.systems)

    def test_score(self):
        s = Solution([self.systems[k] for k in 'CD'])
        self.assertEqual(s.score, 4.5)

        s = Solution([self.systems[k] for k in 'HIJ'])
        self.assertEqual(s.score, 4.5)

    def test_average_wholeness(self):
        s = Solution([self.systems[k] for k in 'CD'])
        self.assertEqual(s.average_wholeness, 0.9)

        s = Solution([self.systems[k] for k in 'HIJ'])
        self.assertEqual(s.average_wholeness, 0.7222222222222222)

    def test_hits_number(self):
        s = Solution([self.systems[k] for k in 'CD'])
        self.assertEqual(s.hits_number, 6)

        s = Solution([self.systems[k] for k in 'AB'])
        self.assertEqual(s.hits_number, 8)


    def test_hits_positions(self):
        s = Solution([self.systems[k] for k in 'CD'])
        self.assertEqual(s.hits_positions,
                         [2, 9, 3, 6, 7, 8])

        s = Solution([self.systems[k] for k in 'AB'])
        self.assertEqual(s.hits_positions,
                         [2, 2, 3, 5, 5, 6, 7, 8])

    def test_iteration(self):
        # be careful the systems are ordered in Solution
        systems = [self.systems[k] for k in 'HIJ']
        s = Solution(systems)
        got = [syst for syst in s]
        self.assertEqual(systems,
                         got)


class SolutionExplorerTest(MacsyTest):

    @classmethod
    def setUpClass(cls) -> None:
        # to turn on debugging
        # uncomment the 3 following lines
        # import macsypy
        # macsypy.init_logger()
        # macsypy.logger_set_level('DEBUG')

        pass

    def setUp(self) -> None:
        args = argparse.Namespace()
        args.sequence_db = self.find_data("base", "test_1.fasta")
        args.db_type = 'gembase'
        args.models_dir = self.find_data('models')
        self.cfg = Config(MacsyDefaults(), args)
        # we need to reset the ProfileFactory
        # because it's a like a singleton
        # so other tests are influenced by ProfileFactory and it's configuration
        # for instance search_genes get profile without hmmer_exe
        self.profile_factory = ProfileFactory(self.cfg)
        self.models, self.clusters = _build_clusters(self.cfg, self.profile_factory)
        self.systems = _build_systems(self.models, self.clusters, self.cfg)


    def test_find_best_solution(self):
        systems = [self.systems[k] for k in 'ABCD']
        sorted_syst = sorted(systems, key=lambda s: (- s.score, s.id))
        # sorted_syst = [('replicon_id_C', 3.0), ('replicon_id_B', 2.0), ('replicon_id_A', 1.5), ('replicon_id_D', 1.5)]
        # replicon_id_C ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_B ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB']
        # replicon_id_A ['hit_sctj', 'hit_sctn', 'hit_gspd', 'hit_sctj', 'hit_sctn']
        # replicon_id_D ['hit_abc', 'hit_sctn']
        # C and D are compatible 4.5
        # B and A are compatible 3.5
        # B and D are compatible 3.5
        # So the best Solution expected is C D 4.5
        best_sol, score = find_best_solutions(sorted_syst)
        expected_sol = [Solution([self.systems[k] for k in 'CD'])]
        # The order of solutions are not relevant
        # The order of systems in each solutions are not relevant
        # transform list in set to compare them

        self.assertEqual(score, 4.5)
        self.assertEqual(best_sol, expected_sol)

        systems = [self.systems[k] for k in 'ABC']
        sorted_syst = sorted(systems, key=lambda s: (- s.score, s.id))
        # sorted_syst = [('replicon_id_C', 3.0), ('replicon_id_B', 2.0), ('replicon_id_A', 1.5)]
        # replicon_id_C ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_B ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB']
        # replicon_id_A ['hit_sctj', 'hit_sctn', 'hit_gspd', 'hit_sctj', 'hit_sctn']
        # C is alone 3.0
        # B and A are compatible 3.5
        # So the best Solution expected is B and A
        best_sol, score = find_best_solutions(sorted_syst)
        expected_sol = [Solution([self.systems[k] for k in 'BA'])]
        self.assertEqual(score, 3.5)
        self.assertEqual(best_sol, expected_sol)

        systems = [self.systems[k] for k in 'BCE']
        sorted_syst = sorted(systems, key=lambda s: (- s.score, s.id))
        # sorted_syst = [('replicon_id_C', 3.0), ('replicon_id_B', 2.0), ('replicon_id_E', 0.5)]
        # replicon_id_C ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_B ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB']
        # replicon_id_E ['hit_gspd']
        # C is alone 3.0
        # B and E are compatible 2.5
        # So the best Solution expected is C
        best_sol, score = find_best_solutions(sorted_syst)
        expected_sol = [Solution([self.systems[k] for k in 'C'])]
        self.assertEqual(score, 3.0)
        self.assertEqual(best_sol, expected_sol)

        # systems = [('replicon_id_C', 3.0), ('replicon_id_B', 2.0), ('replicon_id_A', 1.5),
        #            ('replicon_id_D', 1.5), ('replicon_id_E', 0.5)]
        # replicon_id_C ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_B ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB']
        # replicon_id_A ['hit_sctj', 'hit_sctn', 'hit_gspd', 'hit_sctj', 'hit_sctn']
        # replicon_id_D ['hit_abc', 'hit_sctn']
        # replicon_id_E ['hit_gspd']
        # C and D are compatible 4.5
        # B and A are compatible 3.5
        # B and E are compatible 2.5
        # D and E are compatible 2.0
        systems = [self.systems[k] for k in 'ABCDE']
        sorted_syst = sorted(systems, key=lambda s: (- s.score, s.id))
        best_sol, score = find_best_solutions(sorted_syst)
        expected_sol = [Solution([self.systems[k] for k in 'CD'])]
        self.assertEqual(score, 4.5)
        self.assertEqual(best_sol, expected_sol)

        # systems = [('replicon_id_C', 3.0), ('replicon_id_B', 2.0), ('replicon_id_A', 1.5),
        #            ('replicon_id_D', 1.5), ('replicon_id_E', 0.5), ('replicon_id_F', 1.0)]
        # replicon_id_C ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_B ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB']
        # replicon_id_A ['hit_sctj', 'hit_sctn', 'hit_gspd', 'hit_sctj', 'hit_sctn']
        # replicon_id_D ['hit_abc', 'hit_sctn']
        # replicon_id_E ['hit_gspd']
        # replicon_id_F ['hit_abc']
        # C and D are compatible 4.5
        # C and F are compatible 4.0
        # B and A and F are compatible 4.5
        # B and D and E are compatible 4.0
        # B and E and F are compatible 3.5
        # So the best Solution expected are C D / B A F
        systems = [self.systems[k] for k in 'ABCDEF']
        sorted_syst = sorted(systems, key=lambda s: (- s.score, s.id))
        best_sol, score = find_best_solutions(sorted_syst)
        expected_sol = [Solution([self.systems[k] for k in 'BAF']),  # 3 + 5 + 1 = 9 hits
                        Solution([self.systems[k] for k in 'CD'])]   # 4 + 2 = 7 hits

        self.assertEqual(best_sol[0].score, 4.5)
        self.assertEqual(best_sol[1].score, 4.5)
        # test if the composition is right
        self.assertEqual(best_sol, expected_sol)
        # test if solution order is right

        systems = [self.systems[k] for k in 'ABCDGH']
        sorted_syst = sorted(systems, key=lambda s: (- s.score, s.id))
        # sorted_syst = [('replicon_id_C', 3.0), ('replicon_id_B', 2.0), ('replicon_id_A', 1.5), ('replicon_id_D', 1.5)
        #                ('replicon_id_G', 3.0), ('replicon_id_H', 1.5)]
        # replicon_id_A ['hit_sctj', 'hit_sctn', 'hit_gspd', 'hit_sctj', 'hit_sctn']
        # replicon_id_B ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB']
        # replicon_id_C ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_D ['hit_abc', 'hit_sctn']
        # replicon_id_G ['hit_sctj_flg', 'hit_tadZ', 'hit_flgB', 'hit_gspd']
        # replicon_id_H ['hit_abc', 'hit_sctn']
        # C and D are compatible 4.5   wholeness = 0.8  + 1.0 = 1.8
        # C and H are compatible 4.5               0.8  + 1.0 = 1.8
        # G and D are compatible 4.5               0.8  + 1.0 = 1.8
        # G and H are compatible 4.5               0.8  + 1.0 = 1.8
        # B and A are compatible 3.5               0.75 + 1.0 = 1.75
        # So the best Solution expected are C D / C H / G D / G H with score 4.5

        best_sol, score = find_best_solutions(sorted_syst)
        expected_sol = [Solution([self.systems[k] for k in 'CD']),  # 4 + 2 hits
                        Solution([self.systems[k] for k in 'CH']),  # 4 + 2
                        Solution([self.systems[k] for k in 'GD']),  # 4 + 2
                        Solution([self.systems[k] for k in 'GH'])]  # 4 + 2
        self.assertEqual(score, 4.5)
        self.assertEqual(best_sol, expected_sol)

        systems = [self.systems[k] for k in 'HJKI']
        sorted_syst = sorted(systems, key=lambda s: (- s.score, s.id))
        best_sol, score = find_best_solutions(sorted_syst)

        # check if solution is ordered by wholeness average (3rd criterion)
        # first criterion nb of hits
        # second citerion nb of systems
        # replicon_id_H ['hit_abc', 'hit_sctn']
        # replicon_id_I ['hit_flgB', 'hit_tadZ']
        # replicon_id_J ['hit_abc', 'hit_tadZ']
        # replicon_id_K ['hit_flgB', 'hit_sctn']
        #                                                             score  Nb hits  nb sys wholeness
        expected_sol = [Solution([self.systems[k] for k in 'HI']),  # 1.5 + 1.5 = 3.0    4        2       1.0
                        Solution([self.systems[k] for k in 'JK'])]  # 1.5 + 1.5 = 3.0    4        2       0.5
        self.assertEqual(score, 3.0)
        self.assertEqual(best_sol, expected_sol)


    def test_combine_clusters(self):
        ##################################################
        # with 3 regular clusters 0 loner 0 multisystyem
        ##################################################
        combinations = combine_clusters([self.clusters['c1'], self.clusters['c2'],  self.clusters['c3']],
                                       {},
                                       multi_loci=False)
        self.assertEqual(combinations,
                         [
                             (self.clusters['c1'],),
                             (self.clusters['c2'],),
                             (self.clusters['c3'],)
                         ])
        # the same in multi loci
        combinations = combine_clusters([self.clusters['c1'], self.clusters['c2'], self.clusters['c3']],
                                        {},
                                        multi_loci=True)
        exp_combs = [
                    (self.clusters['c1'],),
                    (self.clusters['c2'],),
                    (self.clusters['c3'],),
                    (self.clusters['c1'], self.clusters['c2']),
                    (self.clusters['c1'], self.clusters['c3']),
                    (self.clusters['c2'], self.clusters['c3']),
                    (self.clusters['c1'], self.clusters['c2'], self.clusters['c3'])
                 ]
        self.assertEqual(combinations, exp_combs)

        ###########################################
        # with 2 RC + 1 L not included in cluster
        ###########################################
        combinations = combine_clusters([self.clusters['c11'], self.clusters['c12']],
                                        {},
                                        multi_loci=False)
        exp_combs = [
                    (self.clusters['c11'],),
                    (self.clusters['c12'],)
                   ]
        self.assertEqual(combinations, exp_combs)
        # the same in multi loci
        combinations = combine_clusters([self.clusters['c11'], self.clusters['c12']],
                                       {},
                                       multi_loci=True)
        exp_combs = [
                    (self.clusters['c11'],),
                    (self.clusters['c12'],),
                    (self.clusters['c11'], self.clusters['c12'])
                   ]
        self.assertEqual(combinations, exp_combs)

        ##################################
        # with 2 RC + 1 L already in RC 2
        ##################################
        # c11 = flgB, sctn_flg
        # c12 = sctj_flg, sctn
        # c13 = Loner sctn
        combinations = combine_clusters([self.clusters['c11'], self.clusters['c12']],
                                        {'sctN': self.clusters['c13']},
                                        multi_loci=False)

        exp_combs = [
            (self.clusters['c11'],),
            (self.clusters['c12'],),
            (self.clusters['c11'], self.clusters['c13']),
            (self.clusters['c13'],),
        ]

        self.assertEqual(combinations, exp_combs)

        # the same in multi loci
        combinations = combine_clusters([self.clusters['c11'], self.clusters['c12']],
                                        {'sctN': self.clusters['c13']},
                                        multi_loci=True)

        exp_combs = [
            (self.clusters['c11'],),
            (self.clusters['c12'],),
            (self.clusters['c11'], self.clusters['c12']),
            (self.clusters['c11'], self.clusters['c13']),
            (self.clusters['c13'],),
        ]
        self.assertEqual(combinations, exp_combs)

        ###########################################
        # with 2 RC  with one containing a MS     #
        ###########################################
        combinations = combine_clusters([self.clusters['c14'], self.clusters['c15']],
                                        {'sctN': self.clusters['c16']},
                                        multi_loci=True)
        # c14 contains a MS
        # c15 do not contains MS
        # c16 is the artificial cluster with only sctn
        exp_combs = [
            (self.clusters['c14'],),
            (self.clusters['c15'],),
            (self.clusters['c14'], self.clusters['c15']),
            (self.clusters['c15'], self.clusters['c16']),
            (self.clusters['c16'],)]
        self.assertEqual(combinations, exp_combs)

        ###########################################
        # with 2 RC + 2 L not included in cluster
        ###########################################
        combinations = combine_clusters([self.clusters['c17'], self.clusters['c18']],
                                        {'sctn': self.clusters['c19'],
                                         'tadZ': self.clusters['c20']
                                         },
                                        multi_loci=False)
        exp_combs = [
            (self.clusters['c17'],),
            (self.clusters['c18'],),
            (self.clusters['c17'], self.clusters['c19']),
            (self.clusters['c18'], self.clusters['c19']),
            (self.clusters['c19'],),
            (self.clusters['c17'], self.clusters['c20']),
            (self.clusters['c18'], self.clusters['c20']),
            (self.clusters['c20'],),
            (self.clusters['c17'], self.clusters['c19'], self.clusters['c20']),
            (self.clusters['c18'], self.clusters['c19'], self.clusters['c20']),
            (self.clusters['c19'], self.clusters['c20'])
        ]

        self.assertEqual(combinations, exp_combs)
        # the same in multi loci
        combinations = combine_clusters([self.clusters['c17'], self.clusters['c18']],
                                        {'sctn': self.clusters['c19'],
                                         'tadZ': self.clusters['c20']
                                         },
                                        multi_loci=True)
        exp_combs = [
            (self.clusters['c17'],),
            (self.clusters['c18'],),
            (self.clusters['c17'], self.clusters['c18']),
            (self.clusters['c17'], self.clusters['c19']),
            (self.clusters['c18'], self.clusters['c19']),
            (self.clusters['c17'],self.clusters['c18'], self.clusters['c19']),
            (self.clusters['c19'],),
            (self.clusters['c17'], self.clusters['c20']),
            (self.clusters['c18'], self.clusters['c20']),
            (self.clusters['c17'],self.clusters['c18'], self.clusters['c20']),
            (self.clusters['c20'],),
            (self.clusters['c17'], self.clusters['c19'], self.clusters['c20']),
            (self.clusters['c18'], self.clusters['c19'], self.clusters['c20']),
            (self.clusters['c17'], self.clusters['c18'], self.clusters['c19'], self.clusters['c20']),
            (self.clusters['c19'], self.clusters['c20'])
        ]
        # print("\n####################################################################")
        # for comb in exp_combs:
        #     print([c.id for c in comb])
        # print("==============================")
        # for comb in combinations:
        #     print([c.id for c in comb])
        self.assertEqual(combinations, exp_combs)


    def test_combine_multisystems(self):
        # key    id     hits
        # c21    c47    sctj abc tadz
        # c22    c48    sctn_flg gspd tadz
        # c23    c49    Loner gspd
        # c24    c50    MS gspd
        # c25    c51    MS sctn
        # c26    c52    MS sctj_flg

        rejected_clst = RejectedCandidate(self.models['O'],
                                         [self.clusters['c21']],
                                         'fake_reason')
        combinations = combine_multisystems(rejected_clst,
                                            [self.clusters['c24'], self.clusters['c25'], self.clusters['c26']]
                                            )
        exp_combs = [
            (self.clusters['c21'], self.clusters['c24']),
            (self.clusters['c21'], self.clusters['c25']),
            (self.clusters['c21'], self.clusters['c24'], self.clusters['c25'])
            ]
        self.assertEqual(combinations, exp_combs)

        rejected_clst = RejectedCandidate(self.models['O'],
                                         [self.clusters['c22']],
                                         'fake_reason')
        combinations = combine_multisystems([rejected_clst],
                                            [self.clusters['c24'], self.clusters['c25'], self.clusters['c26']]
                                            )
        exp_combs = [
            (self.clusters['c22'], self.clusters['c26'])
            ]
        self.assertEqual(combinations, exp_combs)

        rejected_clst = RejectedCandidate(self.models['O'],
                                         [self.clusters['c21'], self.clusters['c23']],
                                         'fake_reason')
        combinations = combine_multisystems([rejected_clst],
                                            [self.clusters['c24'], self.clusters['c25'], self.clusters['c26']]
                                            )
        exp_combs = [
            (self.clusters['c21'], self.clusters['c23'], self.clusters['c25'])
            ]
        self.assertEqual(combinations, exp_combs)
        # print("\n####################################################################")
        # for comb in exp_combs:
        #     print([c.id for c in comb])
        # print("==============================")
        # for comb in combinations:
        #     print([c.id for c in comb])