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

import os
import logging
from time import strftime

from macsypy.config import MacsyDefaults

from tests import MacsyTest


class TestMacsyDefaults(MacsyTest):

    def setUp(self):
        self.defaults = {'cfg_file': None,
                         'coverage_profile': 0.5,
                         'e_value_search': 0.1,
                         'no_cut_ga': False,
                         'db_type': None,
                         'hmmer': 'hmmsearch',
                         'i_evalue_sel': 0.001,
                         'idx': False,
                         'inter_gene_max_space': None,
                         'log_file': 'macsyfinder.log',
                         'log_level': logging.INFO,
                         'max_nb_genes': None,
                         'min_genes_required': None,
                         'min_mandatory_genes_required': None,
                         'models': [],
                         'models_dir': os.path.normpath(os.path.join(os.path.dirname(__file__), 'data', 'models')),
                         'multi_loci': set(),
                         'mute': False,
                         'out_dir': None,
                         'previous_run': None,
                         'profile_suffix': '.hmm',
                         'quiet': 0,
                         'relative_path': False,
                         'replicon_topology': 'circular',
                         'res_extract_suffix': '.res_hmm_extract',
                         'res_search_dir': os.getcwd(),
                         'res_search_suffix': '.search_hmm.out',
                         'sequence_db': None,
                         'topology_file': None,
                         'verbosity': 0,
                         'worker': 1,
                         'mandatory_weight': 1.0,
                         'accessory_weight': .5,
                         'neutral_weight': 0.0,
                         'exchangeable_weight': .8,
                         'itself_weight': 1.0,
                         'redundancy_penalty': 1.5,
                         'loner_multi_system_weight': 0.7
                         }

    def test_MacsyDefaults(self):
        defaults = MacsyDefaults()
        self.assertDictEqual(defaults, self.defaults)

        new_defaults = {k: v for k, v in self.defaults.items()}
        new_defaults['previous_run'] = True
        new_defaults['worker'] = 5
        defaults = MacsyDefaults(previous_run=True, worker=5)
        self.assertDictEqual(defaults, new_defaults)


    def test_MacsyDefaults_with_MACSY_DATA(self):
        import macsypy.config
        macsydata = macsypy.config.__MACSY_DATA__
        macsypy.config.__MACSY_DATA__ = 'niportnaoik'
        self.defaults['models_dir'] = 'niportnaoik/data/models'
        try:
            defaults = MacsyDefaults()
            self.assertDictEqual(defaults, self.defaults)
        finally:
            macsypy.config.__MACSY_DATA__ = macsydata
