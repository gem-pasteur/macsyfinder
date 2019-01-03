# -*- coding: utf-8 -*-

################################################################################
# MacSyFinder - Detection of macromolecular systems in protein datasets        #
#               using systems modelling and similarity search.                 #
# Authors: Sophie Abby, Bertrand Néron                                         #
# Copyright © 2014  Institut Pasteur, Paris.                                   #
# See the COPYRIGHT file for details                                           #
#                                                                              #
# MacsyFinder is distributed under the terms of the GNU General Public License #
# (GPLv3). See the COPYING file for details.                                   #
################################################################################

import os
import logging
from time import strftime

from macsypy.config import MacsyDefaults

from tests import MacsyTest


class TestMacsyDefaults(MacsyTest):

    def setUp(self):
        self.defaults = {'cfg_file': None,
                         'coverage_profile': 0.5,
                         'e_value_search': 1.0,
                         'db_type': None,
                         'hmmer': 'hmmsearch',
                         'i_evalue_sel': 0.001,
                         'idx': False,
                         'inter_gene_max_space': None,
                         'log_file': 'macsyfinder.log',
                         'log_level': logging.WARNING,
                         'max_nb_genes': None,
                         'min_genes_required': None,
                         'min_mandatory_genes_required': None,
                         'models': [],
                         'models_dir': os.path.normpath(os.path.join(os.path.dirname(__file__), '..', 'data', 'models')),
                         'multi_loci': set(),
                         'out_dir': None,
                         'previous_run': False,
                         'profile_suffix': '.hmm',
                         'relative_path': False,
                         'replicon_topology': 'circular',
                         'res_extract_suffix': '.res_hmm_extract',
                         'res_search_dir': os.getcwd(),
                         'res_search_suffix': '.search_hmm.out',
                         'sequence_db': None,
                         'topology_file': None,
                         'verbosity': 0,
                         'worker': 1,
                         }

    def test_MacsyDefaults(self):
        defaults = MacsyDefaults()
        self.assertDictEqual(defaults, self.defaults)

        new_defaults = {k: v for k, v in self.defaults.items()}
        new_defaults['previous_run'] = True
        new_defaults['worker'] = 5
        defaults = MacsyDefaults(previous_run=True,
                                 worker=5)
        self.assertDictEqual(defaults, new_defaults)


    def test_MacsyDefaults_with_MACSY_DATA(self):
        import macsypy.config
        macsydata = macsypy.config.__MACSY_DATA__
        macsypy.config.__MACSY_DATA__ = os.path.normpath(os.path.join(os.path.dirname(__file__), '..'))
        try:
            defaults = MacsyDefaults()
            self.assertDictEqual(defaults, self.defaults)
        finally:
            macsypy.config.__MACSY_DATA__ = macsydata
