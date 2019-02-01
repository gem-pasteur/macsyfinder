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


from macsypy.search_systems import SystemNameGenerator
from tests import MacsyTest


class Test(MacsyTest):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_getSystemName(self):
        sng = SystemNameGenerator()
        name = sng.get_system_name('foo', 'bar')
        self.assertEqual(name, 'foo_bar_1')

        name = sng.get_system_name('foo', 'bar')
        self.assertEqual(name, 'foo_bar_2')
