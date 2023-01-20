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
import logging
import colorlog
import tempfile

import macsypy
from tests import MacsyTest


class Test(MacsyTest):

    def tearDown(self) -> None:
        logger = logging.getLogger('macsypy')
        del logger.manager.loggerDict['macsypy']

    def test_init_logger_default(self):
        handlers = macsypy.init_logger()
        self.assertEqual(len(handlers), 1)
        self.assertTrue(isinstance(handlers[0], logging.StreamHandler))
        logger = logging.getLogger('macsypy')
        self.assertEqual(logger.getEffectiveLevel(),
                         logging.WARNING)

    def test_init_logger_no_out(self):
        handlers = macsypy.init_logger(out=False)
        self.assertEqual(len(handlers), 1)
        self.assertTrue(isinstance(handlers[0], logging.NullHandler))

    def test_init_logger_logfile(self):
        with tempfile.NamedTemporaryFile() as fp:
            handlers = macsypy.init_logger(log_file=fp.name)
            self.assertEqual(len(handlers), 2)
            self.assertTrue(isinstance(handlers[1], logging.FileHandler))

    def test_logger_set_level_default(self):
        macsypy.init_logger()
        macsypy.logger_set_level()
        logger = logging.getLogger('macsypy')
        self.assertEqual(logger.getEffectiveLevel(),
                         logging.INFO)

    def test_logger_set_level_error(self):
        macsypy.init_logger()
        macsypy.logger_set_level(level='ERROR')
        logger = logging.getLogger('macsypy')
        self.assertEqual(logger.getEffectiveLevel(),
                         logging.ERROR)

    def test_logger_set_level_bad_level(self):
        macsypy.init_logger()
        with self.assertRaises(ValueError) as ctx:
            macsypy.logger_set_level(level=-1)
        self.assertEqual(str(ctx.exception),
                         'Level must be NOTSET, DEBUG, INFO, WARNING, ERROR, CRITICAL or a positive integer')

    def test_logger_set_level_handlers(self):
        macsypy.init_logger()
        macsypy.logger_set_level(level='DEBUG')
        logger = colorlog.getLogger('macsypy')
        self.assertEqual(logger.handlers[0].formatter.log_colors,
            {'DEBUG': 'cyan', 'INFO': 'green', 'WARNING': 'yellow', 'ERROR': 'red', 'CRITICAL': 'bold_red'}
                         )