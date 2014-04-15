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
import glob

_prefix_data = '$PREFIXDATA'
if 'MACSY_HOME' in os.environ and os.environ['MACSY_HOME']:
    _prefix_data = os.path.join(os.environ['MACSY_HOME'], 'data')


class ProfilesRegistry(object):
    """
    classdocs
    """


    def __init__(self, cfg):
        """    
        Constructor
        """
        self._register = {}
        global_path = os.path.join(_prefix_data, 'macsyfinder' , 'profiles')
        self._fill_profile(global_path , cfg)
        local_path = cfg.profile_dir
        if local_path:
            self._fill_profile(local_path, cfg)
        
    def _fill_profile(self, dir_path, cfg):
        for path in glob.glob(os.path.join(dir_path, '*' + cfg.profile_suffix)):
            name = os.path.basename(path)
            name = name[:-1 * len(cfg.profile_suffix)]
            self._register[name] = path
            
    def __getattr__(self, name):
        return getattr(self._register, name)


class DefinitionsRegistry(object):
    """
    classdocs
    """


    def __init__(self, cfg):
        """    
        Constructor
        """
        self._register = {}
        global_path = os.path.join(_prefix_data, 'macsyfinder' , 'DEF')
        self._fill_def(global_path , cfg)
        local_path = cfg.def_dir
        if local_path:
            self._fill_def(local_path, cfg)
        
    def _fill_def(self, dir_path, cfg):
        for path in glob.glob(os.path.join(dir_path, '*.xml')):
            name = os.path.basename(path)
            name = os.path.splitext(name)[0]
            self._register[name] = path
            
    def __getattr__(self, name):
        return getattr(self._register, name)
        
        
