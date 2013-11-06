# -*- coding: utf-8 -*-

#=========================
# Created on Nov 5, 2013
#
# @author: bneron
# @contact: bneron <at> pasteur <dot> fr
# @organization: Institut Pasteur
# @license: license
#=========================

import os
import glob

_prefix_path = '$PREFIX'
if os.environ['TXSSCAN_HOME']:
    _prefix_path = os.environ['TXSSCAN_HOME']


class ProfilesRegistry(object):
    """
    classdocs
    """


    def __init__(self, cfg):
        """    
        Constructor
        """
        self._register = {}
        global_path = os.path.join(_prefix_path, 'txsscan' , 'profiles')
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
        global_path = os.path.join(_prefix_path, 'txsscan' , 'DEF')
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
        
        
