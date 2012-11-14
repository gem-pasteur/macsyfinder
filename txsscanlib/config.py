# -*- coding: utf-8 -*-


#===============================================================================
# Created on Nov 12, 2012
# 
# @author: bertrand Néron
# @contact: bneron <at> pasteur <dot> fr
# @organization: organization_name
# @license: license
#===============================================================================

import os
import inspect
from ConfigParser import SafeConfigParser

_prefix_path = '$PREFIX'

class Config(object):
    """
    parse configuration files and handle the configuration according the file location precedence
    /etc/txsscan.conf < ~/.txsscan.conf < .txsscan.conf
    if a configuration file is given on the command line, only this file is used.
    in fine the arguments passed have the highest priority
    """
    
    options = ('hmmer_exe' , 'e_value_res', 'res_search_dir', 'res_search_suffix', 'profile_dir', 'profile_suffix', 'res_extract_suffix', 
               'log_level')

    def __init__(self, cfg_file = "", 
                hmmer_exe = None,
                e_value_res = None,
                res_search_dir = None,
                res_search_suffix = None,
                profile_dir = None,
                profile_suffix = None,
                res_extract_suffix = None,
                log_level = None
                ):
        """
        :param hmmer_exe: the hmmsearch executabe
        :type hmmer_exe: string
        :param e_value_res: à déterminer
        :type  e_value_res: float
        :param res_search_dir: à déterminer
        :type  res_search_dir: string
        :param res_search_suffix: à déterminer
        :type  res_search_suffix: string
        :param profile_dir: à déterminer
        :type  profile_dir: string
        :param profile_suffix: à déterminer
        :type  profile_suffix: string
        :param res_extract_suffix: à déterminer
        :type  res_extract_suffix: string
        :param log_level: the level of log output
        :type log_level: int
        """
        
        config_files = [cfg_file] if cfg_file else [ os.path.join( _prefix_path , '/etc/txsscan/txsscan.conf'),
                                                     os.path.expanduser('~/.txsscan/txsscan.conf'),
                                                      '.txsscan.conf']
        print config_files
        self.parser = SafeConfigParser(defaults={'hmmer_exe' : 'hmmsearch',
                                            'e_value_res' : 1,
                                            'res_search_dir' : './datatest/res_search',
                                            'res_search_suffix' : '.search_hmm.out',
                                            'profile_dir' : './profiles',
                                            'profile_suffix' : '.fasta-aln_edit.hmm',
                                            'res_extract_suffix' : '.res_hmm_extract'
                                            },
                                  )
        used_files = self.parser.read(config_files)
        print "used_files = ",used_files
        frame = inspect.currentframe()
        args, _, _, values = inspect.getargvalues(frame)

        self.cmde_line_opt = {}
        for arg in args:
            if arg in self.options and values[arg]:
                self.cmde_line_opt[arg] = values[arg]
        
    @property
    def hmmer_exe(self):
        return self.parser.get('hmmer', 'hmmer_exe', vars = self.cmde_line_opt)

    @property
    def e_value_res(self):
        try:
            return float( self.parser.get('hmmer', 'e_value_res', vars = self.cmde_line_opt))
        except ValueError:
            msg =  "Invalid value for hmmer e_value_res :{0}: (float expected)".format(self.parser.get('hmmer', 'e_value_res', vars = self.cmde_line_opt))
            self.log.error( msg )
            raise ValueError( msg )
        
    @property
    def res_search_dir(self):
        return self.parser.get('directories', 'res_search_dir', vars = self.cmde_line_opt)

    @property
    def res_search_suffix(self):
        return self.parser.get('directories', 'res_search_suffix', vars = self.cmde_line_opt)

    @property
    def profile_dir(self):
        return self.parser.get('directories', 'profile_dir', vars = self.cmde_line_opt)
    
    @property
    def profile_suffix(self):
        return self.parser.get('directories', 'profile_suffix', vars = self.cmde_line_opt)

    @property
    def res_extract_suffix(self):
        return self.parser.get('directories', 'res_extract_suffix', vars = self.cmde_line_opt)

    @property
    def log_level(self):
        try:
            level = int( self.parser.get('general', 'log_level', vars = self.cmde_line_opt))
        except ValueError:
            level = self.parser.get('general', 'log_level', vars = self.cmde_line_opt)
            import logging
            try:
                level = getattr( logging, level.upper() )
            except AttributeError:
                level = logging.WARNING
                
        return  level
