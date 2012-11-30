#! /usr/bin/env python
# -*- coding: utf-8 -*-

#===============================================================================
# Created on Nov 30, 2012
# 
# @author: bneron
# @contact: user_email
# @organization: organization_name
# @license: license
#===============================================================================

import os
import sys

TXSSCAN_HOME = os.path.abspath('..')
if not TXSSCAN_HOME in sys.path: 
    sys.path.append(os.path.abspath( '..'))

import unittest
     
from argparse import ArgumentParser    
parser = ArgumentParser()
parser.add_argument("tests",
                    nargs = '*',
                    help = "name of test to execute")

parser.add_argument("-v", "--verbose" , 
                    dest= "verbosity" , 
                    action="count" , 
                    help= "set the verbosity level of output",
                    default = 0
                    )

args = parser.parse_args()

if not args.tests:
    suite = unittest.TestLoader().discover('.') 
else:
    suite = unittest.TestSuite()
    for test in args.tests: 
        if os.path.exists(test):
            if os.path.isfile(test):
                fpath, fname =  os.path.split( test )
                suite.addTests(unittest.TestLoader().discover(fpath , pattern = fname )) 
            elif os.path.isdir(test ):  
                suite.addTests(unittest.TestLoader().discover(test )) 
        else:
            sys.stderr.write(  "%s : no such file or directory\n" % test) 

unittest.TextTestRunner(verbosity = args.verbosity).run(suite)
