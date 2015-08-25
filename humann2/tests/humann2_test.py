#!/usr/bin/env python

"""
HUMAnN2 : HMP Unified Metabolic Analysis Network 2

HUMAnN2 is a pipeline for efficiently and accurately determining
the coverage and abundance of microbial pathways in a community
from metagenomic data. Sequencing a metagenome typically produces millions
of short DNA/RNA reads.

This software is used to test the HUMAnN2 pipeline.

Dependencies: HUMAnN2 (and all HUMAnN2 dependencies)

"""

import os
import sys
import unittest

# Try to load one of the humann2 modules to check the installation
try:
    from humann2 import check
except ImportError:
    sys.exit("CRITICAL ERROR: Unable to find the HUMAnN2 python package." +
        " Please check your install.") 

# Check the python version
check.python_version()

import argparse

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "HUMAnN2 Test\n",
        formatter_class=argparse.RawTextHelpFormatter,
        prog="humann2_test")
    parser.add_argument(
        "-r","--run-functional-tests", 
        help="also run the functional tests\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "-b","--bypass-unit-tests", 
        help="do not run the unit tests\n", 
        action="store_true",
        default=False)
    
    return parser.parse_args()

def get_testdirectory():
    """ Return the location of all of the tests """
    
    return os.path.dirname(os.path.abspath(__file__))

def get_funtionaltests():
    """ Return all of the functional tests """
    
    directory_of_tests=get_testdirectory()
    
    functional_suite = unittest.TestLoader().discover(directory_of_tests,pattern='functional_tests_*.py')
    
    return functional_suite

def get_unittests():
    """ Return all of the unit tests """
    
    directory_of_tests=get_testdirectory()
    
    basic_suite = unittest.TestLoader().discover(directory_of_tests,pattern='basic_tests_*.py')
    advanced_suite = unittest.TestLoader().discover(directory_of_tests, pattern='advanced_tests_*.py')    
    
    return [basic_suite, advanced_suite]

def unittests_suite_only():
    """ Return a TestSuite of just the unit tests """
    
    return unittest.TestSuite(get_unittests())

def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    # Get the unittests
    test_suites=[]
    if not args.bypass_unit_tests:
        test_suites=get_unittests()
    
    # Get the functional tests if requested
    if args.run_functional_tests:
        test_suites+=get_funtionaltests()

    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(test_suites))

