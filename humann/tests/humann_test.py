#!/usr/bin/env python

"""
HUMAnN : HMP Unified Metabolic Analysis Network

HUMAnN is a pipeline for efficiently and accurately determining
the coverage and abundance of microbial pathways in a community
from metagenomic data. Sequencing a metagenome typically produces millions
of short DNA/RNA reads.

This software is used to test the HUMAnN pipeline.

Dependencies: HUMAnN (and all HUMAnN dependencies)

"""

import os
import sys
import unittest

# Try to load one of the humann modules to check the installation
try:
    from humann import check
except ImportError:
    sys.exit("CRITICAL ERROR: Unable to find the HUMAnN python package." +
        " Please check your install.") 

# Check the python version
check.python_version()

# check for the biom install
biom_installed=True
try:
    import biom
    import h5py
except ImportError:
    biom_installed=False

import argparse

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "HUMAnN Test\n",
        formatter_class=argparse.RawTextHelpFormatter,
        prog="humann_test")
    parser.add_argument(
        "--run-functional-tests-tools", 
        help="run the functional tests for tools\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "--run-functional-tests-end-to-end", 
        help="run the humann end to end functional tests\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "--bypass-unit-tests", 
        help="do not run the unit tests\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "--run-all-tests", 
        help="run all tests\n", 
        action="store_true",
        default=False)
    
    return parser.parse_args()

def get_testdirectory():
    """ Return the location of all of the tests """
    
    return os.path.dirname(os.path.abspath(__file__))

def get_funtionaltests_tools(run_all_tests=None):
    """ Return all of the functional tests for tools"""
    
    directory_of_tests=get_testdirectory()
    
    functional_suite = [unittest.TestLoader().discover(directory_of_tests,pattern='functional_tests_tools*.py')]
    
    # if biom is installed, add the functional tools tests with biom
    if biom_installed or run_all_tests:
        functional_suite += [unittest.TestLoader().discover(directory_of_tests,pattern='functional_tests_biom_tools*.py')]
    
    return functional_suite

def get_funtionaltests_other(run_all_tests=None):
    """ Return all of the other functional tests """
    
    directory_of_tests=get_testdirectory()
    
    functional_suite = [unittest.TestLoader().discover(directory_of_tests,pattern='functional_tests_humann*.py')]
    
    # if biom is installed, add the functional end to end tests with biom
    if biom_installed or run_all_tests:
        functional_suite += [unittest.TestLoader().discover(directory_of_tests,pattern='functional_tests_biom_humann*.py')]
    
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
    if args.run_functional_tests_tools or args.run_all_tests:
        test_suites+=get_funtionaltests_tools(args.run_all_tests)
    if args.run_functional_tests_end_to_end or args.run_all_tests:
        print("\n\nPlease note running functional end to end tests requires all dependencies of HUMAnN.\n")
        test_suites+=get_funtionaltests_other(args.run_all_tests)

    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(test_suites))

