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
import importlib

def main():
    directory_of_tests=os.path.dirname(os.path.realpath(__file__))
    
    basic_suite = unittest.TestLoader().discover(directory_of_tests,pattern='basic_tests_*.py')
    advanced_suite = unittest.TestLoader().discover(directory_of_tests, pattern='advanced_tests_*.py')
    full_suite = unittest.TestSuite([basic_suite,advanced_suite])
   
    return full_suite 
