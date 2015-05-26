import unittest
import logging

import cfg
import utils

from humann2.quantify import modules

class TestHumann2QuantifyModulesFunctions(unittest.TestCase):
    """
    Test the functions found in humann2.quantify.modules
    """
    
    def setUp(self):
        # set up nullhandler for logger
        logging.getLogger('humann2.quantify.modules').addHandler(logging.NullHandler())
        
    def test_harmonic_mean(self):
        """
        Test the harmonmic mean function
        """
        
        values=[1,2,3,4]
        result=modules.harmonic_mean(values)
        expect_result=len(values)/sum(1.0/v for v in values)
        
        self.assertEqual(result, expect_result)
        
    def test_harmonic_mean_empty_list(self):
        """
        Test the harmonmic mean function with an empty list
        """
        
        result=modules.harmonic_mean([])
        expect_result=0
        
        self.assertEqual(result, expect_result)
        
    def test_harmonic_mean_one_zero(self):
        """
        Test the harmonmic mean function with a set of values with on zero
        """
        
        result=modules.harmonic_mean([1,2,3,0,4])
        expect_result=0
        
        self.assertEqual(result, expect_result)
        
