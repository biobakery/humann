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
        
    def test_gap_fill_less_than_threshold(self):
        """
        Test the gap fill function, with a set of scores that are less than the
        threshold of 75% to apply gap filling
        """
        
        key_reactions=["A","B","C"]
        reaction_scores={ "A": 1, "B": 2 }
        
        gap_filled_reaction_scores=modules.gap_fill(key_reactions, reaction_scores)
        
        self.assertDictEqual(gap_filled_reaction_scores, reaction_scores)
        
    def test_gap_fill_greater_than_threshold(self):
        """
        Test the gap fill function, with a set of scores that are greater than the
        threshold of 75% to apply gap filling
        """
        
        key_reactions=["A","B","C","D","E"]
        reaction_scores={ "A": 1, "B": 2 , "C": 2, "D": 1}
        
        expected_result={ "A": 1, "B": 2 , "C": 2, "D": 1, "E": 1}
        
        gap_filled_reaction_scores=modules.gap_fill(key_reactions, reaction_scores)
        
        self.assertDictEqual(gap_filled_reaction_scores, expected_result)

    def test_gap_fill_equal_threshold(self):
        """
        Test the gap fill function, with a set of scores that are equal to the
        threshold of 75% to apply gap filling
        """
        
        key_reactions=["A","B","C","D"]
        reaction_scores={ "A": 1, "B": 2 , "C": 2}
        
        expected_result={ "A": 1, "B": 2 , "C": 2, "D": 1}
        
        gap_filled_reaction_scores=modules.gap_fill(key_reactions, reaction_scores)
        
        self.assertDictEqual(gap_filled_reaction_scores, expected_result)
        
    def test_gap_fill_optional_reactions_less_than_threshold(self):
        """
        Test the gap fill function, with a set of scores that include optional reactions
        where just considering the required reactions it is less than the threshold
        for gap filling
        """
        
        key_reactions=["A","B","C","D"]
        reaction_scores={ "A": 1, "B": 2 , "E": 0.1}
        
        gap_filled_reaction_scores=modules.gap_fill(key_reactions, reaction_scores)
        
        self.assertDictEqual(gap_filled_reaction_scores, reaction_scores)
        
    def test_gap_fill_optional_reactions_greater_than_threshold(self):
        """
        Test the gap fill function, with a set of scores that include optional reactions
        where just considering the required reactions it is greater than the threshold
        for gap filling
        Test with a minimum score lower for all reactions that the required reactions
        """
        
        key_reactions=["A","B","C","D","E"]
        reaction_scores={ "A": 1, "B": 2 , "C": 2, "D": 1, "F": 0.1}
        
        expected_result={ "A": 1, "B": 2 , "C": 2, "D": 1, "E": 1, "F": 0.1}
        
        gap_filled_reaction_scores=modules.gap_fill(key_reactions, reaction_scores)
        
        self.assertDictEqual(gap_filled_reaction_scores, expected_result)
        
        
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
        
