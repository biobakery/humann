import unittest
import logging

import cfg
import utils

from humann.quantify import modules
from humann import config

class TestHumannQuantifyModulesFunctions(unittest.TestCase):
    """
    Test the functions found in humann.quantify.modules
    """
    
    def setUp(self):
        # set gap fill on
        config.gap_fill_toggle="on"
        
        # set up nullhandler for logger
        logging.getLogger('humann.quantify.modules').addHandler(logging.NullHandler())
        
    def test_gap_fill_zero_gaps(self):
        """
        Test the gap fill function, with a set of scores that do not have gaps
        Test for boost of lowest score
        """
        
        key_reactions=["A","B"]
        reaction_scores={ "A": 1, "B": 2 }
        
        gap_filled_reaction_scores=modules.gap_fill(key_reactions, reaction_scores)
        expected_gap_filled_reaction_scores={ "A": 2, "B": 2 }
        
        self.assertDictEqual(gap_filled_reaction_scores, expected_gap_filled_reaction_scores)
        
    def test_gap_fill_greater_than_threshold(self):
        """
        Test the gap fill function, with a set of scores where the gaps are
        greater than threshold to apply gap filling
        """
        
        key_reactions=["A","B","C","D","E","G"]
        reaction_scores={ "A": 1, "B": 2 , "C": 2, "D": 1}
        
        gap_filled_reaction_scores=modules.gap_fill(key_reactions, reaction_scores)
        
        self.assertDictEqual(gap_filled_reaction_scores, reaction_scores)
        
    def test_gap_fill_equal_threshold(self):
        """
        Test the gap fill function, with a set of scores where the gaps equal the
        threshold to apply gap filling
        """
        
        key_reactions=["A","B","C","D","E"]
        reaction_scores={ "A": 1, "B": 2 , "C": 2, "D": 1}
        
        expected_result={ "A": 1, "B": 2 , "C": 2, "D": 1, "E": 1}
        
        gap_filled_reaction_scores=modules.gap_fill(key_reactions, reaction_scores)
        
        self.assertDictEqual(gap_filled_reaction_scores, expected_result)
        
    def test_gap_fill_optional_reactions_zero_gaps(self):
        """
        Test the gap fill function, with a set of scores that include optional reactions
        where just considering the required reactions it does not require gap filling
        Test boost lowest abundance score of key reactions
        """
        
        key_reactions=["A","B"]
        reaction_scores={ "A": 1, "B": 2 , "E": 0.1}
        
        gap_filled_reaction_scores=modules.gap_fill(key_reactions, reaction_scores)
        expected_gap_filled_reaction_scores={ "A": 2, "B": 2 , "E": 0.1}
        
        self.assertDictEqual(gap_filled_reaction_scores, expected_gap_filled_reaction_scores)
        
    def test_gap_fill_optional_reactions_equal_threshold(self):
        """
        Test the gap fill function, with a set of scores that include optional reactions
        where just considering the required reactions the gaps equal the threshold
        for gap filling
        Test with a minimum score lower for all reactions that the required reactions
        """
        
        key_reactions=["A","B","C","D","E"]
        reaction_scores={ "A": 1, "B": 2 , "C": 2, "D": 1, "F": 0.1}
        
        expected_result={ "A": 1, "B": 2 , "C": 2, "D": 1, "E": 1, "F": 0.1}
        
        gap_filled_reaction_scores=modules.gap_fill(key_reactions, reaction_scores)
        
        self.assertDictEqual(gap_filled_reaction_scores, expected_result)
        
    def test_gap_fill_optional_reactions_greater_than_threshold(self):
        """
        Test the gap fill function, with a set of scores that include optional reactions
        where just considering the required reactions the gaps are greater than the threshold
        for gap filling
        Test with a minimum score lower for all reactions that the required reactions
        """
        
        key_reactions=["A","B","C","D","E","G"]
        reaction_scores={ "A": 1, "B": 2 , "C": 2, "D": 1, "F": 0.1}
        
        gap_filled_reaction_scores=modules.gap_fill(key_reactions, reaction_scores)
        
        self.assertDictEqual(gap_filled_reaction_scores, reaction_scores)
        
    def test_gap_fill_all_required_reactions(self):
        """
        Test the gap fill function, with a set of scores of all required reactions
        Test the lowest score is boosted
        """
        
        key_reactions=["A","B","C","D","E"]
        reaction_scores={ "A": 1, "B": 2 , "C": 2, "D": 1, "E": 0.1}
        
        gap_filled_reaction_scores=modules.gap_fill(key_reactions, reaction_scores)
        
        expected_reaction_scores={ "A": 1, "B": 2 , "C": 2, "D": 1, "E": 1}
        
        self.assertDictEqual(gap_filled_reaction_scores, expected_reaction_scores)
        
    def test_gap_fill_all_required_reactions_one_optional(self):
        """
        Test the gap fill function, with a set of scores of all required reactions
        Test the lowest score is boosted
        Test the optional reaction is not boosted
        """
        
        key_reactions=["A","B","C","D","E"]
        reaction_scores={ "A": 1, "B": 2 , "C": 2, "D": 1, "E": 0.1, "F": 0.1}
        
        gap_filled_reaction_scores=modules.gap_fill(key_reactions, reaction_scores)
        
        expected_reaction_scores={ "A": 1, "B": 2 , "C": 2, "D": 1, "E": 1, "F": 0.1}
        
        self.assertDictEqual(gap_filled_reaction_scores, expected_reaction_scores)
        
    def test_gap_fill_all_required_two_lowest_scores(self):
        """
        Test the gap fill function, with a set of scores of all required reactions
        Test the lowest score is boosted
        Test the two lowest scores are unchanged
        """
        
        key_reactions=["A","B","C","D","E"]
        reaction_scores={ "A": 1, "B": 2 , "C": 2, "D": 0.1, "E": 0.1}
        
        gap_filled_reaction_scores=modules.gap_fill(key_reactions, reaction_scores)
        
        self.assertDictEqual(gap_filled_reaction_scores, reaction_scores)
        
        
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
        
