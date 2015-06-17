import unittest
import logging

import cfg
import utils

from humann2.quantify import modules
from humann2 import store
from humann2 import config
from humann2.quantify import chi2cdf

class TestHumann2QuantifyModulesFunctions(unittest.TestCase):
    """
    Test the functions found in humann2.quantify.modules
    """
    
    def setUp(self):
        # set up nullhandler for logger
        logging.getLogger('humann2.quantify.modules').addHandler(logging.NullHandler())
        
    def test_compute_structured_pathway_abundance_or_coverage_test_abundance(self):
        """
        Test the compute_structured_pathway_abundance_or_coverage function for a simple structure with abundance
        Test the PathwaysDatabase add and get pathway structure along with key reactions
        """
        
        # Create the database structure
        pathways_database_store=store.PathwaysDatabase()
        structure_string=" A B C "
        pathways_database_store.add_pathway_structure("pathway1",structure_string)
        
        reaction_scores={ "A": 1, "B": 2, "C": 3}
        structure=pathways_database_store.get_structure_for_pathway("pathway1")
        key_reactions=pathways_database_store.get_key_reactions_for_pathway("pathway1")
        
        # Compute the abundance
        abundance=modules.compute_structured_pathway_abundance_or_coverage(structure, key_reactions, reaction_scores, 
            False, 0)
        
        # Compute the expected abundance which is the harmonic mean of the values
        expected_abundance=len(reaction_scores.values())/sum(1.0/v for v in reaction_scores.values())
        
        self.assertEqual(abundance, expected_abundance)
        
    def test_compute_structured_pathway_abundance_or_coverage_test_abundance_with_OR(self):
        """
        Test the compute_structured_pathway_abundance_or_coverage function for abundance
        Test the PathwaysDatabase add and get pathway structure along with key reactions
        Test with an OR structure
        """
        
        # Create the database structure
        pathways_database_store=store.PathwaysDatabase()
        structure_string=" A B C ( E , F )"
        pathways_database_store.add_pathway_structure("pathway1",structure_string)
        
        reaction_scores={ "A": 1, "B": 2, "C": 3, "E": 4, "F": 5}
        structure=pathways_database_store.get_structure_for_pathway("pathway1")
        key_reactions=pathways_database_store.get_key_reactions_for_pathway("pathway1")
        
        # Compute the abundance
        abundance=modules.compute_structured_pathway_abundance_or_coverage(structure, key_reactions, reaction_scores, 
            False, 0)
        
        # Compute the expected abundance which is the harmonic mean of the values with the max for the OR
        or_abundance=max([reaction_scores["E"]]+[reaction_scores["F"]])
        del reaction_scores["E"]
        del reaction_scores["F"]
        reaction_scores["E_or_F"]=or_abundance
        expected_abundance=len(reaction_scores.values())/sum(1.0/v for v in reaction_scores.values())
        
        self.assertEqual(abundance, expected_abundance)
        
    def test_compute_structured_pathway_abundance_or_coverage_test_abundance_with_OR_embedded(self):
        """
        Test the compute_structured_pathway_abundance_or_coverage function for abundance
        Test the PathwaysDatabase add and get pathway structure along with key reactions
        Test with an OR structure embedded
        """
        
        # Create the database structure
        pathways_database_store=store.PathwaysDatabase()
        structure_string="( ( A B ) , ( C D ) ) E"
        pathways_database_store.add_pathway_structure("pathway1",structure_string)
        
        reaction_scores={ "A": 1, "B": 2, "C": 3, "D": 4, "E": 5}
        structure=pathways_database_store.get_structure_for_pathway("pathway1")
        key_reactions=pathways_database_store.get_key_reactions_for_pathway("pathway1")
        
        # Compute the abundance
        abundance=modules.compute_structured_pathway_abundance_or_coverage(structure, key_reactions, reaction_scores, 
            False, 0)
        
        # Compute the expected abundance which is the harmonic mean of the values with the max for the ORs
        set_1=[reaction_scores["A"]]+[reaction_scores["B"]]
        del reaction_scores["A"]
        del reaction_scores["B"]
        reaction_scores["AB"]=len(set_1)/sum(1.0/v for v in set_1)
        
        set_2=[reaction_scores["C"]]+[reaction_scores["D"]]
        del reaction_scores["C"]
        del reaction_scores["D"]
        reaction_scores["CD"]=len(set_2)/sum(1.0/v for v in set_2)
        
        reaction_scores["ABCD"]=max([reaction_scores["AB"]]+[reaction_scores["CD"]])
        del reaction_scores["AB"]
        del reaction_scores["CD"]
        
        expected_abundance=len(reaction_scores.values())/sum(1.0/v for v in reaction_scores.values())
        
        self.assertEqual(abundance, expected_abundance)
        
    def test_compute_structured_pathway_abundance_or_coverage_test_abundance_with_AND(self):
        """
        Test the compute_structured_pathway_abundance_or_coverage function for abundance
        Test the PathwaysDatabase add and get pathway structure along with key reactions
        Test with an AND structure
        """
        
        # Create the database structure
        pathways_database_store=store.PathwaysDatabase()
        structure_string=" A B C ( E + F )"
        pathways_database_store.add_pathway_structure("pathway1",structure_string)
        
        reaction_scores={ "A": 1, "B": 2, "C": 3, "E": 4, "F": 5}
        structure=pathways_database_store.get_structure_for_pathway("pathway1")
        key_reactions=pathways_database_store.get_key_reactions_for_pathway("pathway1")
        
        # Compute the abundance
        abundance=modules.compute_structured_pathway_abundance_or_coverage(structure, key_reactions, reaction_scores, 
            False, 0)
        
        # Compute the expected abundance which is the harmonic mean of the values with the harmonic mean for the AND
        or_abundance=2/sum(1.0/v for v in [reaction_scores["E"]]+[reaction_scores["F"]])
        del reaction_scores["E"]
        del reaction_scores["F"]
        reaction_scores["E_or_F"]=or_abundance
        expected_abundance=len(reaction_scores.values())/sum(1.0/v for v in reaction_scores.values())
        
        self.assertEqual(abundance, expected_abundance)
        
    def test_compute_structured_pathway_abundance_or_coverage_test_abundance_missing_required_reaction(self):
        """
        Test the compute_structured_pathway_abundance_or_coverage function for a simple structure with abundance
        Test the PathwaysDatabase add and get pathway structure along with key reactions
        Test with a required reaction missing
        """
        
        # Create the database structure
        pathways_database_store=store.PathwaysDatabase()
        structure_string=" A B C "
        pathways_database_store.add_pathway_structure("pathway1",structure_string)
        
        reaction_scores={ "A": 1, "B": 0, "C": 3}
        structure=pathways_database_store.get_structure_for_pathway("pathway1")
        key_reactions=pathways_database_store.get_key_reactions_for_pathway("pathway1")
        
        # Compute the abundance
        abundance=modules.compute_structured_pathway_abundance_or_coverage(structure, key_reactions, reaction_scores, 
            False, 0)
        
        # Compute the expected abundance which is the harmonic mean of the values that is 0 in the case of a missing reaction
        expected_abundance=0
        
        self.assertEqual(abundance, expected_abundance)

    def test_compute_structured_pathway_abundance_or_coverage_test_abundance_missing_optional_reaction(self):
        """
        Test the compute_structured_pathway_abundance_or_coverage function for a simple structure with abundance
        Test the PathwaysDatabase add and get pathway structure along with key reactions
        Test with a optional reaction missing
        """
        
        # Create the database structure
        pathways_database_store=store.PathwaysDatabase()
        structure_string=" A -B C "
        pathways_database_store.add_pathway_structure("pathway1",structure_string)
        
        reaction_scores={ "A": 1, "B": 0, "C": 3}
        structure=pathways_database_store.get_structure_for_pathway("pathway1")
        key_reactions=pathways_database_store.get_key_reactions_for_pathway("pathway1")
        
        # Compute the abundance
        abundance=modules.compute_structured_pathway_abundance_or_coverage(structure, key_reactions, reaction_scores, 
            False, 0)
        
        # Compute the expected abundance which is the harmonic mean of the values 
        # from the required reactions
        del reaction_scores["B"]
        expected_abundance=len(reaction_scores.values())/sum(1.0/v for v in reaction_scores.values())
        
        self.assertEqual(abundance, expected_abundance)
        
    def test_compute_structured_pathway_abundance_or_coverage_test_coverage(self):
        """
        Test the compute_structured_pathway_abundance_or_coverage function for a simple structure with coverage
        Test the PathwaysDatabase add and get pathway structure along with key reactions
        """
        
        # Create the database structure
        pathways_database_store=store.PathwaysDatabase()
        structure_string=" A B C "
        pathways_database_store.add_pathway_structure("pathway1",structure_string)
        
        reaction_scores={ "A": 1, "B": 2, "C": 3}
        structure=pathways_database_store.get_structure_for_pathway("pathway1")
        key_reactions=pathways_database_store.get_key_reactions_for_pathway("pathway1")
        median=2
        
        # Compute the coverage
        coverage=modules.compute_structured_pathway_abundance_or_coverage(structure, key_reactions, reaction_scores, 
            True, median)
        
        # Compute the expected coverage which is the harmonic mean of the chi2cdf values
        expected_coverage=modules.harmonic_mean([chi2cdf.chi2cdf(v,median) for v in reaction_scores.values()])
        
        self.assertEqual(coverage, expected_coverage)
        
    def test_compute_structured_pathway_abundance_or_coverage_test_coverage_missing_required_reaction(self):
        """
        Test the compute_structured_pathway_abundance_or_coverage function for a simple structure with coverage
        Test the PathwaysDatabase add and get pathway structure along with key reactions
        Test with a required reaction missing
        """
        
        # Create the database structure
        pathways_database_store=store.PathwaysDatabase()
        structure_string=" A B C "
        pathways_database_store.add_pathway_structure("pathway1",structure_string)
        
        reaction_scores={ "A": 1, "B": 0, "C": 3}
        structure=pathways_database_store.get_structure_for_pathway("pathway1")
        key_reactions=pathways_database_store.get_key_reactions_for_pathway("pathway1")
        median=1
        
        # Compute the coverage
        coverage=modules.compute_structured_pathway_abundance_or_coverage(structure, key_reactions, reaction_scores, 
            True, median)
        
        # Compute the expected coverage which is the harmonic mean of the chi2cdf
        # This is zero since one required reaction is missing
        expected_coverage=0
        
        self.assertEqual(coverage, expected_coverage)
        
    def test_compute_structured_pathway_abundance_or_coverage_test_coverage_missing_optional_reaction(self):
        """
        Test the compute_structured_pathway_abundance_or_coverage function for a simple structure with coverage
        Test the PathwaysDatabase add and get pathway structure along with key reactions
        Test with an optional reaction missing
        """
        
        # Create the database structure
        pathways_database_store=store.PathwaysDatabase()
        structure_string=" A -B C "
        pathways_database_store.add_pathway_structure("pathway1",structure_string)
        
        reaction_scores={ "A": 1, "B": 0, "C": 3}
        structure=pathways_database_store.get_structure_for_pathway("pathway1")
        key_reactions=pathways_database_store.get_key_reactions_for_pathway("pathway1")
        median=1
        
        # Compute the coverage
        coverage=modules.compute_structured_pathway_abundance_or_coverage(structure, key_reactions, reaction_scores, 
            True, median)
        
        # Compute the expected coverage which is the harmonic mean of the chi2cdf for the required reactions
        del reaction_scores["B"]
        expected_coverage=modules.harmonic_mean([chi2cdf.chi2cdf(v,median) for v in reaction_scores.values()])
        
        self.assertEqual(coverage, expected_coverage)
        
    def test_compute_pathways_abundance_unstructured(self):
        """
        Test the compute_pathways_abundance function
        Test PathwaysDatabase add
        Test PathwaysAndReactions store
        Test Pathways store
        Test with unstructured pathways
        """
      
        # Create the database structure
        pathways_database_store=store.PathwaysDatabase()
        pathways_database_store.add_pathway("pathway1",["A","B","C","D"])
        pathways_database_store.add_pathway("pathway2",["A","B","C","D","E","F"])
        
        # Have all test data be from the same bug
        bug="bug"
        pathways_and_reactions_store=store.PathwaysAndReactions()
        pathways_and_reactions_store.add(bug, "A", "pathway1", 1)
        # Note B is not recored for pathway1 which will result in a zero value
        pathways_and_reactions_store.add(bug, "C", "pathway1", 3)
        pathways_and_reactions_store.add(bug, "D", "pathway1", 4)
        pathways_and_reactions_store.add(bug, "A", "pathway2", 10)
        pathways_and_reactions_store.add(bug, "B", "pathway2", 20)
        # Note C is not recored for pathway2 which will result in a zero value
        pathways_and_reactions_store.add(bug, "D", "pathway2", 40)
        pathways_and_reactions_store.add(bug, "E", "pathway2", 50)
        # Note F is not recored for pathway2 which will result in a zero value
        
        # The abundance for each pathway is the average of the largest half of the reaction values
        # For unstructured pathways, if the reaction is not included it does not result in a zero abundance
        pathway1_values=[0,1,3,4]
        pathway1_abundance_set=pathway1_values[(len(pathway1_values)/2):]
        pathway1_abundance=sum(pathway1_abundance_set)/len(pathway1_abundance_set)
        
        pathway2_values=[0,0,10,20,40,50]
        pathway2_abundance_set=pathway2_values[(len(pathway2_values)/2):]
        pathway2_abundance=sum(pathway2_abundance_set)/len(pathway2_abundance_set)      
        
        pathways_abundance_store_result=modules.compute_pathways_abundance(pathways_and_reactions_store, pathways_database_store)
        
        # Test the pathways abundance match those expected
        self.assertEqual(pathways_abundance_store_result.get_score_for_bug(bug,"pathway1"), pathway1_abundance)
        self.assertEqual(pathways_abundance_store_result.get_score_for_bug(bug,"pathway2"), pathway2_abundance)
        
    def test_compute_pathways_abundance_structured(self):
        """
        Test the compute_pathways_abundance function
        Test PathwaysDatabase add
        Test PathwaysAndReactions store
        Test Pathways store
        Test with structured pathways
        """
      
        # Create the database structure
        pathways_database_store=store.PathwaysDatabase()
        pathways_database_store.add_pathway_structure("pathway1"," A B C D ")
        pathways_database_store.add_pathway_structure("pathway2"," A B C D E F ")
        
        # Have all test data be from the same bug
        bug="bug"
        pathways_and_reactions_store=store.PathwaysAndReactions()
        # Just a note that a value of 1 has a chi2cdf value of 0
        # Also values ~10 or less have small chi2cdf values
        pathways_and_reactions_store.add(bug, "A", "pathway1", 11)
        pathways_and_reactions_store.add(bug, "B", "pathway1", 12)
        pathways_and_reactions_store.add(bug, "C", "pathway1", 13)
        pathways_and_reactions_store.add(bug, "D", "pathway1", 14)
        pathways_and_reactions_store.add(bug, "A", "pathway2", 19)
        pathways_and_reactions_store.add(bug, "B", "pathway2", 20)
        pathways_and_reactions_store.add(bug, "C", "pathway2", 30)
        pathways_and_reactions_store.add(bug, "D", "pathway2", 40)
        pathways_and_reactions_store.add(bug, "E", "pathway2", 50)
        pathways_and_reactions_store.add(bug, "F", "pathway2", 60)
        
        # The abundance for each pathway is the harmonic mean of the values
        pathway1_values=[11,12,13,14]
        pathway1_abundance=len(pathway1_values)/sum(1.0/v for v in pathway1_values)
        
        pathway2_values=[19,20,30,40,50,60]
        pathway2_abundance=len(pathway2_values)/sum(1.0/v for v in pathway2_values)    
        
        # Find the actual result
        pathways_abundance_store_result=modules.compute_pathways_abundance(pathways_and_reactions_store, pathways_database_store)
        
        # Test the pathways abundance match those expected
        self.assertEqual(pathways_abundance_store_result.get_score_for_bug(bug,"pathway1"), pathway1_abundance)
        self.assertEqual(pathways_abundance_store_result.get_score_for_bug(bug,"pathway2"), pathway2_abundance)
        
    def test_compute_pathways_coverage_unstructured(self):
        """
        Test the compute_pathways_coaverage function
        Test PathwaysDatabase add
        Test PathwaysAndReactions store
        Test Pathways store
        Test with unstructured pathways
        """
        
        # Set xipe to off
        config.xipe_toggle = "off"
      
        # Create the database structure
        pathways_database_store=store.PathwaysDatabase()
        pathways_database_store.add_pathway("pathway1",["A","B","C","D"])
        pathways_database_store.add_pathway("pathway2",["A","B","C","D","E","F"])
        
        # Have all test data be from the same bug
        bug="bug"
        pathways_and_reactions_store=store.PathwaysAndReactions()
        pathways_and_reactions_store.add(bug, "A", "pathway1", 1)
        # Note B is not recored for pathway1 which will result in a zero value
        pathways_and_reactions_store.add(bug, "C", "pathway1", 3)
        pathways_and_reactions_store.add(bug, "D", "pathway1", 40)
        pathways_and_reactions_store.add(bug, "A", "pathway2", 10)
        pathways_and_reactions_store.add(bug, "B", "pathway2", 20)
        # Note C is not recored for pathway2 which will result in a zero value
        pathways_and_reactions_store.add(bug, "D", "pathway2", 40)
        pathways_and_reactions_store.add(bug, "E", "pathway2", 50)
        # Note F is not recored for pathway2 which will result in a zero value
        
        # The coverage for each pathway if the number of reactions greater than the median
        # divided by the total reactions in the pathway
        # The median for this set is 20
        pathway1_values=[0,1,3,40]
        count_greater_than_median=1
        pathway1_coverage=count_greater_than_median/float(len(pathway1_values))

        
        pathway2_values=[0,0,10,20,40,50]
        count_greater_than_median=2
        pathway2_coverage=count_greater_than_median/float(len(pathway2_values))    
        
        pathways_coverage_store_result=modules.compute_pathways_coverage(pathways_and_reactions_store, pathways_database_store)
        
        # Test the pathways abundance match those expected
        self.assertEqual(pathways_coverage_store_result.get_score_for_bug(bug,"pathway1"), pathway1_coverage)
        self.assertEqual(pathways_coverage_store_result.get_score_for_bug(bug,"pathway2"), pathway2_coverage)
        
    def test_compute_pathways_coverage_structured(self):
        """
        Test the compute_pathways_coverage function
        Test PathwaysDatabase add
        Test PathwaysAndReactions store
        Test Pathways store
        Test with structured pathways
        """
      
        # Set xipe to off
        config.xipe_toggle = "off"
      
        # Create the database structure
        pathways_database_store=store.PathwaysDatabase()
        pathways_database_store.add_pathway_structure("pathway1"," A B C D ")
        pathways_database_store.add_pathway_structure("pathway2"," A B C D E F ")
        
        # Have all test data be from the same bug
        bug="bug"
        pathways_and_reactions_store=store.PathwaysAndReactions()
        # Just a note that a value of 1 has a chi2cdf value of 0
        # Also values ~10 or less have small chi2cdf values
        pathways_and_reactions_store.add(bug, "A", "pathway1", 11)
        pathways_and_reactions_store.add(bug, "B", "pathway1", 12)
        pathways_and_reactions_store.add(bug, "C", "pathway1", 13)
        pathways_and_reactions_store.add(bug, "D", "pathway1", 14)
        pathways_and_reactions_store.add(bug, "A", "pathway2", 19)
        pathways_and_reactions_store.add(bug, "B", "pathway2", 20)
        pathways_and_reactions_store.add(bug, "C", "pathway2", 30)
        pathways_and_reactions_store.add(bug, "D", "pathway2", 40)
        pathways_and_reactions_store.add(bug, "E", "pathway2", 50)
        pathways_and_reactions_store.add(bug, "F", "pathway2", 60)    
        
        # Get the coverage result
        # The median is the median of all of the reactions of all of the pathways for this bug
        median_score_value=19.5
        pathway1_values=[11,12,13,14]
        coverage_pathway1=len(pathway1_values)/sum(1.0/chi2cdf.chi2cdf(v,median_score_value) for v in pathway1_values)

        pathway2_values=[19,20,30,40,50,60]
        coverage_pathway2=len(pathway2_values)/sum(1.0/chi2cdf.chi2cdf(v,median_score_value) for v in pathway2_values)
        
        # Find the actual result
        pathways_abundance_store_result=modules.compute_pathways_coverage(pathways_and_reactions_store, pathways_database_store)
        
        # Test the pathways abundance match those expected
        self.assertEqual(pathways_abundance_store_result.get_score_for_bug(bug,"pathway1"), coverage_pathway1)
        self.assertEqual(pathways_abundance_store_result.get_score_for_bug(bug,"pathway2"), coverage_pathway2)
        
