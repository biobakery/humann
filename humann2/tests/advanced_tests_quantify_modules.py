import unittest
import logging
import tempfile
import filecmp

import cfg
import utils
import os

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
        
        config.unnamed_temp_dir=tempfile.gettempdir()
        
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
        
        self.assertAlmostEqual(abundance, expected_abundance)
        
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
        pathway1_abundance_set=pathway1_values[int(len(pathway1_values)/2):]
        pathway1_abundance=sum(pathway1_abundance_set)/len(pathway1_abundance_set)
        
        pathway2_values=[0,0,10,20,40,50]
        pathway2_abundance_set=pathway2_values[int(len(pathway2_values)/2):]
        pathway2_abundance=sum(pathway2_abundance_set)/len(pathway2_abundance_set)      
        
        pathways_abundance_store_result, reactions_in_pathways_present=modules.compute_pathways_abundance(pathways_and_reactions_store, pathways_database_store)
        
        # Test the pathways abundance match those expected
        self.assertEqual(pathways_abundance_store_result.get_score_for_bug(bug,"pathway1"), pathway1_abundance)
        self.assertEqual(pathways_abundance_store_result.get_score_for_bug(bug,"pathway2"), pathway2_abundance)
        
    def test_compute_pathways_abundance_unstructured_reactions_list(self):
        """
        Test the compute_pathways_abundance function
        Test PathwaysDatabase add
        Test PathwaysAndReactions store
        Test Pathways store
        Test with unstructured pathways
        Test the resulting list of reactions included in pathways
        """
      
        # Create the database structure
        pathways_database_store=store.PathwaysDatabase()
        pathways_database_store.add_pathway("pathway1",["A","B","C","D"])
        pathways_database_store.add_pathway("pathway2",["A","B","C","D","E","F"])
        pathways_database_store.add_pathway("pathway3",["A","B","G"])
        
        # Have all test data be from two bugs
        bug="bug1"
        pathways_and_reactions_store=store.PathwaysAndReactions()
        pathways_and_reactions_store.add(bug, "A", "pathway1", 1)
        pathways_and_reactions_store.add(bug, "B", "pathway1", 2)
        pathways_and_reactions_store.add(bug, "B", "pathway2", 2)
        pathways_and_reactions_store.add(bug, "C", "pathway2", 3)

        expected_reactions_in_pathways_present={}
        expected_reactions_in_pathways_present[bug]=["A","B","C"]
       
        bug="bug2"
        pathways_and_reactions_store.add(bug, "A", "pathway1", 1)
        pathways_and_reactions_store.add(bug, "D", "pathway1", 3)
        pathways_and_reactions_store.add(bug, "B", "pathway1", 2)
        pathways_and_reactions_store.add(bug, "B", "pathway2", 2)
       
        expected_reactions_in_pathways_present[bug]=["A","B","D"]
        
        pathways_abundance_store_result, reactions_in_pathways_present=modules.compute_pathways_abundance(
            pathways_and_reactions_store, pathways_database_store)
        
        # Test the reactions match those expected
        self.assertEqual(sorted(expected_reactions_in_pathways_present["bug1"]),
                         sorted(list(reactions_in_pathways_present["bug1"])))
        self.assertEqual(sorted(expected_reactions_in_pathways_present["bug2"]),
                         sorted(list(reactions_in_pathways_present["bug2"])))
        
    def test_compute_pathways_abundance_structured_reactions_list(self):
        """
        Test the compute_pathways_abundance function
        Test PathwaysDatabase add
        Test PathwaysAndReactions store
        Test Pathways store
        Test with structured pathways
        Test gap fill
        Test the resulting list of reactions included in pathways
        """
      
        # Create the database structure
        pathways_database_store=store.PathwaysDatabase()
        pathways_database_store.add_pathway_structure("pathway1"," A B C D ")
        pathways_database_store.add_pathway_structure("pathway2"," A B C D E F ")
        pathways_database_store.add_pathway_structure("pathway3"," A B G")
        
        # Have all test data be from three bugs
        bug="bug1"
        pathways_and_reactions_store=store.PathwaysAndReactions()
        pathways_and_reactions_store.add(bug, "A", "pathway1", 1)
        pathways_and_reactions_store.add(bug, "B", "pathway1", 2)
        pathways_and_reactions_store.add(bug, "C", "pathway1", 3)

        expected_reactions_in_pathways_present={}
        # This pathway is present because D is filled in
        # Though D does not have abundance so it is not included in the list
        expected_reactions_in_pathways_present[bug]=["A","B","C"]
       
        bug="bug2"
        pathways_and_reactions_store.add(bug, "A", "pathway1", 1)
        pathways_and_reactions_store.add(bug, "B", "pathway1", 3)
        pathways_and_reactions_store.add(bug, "D", "pathway2", 2)
        pathways_and_reactions_store.add(bug, "B", "pathway2", 2)
       
        # The pathways for this bug are missing too many reactions to have abundance
        expected_reactions_in_pathways_present[bug]=[]
        
        bug="bug3"
        pathways_and_reactions_store.add(bug, "A", "pathway3", 1)
        pathways_and_reactions_store.add(bug, "B", "pathway3", 3)
        pathways_and_reactions_store.add(bug, "G", "pathway3", 2)
        pathways_and_reactions_store.add(bug, "B", "pathway2", 2)
       
        # One pathway for this bug includes all reactions
        expected_reactions_in_pathways_present[bug]=["A","B","G"]
        
        pathways_abundance_store_result, reactions_in_pathways_present=modules.compute_pathways_abundance(
            pathways_and_reactions_store, pathways_database_store)
        
        # Test the reactions match those expected
        self.assertEqual(sorted(expected_reactions_in_pathways_present["bug1"]),
                         sorted(list(reactions_in_pathways_present["bug1"])))
        self.assertEqual(sorted(expected_reactions_in_pathways_present["bug2"]),
                         sorted(list(reactions_in_pathways_present["bug2"])))
        self.assertEqual(sorted(expected_reactions_in_pathways_present["bug3"]),
                         sorted(list(reactions_in_pathways_present["bug3"])))
                
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
        # boost the lowest value in the pathway
        pathway1_values_boosted=[12,12,13,14]
        pathway1_abundance=len(pathway1_values_boosted)/sum(1.0/v for v in pathway1_values_boosted)
        
        pathway2_values_boosted=[20,20,30,40,50,60]
        pathway2_abundance=len(pathway2_values_boosted)/sum(1.0/v for v in pathway2_values_boosted)    
        
        # Find the actual result
        pathways_abundance_store_result, reactions_in_pathways_present=modules.compute_pathways_abundance(pathways_and_reactions_store, pathways_database_store)
        
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
        # boost the pathway values
        pathway1_values_boosted=[12,12,13,14]
        coverage_pathway1=len(pathway1_values_boosted)/sum(1.0/chi2cdf.chi2cdf(v,median_score_value) for v in pathway1_values_boosted)

        pathway2_values_boosted=[20,20,30,40,50,60]
        coverage_pathway2=len(pathway2_values_boosted)/sum(1.0/chi2cdf.chi2cdf(v,median_score_value) for v in pathway2_values_boosted)
        
        # Find the actual result
        pathways_abundance_store_result=modules.compute_pathways_coverage(pathways_and_reactions_store, pathways_database_store)
        
        # Test the pathways abundance match those expected
        self.assertEqual(pathways_abundance_store_result.get_score_for_bug(bug,"pathway1"), coverage_pathway1)
        self.assertEqual(pathways_abundance_store_result.get_score_for_bug(bug,"pathway2"), coverage_pathway2)
        
    def test_pathways_coverage_with_names(self):
        """
        Test the pathways coverage computation (xipe and minpath are off)
        Test the pathways print function
        Test the pathways mapping to names
        Test the unmapped and unintegrated values are printed
        """
        
        # update the max decimals to allow for rounding
        config.output_max_decimals=7
        
        # Load in the pathways databases
        reactions_database=store.ReactionsDatabase(config.pathways_database_part1)
        pathways_database=store.PathwaysDatabase(config.pathways_database_part2, reactions_database)
        
        # Load in the gene scores from the file
        # This file has the gene names included
        gene_scores=store.GeneScores()
        gene_scores.add_from_file(cfg.larger_gene_families_uniref50_with_names_file)
        
        # Turn off xipe and minpath
        minpath_toggle_original=config.minpath_toggle
        config.minpath_toggle="off"
        xipe_toggle_original=config.xipe_toggle
        config.xipe_toggle="off"
        
        pathways_and_reactions_store=modules.identify_reactions_and_pathways(
        gene_scores, reactions_database, pathways_database)
        
        # set the locations to write as temp files
        file_out, abundance_file=tempfile.mkstemp()
        os.close(file_out)
        config.pathabundance_file=abundance_file
        
        file_out, coverage_file=tempfile.mkstemp()
        os.close(file_out)
        config.pathcoverage_file=coverage_file
        
        unaligned_reads_count=10
        abundance_file, coverage_file=modules.compute_pathways_abundance_and_coverage(
        gene_scores, reactions_database, pathways_and_reactions_store, pathways_database,
        unaligned_reads_count)
        
        # Reset xipe and minpath
        config.minpath_toggle=minpath_toggle_original
        config.xipe_toggle=xipe_toggle_original
        
        # check the output is as expected
        self.assertTrue(filecmp.cmp(coverage_file,
            cfg.demo_pathcoverage_file, shallow=False))
        
        utils.remove_temp_file(abundance_file)
        utils.remove_temp_file(coverage_file)
        
    def test_pathways_abundance_with_names(self):
        """
        Test the pathways abundance computation (xipe and minpath are off)
        Test the pathways print function
        Test the pathways mapping to names
        Test the unmapped and unintegrated values are printed
        """
        
        # update the max decimals to allow for rounding
        config.output_max_decimals=7
        
        # Load in the pathways databases
        reactions_database=store.ReactionsDatabase(config.pathways_database_part1)
        pathways_database=store.PathwaysDatabase(config.pathways_database_part2, reactions_database)
        
        # Load in the gene scores from the file
        # This file has the gene names included
        gene_scores=store.GeneScores()
        gene_scores.add_from_file(cfg.larger_gene_families_uniref50_with_names_file)
        
        # Turn off xipe and minpath
        minpath_toggle_original=config.minpath_toggle
        config.minpath_toggle="off"
        xipe_toggle_original=config.xipe_toggle
        config.xipe_toggle="off"
        
        pathways_and_reactions_store=modules.identify_reactions_and_pathways(
        gene_scores, reactions_database, pathways_database)
        
        # set the locations to write as temp files
        file_out, abundance_file=tempfile.mkstemp()
        os.close(file_out)
        config.pathabundance_file=abundance_file
        
        file_out, coverage_file=tempfile.mkstemp()
        os.close(file_out)
        config.pathcoverage_file=coverage_file
        
        unaligned_reads_count=10
        abundance_file, coverage_file=modules.compute_pathways_abundance_and_coverage(
        gene_scores, reactions_database, pathways_and_reactions_store, pathways_database,
        unaligned_reads_count)
        
        # Reset xipe and minpath
        config.minpath_toggle=minpath_toggle_original
        config.xipe_toggle=xipe_toggle_original
        
        # check the output is as expected
        self.assertTrue(filecmp.cmp(abundance_file,
            cfg.demo_pathabundance_file, shallow=False))
        
        utils.remove_temp_file(abundance_file)
        utils.remove_temp_file(coverage_file)
        
    def test_compute_gene_abundance_in_pathways_without_reactions_database(self):
        """
        Test the compute gene abundance function
        Test the GeneScores add function
        Test without a reactions database (the pathways database is composed of genes)
        """
        
        gene_scores=store.GeneScores()
        # Add gene scores for two bugs
        reactions_in_pathways_present={}
        bug="bug1"
        gene_scores.add_single_score(bug, "gene1", 1)
        gene_scores.add_single_score(bug, "gene2", 2)
        gene_scores.add_single_score(bug, "gene4", 4)
        reactions_in_pathways_present[bug]=["gene1","gene2"]
        
        bug="bug2"
        # Test with different values of gene1 for each bug
        gene_scores.add_single_score(bug, "gene1", 1.1)
        gene_scores.add_single_score(bug, "gene7", 7)
        gene_scores.add_single_score(bug, "gene6", 6)
        reactions_in_pathways_present[bug]=["gene6"]
        
        reactions_database=None
        gene_abundance_in_pathways, remaining_gene_abundance=modules.compute_gene_abundance_in_pathways(
            gene_scores, reactions_database, reactions_in_pathways_present)
        
        # Check the gene abundances in pathways are correct
        self.assertEqual(gene_abundance_in_pathways["bug1"], 3)
        self.assertEqual(gene_abundance_in_pathways["bug2"], 6)
        
        # Check the gene abundances not in pathways are correct
        self.assertEqual(remaining_gene_abundance["bug1"], 4)
        self.assertAlmostEqual(remaining_gene_abundance["bug2"], 8.1)
        
    def test_compute_gene_abundance_in_pathways_with_reactions_database(self):
        """
        Test the compute gene abundance function
        Test the GeneScores add function
        Test the ReactionsDatabase add function
        Test with a reactions database (the pathways database is composed of reactions
        and these reactions map to genes, with some genes mapping to multiple reactions)
        """
        
        gene_scores=store.GeneScores()
        # Add gene scores for two bugs
        reactions_in_pathways_present={}
        bug="bug1"
        gene_scores.add_single_score(bug, "gene1", 1)
        gene_scores.add_single_score(bug, "gene2", 2)
        gene_scores.add_single_score(bug, "gene4", 4)
        
        bug="bug2"
        # Test with different values of gene1 for each bug
        gene_scores.add_single_score(bug, "gene1", 1.1)
        gene_scores.add_single_score(bug, "gene7", 7)
        gene_scores.add_single_score(bug, "gene6", 6)
        gene_scores.add_single_score(bug, "gene8", 0.2)
        
        reactions_database=store.ReactionsDatabase()
        reactions={"reaction1":["gene1","gene6"], "reaction2":["gene1","gene2"],
                   "reaction3":["gene4","gene7"],"reaction4":["gene8"]}
        reactions_database.add_reactions(reactions)
        
        # Test one bug with two reactions and one bug with a single reaction
        # For the bug with two reactions, test with both reactions including
        # The same gene (to test this value is not added twice in the abundance result)
        reactions_in_pathways_present["bug1"]=["reaction1","reaction2"]
        reactions_in_pathways_present["bug2"]=["reaction1"]
        
        gene_abundance_in_pathways, remaining_gene_abundance=modules.compute_gene_abundance_in_pathways(
            gene_scores, reactions_database, reactions_in_pathways_present)
        
        # Check the gene abundances in pathways are correct
        self.assertEqual(gene_abundance_in_pathways["bug1"], 3)
        self.assertAlmostEqual(gene_abundance_in_pathways["bug2"], 7.1)
        
        # Check the gene abundances not in pathways are correct
        self.assertEqual(remaining_gene_abundance["bug1"], 4)
        self.assertAlmostEqual(remaining_gene_abundance["bug2"], 7.2)
        
    def test_compute_unmapped_and_unintegrated(self):
        """
        Test the unmapped and unintegrated function
        Test Pathways add
        """
        
        # Add pathways for two bugs
        pathways_abundance=store.Pathways()
        pathways_abundance.add("bug1","pathway1",1)
        pathways_abundance.add("bug1","pathway2",2)
        pathways_abundance.add("bug1","pathway3",3)
        
        pathways_abundance.add("bug2","pathway1",10)
        pathways_abundance.add("bug2","pathway2",20)
        pathways_abundance.add("bug2","pathway4",4)
        
        pathways_abundance.add("all","pathway1",100)
        pathways_abundance.add("all","pathway2",200)
        pathways_abundance.add("all","pathway3",300)
        pathways_abundance.add("all","pathway4",40)
        
        gene_abundance_in_pathways={}
        gene_abundance_in_pathways["bug1"]=5
        gene_abundance_in_pathways["bug2"]=15
        gene_abundance_in_pathways["all"]=155
        
        remaining_gene_abundance={}
        remaining_gene_abundance["bug1"]=0.1
        remaining_gene_abundance["bug2"]=0.5
        remaining_gene_abundance["all"]=0.15

        unaligned_reads_count=100
        
        unmapped_all, unintegrated_all, unintegrated_per_bug=modules.compute_unmapped_and_unintegrated(
            gene_abundance_in_pathways, remaining_gene_abundance, unaligned_reads_count, pathways_abundance)
        
        # Compression constant = total abundance of all pathways / total abundance of genes in pathways 
        # Unmapped = Compression constant * number of unmapped reads
        expected_unmapped_all=((100+200+300+40)/(155*1.0)) * unaligned_reads_count
        
        # Unintegrated = Compression constant * total abundance of genes NOT in pathways
        expected_unintegrated_all=((100+200+300+40)/(155*1.0)) * 0.15
        expected_unintegrated_bug1=((1+2+3)/(5*1.0)) * 0.1
        expected_unintegrated_bug2=((10+20+4)/(15*1.0)) * 0.5
        
        # Test the computed values are almost equal to the expected
        self.assertAlmostEqual(unmapped_all, expected_unmapped_all)
        self.assertAlmostEqual(unintegrated_all, expected_unintegrated_all)   
        self.assertAlmostEqual(unintegrated_per_bug["bug1"], expected_unintegrated_bug1)
        self.assertAlmostEqual(unintegrated_per_bug["bug2"], expected_unintegrated_bug2)
        
