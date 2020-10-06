import unittest
import logging
import os
import sys
import math
import tempfile
import logging

import cfg

from humann import store
from humann import config

class TestAdvancedHumannStoreFunctions(unittest.TestCase):
    """
    Test the functions found in humann/src/store.py
    """
    
    def setUp(self):
        config.unnamed_temp_dir=tempfile.gettempdir()
        
        # set up nullhandler for logger
        logging.getLogger('humann.store').addHandler(logging.NullHandler())
        
    def test_Alignments_compute_gene_scores_single_gene_single_query(self):
        """
        Test the compute_gene_scores function
        Test one hit for gene with one hit for query
        """
        
        # create a set of hits
        matches1=41.0
        matches2=57.1
        matches3=61.0
        matches4=72.1
        
        gene1_length=2
        gene2_length=3
        gene3_length=4
        
        # Create a set of alignments
        alignments_store=store.Alignments()
        alignments_store.add("gene1",gene1_length,"query1",matches1,"bug1")
        alignments_store.add("gene2",gene2_length,"query1",matches2,"bug1")
        alignments_store.add("gene2",gene2_length,"query2",matches3,"bug1")
        alignments_store.add("gene3",gene3_length,"query3",matches4,"bug1")
        
        gene_scores_store=store.GeneScores()
        
        # compute gene scores
        alignments_store.convert_alignments_to_gene_scores(gene_scores_store)
        
        # convert lengths to per kb
        gene3_length=gene3_length/1000.0
        
        # gene3
        hit4_score=math.pow(matches4, config.match_power)
        query3_sum=hit4_score
        gene_score=hit4_score/query3_sum/gene3_length

        self.assertEqual(gene_scores_store.get_score("bug1","gene3"),gene_score)

    def test_Alignments_compute_gene_scores_single_gene_double_query(self):
        """
        Test the compute_gene_scores function
        Test one hit for gene with more than one hit per query
        """
        
        # create a set of hits
        # bug, reference, reference_length, query, matches = hit
        
        matches1=41.0
        matches2=57.1
        matches3=61.0
        matches4=72.1
        
        gene1_length=2
        gene2_length=3
        gene3_length=4
        
        # Create a set of alignments
        alignments_store=store.Alignments()
        alignments_store.add("gene1",gene1_length,"query1",matches1,"bug1")
        alignments_store.add("gene2",gene2_length,"query1",matches2,"bug1")
        alignments_store.add("gene2",gene2_length,"query2",matches3,"bug1")
        alignments_store.add("gene3",gene3_length,"query3",matches4,"bug1")
        
        gene_scores_store=store.GeneScores()
        
        # compute gene scores
        alignments_store.convert_alignments_to_gene_scores(gene_scores_store)
        
        # convert lengths to per kb
        gene1_length=gene1_length/1000.0
        
        # gene1
        hit1_score=math.pow(matches1, config.match_power)
        hit2_score=math.pow(matches2, config.match_power)
        query1_sum=hit1_score+hit2_score
        gene_score=hit1_score/query1_sum/gene1_length

        self.assertEqual(gene_scores_store.get_score("bug1","gene1"),gene_score)

    def test_Alignments_compute_gene_scores_double_gene_double_query(self):
        """
        Test the compute_gene_scores function
        Test two hits to gene with more than one hit per query
        """
        
        # create a set of hits
        # bug, reference, reference_length, query, matches = hit
        
        matches1=41.0
        matches2=57.1
        matches3=61.0
        matches4=72.1
        
        gene1_length=2
        gene2_length=3
        gene3_length=4
        
        # Create a set of alignments
        alignments_store=store.Alignments()
        alignments_store.add("gene1",gene1_length,"query1",matches1,"bug1")
        alignments_store.add("gene2",gene2_length,"query1",matches2,"bug1")
        alignments_store.add("gene2",gene2_length,"query2",matches3,"bug1")
        alignments_store.add("gene3",gene3_length,"query3",matches4,"bug1")
        
        gene_scores_store=store.GeneScores()
        
        # compute gene scores
        alignments_store.convert_alignments_to_gene_scores(gene_scores_store)
        
        # gene1
        hit1_score=math.pow(matches1, config.match_power)
        hit2_score=math.pow(matches2, config.match_power)
        query1_sum=hit1_score+hit2_score
        
        # convert lengths to per kb
        gene2_length=gene2_length/1000.0
        
        # gene2
        hit3_score=math.pow(matches3, config.match_power)
        query2_sum=hit3_score
        gene_score=hit3_score/query2_sum/gene2_length + hit2_score/query1_sum/gene2_length

        self.assertAlmostEqual(gene_scores_store.get_score("bug1","gene2"),gene_score,places=7)
        
    def test_Alignments_id_mapping_all_gene_list(self):
        """
        Test the store_id_mapping function
        Test the add_annotated and process_reference_annotation with id mapping
        Test the genes are mapped correctly
        """
        
        alignments_store=store.Alignments()
        
        # load in the id_mapping file
        alignments_store.process_id_mapping(cfg.id_mapping_file)
        
        # store some alignments
        alignments_store.add_annotated("query1",1,"ref1")
        alignments_store.add_annotated("query2",1,"ref2")
        alignments_store.add_annotated("query3",1,"ref3")
        
        # test the genes are correct
        self.assertEqual(sorted(alignments_store.gene_list()),sorted(["gene1","gene2","gene3"]))
        
    def test_Alignments_id_mapping_all_bug_list(self):
        """
        Test the store_id_mapping function
        Test the add_annotated and process_reference_annotation with id mapping
        Test the bugs are mapped correctly
        """
        
        alignments_store=store.Alignments()
        
        # load in the id_mapping file
        alignments_store.process_id_mapping(cfg.id_mapping_file)
        
        # store some alignments
        alignments_store.add_annotated("query1",1,"ref1")
        alignments_store.add_annotated("query2",1,"ref2")
        alignments_store.add_annotated("query3",1,"ref3")
        
        # test the bugs are correct
        self.assertEqual(sorted(alignments_store.bug_list()),sorted(["bug3","unclassified"]))
        
    def test_Alignments_id_mapping_all_hits(self):
        """
        Test the store_id_mapping function
        Test the add_annotated and process_reference_annotation with id mapping
        Test the lengths are mapped correctly
        """
        
        alignments_store=store.Alignments()
        
        # load in the id_mapping file
        alignments_store.process_id_mapping(cfg.id_mapping_file)
        
        # store some alignments
        alignments_store.add_annotated("query1",1,"ref1")
        alignments_store.add_annotated("query2",1,"ref2")
        alignments_store.add_annotated("query3",1,"ref3")
        
        # test the lengths are correct
        stored_lengths=[item[-1] for item in alignments_store.get_hit_list()]
        self.assertEqual(sorted(stored_lengths),sorted([1/1000.0,10/1000.0,1000/1000.0]))
        
    def test_Alignments_id_mapping_half_hits(self):
        """
        Test the store_id_mapping function
        Test the add_annotated and process_reference_annotation with id mapping
        Test the lengths are mapped correctly with only some references included
        in those provided for id mapping
        """
        
        alignments_store=store.Alignments()
        
        # load in the id_mapping file
        alignments_store.process_id_mapping(cfg.id_mapping_file)
        
        # store some alignments
        alignments_store.add_annotated("query1",1,"ref1")
        alignments_store.add_annotated("query2",1,"ref2")
        alignments_store.add_annotated("query3",1,"ref1|100")
        alignments_store.add_annotated("query3",1,"200|ref2")
        
        # test the lengths are correct
        stored_lengths=[item[-1] for item in alignments_store.get_hit_list()]
        self.assertEqual(sorted(stored_lengths),sorted([1/1000.0,100/1000.0,
            200/1000.0,1000/1000.0]))
        
    def test_Alignments_compute_gene_scores_single_gene_single_query_with_temp_alignment_file(self):
        """
        Test the compute_gene_scores function
        Test one hit for gene with one hit for query
        Test with the temp alignment file
        """
        
        # create a set of hits
        matches1=41.0
        matches2=57.1
        matches3=61.0
        matches4=72.1
        
        gene1_length=2
        gene2_length=3
        gene3_length=4
        
        # Create a set of alignments
        alignments_store=store.Alignments(minimize_memory_use=True)
        alignments_store.add("gene1",gene1_length,"query1",matches1,"bug1")
        alignments_store.add("gene2",gene2_length,"query1",matches2,"bug1")
        alignments_store.add("gene2",gene2_length,"query2",matches3,"bug1")
        alignments_store.add("gene3",gene3_length,"query3",matches4,"bug1")
        
        gene_scores_store=store.GeneScores()
        
        # compute gene scores
        alignments_store.convert_alignments_to_gene_scores(gene_scores_store)
        
        # convert lengths to per kb
        gene3_length=gene3_length/1000.0
        
        # gene3
        hit4_score=math.pow(matches4, config.match_power)
        query3_sum=hit4_score
        expected_gene_score=hit4_score/query3_sum/gene3_length
        
        actual_gene_score=gene_scores_store.get_score("bug1","gene3")
        
        # delete the temp alignment file
        alignments_store.delete_temp_alignments_file()  
        
        self.assertEqual(actual_gene_score,expected_gene_score)

    def test_Alignments_compute_gene_scores_single_gene_double_query_with_temp_alignment_file(self):
        """
        Test the compute_gene_scores function
        Test one hit for gene with more than one hit per query
        Test with the temp alignment file
        """
        
        # create a set of hits
        # bug, reference, reference_length, query, matches = hit
        
        matches1=41.0
        matches2=57.1
        matches3=61.0
        matches4=72.1
        
        gene1_length=2
        gene2_length=3
        gene3_length=4
        
        # Create a set of alignments
        alignments_store=store.Alignments(minimize_memory_use=True)
        alignments_store.add("gene1",gene1_length,"query1",matches1,"bug1")
        alignments_store.add("gene2",gene2_length,"query1",matches2,"bug1")
        alignments_store.add("gene2",gene2_length,"query2",matches3,"bug1")
        alignments_store.add("gene3",gene3_length,"query3",matches4,"bug1")
        
        gene_scores_store=store.GeneScores()
        
        # compute gene scores
        alignments_store.convert_alignments_to_gene_scores(gene_scores_store)
        
        # convert lengths to per kb
        gene1_length=gene1_length/1000.0
        
        # gene1
        hit1_score=math.pow(matches1, config.match_power)
        hit2_score=math.pow(matches2, config.match_power)
        query1_sum=hit1_score+hit2_score
        expected_gene_score=hit1_score/query1_sum/gene1_length   
        
        actual_gene_score=gene_scores_store.get_score("bug1","gene1")
        
        # delete the temp alignment file
        alignments_store.delete_temp_alignments_file()  

        self.assertAlmostEqual(actual_gene_score,expected_gene_score)

    def test_Alignments_compute_gene_scores_double_gene_double_query_with_temp_alignment_file(self):
        """
        Test the compute_gene_scores function
        Test two hits to gene with more than one hit per query
        Test with the temp alignment file
        """
        
        # create a set of hits
        # bug, reference, reference_length, query, matches = hit
        
        matches1=41.0
        matches2=57.1
        matches3=61.0
        matches4=72.1
        
        gene1_length=2
        gene2_length=3
        gene3_length=4
        
        # Create a set of alignments
        alignments_store=store.Alignments(minimize_memory_use=True)
        alignments_store.add("gene1",gene1_length,"query1",matches1,"bug1")
        alignments_store.add("gene2",gene2_length,"query1",matches2,"bug1")
        alignments_store.add("gene2",gene2_length,"query2",matches3,"bug1")
        alignments_store.add("gene3",gene3_length,"query3",matches4,"bug1")
        
        gene_scores_store=store.GeneScores()
        
        # compute gene scores
        alignments_store.convert_alignments_to_gene_scores(gene_scores_store)
        
        # gene1
        hit1_score=math.pow(matches1, config.match_power)
        hit2_score=math.pow(matches2, config.match_power)
        query1_sum=hit1_score+hit2_score
        
        # convert lengths to per kb
        gene2_length=gene2_length/1000.0
        
        # gene2
        hit3_score=math.pow(matches3, config.match_power)
        query2_sum=hit3_score
        expected_gene_score=hit3_score/query2_sum/gene2_length + hit2_score/query1_sum/gene2_length

        actual_gene_score=gene_scores_store.get_score("bug1","gene2")

        # delete the temp alignment file
        alignments_store.delete_temp_alignments_file()  

        self.assertAlmostEqual(actual_gene_score,expected_gene_score,places=7)
        
    def test_Alignments_id_mapping_all_gene_list_with_temp_alignment_file(self):
        """
        Test the store_id_mapping function
        Test the add_annotated and process_reference_annotation with id mapping
        Test the genes are mapped correctly
        Test with the temp alignment file
        """
        
        alignments_store=store.Alignments(minimize_memory_use=True)
        
        # load in the id_mapping file
        alignments_store.process_id_mapping(cfg.id_mapping_file)
        
        # store some alignments
        alignments_store.add_annotated("query1",1,"ref1")
        alignments_store.add_annotated("query2",1,"ref2")
        alignments_store.add_annotated("query3",1,"ref3")
        
        gene_list=alignments_store.gene_list()
        
        # delete the temp alignment file
        alignments_store.delete_temp_alignments_file() 
        
        # test the genes are correct
        self.assertEqual(sorted(gene_list),sorted(["gene1","gene2","gene3"]))
        
    def test_Alignments_id_mapping_all_bug_list_with_temp_alignment_file(self):
        """
        Test the store_id_mapping function
        Test the add_annotated and process_reference_annotation with id mapping
        Test the bugs are mapped correctly
        Test with the temp alignment file
        """
        
        alignments_store=store.Alignments(minimize_memory_use=True)
        
        # load in the id_mapping file
        alignments_store.process_id_mapping(cfg.id_mapping_file)
        
        # store some alignments
        alignments_store.add_annotated("query1",1,"ref1")
        alignments_store.add_annotated("query2",1,"ref2")
        alignments_store.add_annotated("query3",1,"ref3")
        
        bug_list=alignments_store.bug_list()
        
        # delete the temp alignment file
        alignments_store.delete_temp_alignments_file() 
        
        # test the bugs are correct
        self.assertEqual(sorted(bug_list),sorted(["bug3","unclassified"]))
        
    def test_Alignments_id_mapping_all_hits_with_temp_alignment_file(self):
        """
        Test the store_id_mapping function
        Test the add_annotated and process_reference_annotation with id mapping
        Test the lengths are mapped correctly
        Test with the temp alignment file
        """
        
        alignments_store=store.Alignments(minimize_memory_use=True)
        
        # load in the id_mapping file
        alignments_store.process_id_mapping(cfg.id_mapping_file)
        
        # store some alignments
        alignments_store.add_annotated("query1",1,"ref1")
        alignments_store.add_annotated("query2",1,"ref2")
        alignments_store.add_annotated("query3",1,"ref3")
        
        hit_list=alignments_store.get_hit_list()
        
        # delete the temp alignment file
        alignments_store.delete_temp_alignments_file() 
        
        # test the lengths are correct
        stored_lengths=[item[-1] for item in hit_list]
        self.assertEqual(sorted(stored_lengths),sorted([1/1000.0,10/1000.0,1000/1000.0]))
        
    def test_Alignments_id_mapping_half_hits_with_temp_alignment_file(self):
        """
        Test the store_id_mapping function
        Test the add_annotated and process_reference_annotation with id mapping
        Test the lengths are mapped correctly with only some references included
        in those provided for id mapping
        Test with the temp alignment file
        """
        
        alignments_store=store.Alignments(minimize_memory_use=True)
        
        # load in the id_mapping file
        alignments_store.process_id_mapping(cfg.id_mapping_file)
        
        # store some alignments
        alignments_store.add_annotated("query1",1,"ref1")
        alignments_store.add_annotated("query2",1,"ref2")
        alignments_store.add_annotated("query3",1,"ref1|100")
        alignments_store.add_annotated("query3",1,"200|ref2")
        
        hit_list=alignments_store.get_hit_list()
        
        # delete the temp alignment file
        alignments_store.delete_temp_alignments_file() 
        
        # test the lengths are correct
        stored_lengths=[item[-1] for item in hit_list]
        self.assertEqual(sorted(stored_lengths),sorted([1/1000.0,100/1000.0,
            200/1000.0,1000/1000.0]))
        
    def test_GeneScores_add_from_file_id_mapping_bug_list(self):
        """
        GeneScores class: Test add_from_file bug list with id mapping
        """
        
        gene_scores=store.GeneScores()
        
        gene_scores.add_from_file(cfg.genetable_file, 
            id_mapping_file=cfg.id_mapping_gene_table_file)
        
        # Test the bug list is as expected
        self.assertEqual(sorted(cfg.genetable_file_bug_scores_id_mapping.keys()),
            sorted(gene_scores.bug_list()))
        
    def test_GeneScores_add_from_file_id_mapping_gene_list(self):
        """
        GeneScores class: Test add_from_file gene list with id mapping
        """
        
        gene_scores=store.GeneScores()
        
        gene_scores.add_from_file(cfg.genetable_file,
            id_mapping_file=cfg.id_mapping_gene_table_file)
        
        # Create a list of all of the genes in the table
        genes={}
        for bug in cfg.genetable_file_bug_scores_id_mapping:
            for gene in cfg.genetable_file_bug_scores_id_mapping[bug]:
                genes[gene]=1
        
        # Test the gene list is as expected
        self.assertEqual(sorted(genes.keys()),sorted(gene_scores.gene_list()))
        
    def test_PathwaysDatabase_is_structured_structure(self):
        """
        Pathways database class: Test the storing of a structured set of pathways
        Test this file is identified as structured
        """
        
        pathways_database_store=store.PathwaysDatabase(cfg.pathways_file)
        
        self.assertTrue(pathways_database_store.is_structured())
        
    def test_PathwaysDatabase_is_structured_unstructure(self):
        """
        Pathways database class: Test the storing of a unstructured set of pathways
        Test this file is identified as unstructured
        """
        
        pathways_database_flat_store=store.PathwaysDatabase(cfg.pathways_flat_file)
        
        self.assertTrue(not pathways_database_flat_store.is_structured())
        
    def test_PathwaysDatabase_add_pathway_structure_test_key_reactions_not_included_in_reactions_database(self):
        """
        Pathways database class: Test the add pathway structure
        Test the function with a structure with two starting points that contract
        Test the key reactions are correct for reactions that are not included in the reactions database
        Test that key reactions included are correct for reactions with 1 and 3 optional indicators
        """
        
        # Create a reactions database of a subset of the reactions in the pathways
        reactions_database_store=store.ReactionsDatabase()
        
        reactions={ "A": ["gene1","gene2"], "---B": ["gene3"], "--Z":["gene4"],"-F":["gene5"]}
        reactions_database_store.add_reactions(reactions)
        
        pathways_database_store=store.PathwaysDatabase()
        
        structure_string="( (  L A ---B ) , ( --Z ---C D ) )  -E -F"
        
        pathways_database_store.add_pathway_structure("pathway1",structure_string,reactions_database_store)
        
        expected_key_reactions=["L","A","---B","--Z","D","-F"]
        
        self.assertEqual(expected_key_reactions,pathways_database_store.get_key_reactions_for_pathway("pathway1"))
        
        
