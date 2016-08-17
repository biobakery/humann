import unittest
import importlib
import re
import tempfile
import os
import logging

import cfg
import utils

from humann2 import store
from humann2 import config

class TestHumann2StoreFunctions(unittest.TestCase):
    """
    Test the functions found in humann2.store
    """
    
    def setUp(self):
        # set up nullhandler for logger
        logging.getLogger('humann2.store').addHandler(logging.NullHandler())
        
        config.unnamed_temp_dir=tempfile.gettempdir()
        
    def test_Read_add_count_reads_duplicate(self):
        """
        Read class: Test the adding of a set of reads
        Test the count is correct with duplicate ids
        """
        
        reads_store=store.Reads()
        
        reads_store.add("id1","ATCG")
        reads_store.add("id2","ATTG")
        reads_store.add("id1","ATCG")
        
        self.assertEqual(reads_store.count_reads(), 2)
        
    def test_Read_add_count_reads_duplicate_minimize_memory_use(self):
        """
        Read class: Test the adding of a set of reads
        Test the count is correct with duplicate ids
        Test with minimize memory use
        """

        reads_store=store.Reads(minimize_memory_use=True)
        
        reads_store.add("id1","ATCG")
        reads_store.add("id2","ATTG")
        reads_store.add("id1","ATCG")
        
        self.assertEqual(reads_store.count_reads(), 2)

    def test_Read_print_fasta_id_count(self):
        """
        Read class: Test the loading of a full fasta file
        Test the total number of expected ids are loaded
        """
        
        reads_store=store.Reads(cfg.small_fasta_file)
        
        # Check that the total number of expected reads are loaded
        self.assertEqual(len(reads_store.id_list()), cfg.small_fasta_file_total_sequences)
        
    def test_Read_print_fasta_count_reads(self):
        """
        Read class: Test the loading of a full fasta file
        Test the total number of expected reads counted
        """
        
        reads_store=store.Reads(cfg.small_fasta_file)
        
        # Check that the total number of expected reads are counted
        self.assertEqual(reads_store.count_reads(), cfg.small_fasta_file_total_sequences)
            
    def test_Read_print_fasta_count_reads_minimize_memory_use(self):
        """
        Read class: Test the loading of a full fasta file
        Test the total number of expected reads counted
        Test with minimize memory use
        """
        
        reads_store=store.Reads(cfg.small_fasta_file, minimize_memory_use=True)
        
        # Check that the total number of expected reads are counted
        self.assertEqual(reads_store.count_reads(), cfg.small_fasta_file_total_sequences)            


    def test_Read_print_fasta_id_list(self):
        """
        Read class: Test the loading of a full fasta file
        Test the expected ids are loaded
        """
        
        reads_store=store.Reads(cfg.small_fasta_file)
        
        # Check the reads are printed correctly
        stored_fasta=[]
        for line in reads_store.get_fasta():
            stored_fasta.append(line)
            
        printed_stored_fasta="\n".join(stored_fasta)
        
        compare_fasta={}
        # organize the fasta from the read class and the 
        # file of correct fasta output
        file_handle=open(cfg.small_fasta_file_single_line_sequences)
        for input in [printed_stored_fasta.split("\n"), file_handle]:
            id=""
            seq=""
            for line in input:
                if re.search(">",line):
                    # store prior id
                    if id and seq:
                        compare_fasta[id]=compare_fasta.get(id,[])+[seq]
                    id=line.strip()
                    seq=""
                else:
                    seq=line.strip()
                    
            # store the last sequence found
            if id and seq:
                compare_fasta[id]=compare_fasta.get(id,[])+[seq]
        
        file_handle.close()
        
        # check there are still the same number of ids
        self.assertEqual(len(compare_fasta.keys()),cfg.small_fasta_file_total_sequences)
        
        # check the sequences match
        for id, sequences in compare_fasta.items():
            self.assertTrue(len(sequences)==2)
            self.assertEqual(sequences[0], sequences[1])
            
    def test_Read_print_fasta_id_list_minimize_memory_use(self):
        """
        Read class: Test the loading of a full fasta file
        Test the expected ids are loaded
        Test with minimize memory use
        """
        
        reads_store=store.Reads(cfg.small_fasta_file, minimize_memory_use=True)
        
        # Check the reads are printed correctly
        stored_fasta=[]
        for line in reads_store.get_fasta():
            stored_fasta.append(line)
            
        printed_stored_fasta="\n".join(stored_fasta)
        
        compare_fasta={}
        # organize the fasta from the read class and the 
        # file of correct fasta output
        file_handle=open(cfg.small_fasta_file_single_line_sequences)
        for input in [printed_stored_fasta.split("\n"), file_handle]:
            id=""
            seq=""
            for line in input:
                if re.search(">",line):
                    # store prior id
                    if id and seq:
                        compare_fasta[id]=compare_fasta.get(id,[])+[seq]
                    id=line.strip()
                    seq=""
                else:
                    seq=line.strip()
                    
            # store the last sequence found
            if id and seq:
                compare_fasta[id]=compare_fasta.get(id,[])+[seq]
        
        file_handle.close()
        
        # check there are still the same number of ids
        self.assertEqual(len(compare_fasta.keys()),cfg.small_fasta_file_total_sequences)
        
        # check the sequences match
        for id, sequences in compare_fasta.items():
            self.assertTrue(len(sequences)==2)
            self.assertEqual(sequences[0], sequences[1])
            
    def test_Read_print_fasta_sequence_list(self):
        """
        Read class: Test the loading of a full fasta file
        Test the sequences are loaded
        """
        
        reads_store=store.Reads(cfg.small_fasta_file)
        
        # Check the reads are printed correctly
        stored_fasta=[]
        for line in reads_store.get_fasta():
            stored_fasta.append(line)
            
        printed_stored_fasta="\n".join(stored_fasta)
        
        compare_fasta={}
        # organize the fasta from the read class and the 
        # file of correct fasta output
        file_handle=open(cfg.small_fasta_file_single_line_sequences)
        for input in [printed_stored_fasta.split("\n"), file_handle]:
            id=""
            seq=""
            for line in input:
                if re.search(">",line):
                    # store prior id
                    if id and seq:
                        compare_fasta[id]=compare_fasta.get(id,[])+[seq]
                    id=line.strip()
                    seq=""
                else:
                    seq=line.strip()
                    
            # store the last sequence found
            if id and seq:
                compare_fasta[id]=compare_fasta.get(id,[])+[seq]
        
        file_handle.close()
        
        # check the sequences match
        for id, sequences in compare_fasta.items():
            self.assertTrue(len(sequences)==2)
            self.assertEqual(sequences[0], sequences[1])
            
    def test_Read_print_fasta_sequence_list_minimize_memory_use(self):
        """
        Read class: Test the loading of a full fasta file
        Test the sequences are loaded
        Test with minimize memory use
        """
        
        reads_store=store.Reads(cfg.small_fasta_file, minimize_memory_use=True)
        
        # Check the reads are printed correctly
        stored_fasta=[]
        for line in reads_store.get_fasta():
            stored_fasta.append(line)
            
        printed_stored_fasta="\n".join(stored_fasta)
        
        compare_fasta={}
        # organize the fasta from the read class and the 
        # file of correct fasta output
        file_handle=open(cfg.small_fasta_file_single_line_sequences)
        for input in [printed_stored_fasta.split("\n"), file_handle]:
            id=""
            seq=""
            for line in input:
                if re.search(">",line):
                    # store prior id
                    if id and seq:
                        compare_fasta[id]=compare_fasta.get(id,[])+[seq]
                    id=line.strip()
                    seq=""
                else:
                    seq=line.strip()
                    
            # store the last sequence found
            if id and seq:
                compare_fasta[id]=compare_fasta.get(id,[])+[seq]
        
        file_handle.close()
        
        # check the sequences match
        for id, sequences in compare_fasta.items():
            self.assertTrue(len(sequences)==2)
            self.assertEqual(sequences[0], sequences[1])
        
    def test_Read_delete_id(self):
        """
        Read class: Test the deleting of ids
        """
        
        reads_store=store.Reads(cfg.small_fasta_file)
        
        # delete all but one of the reads and check structure is empty
        id_list=reads_store.id_list()
        keep_id=id_list.pop()
        
        for id in id_list:
            reads_store.remove_id(id)
            
        self.assertEqual(reads_store.id_list(), [keep_id])
        
    def test_Read_delete_id_minimize_memory_use(self):
        """
        Read class: Test the deleting of ids
        Test with minimial memory use
        """
        
        reads_store=store.Reads(cfg.small_fasta_file, minimize_memory_use=True)
        
        # delete all but one of the reads and check structure is empty
        id_list=reads_store.id_list()
        keep_id=id_list.pop()
        
        for id in id_list:
            reads_store.remove_id(id)
            
        self.assertEqual(reads_store.id_list(), [keep_id])
        
    def test_PathwaysDatabase_read_pathways_count(self):
        """
        Pathways database class: Test the storing of a structured set of pathways
        Test for the number of pathways
        """
        
        pathways_database_store=store.PathwaysDatabase(cfg.pathways_file)
        pathways_database_flat_store=store.PathwaysDatabase(cfg.pathways_flat_file)
        
        # check for the same number of pathways
        pathway_list=pathways_database_store.pathway_list()
        pathway_list_flat=pathways_database_flat_store.pathway_list()
        
        self.assertEqual(len(pathway_list),len(pathway_list_flat))
            
    def test_PathwaysDatabase_read_pathways_ids(self):
        """
        Pathways database class: Test the storing of a structured set of pathways
        Test for the pathways ids
        """
        
        pathways_database_store=store.PathwaysDatabase(cfg.pathways_file)
        pathways_database_flat_store=store.PathwaysDatabase(cfg.pathways_flat_file)
        
        # check for the same number of pathways
        pathway_list=pathways_database_store.pathway_list()
        pathway_list_flat=pathways_database_flat_store.pathway_list()
        
        # check that the pathway ids are identical
        for pathway in pathway_list:
            self.assertTrue(pathway in pathway_list_flat)

    def test_PathwaysDatabase_read_reactions_list(self):
        """
        Pathways database class: Test the storing of a structured set of pathways
        Test for the reactions list
        """
        
        pathways_database_store=store.PathwaysDatabase(cfg.pathways_file)
        pathways_database_flat_store=store.PathwaysDatabase(cfg.pathways_flat_file)
        
        # check for the same number of pathways
        pathway_list=pathways_database_store.pathway_list()
        pathway_list_flat=pathways_database_flat_store.pathway_list()
            
        # check that the reactions list for each pathway is identical
        for pathway in pathway_list:
            self.assertEqual(sorted(pathways_database_store.find_reactions(pathway)),
                sorted(pathways_database_flat_store.find_reactions(pathway)))
            
    def test_PathwaysDatabase_read_reactions_count(self):
        """
        Pathways database class: Test the storing of a structured set of pathways
        Test for the number of reactions 
        """
        
        pathways_database_store=store.PathwaysDatabase(cfg.pathways_file)
        pathways_database_flat_store=store.PathwaysDatabase(cfg.pathways_flat_file)
        
        # check for the same number of pathways
        pathway_list=pathways_database_store.pathway_list()
        pathway_list_flat=pathways_database_flat_store.pathway_list()
        
        # check for the same number of reactions
        self.assertEqual(len(pathways_database_store.reaction_list()), 
            len(pathways_database_flat_store.reaction_list()))
        
        # check that the reactions are identical
        reaction_list=pathways_database_store.reaction_list()
        reaction_list_flat=pathways_database_flat_store.reaction_list()
        for reaction in reaction_list:
            self.assertTrue(reaction in reaction_list_flat)

    def test_PathwaysDatabase_read_reactions_ids(self):
        """
        Pathways database class: Test the storing of a structured set of pathways
        Test for the reactions ids 
        """
        
        pathways_database_store=store.PathwaysDatabase(cfg.pathways_file)
        pathways_database_flat_store=store.PathwaysDatabase(cfg.pathways_flat_file)
        
        # check for the same number of pathways
        pathway_list=pathways_database_store.pathway_list()
        pathway_list_flat=pathways_database_flat_store.pathway_list()
        
        # check that the reactions are identical
        reaction_list=pathways_database_store.reaction_list()
        reaction_list_flat=pathways_database_flat_store.reaction_list()
        for reaction in reaction_list:
            self.assertTrue(reaction in reaction_list_flat)
            
    def test_PathwaysDatabase_print_flat_file_pathways_count(self):
        """
        Pathways database class: Test the printing of a flat file from a structured file
        Test for the total number of pathways
        """
 
        pathways_database_store=store.PathwaysDatabase(cfg.pathways_file)
        pathways_database_flat_store=store.PathwaysDatabase(cfg.pathways_flat_file)       
        
        # write the flat file created from a structured file to a temp file
        file_out, new_file=tempfile.mkstemp()
        os.close(file_out)
        with open(new_file, "w") as file_handle:
            file_handle.write(pathways_database_store.get_database())
        
        # load in the flat file and compare with the correct flat file
        pathways_database_flat_store_write=store.PathwaysDatabase(new_file)
        
        # remove the temp file
        utils.remove_temp_file(new_file)
        
        # check for the same number of pathways
        pathway_list=pathways_database_flat_store_write.pathway_list()
        pathway_list_flat=pathways_database_flat_store.pathway_list()
        
        self.assertEqual(len(pathway_list),len(pathway_list_flat))
            
    def test_PathwaysDatabase_print_flat_file_pathways_list(self):
        """
        Pathways database class: Test the printing of a flat file from a structured file
        Test for the pathways list
        """
 
        pathways_database_store=store.PathwaysDatabase(cfg.pathways_file)
        pathways_database_flat_store=store.PathwaysDatabase(cfg.pathways_flat_file)       
        
        # write the flat file created from a structured file to a temp file
        file_out, new_file=tempfile.mkstemp()
        os.close(file_out)
        with open(new_file, "w") as file_handle:
            file_handle.write(pathways_database_store.get_database())      
  
        # load in the flat file and compare with the correct flat file
        pathways_database_flat_store_write=store.PathwaysDatabase(new_file)
        
        # remove the temp file
        utils.remove_temp_file(new_file)
        
        # check for the same number of pathways
        pathway_list=pathways_database_flat_store_write.pathway_list()
        pathway_list_flat=pathways_database_flat_store.pathway_list()
        
        # check that the pathway ids are identical
        for pathway in pathway_list:
            self.assertTrue(pathway in pathway_list_flat)
            
    def test_PathwaysDatabase_print_flat_file_reactions_list(self):
        """
        Pathways database class: Test the printing of a flat file from a structured file
        Test the reactions list
        """
 
        pathways_database_store=store.PathwaysDatabase(cfg.pathways_file)
        pathways_database_flat_store=store.PathwaysDatabase(cfg.pathways_flat_file)       
        
        # write the flat file created from a structured file to a temp file
        file_out, new_file=tempfile.mkstemp()
        os.close(file_out)
        with open(new_file, "w") as file_handle:
            file_handle.write(pathways_database_store.get_database())       
 
        # load in the flat file and compare with the correct flat file
        pathways_database_flat_store_write=store.PathwaysDatabase(new_file)
        
        # remove the temp file
        utils.remove_temp_file(new_file)
        
        # check for the same number of pathways
        pathway_list=pathways_database_flat_store_write.pathway_list()
        pathway_list_flat=pathways_database_flat_store.pathway_list()
        
        # check that the reactions list for each pathway is identical
        for pathway in pathway_list:
            self.assertEqual(sorted(pathways_database_flat_store_write.find_reactions(pathway)),
                sorted(pathways_database_flat_store.find_reactions(pathway)))
            
    def test_PathwaysDatabase_add_pathway_structure_test_structure_expansion(self):
        """
        Pathways database class: Test the add pathway structure
        Test the function with a structure with one starting point and an expansion
        Test the structure is being stored correctly
        """
        
        pathways_database_store=store.PathwaysDatabase()
        
        structure_string="A ( ( B C ) , ( D E ) )"
        
        pathways_database_store.add_pathway_structure("pathway1",structure_string)
        
        expected_structure=[" ","A",[",",[" ","B","C"],[" ","D","E"]]]
        
        self.assertEqual(expected_structure,pathways_database_store.get_structure_for_pathway("pathway1"))
        
    def test_PathwaysDatabase_add_pathway_structure_test_structure_contraction(self):
        """
        Pathways database class: Test the add pathway structure
        Test the function with a structure with two starting points that contract
        Test the structure is being stored correctly
        """
        
        pathways_database_store=store.PathwaysDatabase()
        
        structure_string="( (  L A B ) , ( Z C D ) )  E F"
        
        pathways_database_store.add_pathway_structure("pathway1",structure_string)
        
        expected_structure=[" ",[",",[" ","L","A","B"],[" ","Z","C","D"]],"E","F"]
        
        self.assertEqual(expected_structure,pathways_database_store.get_structure_for_pathway("pathway1"))
        
    def test_PathwaysDatabase_add_pathway_structure_test_reaction_list_expansion(self):
        """
        Pathways database class: Test the add pathway structure
        Test the function with a structure with one starting point and an expansion
        Test the reaction list is correct
        """
        
        pathways_database_store=store.PathwaysDatabase()
        
        structure_string="A ( ( B -C ) , ( -D E ) )"
        
        pathways_database_store.add_pathway_structure("pathway1",structure_string)
        
        expected_reaction_list=["A","B","C","D","E"]
        
        self.assertEqual(expected_reaction_list,pathways_database_store.find_reactions("pathway1"))
        
    def test_PathwaysDatabase_add_pathway_structure_test_reaction_list_contraction(self):
        """
        Pathways database class: Test the add pathway structure
        Test the function with a structure with two starting points that contract
        Test the reaction list is correct
        """
        
        pathways_database_store=store.PathwaysDatabase()
        
        structure_string="( (  L A -B ) , ( Z -C D ) )  -E F"
        
        pathways_database_store.add_pathway_structure("pathway1",structure_string)
        
        expected_reaction_list=["L","A","B","Z","C","D","E","F"]
        
        self.assertEqual(expected_reaction_list,pathways_database_store.find_reactions("pathway1"))
        
    def test_PathwaysDatabase_add_pathway_structure_test_key_reactions_expansion(self):
        """
        Pathways database class: Test the add pathway structure
        Test the function with a structure with one starting point and an expansion
        Test the key reactions are correct
        """
        
        pathways_database_store=store.PathwaysDatabase()
        
        structure_string="A ( ( B -C ) , ( -D E ) )"
        
        pathways_database_store.add_pathway_structure("pathway1",structure_string)
        
        expected_key_reactions=["A","B","E"]
        
        self.assertEqual(expected_key_reactions,pathways_database_store.get_key_reactions_for_pathway("pathway1"))
        
    def test_PathwaysDatabase_add_pathway_structure_test_key_reactions_contraction(self):
        """
        Pathways database class: Test the add pathway structure
        Test the function with a structure with two starting points that contract
        Test the key reactions are correct
        """
        
        pathways_database_store=store.PathwaysDatabase()
        
        structure_string="( (  L A -B ) , ( Z -C D ) )  -E F"
        
        pathways_database_store.add_pathway_structure("pathway1",structure_string)
        
        expected_key_reactions=["L","A","Z","D","F"]
        
        self.assertEqual(expected_key_reactions,pathways_database_store.get_key_reactions_for_pathway("pathway1"))

    def test_PathwaysDatabase_add_pathway_structure_test_key_reactions_with_optional_indicator(self):
        """
        Pathways database class: Test the add pathway structure
        Test the function with a structure with two starting points that contract
        Test the key reactions are correct for reactions with "--" at the beginning of their name
        """
        
        pathways_database_store=store.PathwaysDatabase()
        
        structure_string="( (  L A -B ) , ( --Z ---C D ) )  -E F"
        
        pathways_database_store.add_pathway_structure("pathway1",structure_string)
        
        expected_key_reactions=["L","A","--Z","D","F"]
        
        self.assertEqual(expected_key_reactions,pathways_database_store.get_key_reactions_for_pathway("pathway1"))
            
    def test_ReactionsDatabase_read_reactions_count(self):
        """
        Reactions Database class: Test the storing of reactions
        Test the total number of reactions
        """
        
        reactions_database_store=store.ReactionsDatabase(cfg.reactions_file)
        
        # read in the reactions directly from the file
        file_handle=open(cfg.reactions_file)
        
        reactions={}
        for line in file_handle:
            data=line.strip().split("\t")
            reactions[data[0]]=data[2:]
        file_handle.close()
        
        # test for the same number of reactions
        self.assertEqual(len(reactions.keys()), len(reactions_database_store.reaction_list()))
            
    def test_ReactionsDatabase_read_gene_list(self):
        """
        Reactions Database class: Test the storing of reactions
        Test for the gene list
        """
        
        reactions_database_store=store.ReactionsDatabase(cfg.reactions_file)
        
        # read in the reactions directly from the file
        file_handle=open(cfg.reactions_file)
        
        reactions={}
        for line in file_handle:
            data=line.strip().split("\t")
            reactions[data[0]]=data[2:]
        file_handle.close()
        
        # test for the same reactions and genes
        for rxn in reactions:
            self.assertEqual(reactions[rxn],reactions_database_store.find_genes(rxn))
            
    def test_PathwaysAndReactions_median_score_odd_number_vary_reactions(self):
        """
        Pathways and Reactions class: Test add and median score
        Test and odd number of values and different reactions
        """
        
        pathways_and_reactions=store.PathwaysAndReactions()
        
        # add scores all for same pathway and different reaction
        pathways_and_reactions.add("bug","R1","P1",1)
        pathways_and_reactions.add("bug","R2","P1",2)        
        pathways_and_reactions.add("bug","R3","P1",3)
        pathways_and_reactions.add("bug","R4","P1",4)
        pathways_and_reactions.add("bug","R5","P1",5)  
        
        # test median score for odd number of values
        self.assertEqual(pathways_and_reactions.median_score("bug"),3)
        
    def test_PathwaysAndReactions_median_score_even_number_vary_reactions(self):
        """
        Pathways and Reactions class: Test add and median score
        Test an even number of values and different reactions
        """
        
        pathways_and_reactions=store.PathwaysAndReactions()
        
        # add scores all for same pathway and  different reaction
        pathways_and_reactions.add("bug","R1","P1",1)
        pathways_and_reactions.add("bug","R2","P1",2)        
        pathways_and_reactions.add("bug","R3","P1",3)
        pathways_and_reactions.add("bug","R4","P1",4) 
        
        # test median score for an even number of values
        self.assertEqual(pathways_and_reactions.median_score("bug"),2.5)    
        
    def test_PathwaysAndReactions_median_score_odd_number_vary_pathways(self):
        """
        Pathways and Reactions class: Test add and median score
        Test an odd number of values and different pathways
        """
        
        pathways_and_reactions=store.PathwaysAndReactions()        
        
        # add scores all for same pathway and reaction
        pathways_and_reactions.add("bug","R1","P1",1)
        pathways_and_reactions.add("bug","R1","P2",2)        
        pathways_and_reactions.add("bug","R1","P3",3)
        pathways_and_reactions.add("bug","R1","P4",4)
        pathways_and_reactions.add("bug","R1","P5",5)  
        
        # test median score for odd number of values
        self.assertEqual(pathways_and_reactions.median_score("bug"),3)
        
    def test_PathwaysAndReactions_median_score_even_number_vary_pathways(self):
        """
        Pathways and Reactions class: Test add and median score
        Test an even number of values and different pathways
        """
        
        pathways_and_reactions=store.PathwaysAndReactions()
        
        # add scores all for same pathway and reaction
        pathways_and_reactions.add("bug","R1","P1",1)
        pathways_and_reactions.add("bug","R1","P2",2)        
        pathways_and_reactions.add("bug","R1","P3",3)
        pathways_and_reactions.add("bug","R1","P4",4) 
        
        # test median score for an even number of values
        self.assertEqual(pathways_and_reactions.median_score("bug"),2.5)     
        
    def test_PathwaysAndReactions_max_median_score_odd_number_vary_reactions(self):
        """
        Pathways and Reactions class: Test add and max median score
        Test and odd number of values and different reactions
        """
        
        pathways_and_reactions=store.PathwaysAndReactions()
        
        # add scores all for same pathway and different reaction
        pathways_and_reactions.add("bug","R1","P1",1)
        pathways_and_reactions.add("bug","R2","P1",2)        
        pathways_and_reactions.add("bug","R3","P1",3)
        pathways_and_reactions.add("bug","R4","P1",4)
        pathways_and_reactions.add("bug","R5","P1",5)  
        
        # test median score for odd number of values
        self.assertEqual(pathways_and_reactions.max_median_score("bug"),5)
        
    def test_PathwaysAndReactions_max_median_score_even_number_vary_reactions(self):
        """
        Pathways and Reactions class: Test add and max median score
        Test an even number of values and different reactions
        """
        
        pathways_and_reactions=store.PathwaysAndReactions()
        
        # add scores all for same pathway and  different reaction
        pathways_and_reactions.add("bug","R1","P1",1)
        pathways_and_reactions.add("bug","R2","P1",2)        
        pathways_and_reactions.add("bug","R3","P1",3)
        pathways_and_reactions.add("bug","R4","P1",4) 
        
        # test median score for an even number of values
        self.assertEqual(pathways_and_reactions.max_median_score("bug"),4)    
        
    def test_PathwaysAndReactions_max_median_score_odd_number_vary_pathways(self):
        """
        Pathways and Reactions class: Test add and max median score
        Test an odd number of values and different pathways
        """
        
        pathways_and_reactions=store.PathwaysAndReactions()        
        
        # add scores all for same pathway and reaction
        pathways_and_reactions.add("bug","R1","P1",0.9)
        pathways_and_reactions.add("bug","R1","P1",1)
        pathways_and_reactions.add("bug","R1","P2",2)        
        pathways_and_reactions.add("bug","R1","P3",3)
        pathways_and_reactions.add("bug","R1","P4",4)
        pathways_and_reactions.add("bug","R1","P5",5)  
        
        # test median score for odd number of values
        self.assertEqual(pathways_and_reactions.max_median_score("bug"),3)
        
    def test_PathwaysAndReactions_max_median_score_even_number_vary_pathways(self):
        """
        Pathways and Reactions class: Test add and max median score
        Test an even number of values and different pathways
        """
        
        pathways_and_reactions=store.PathwaysAndReactions()
        
        # add scores all for same pathway and reaction
        pathways_and_reactions.add("bug","R1","P1",0.9)
        pathways_and_reactions.add("bug","R1","P1",1)
        pathways_and_reactions.add("bug","R1","P2",2)        
        pathways_and_reactions.add("bug","R1","P3",3)
        pathways_and_reactions.add("bug","R1","P4",4) 
        
        # test median score for an even number of values
        self.assertEqual(pathways_and_reactions.max_median_score("bug"),2.5)       
        
    def test_Alignments_add_bug_count(self):
        """
        Alignments class: Test add function
        Test the total bugs
        """             
        
        alignments_store=store.Alignments()
        
        
        alignments_store.add("gene2", 1, "Q3", 0.01, "bug1",1)
        alignments_store.add("gene1", 1, "Q1", 0.01, "bug2",1)
        alignments_store.add("gene3", 1, "Q2", 0.01, "bug3",1)
        alignments_store.add("gene1", 1, "Q1", 0.01, "bug1",1)
        
        # check the total bugs
        self.assertEqual(alignments_store.count_bugs(),3)
        
    def test_Alignments_add_gene_count(self):
        """
        Alignments class: Test add function
        Test the total genes
        """             
        
        alignments_store=store.Alignments()
        
        alignments_store.add("gene2", 1, "Q3", 0.01, "bug1",1)
        alignments_store.add("gene1", 1, "Q1", 0.01, "bug2",1)
        alignments_store.add("gene3", 1, "Q2", 0.01, "bug3",1)
        alignments_store.add("gene1", 1, "Q1", 0.01, "bug1",1)
        
        # check the total genes
        self.assertEqual(alignments_store.count_genes(),3)
        
    def test_Alignments_add_bug_list(self):
        """
        Alignments class: Test add function
        Test the bug list
        """             
        
        alignments_store=store.Alignments()
        
        alignments_store.add("gene2", 1, "Q3", 0.01, "bug1",1)
        alignments_store.add("gene1", 1, "Q1", 0.01, "bug2",1)
        alignments_store.add("gene3", 1, "Q2", 0.01, "bug3",1)
        alignments_store.add("gene1", 1, "Q1", 0.01, "bug1",1)
        
        # check bug list
        self.assertEqual(sorted(alignments_store.bug_list()),["bug1","bug2","bug3"])
        
    def test_Alignments_add_gene_list(self):
        """
        Alignments class: Test add function
        Test the gene list
        """             
        
        alignments_store=store.Alignments()
        
        alignments_store.add("gene2", 1, "Q3", 0.01, "bug1",1)
        alignments_store.add("gene1", 1, "Q1", 0.01, "bug2",1)
        alignments_store.add("gene3", 1, "Q2", 0.01, "bug3",1)
        alignments_store.add("gene1", 1, "Q1", 0.01, "bug1",1)
        
        # check gene list
        self.assertEqual(sorted(alignments_store.gene_list()),["gene1","gene2","gene3"])     
        
    def test_Alignments_add_gene_lengths(self):
        """
        Alignments class: Test add function
        Test the gene lengths
        """             
        
        alignments_store=store.Alignments()
        
        alignments_store.add("gene2", 10, "Q3", 0.01, "bug1",1)
        alignments_store.add("gene1", 100, "Q1", 0.01, "bug2",1)
        alignments_store.add("gene3", 1000, "Q2", 0.01, "bug3",1)
        alignments_store.add("gene1", 0, "Q1", 0.01, "bug1",1)
        
        # test the lengths are correct
        stored_lengths=[item[-1] for item in alignments_store.get_hit_list()]
        self.assertEqual(sorted(stored_lengths),sorted([10/1000.0,100/1000.0,1000/1000.0,1000/1000.0])) 
        
    def test_Alignments_add_gene_lengths_with_read_length_normalization(self):
        """
        Alignments class: Test add function
        Test the gene lengths with read length normalization
        Test setting the average read length
        """             
        
        alignments_store=store.Alignments()
        
        # set the average read length
        average_read_length=100
        
        alignments_store.add("gene2", 10, "Q3", 0.01, "bug1",average_read_length)
        alignments_store.add("gene1", 100, "Q1", 0.01, "bug2",average_read_length)
        alignments_store.add("gene3", 1000, "Q2", 0.01, "bug3",average_read_length)
        alignments_store.add("gene1", 0, "Q1", 0.01, "bug1",average_read_length)
        
        # test the lengths are correct
        stored_lengths=[item[-1] for item in alignments_store.get_hit_list()]
        self.assertEqual(sorted(stored_lengths),sorted([1/1000.0,91/1000.0,901/1000.0,901/1000.0]))   
              
    def test_Alignments_add_gene_lengths_with_temp_alignment_file(self):
        """
        Alignments class: Test add function
        Test the gene lengths
        Test using the temp alignment file instead of storing data in memory
        """             
        
        alignments_store=store.Alignments(minimize_memory_use=True)
        
        alignments_store.add("gene2", 10, "Q3", 0.01, "bug1",1)
        alignments_store.add("gene1", 100, "Q1", 0.01, "bug2",1)
        alignments_store.add("gene3", 1000, "Q2", 0.01, "bug3",1)
        alignments_store.add("gene1", 0, "Q1", 0.01, "bug1",1)
        
        # test the lengths are correct
        stored_lengths=[item[-1] for item in alignments_store.get_hit_list()]
        
        # delete the temp alignment file
        alignments_store.delete_temp_alignments_file()
        
        self.assertEqual(sorted(stored_lengths),sorted([10/1000.0,100/1000.0,1000/1000.0,1000/1000.0])) 
        
    def test_Alignments_add_gene_lengths_with_read_length_normalization(self):
        """
        Alignments class: Test add function
        Test the gene lengths with read length normalization
        """             
        
        alignments_store=store.Alignments(minimize_memory_use=True)
        
        # set the average read length
        average_read_length=100
        
        alignments_store.add("gene2", 10, "Q3", 0.01, "bug1",average_read_length)
        alignments_store.add("gene1", 100, "Q1", 0.01, "bug2",average_read_length)
        alignments_store.add("gene3", 1000, "Q2", 0.01, "bug3",average_read_length)
        alignments_store.add("gene1", 0, "Q1", 0.01, "bug1",average_read_length)
        
        # test the lengths are correct
        stored_lengths=[item[-1] for item in alignments_store.get_hit_list()]
        
        # delete the temp alignment file
        alignments_store.delete_temp_alignments_file()
        
        self.assertEqual(sorted(stored_lengths),sorted([1/1000.0,91/1000.0,901/1000.0,901/1000.0]))    
        
    def test_Alignments_process_chocophlan_length(self):
        """
        Test the process_chocophlan_length with standard length format
        """
        
        alignments_store=store.Alignments()
        
        length=alignments_store.process_chocophlan_length("1-100","gene")
        
        self.assertEqual(length, 100)
        
    def test_Alignments_process_chocophlan_length_multiple(self):
        """
        Test the process_chocophlan_length with multiple lengths
        Test with one length on the reverse strand
        """
        
        alignments_store=store.Alignments()
        
        length=alignments_store.process_chocophlan_length("c:100-1,1-100","gene")
        
        self.assertEqual(length, 200) 
        
    def test_Alignments_process_reference_annotation_gene_length(self):
        """
        Test the process reference annotation function with a gene and length
        """
        
        alignments_store=store.Alignments()
        
        output=alignments_store.process_reference_annotation("gene|3000")
        
        expected_output=["gene",3000,"unclassified"]
        
        self.assertEqual(expected_output,output) 
        
    def test_Alignments_process_reference_annotation_gene_length_with_bug(self):
        """
        Test the process reference annotation function with a gene and length and bug
        """
        
        alignments_store=store.Alignments()
        
        output=alignments_store.process_reference_annotation("gene|3000|bug")
        
        expected_output=["gene",3000,"bug"]
        
        self.assertEqual(expected_output,output) 
        
    def test_Alignments_process_reference_annotation_gene_length_reversed(self):
        """
        Test the process reference annotation function with a gene and length reversed
        """
        
        alignments_store=store.Alignments()
        
        output=alignments_store.process_reference_annotation("3000|gene")
        
        expected_output=["gene",3000,"unclassified"]
        
        self.assertEqual(expected_output,output) 
        
    def test_Alignments_process_reference_annotation_numerical_gene_length(self):
        """
        Test the process reference annotation function with gene (as number) and length
        """
        
        alignments_store=store.Alignments()
        
        output=alignments_store.process_reference_annotation("59787|5000")
        
        expected_output=["59787",5000,"unclassified"]
        
        self.assertEqual(expected_output,output) 
        
    def test_Alignments_process_reference_annotation_original_chocophlan_annotations(self):
        """
        Test the process reference annotation function with the original chocophlan annotations
        """
        
        alignments_store=store.Alignments()
        
        output=alignments_store.process_reference_annotation(
            "gi|554771211|gb|ACIN03000006.1|:c1189-5|46125|g__Abiotrophia.s__Abiotrophia_defectiva|UniRef90_W1Q3F0|UniRef50_P59787")
        
        expected_output=["UniRef50_P59787",(1189-5+1),"g__Abiotrophia.s__Abiotrophia_defectiva"]
        
        self.assertEqual(expected_output,output) 
        
    def test_Alignments_process_reference_annotation_new_chocophlan_annotations(self):
        """
        Test the process reference annotation function with the new chocophlan annotations
        """
        
        alignments_store=store.Alignments()
        
        output=alignments_store.process_reference_annotation(
            "gi|554771211|gb|ACIN03000006.1|:c1189-5|46125|g__Abiotrophia.s__Abiotrophia_defectiva|UniRef90_W1Q3F0|UniRef50_P59787|5000")
        
        expected_output=["UniRef50_P59787",5000,"g__Abiotrophia.s__Abiotrophia_defectiva"]
        
        self.assertEqual(expected_output,output) 
        
    def test_Alignments_process_reference_annotation_unknown_annotations_three_items_bug_int(self):
        """
        Test the process reference annotation function with unknown annotations (three items) with int as bug
        """
        
        alignments_store=store.Alignments()
        
        output=alignments_store.process_reference_annotation("UniRef90_W1Q3F0|5000|5000")
        
        expected_output=["UniRef90_W1Q3F0|5000|5000",0,"unclassified"]
        
        self.assertEqual(expected_output,output) 
        
    def test_Alignments_process_reference_annotation_unknown_annotations_three_items_length_string(self):
        """
        Test the process reference annotation function with unknown annotations (three items) with string for length
        """
        
        alignments_store=store.Alignments()
        
        output=alignments_store.process_reference_annotation("UniRef90_W1Q3F0|UniRef50_P59787|5000")
        
        expected_output=["UniRef90_W1Q3F0|UniRef50_P59787|5000",0,"unclassified"]
        
        self.assertEqual(expected_output,output) 
        
    def test_Alignments_process_reference_annotation_unknown_annotations_four_items(self):
        """
        Test the process reference annotation function with unknown annotations (four items)
        """
        
        alignments_store=store.Alignments()
        
        output=alignments_store.process_reference_annotation("UniRef90_W1Q3F0|UniRef50_P59787|5000|bug")
        
        expected_output=["UniRef90_W1Q3F0|UniRef50_P59787|5000|bug",0,"unclassified"]
        
        self.assertEqual(expected_output,output) 

    def test_GeneScores_add(self):
        """
        GeneScores class: Test add function
        """
        
        gene_scores=store.GeneScores()
        
        bug1_scores={"gene1":1,"gene2":2}
        gene_scores.add(bug1_scores,"bug1")
        
        bug2_scores={"gene1":1,"gene2":2}
        gene_scores.add(bug2_scores,"bug2")
        
        self.assertEqual(gene_scores.count_genes_for_bug("bug1"),2)
        
    def test_GeneScores_add_second_set(self):
        """
        GeneScores class: Test add function
        Test adding a second set of scores to bug set
        """
        
        gene_scores=store.GeneScores()
        
        bug1_scores={"gene1":1,"gene2":2}
        gene_scores.add(bug1_scores,"bug1")
        
        bug1_scores_2={"gene3":1,"gene4":2, "gene2":22}
        gene_scores.add(bug1_scores_2,"bug1")
        
        self.assertEqual(gene_scores.count_genes_for_bug("bug1"),4)
        
    def test_GeneScores_get_score(self):
        """
        GeneScores class: Test get_score function
        """
        
        gene_scores=store.GeneScores()
        
        bug1_scores={"gene1":1,"gene2":2}
        gene_scores.add(bug1_scores,"bug1")
        
        bug2_scores={"gene1":1,"gene2":2}
        gene_scores.add(bug2_scores,"bug2")
        
        self.assertEqual(gene_scores.get_score("bug1","gene2"),2)
        
    def test_GeneScores_get_score_second_set(self):
        """
        GeneScores class: Test get_score function
        Test getting the score for a second set of scores added to bug set
        """
        
        gene_scores=store.GeneScores()
        
        bug1_scores={"gene1":1,"gene2":2}
        gene_scores.add(bug1_scores,"bug1")
        
        bug1_scores_2={"gene3":1,"gene4":2, "gene2":22}
        gene_scores.add(bug1_scores_2,"bug1")
        
        # Test that the most recent score for gene2 is returned
        self.assertEqual(gene_scores.get_score("bug1","gene2"),22)
        
    def test_GeneScores_scores_for_bug(self):
        """
        GeneScores class: Test scores_for_bug
        """
  
        gene_scores=store.GeneScores()
        
        bug1_scores={"gene1":1,"gene2":2}
        gene_scores.add(bug1_scores,"bug1")
        
        bug1_scores_2={"gene3":1,"gene4":2, "gene2":22}
        gene_scores.add(bug1_scores_2,"bug1")
        
        # Test that the most recent score for gene2 is returned
        self.assertDictEqual(gene_scores.scores_for_bug("bug1"),
            {"gene1":1,"gene2":22,"gene3":1,"gene4":2})      
    
    def test_GeneScores_add_from_file_bug_list(self):
        """
        GeneScores class: Test add_from_file bug list
        """
        
        gene_scores=store.GeneScores()
        
        gene_scores.add_from_file(cfg.genetable_file)
        
        # Test the bug list is as expected
        self.assertEqual(sorted(cfg.genetable_file_bug_scores.keys()),sorted(gene_scores.bug_list()))
        
    def test_GeneScores_add_from_file_gene_list(self):
        """
        GeneScores class: Test add_from_file gene list
        """
        
        gene_scores=store.GeneScores()
        
        gene_scores.add_from_file(cfg.genetable_file)
        
        # Create a list of all of the genes in the table
        genes={}
        for bug in cfg.genetable_file_bug_scores:
            for gene in cfg.genetable_file_bug_scores[bug]:
                genes[gene]=1
        
        # Test the gene list is as expected
        self.assertEqual(sorted(genes.keys()),sorted(gene_scores.gene_list()))
        
    def test_GeneScores_add_from_file_scores(self):
        """
        GeneScores class: Test add_from_file scores
        """
        
        gene_scores=store.GeneScores()
        
        gene_scores.add_from_file(cfg.genetable_file)
        
        # Test the scores for all bugs and genes
        for bug in cfg.genetable_file_bug_scores:
            self.assertDictEqual(cfg.genetable_file_bug_scores[bug],gene_scores.scores_for_bug(bug))
            
        
