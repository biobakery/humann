import unittest
import re
import tempfile
import os
import filecmp
import logging

import cfg
import utils

from humann2.quantify import families
from humann2 import store
from humann2 import config

class TestAdvancedHumann2QuantifyFamiliesFunctions(unittest.TestCase):
    """
    Test the functions found in humann2.quantify.families
    """
    
    def setUp(self):
        # set up nullhandler for logger
        logging.getLogger('humann2.quantify.gene_families').addHandler(logging.NullHandler())
        logging.getLogger('humann2.store').addHandler(logging.NullHandler())

    def test_gene_families_bug_list(self):
        """
        Test the gene families function and the blast config indexes
        Test UniRef50_unknown is read in and used for gene scores but not printed
        Test the bug list
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        
        # load the usearch output
        file_handle=open(cfg.usearch_file)
        
        for line in file_handle:
            if not re.search("^#",line):
                data=line.strip().split(config.blast_delimiter)
                
                referenceids=data[config.blast_reference_index].split("|")
                queryid=data[config.blast_query_index]
                identity=float(data[config.blast_identity_index])
            
                alignments.add(referenceids[1], 1, queryid, identity,referenceids[0])
            
        file_handle.close()
        
        # check the bugs were loaded correctly
        self.assertEqual(sorted(cfg.usearch_file_bug_list),sorted(alignments.bug_list()))

    def test_gene_families_gene_list(self):
        """
        Test the gene families function and the blast config indexes
        Test UniRef50_unknown is read in and used for gene scores but not printed
        Test the gene list
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        
        # load the usearch output
        file_handle=open(cfg.usearch_file)
        
        for line in file_handle:
            if not re.search("^#",line):
                data=line.strip().split(config.blast_delimiter)
                
                referenceids=data[config.blast_reference_index].split("|")
                queryid=data[config.blast_query_index]
                identity=float(data[config.blast_identity_index])
            
                alignments.add(referenceids[1], 1, queryid, identity,referenceids[0])
            
        file_handle.close()
        
        # check the genes were loaded correctly
        self.assertEqual(sorted(cfg.usearch_file_gene_list),sorted(alignments.gene_list()))
        
    def test_gene_families_tsv_output(self):
        """
        Test the gene families function and the blast config indexes
        Test UniRef50_unknown is read in and used for gene scores but not printed
        Test the tsv output
        """
        
        # update the max decimals to allow for rounding
        config.output_max_decimals=7
        
        # do not use gene name mapping
        original_gene_family_mapping_file=config.gene_family_name_mapping_file
        config.gene_family_name_mapping_file=None
        
        # create a set of alignments
        alignments=store.Alignments()
        
        # load the usearch output
        file_handle=open(cfg.usearch_file)
        
        for line in file_handle:
            if not re.search("^#",line):
                data=line.strip().split(config.blast_delimiter)
                
                referenceids=data[config.blast_reference_index].split("|")
                queryid=data[config.blast_query_index]
                identity=float(data[config.blast_identity_index])
            
                alignments.add(referenceids[1], 1, queryid, identity,referenceids[0])
            
        file_handle.close()
        
        # set the output format
        config.output_format="tsv"
        
        # set the location of the file to write to as a temp file
        file_out, gene_families_file=tempfile.mkstemp()
        os.close(file_out)
        config.genefamilies_file=gene_families_file
        
        # create gene_scores instance
        gene_scores=store.GeneScores()
        
        # obtain the gene families
        gene_families_file=families.gene_families(alignments,gene_scores,0)
        
        # check the gene families output is as expected
        self.assertTrue(filecmp.cmp(gene_families_file,
            cfg.gene_familes_file, shallow=False))
        
        # reset the mapping file
        config.gene_family_name_mapping_file=original_gene_family_mapping_file
        
        # delete the temp file
        utils.remove_temp_file(gene_families_file)
        
    def test_gene_families_tsv_output_with_names(self):
        """
        Test the gene families function and the blast config indexes
        Test UniRef50_unknown is read in and used for gene scores but not printed
        Test the tsv output
        Test that gene families have names applied to them
        Test unmapped reads total is written with the same precision as other lines
        """
        
        # update the max decimals to allow for rounding
        config.output_max_decimals=7
        
        # set to a smaller mapping file
        original_gene_family_mapping_file=config.gene_family_name_mapping_file
        config.gene_family_name_mapping_file=cfg.gene_families_to_names_file
        
        # create a set of alignments
        alignments=store.Alignments()
        
        # load the usearch output
        file_handle=open(cfg.usearch_uniref50_file)
        
        for line in file_handle:
            if not re.search("^#",line):
                data=line.strip().split(config.blast_delimiter)
                
                referenceids=data[config.blast_reference_index].split("|")
                queryid=data[config.blast_query_index]
                identity=float(data[config.blast_identity_index])
            
                alignments.add(referenceids[1], 1, queryid, identity,referenceids[0])
            
        file_handle.close()
        
        # set the output format
        config.output_format="tsv"
        
        # set the location of the file to write to as a temp file
        file_out, gene_families_file=tempfile.mkstemp()
        os.close(file_out)
        config.genefamilies_file=gene_families_file
        
        # create gene_scores instance
        gene_scores=store.GeneScores()
        
        # obtain the gene families
        gene_families_file=families.gene_families(alignments,gene_scores,1)
        
        # check the gene families output is as expected
        self.assertTrue(filecmp.cmp(gene_families_file,
            cfg.gene_familes_uniref50_with_names_file, shallow=False))
        
        # reset the mapping file
        config.gene_family_name_mapping_file=original_gene_family_mapping_file
        
        # delete the temp file
        utils.remove_temp_file(gene_families_file)
