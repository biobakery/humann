import unittest
import re
import tempfile
import os
import filecmp

import cfg
import utils

from humann2.quantify import families
from humann2 import store
from humann2 import config

class TestAdvancedHumann2QuantifyFamiliesFunctions(unittest.TestCase):
    """
    Test the functions found in humann2.quantify.families
    """

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
                evalue=float(data[config.blast_evalue_index])
            
                alignments.add(referenceids[1], 1, queryid, evalue,referenceids[0])
            
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
                evalue=float(data[config.blast_evalue_index])
            
                alignments.add(referenceids[1], 1, queryid, evalue,referenceids[0])
            
        file_handle.close()
        
        # check the genes were loaded correctly
        self.assertEqual(sorted(cfg.usearch_file_gene_list),sorted(alignments.gene_list()))
        
    def test_gene_families_tsv_output(self):
        """
        Test the gene families function and the blast config indexes
        Test UniRef50_unknown is read in and used for gene scores but not printed
        Test the tsv output
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
                evalue=float(data[config.blast_evalue_index])
            
                alignments.add(referenceids[1], 1, queryid, evalue,referenceids[0])
            
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
        gene_families_file=families.gene_families(alignments,gene_scores)
        
        # check the gene families output is as expected
        self.assertTrue(filecmp.cmp(gene_families_file,
            cfg.gene_familes_file, shallow=False))
        
        # delete the temp file
        utils.remove_temp_file(gene_families_file)
