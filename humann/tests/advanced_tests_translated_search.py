import unittest
import logging
import re
import math

import cfg
import utils

from humann.search import translated
from humann.search import blastx_coverage
from humann import store
from humann import config
from humann import utilities

class TestAdvancedHumannTranslatedSearchFunctions(unittest.TestCase):
    """
    Test the functions found in humann.search.translated
    """
    
    def setUp(self):
        config.unnamed_temp_dir="/tmp/"

        # set default identity threshold
        config.identity_threshold = 50.0
        
        # set up nullhandler for logger
        logging.getLogger('humann.search.translated').addHandler(logging.NullHandler())
        logging.getLogger('humann.search.blastx_coverage').addHandler(logging.NullHandler())

    def test_translated_search_unaligned_reads_blastm8(self):
        """
        Test the unaligned reads and the store alignments
        Test with a blastm8-like output file
        Test with empty reads structure
        Test that function does not require gene lengths in reference id
        Test without the coverage filter
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        
        # set the coverage threshold to zero so as to not test with filter on
        current_coverage_threshold=config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold=0
        
        # load the blastm8-like output
        file_handle=open(cfg.rapsearch2_output_file_without_header)
        
        for line in file_handle:
            if not re.search("^#",line):
                data=line.strip().split(config.blast_delimiter)
                
                referenceid=data[config.blast_reference_index]
                queryid=data[config.blast_query_index]
                identity=float(data[config.blast_identity_index])
                alignment_length=float(data[config.blast_aligned_length_index])
            
                alignments.add(referenceid, 0, queryid, identity/100.0*alignment_length,"unclassified",alignment_length)
            
        file_handle.close()
        
        alignments_test=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # load the blastm8-like output with the unaligned reads function
        unaligned_file_fasta=translated.unaligned_reads(unaligned_reads_store, 
            cfg.rapsearch2_output_file_without_header, alignments_test)
        
        # remove temp file
        utils.remove_temp_file(unaligned_file_fasta)
        
        # reset the coverage threshold
        config.translated_subject_coverage_threshold=current_coverage_threshold
        
        # check the values are unchanged
        self.assertEqual(sorted(alignments.get_hit_list()), sorted(alignments_test.get_hit_list()))
        
    def test_translated_search_unaligned_reads_blastm8_coverage_filter(self):
        """
        Test the unaligned reads and the store alignments
        Test with a blastm8-like output file
        Test with empty reads structure
        Test that function does not require gene lengths in reference id
        Test with the coverage filter
        Test with query length annotations
        Test that an alignment with query start larger than query end is not filtered
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        
        # set the coverage threshold to a small value so as to have some alignments pass
        current_coverage_threshold=config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold=0.50
        
        # get the set of allowed proteins
        allowed_proteins = blastx_coverage.blastx_coverage(cfg.rapsearch2_output_file_without_header_coverage,
            config.translated_subject_coverage_threshold, alignments, True)
        
        # load the blastm8-like output
        file_handle=open(cfg.rapsearch2_output_file_without_header_coverage)
        
        for line in file_handle:
            if not re.search("^#",line):
                data=line.strip().split(config.blast_delimiter)
                
                referenceid=data[config.blast_reference_index]
                gene, length, bug = alignments.process_reference_annotation(referenceid)
                queryid, query_length=utilities.get_length_annotation(data[config.blast_query_index])
                identity=float(data[config.blast_identity_index])
                alignment_length=float(data[config.blast_aligned_length_index])
            
                if gene in allowed_proteins:
                    alignments.add(gene, length, queryid, identity/100.0*alignment_length,bug,alignment_length)
            
        file_handle.close()
        
        alignments_test=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # load the blastm8-like output with the unaligned reads function
        unaligned_file_fasta=translated.unaligned_reads(unaligned_reads_store, 
            cfg.rapsearch2_output_file_without_header_coverage, alignments_test)
        
        # remove temp file
        utils.remove_temp_file(unaligned_file_fasta)
        
        # reset the coverage threshold
        config.translated_subject_coverage_threshold=current_coverage_threshold
        
        # check the values are unchanged
        self.assertEqual(sorted(alignments.get_hit_list()), sorted(alignments_test.get_hit_list()))
        
    def test_translated_search_unaligned_reads_rapsearch_log_evalue_threshold(self):
        """
        Test the unaligned reads function
        Test with a rapsearch output file
        Test that log of evalue is taken
        Test the evalue threshold for filtering alignments
        Test without the coverage filter
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        
        # set the coverage threshold to zero so as to not test with filter on
        current_coverage_threshold=config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold=0
        
        # load the rapsearch output
        file_handle=open(cfg.rapsearch2_output_file_with_header)
        
        original_evalue_threshold=config.evalue_threshold
        
        # set a new threshold that will select 3 of the 5 alignments if the log
        # is identified to having been applied
        config.evalue_threshold=1.5e-07
        
        for line in file_handle:
            if not re.search("^#",line):
                data=line.strip().split(config.blast_delimiter)
                
                referenceid=data[config.blast_reference_index]
                queryid=data[config.blast_query_index]
                evalue=float(data[config.blast_evalue_index])
                identity=float(data[config.blast_identity_index])
                alignment_length=float(data[config.blast_aligned_length_index])
                
                # only store those alignments with log evalues that meet threshold
                if math.pow(10.0, evalue) < config.evalue_threshold:
                    alignments.add(referenceid, 0, queryid, identity/100.0*alignment_length,"unclassified",alignment_length)
            
        file_handle.close()
        
        alignments_test=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # load the rapsearch output with the unaligned reads function
        unaligned_file_fasta=translated.unaligned_reads(unaligned_reads_store, 
            cfg.rapsearch2_output_file_with_header, alignments_test)
        
        # remove temp file
        utils.remove_temp_file(unaligned_file_fasta)
        
        # reset the coverage threshold
        config.translated_subject_coverage_threshold=current_coverage_threshold
        
        # set the threshold back to the default
        config.evalue_threshold=original_evalue_threshold
        
        # check the total number of alignments is the same
        self.assertEqual(len(alignments.get_hit_list()),len(alignments_test.get_hit_list()))
        
    def test_translated_search_unaligned_reads_identity_threshold(self):
        """
        Test the unaligned reads function
        Test with a rapsearch output file
        Test the identity threshold filtering
        Test without the coverage filter
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        
        # set the coverage threshold to zero so as to not test with filter on
        current_coverage_threshold=config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold=0
        
        # load the rapsearch output
        file_handle=open(cfg.rapsearch2_output_file_with_header)
        
        original_identity_threshold=config.identity_threshold
        
        # set a new threshold that will select 3 of the 5 alignments
        config.identity_threshold=60.0
        
        for line in file_handle:
            if not re.search("^#",line):
                data=line.strip().split(config.blast_delimiter)
                
                referenceid=data[config.blast_reference_index]
                queryid=data[config.blast_query_index]
                identity=float(data[config.blast_identity_index])
                alignment_length=float(data[config.blast_aligned_length_index])
                
                # only store those alignments with identities that meet threshold
                if identity > config.identity_threshold:
                    alignments.add(referenceid, 0, queryid, identity/100.0*alignment_length,"unclassified",alignment_length)
            
        file_handle.close()
        
        alignments_test=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # load the rapsearch output with the unaligned reads function
        unaligned_file_fasta=translated.unaligned_reads(unaligned_reads_store, 
            cfg.rapsearch2_output_file_with_header, alignments_test)
        
        # remove temp file
        utils.remove_temp_file(unaligned_file_fasta)
        
        # reset the coverage threshold
        config.translated_subject_coverage_threshold=current_coverage_threshold
        
        # set the threshold back to the default
        config.identity_threshold=original_identity_threshold
        
        # check the total number of alignments is the same
        self.assertEqual(len(alignments.get_hit_list()),len(alignments_test.get_hit_list()))

    def test_translated_search_unaligned_reads_rapsearch2_no_log(self):
        """
        Test the unaligned reads function
        Test with a rapsearch2 output file with out the log
        Test that log of evalue is not taken
        Test without the coverage filter
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        
        # set the coverage threshold to zero so as to not test with filter on
        current_coverage_threshold=config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold=0
        
        # load the rapsearch2 output
        file_handle=open(cfg.rapsearch2_output_file_with_header_no_log)
        
        for line in file_handle:
            if not re.search("^#",line):
                data=line.strip().split(config.blast_delimiter)
                
                referenceid=data[config.blast_reference_index]
                queryid=data[config.blast_query_index]
                identity=float(data[config.blast_identity_index])
                alignment_length=float(data[config.blast_aligned_length_index])
            
                alignments.add(referenceid, 0, queryid, identity/100.0*alignment_length,"unclassified",alignment_length)
            
        file_handle.close()
        
        alignments_test=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # load the rapsearch2 output with the unaligned reads function
        unaligned_file_fasta=translated.unaligned_reads(unaligned_reads_store, 
            cfg.rapsearch2_output_file_with_header_no_log, alignments_test)
        
        # remove temp file
        utils.remove_temp_file(unaligned_file_fasta)
        
        # reset the coverage threshold
        config.translated_subject_coverage_threshold=current_coverage_threshold
        
        # check the evalues are unchanged
        self.assertEqual(sorted(alignments.get_hit_list()), sorted(alignments_test.get_hit_list()))
        
        
    def test_translated_search_unaligned_reads_annotations_reference(self):
        """
        Test the unaligned reads and the store alignments
        Test with a rapsearch2 output file
        Test the different annotation formats are recognized for reference
        Test without the coverage filter
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # set the coverage threshold to zero so as to not test with filter on
        current_coverage_threshold=config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold=0
        
        # load the rapsearch2 output with the unaligned reads function
        unaligned_file_fasta=translated.unaligned_reads(unaligned_reads_store, 
            cfg.rapsearch_file_annotations, alignments)
        
        # remove temp file
        utils.remove_temp_file(unaligned_file_fasta)
        
        # reset the coverage threshold
        config.translated_subject_coverage_threshold=current_coverage_threshold
        
        # three of the hits should be for gene "UniRef50"
        hits=alignments.hits_for_gene("UniRef50")
        self.assertEqual(len(hits),3)
        
                
    def test_translated_search_unaligned_reads_annotations_bug(self):
        """
        Test the unaligned reads and the store alignments
        Test with a rapsearch2 output file
        Test the different annotation formats are recognized for bug
        Test without the coverage filter
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # set the coverage threshold to zero so as to not test with filter on
        current_coverage_threshold=config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold=0
        
        # load the rapsearch2 output with the unaligned reads function
        unaligned_file_fasta=translated.unaligned_reads(unaligned_reads_store, 
            cfg.rapsearch_file_annotations, alignments)
        
        # remove temp file
        utils.remove_temp_file(unaligned_file_fasta)
        
        # reset the coverage threshold
        config.translated_subject_coverage_threshold=current_coverage_threshold
        
        # there should be one bug name and the other should be unclassified
        self.assertEqual(sorted(alignments.bug_list()),sorted(["g__Bacteroides.s__Bacteroides_xylanisolvens","unclassified"]))
        
                
    def test_translated_search_unaligned_reads_annotations_gene_length(self):
        """
        Test the unaligned reads and the store alignments
        Test with a rapsearch2 output file
        Test the different annotation formats are recognized for gene length
        Test without the coverage filter
        """
 
         # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # set the coverage threshold to zero so as to not test with filter on
        current_coverage_threshold=config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold=0
        
        # load the rapsearch2 output with the unaligned reads function
        unaligned_file_fasta=translated.unaligned_reads(unaligned_reads_store, 
            cfg.rapsearch_file_annotations, alignments)
        
        # remove temp file
        utils.remove_temp_file(unaligned_file_fasta)   
        
        # reset the coverage threshold
        config.translated_subject_coverage_threshold=current_coverage_threshold    

        # there should be 4 hits identified
        all_hits=alignments.get_hit_list()
        self.assertEqual(len(all_hits),4)
        
        # check for set and default gene lengths
        read_length = 50
        expected_length_uniref50 = (abs(2000 - read_length)+1)/1000.0
        expected_length_other = (abs(1000 - read_length)+1)/1000.0
        
        # check for set and default gene lengths
        for hit in all_hits:
            query, bug, reference, score, length = hit
            if reference == "UniRef50":
                self.assertEqual(length,expected_length_uniref50)
            else:
                self.assertEqual(length,expected_length_other)


