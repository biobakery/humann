import unittest
import logging
import tempfile
import math

import cfg
import utils

from humann.search import nucleotide
from humann import store
from humann import config
from humann import utilities

class TestAdvancedHumannNucleotideSearchFunctions(unittest.TestCase):
    """
    Test the functions found in humann.search.nucleotide
    """
    
    def setUp(self):
        config.unnamed_temp_dir=tempfile.gettempdir()
        config.temp_dir=tempfile.gettempdir()
        config.file_basename="HUMAnN_test"

        # set default identity threshold
        config.identity_threshold = 50.0
        
        # set up nullhandler for logger
        logging.getLogger('humann.search.nucleotide').addHandler(logging.NullHandler())

        # record default config options (for tests that change these)
        self.default_nucleotide_subject_coverage_threshold = config.nucleotide_subject_coverage_threshold
        self.default_nucleotide_query_coverage_threshold = config.nucleotide_query_coverage_threshold

    def test_nucleotide_search_unaligned_reads_output_fasta_format(self):
        """
        Test the unaligned reads and the store alignments
        Test with a bowtie2/sam output file
        Test output file is of fasta format
        Test sam file is not removed
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()
       
        # turn off query/subject filtering
        config.nucleotide_subject_coverage_threshold = 0
        config.nucleotide_query_coverage_threshold = 0
 
        # read in the aligned and unaligned reads
        [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
            cfg.sam_file_unaligned_reads, alignments, unaligned_reads_store, keep_sam=True) 
        
        # reset query/subject filtering
        config.nucleotide_subject_coverage_threshold = self.default_nucleotide_subject_coverage_threshold
        config.nucleotide_query_coverage_threshold = self.default_nucleotide_query_coverage_threshold

        # check for fasta output file format
        file_format=utilities.determine_file_format(unaligned_reads_file_fasta)
        self.assertEqual("fasta",file_format)
        
        # remove temp files
        utils.remove_temp_file(unaligned_reads_file_fasta)
        utils.remove_temp_file(reduced_aligned_reads_file)


    def test_nucleotide_search_unaligned_reads_read_count_aligned(self):
        """
        Test the unaligned reads and the store alignments
        Test with a bowtie2/sam output file
        Test for aligned read counts
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # turn off query/subject filtering
        config.nucleotide_subject_coverage_threshold = 0
        config.nucleotide_query_coverage_threshold = 0
        
        # read in the aligned and unaligned reads
        [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
            cfg.sam_file_unaligned_reads, alignments, unaligned_reads_store, keep_sam=True) 
        
        # reset query/subject filtering
        config.nucleotide_subject_coverage_threshold = self.default_nucleotide_subject_coverage_threshold
        config.nucleotide_query_coverage_threshold = self.default_nucleotide_query_coverage_threshold
        
        # remove temp files
        utils.remove_temp_file(unaligned_reads_file_fasta)
        utils.remove_temp_file(reduced_aligned_reads_file)
        
        # check the aligned reads count
        self.assertEqual(len(alignments.get_hit_list()),cfg.sam_file_unaligned_reads_total_aligned)
        
    def test_nucleotide_search_unaligned_reads_read_count_aligned_subject_coverage(self):
        """
        Test the unaligned reads and the store alignments
        Test with a bowtie2/sam output file
        Test for aligned read counts
        Test with subject coverage filtering
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # turn off subject filtering
        config.nucleotide_query_coverage_threshold = 0
        
        # read in the aligned and unaligned reads
        [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
            cfg.sam_file_unaligned_reads, alignments, unaligned_reads_store, keep_sam=True) 
        
        # reset subject filtering
        config.nucleotide_query_coverage_threshold = self.default_nucleotide_query_coverage_threshold
        
        # remove temp files
        utils.remove_temp_file(unaligned_reads_file_fasta)
        utils.remove_temp_file(reduced_aligned_reads_file)
        
        # check the aligned reads count
        self.assertEqual(len(alignments.get_hit_list()),cfg.sam_file_unaligned_reads_total_aligned_subject_coverage)
        
    def test_nucleotide_search_unaligned_reads_subject_coverage_repeated_gene_family(self):
        """
        Test the unaligned reads and the store alignments
        Test with a bowtie2/sam output file with the same gene family in two species
        Test that subject coverage is computed with the length of each reference sequence
        so the alignment to the poorly covered species is filtered
        """

        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()

        # turn off query filtering and require half of the reference to be covered
        config.nucleotide_query_coverage_threshold = 0
        config.nucleotide_subject_coverage_threshold = 50.0

        # read in the aligned and unaligned reads
        [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
            cfg.sam_file_two_species_same_gene_family, alignments, unaligned_reads_store, keep_sam=True)

        # reset the filtering thresholds
        config.nucleotide_query_coverage_threshold = self.default_nucleotide_query_coverage_threshold
        config.nucleotide_subject_coverage_threshold = self.default_nucleotide_subject_coverage_threshold

        # remove temp files
        utils.remove_temp_file(unaligned_reads_file_fasta)
        utils.remove_temp_file(reduced_aligned_reads_file)

        # only the alignment to the sequence that is covered should be stored
        hits=alignments.get_hit_list()
        self.assertEqual(len(hits),1)
        [query, bug, reference, score, length]=hits[0]
        self.assertEqual(bug,cfg.sam_file_two_species_same_gene_family_covered_bug)

        # the alignment to the other copy of the gene family is treated as unaligned
        self.assertEqual(unaligned_reads_store.count_reads(),1)

    def test_nucleotide_search_unaligned_reads_multiple_alignments_per_read(self):
        """
        Test the unaligned reads and the store alignments
        Test with a sam file where reads have more than one alignment
        Test a read with a single alignment kept is not also reported as unaligned
        Test a read with all of its alignments filtered is only written out once
        """

        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()

        # turn off query filtering and require half of the reference to be covered
        config.nucleotide_query_coverage_threshold = 0
        config.nucleotide_subject_coverage_threshold = 50.0

        # read in the aligned and unaligned reads
        [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
            cfg.sam_file_multiple_alignments_per_read, alignments, unaligned_reads_store, keep_sam=True)

        # reset the filtering thresholds
        config.nucleotide_query_coverage_threshold = self.default_nucleotide_query_coverage_threshold
        config.nucleotide_subject_coverage_threshold = self.default_nucleotide_subject_coverage_threshold

        # count the reads written to the unaligned reads file
        unaligned_ids=[line.rstrip()[1:] for line in open(unaligned_reads_file_fasta)
            if line.startswith(">")]

        # remove temp files
        utils.remove_temp_file(unaligned_reads_file_fasta)
        utils.remove_temp_file(reduced_aligned_reads_file)

        # the read with one alignment to a well covered reference sequence is aligned
        hits=alignments.get_hit_list()
        self.assertEqual(len(hits),cfg.sam_file_multiple_alignments_per_read_total_aligned)

        # the total number of reads is counted once per read and not once per alignment
        self.assertEqual(unaligned_reads_store.get_initial_read_count(),
            cfg.sam_file_multiple_alignments_per_read_total_reads)

        # the aligned read is not also counted as unaligned
        self.assertEqual(unaligned_reads_store.count_reads(),
            cfg.sam_file_multiple_alignments_per_read_total_unaligned)

        # each unaligned read is written to the unaligned reads file exactly once, so
        # the file matches the number of unaligned reads stored
        self.assertEqual(len(unaligned_ids),
            cfg.sam_file_multiple_alignments_per_read_total_unaligned)
        self.assertEqual(len(set(unaligned_ids)),len(unaligned_ids))

        # the read that is aligned is not included in the unaligned reads
        self.assertEqual(sorted(unaligned_reads_store.id_list()),["r2","r3"])

    def test_nucleotide_search_unaligned_reads_duplicate_read_names(self):
        """
        Test the unaligned reads and the store alignments
        Test with a sam file where two reads share a name
        Test a warning is written since the read counts reported will be too small
        """

        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()

        config.nucleotide_query_coverage_threshold = 0
        config.nucleotide_subject_coverage_threshold = 50.0

        with self.assertLogs("humann.search.nucleotide", level="WARNING") as logged:
            [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
                cfg.sam_file_duplicate_read_names, alignments, unaligned_reads_store, keep_sam=True)

        config.nucleotide_query_coverage_threshold = self.default_nucleotide_query_coverage_threshold
        config.nucleotide_subject_coverage_threshold = self.default_nucleotide_subject_coverage_threshold

        utils.remove_temp_file(unaligned_reads_file_fasta)
        utils.remove_temp_file(reduced_aligned_reads_file)

        # the two reads that share a name are tracked as a single read
        self.assertEqual(unaligned_reads_store.get_initial_read_count(),1)

        # the warning names the number of reads that could not be tracked
        self.assertTrue(any("1 read(s) with a name that is not unique" in message
            for message in logged.output))

    def test_nucleotide_search_unaligned_reads_read_count_aligned_query_coverage(self):
        """
        Test the unaligned reads and the store alignments
        Test with a bowtie2/sam output file
        Test for aligned read counts
        Test with query coverage filtering
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # turn off query filtering
        config.nucleotide_subject_coverage_threshold = 0
        
        # read in the aligned and unaligned reads
        [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
            cfg.sam_file_unaligned_reads, alignments, unaligned_reads_store, keep_sam=True) 
        
        # reset query filtering
        config.nucleotide_query_coverage_threshold = self.default_nucleotide_query_coverage_threshold
        
        # remove temp files
        utils.remove_temp_file(unaligned_reads_file_fasta)
        utils.remove_temp_file(reduced_aligned_reads_file)
        
        # check the aligned reads count
        self.assertEqual(len(alignments.get_hit_list()),cfg.sam_file_unaligned_reads_total_aligned_query_coverage)
        
    def test_nucleotide_search_unaligned_reads_read_count_aligned_evalue_threshold(self):
        """
        Test the unaligned reads and the store alignments
        Test with a bowtie2/sam output file
        Test for aligned read counts
        Test the evalue threshold does not filter alignments
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # turn off query/subject filtering
        config.nucleotide_subject_coverage_threshold = 0
        config.nucleotide_query_coverage_threshold = 0

        # update the evalue threshold to a number less than those for the alignment file
        original_evalue_threshold=config.evalue_threshold
        config.evalue_threshold=1e-15
        
        # read in the aligned and unaligned reads
        [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
            cfg.sam_file_unaligned_reads, alignments, unaligned_reads_store, keep_sam=True) 
        
        # reset query/subject filtering
        config.nucleotide_subject_coverage_threshold = self.default_nucleotide_subject_coverage_threshold
        config.nucleotide_query_coverage_threshold = self.default_nucleotide_query_coverage_threshold
        
        # remove temp files
        utils.remove_temp_file(unaligned_reads_file_fasta)
        utils.remove_temp_file(reduced_aligned_reads_file)
        
        # reset the evalue threshold back to the original
        config.evalue_threshold=original_evalue_threshold
        
        # check the aligned reads count (all reads should be aligned even though they do not
        # meet the threshold as the evalue threshold is not applied for this type of alignment)
        self.assertEqual(len(alignments.get_hit_list()),cfg.sam_file_unaligned_reads_total_aligned)
        
    def test_nucleotide_search_unaligned_reads_read_count_aligned_identity_threshold(self):
        """
        Test the unaligned reads and the store alignments
        Test with a bowtie2/sam output file
        Test for aligned read counts
        Test the identity threshold does filter alignments
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # turn off query/subject filtering
        config.nucleotide_subject_coverage_threshold = 0
        config.nucleotide_query_coverage_threshold = 0
        
        # update the identity threshold to a number larger than those in the alignments
        original_identity_threshold=config.identity_threshold
        config.identity_threshold=101.0
        
        # read in the aligned and unaligned reads
        [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
            cfg.sam_file_unaligned_reads, alignments, unaligned_reads_store, keep_sam=True) 
        
        # reset query/subject filtering
        config.nucleotide_subject_coverage_threshold = self.default_nucleotide_subject_coverage_threshold
        config.nucleotide_query_coverage_threshold = self.default_nucleotide_query_coverage_threshold
        
        # remove temp files
        utils.remove_temp_file(unaligned_reads_file_fasta)
        utils.remove_temp_file(reduced_aligned_reads_file)
        
        # reset the identity threshold back to the original
        config.identity_threshold=original_identity_threshold
        
        # check the aligned reads count (it should be two as both should pass the threshold)
        self.assertEqual(len(alignments.get_hit_list()),2)

    def test_nucleotide_search_unaligned_reads_read_count_unaligned(self):
        """
        Test the unaligned reads and the store alignments
        Test with a bowtie2/sam output file
        Test for unaligned read counts
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # turn off query/subject filtering
        config.nucleotide_subject_coverage_threshold = 0
        config.nucleotide_query_coverage_threshold = 0
        
        # read in the aligned and unaligned reads
        [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
            cfg.sam_file_unaligned_reads, alignments, unaligned_reads_store, keep_sam=True) 
        
        # reset query/subject filtering
        config.nucleotide_subject_coverage_threshold = self.default_nucleotide_subject_coverage_threshold
        config.nucleotide_query_coverage_threshold = self.default_nucleotide_query_coverage_threshold
        
        # remove temp files
        utils.remove_temp_file(unaligned_reads_file_fasta)
        utils.remove_temp_file(reduced_aligned_reads_file)
        
        # check the unaligned reads count
        self.assertEqual(unaligned_reads_store.count_reads(),cfg.sam_file_unaligned_reads_total_unaligned)
        
    def test_nucleotide_search_unaligned_reads_read_count_unaligned_minimize_memory_use(self):
        """
        Test the unaligned reads and the store alignments
        Test with a bowtie2/sam output file
        Test for unaligned read counts
        Test with minimize memory use
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads(minimize_memory_use=True)
        
        # turn off query/subject filtering
        config.nucleotide_subject_coverage_threshold = 0
        config.nucleotide_query_coverage_threshold = 0
        
        # read in the aligned and unaligned reads
        [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
            cfg.sam_file_unaligned_reads, alignments, unaligned_reads_store, keep_sam=True) 
        
        # reset query/subject filtering
        config.nucleotide_subject_coverage_threshold = self.default_nucleotide_subject_coverage_threshold
        config.nucleotide_query_coverage_threshold = self.default_nucleotide_query_coverage_threshold
        
        # remove temp files
        utils.remove_temp_file(unaligned_reads_file_fasta)
        utils.remove_temp_file(reduced_aligned_reads_file)
        
        # check the unaligned reads count
        self.assertEqual(unaligned_reads_store.count_reads(),cfg.sam_file_unaligned_reads_total_unaligned)
        
    def test_nucleotide_search_unaligned_reads_annotations_reference(self):
        """
        Test the unaligned reads and the store alignments
        Test with a bowtie2/sam output file
        Test the different annotation formats are recognized for reference
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # turn off query/subject filtering
        config.nucleotide_subject_coverage_threshold = 0
        config.nucleotide_query_coverage_threshold = 0
        
        # read in the aligned and unaligned reads
        [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
            cfg.sam_file_annotations, alignments, unaligned_reads_store, keep_sam=True) 
        
        # reset query/subject filtering
        config.nucleotide_subject_coverage_threshold = self.default_nucleotide_subject_coverage_threshold
        config.nucleotide_query_coverage_threshold = self.default_nucleotide_query_coverage_threshold
        
        # remove temp files
        utils.remove_temp_file(unaligned_reads_file_fasta)
        utils.remove_temp_file(reduced_aligned_reads_file)
        
        # two of the hits should be for gene "UniRef50"
        hits=alignments.hits_for_gene("UniRef50")
        self.assertEqual(len(hits),2)
        
                
    def test_nucleotide_search_unaligned_reads_annotations_bug(self):
        """
        Test the unaligned reads and the store alignments
        Test with a bowtie2/sam output file
        Test the different annotation formats are recognized for bug
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # turn off query/subject filtering
        config.nucleotide_subject_coverage_threshold = 0
        config.nucleotide_query_coverage_threshold = 0
        
        # read in the aligned and unaligned reads
        [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
            cfg.sam_file_annotations, alignments, unaligned_reads_store, keep_sam=True) 
        
        # reset query/subject filtering
        config.nucleotide_subject_coverage_threshold = self.default_nucleotide_subject_coverage_threshold
        config.nucleotide_query_coverage_threshold = self.default_nucleotide_query_coverage_threshold
        
        # remove temp files
        utils.remove_temp_file(unaligned_reads_file_fasta)
        utils.remove_temp_file(reduced_aligned_reads_file)
        
        # there should be one bug which is unclassified
        self.assertEqual(alignments.bug_list(),["unclassified"])
        
    def test_nucleotide_search_unaligned_reads_annotations_gene_length(self):
        """
        Test the unaligned reads and the store alignments
        Test with a bowtie2/sam output file
        Test the different annotation formats are recognized for gene length
        Test the gene length uses the read length from the sam file
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # turn off query/subject filtering
        config.nucleotide_subject_coverage_threshold = 0
        config.nucleotide_query_coverage_threshold = 0
        
        # read in the aligned and unaligned reads
        [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
            cfg.sam_file_annotations, alignments, unaligned_reads_store, keep_sam=True) 
        
        # reset query/subject filtering
        config.nucleotide_subject_coverage_threshold = self.default_nucleotide_subject_coverage_threshold
        config.nucleotide_query_coverage_threshold = self.default_nucleotide_query_coverage_threshold
        
        # remove temp files
        utils.remove_temp_file(unaligned_reads_file_fasta)
        utils.remove_temp_file(reduced_aligned_reads_file)
        
        # there should be 4 hits identified
        all_hits=alignments.get_hit_list()
        self.assertEqual(len(all_hits),4)
        
        # check for set and default gene lengths
        read_length = 151
        expected_length_uniref50 = (abs(2000 - read_length)+1)/1000.0
        expected_length_other = (abs(1000 - read_length)+1)/1000.0
        
        for hit in all_hits:
            query, bug, reference, score, length = hit
            if reference == "UniRef50":
                self.assertEqual(length,expected_length_uniref50)
            else:
                self.assertEqual(length,expected_length_other)
                
    def test_nucleotide_search_unaligned_reads_scores(self):
        """
        Test the unaligned reads and the store alignments
        Test with a bowtie2/sam output file
        Test the scores are based on percent identities
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # turn off query/subject filtering
        config.nucleotide_subject_coverage_threshold = 0
        config.nucleotide_query_coverage_threshold = 0
        
        # read in the aligned and unaligned reads
        [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
            cfg.sam_file_annotations, alignments, unaligned_reads_store, keep_sam=True) 
        
        # reset query/subject filtering
        config.nucleotide_subject_coverage_threshold = self.default_nucleotide_subject_coverage_threshold
        config.nucleotide_query_coverage_threshold = self.default_nucleotide_query_coverage_threshold
        
        # remove temp files
        utils.remove_temp_file(unaligned_reads_file_fasta)
        utils.remove_temp_file(reduced_aligned_reads_file)
        
        # there should be 4 hits identified
        all_hits=alignments.get_hit_list()
        
        # check for set and default gene lengths
        expected_score=math.pow(151.0, config.match_power)
        
        for hit in all_hits:
            query, bug, reference, score, length = hit
            self.assertEqual(score,expected_score)

    def test_nucleotide_search_unaligned_reads_output_blast_format(self):
        """
        Test the unaligned reads and the store alignments
        Test with a bowtie2/sam output file
        Test the aligned reads file created is of the blastm8 format
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        unaligned_reads_store=store.Reads()
        
        # turn off query/subject filtering
        config.nucleotide_subject_coverage_threshold = 0
        config.nucleotide_query_coverage_threshold = 0
        
        config.file_basename="TEST"
        
        # read in the aligned and unaligned reads
        [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
            cfg.sam_file_annotations, alignments, unaligned_reads_store, keep_sam=True) 
        
        # reset query/subject filtering
        config.nucleotide_subject_coverage_threshold = self.default_nucleotide_subject_coverage_threshold
        config.nucleotide_query_coverage_threshold = self.default_nucleotide_query_coverage_threshold
        
        # test file is of the blastm8 format
        file_format=utilities.determine_file_format(reduced_aligned_reads_file)
        
        # remove temp files
        utils.remove_temp_file(unaligned_reads_file_fasta)
        utils.remove_temp_file(reduced_aligned_reads_file)           
        
        self.assertEqual(file_format,"blastm8")
        
        
        
        
        


