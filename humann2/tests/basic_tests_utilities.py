import unittest
import importlib
import os
import filecmp
import sys
import logging
import gzip
import tempfile
import shutil

import cfg
import utils

from humann2 import utilities
from humann2 import config

class TestHumann2UtilitiesFunctions(unittest.TestCase):
    """
    Test the functions found in humann2.utilities
    """
    
    def setUp(self):
        config.unnamed_temp_dir=tempfile.gettempdir()
        
        # set up nullhandler for logger
        logging.getLogger('humann2.utilities').addHandler(logging.NullHandler())

    def test_file_exists_readable(self):
        """
        Test the file_exists_readable function with a config file
        """
        
        utilities.file_exists_readable(cfg.small_fasta_file)
        
    def test_file_exists_readable_raise(self):
        """
        Test the file_exists_readable function with out a file
        """
        
        # Redirect stdout
        sys.stdout=open(os.devnull,"w")

        with self.assertRaises(IOError):
            utilities.file_exists_readable(
            os.path.join(cfg.data_folder,"not_a_file"),raise_IOError=True)

        # Undo stdout redirect
        sys.stdout=sys.__stdout__

    def test_count_reads_fasta(self):
        """
        Test the count_reads function on a fasta file
        """
        
        fasta_seqs=utilities.count_reads(cfg.small_fasta_file)
        self.assertEqual(fasta_seqs, cfg.small_fasta_file_total_sequences)                
        
    def test_count_reads_fastq(self):
        """
        Test the count_reads function on a fastq file
        """               

        fastq_seqs=utilities.count_reads(cfg.small_fastq_file)
        self.assertEqual(fastq_seqs, cfg.small_fastq_file_total_sequences)      

    def test_estimate_unaligned_reads_identical_files(self):
        """
        Test the estimate_unaligned_reads function on identical files
        """ 

        # Test with identical files
        percent_unaligned=utilities.estimate_unaligned_reads(
            cfg.small_fastq_file, cfg.small_fastq_file)

        percent_expected=(cfg.small_fastq_file_total_sequences / 
            float(cfg.small_fastq_file_total_sequences)*100.0)

        self.assertAlmostEqual(float(percent_unaligned), percent_expected)   
        
    def test_estimate_unaligned_reads_fasta_and_fastq_files(self):
        """
        Test the estimate_unaligned_reads function on fasta/fastq file
        """
        
        # Test with a fasta and a fastq file with a different number of reads
        percent_unaligned=utilities.estimate_unaligned_reads(
            cfg.small_fasta_file, cfg.small_fastq_file)

        percent_expected=(cfg.small_fastq_file_total_sequences / 
            float(cfg.small_fasta_file_total_sequences)*100.0)

        self.assertAlmostEqual(float(percent_unaligned), percent_expected)   

    def test_break_up_fasta_file(self):
        """
        Test the break_up_fasta_file function
        """
        
        # Break up the file into smaller files each containing a single read
        new_fasta_files=utilities.break_up_fasta_file(
            cfg.small_fasta_file,1)

        for file in new_fasta_files:
            sequence_count=utilities.count_reads(file)
            self.assertEqual(sequence_count,1)
            utils.remove_temp_file(file) 

    def test_fastq_to_fasta(self):
        """
        Test the fastq_to_fasta function
        """
        
        new_fasta_file=utilities.fastq_to_fasta(
            cfg.convert_fastq_file)
        self.assertTrue(filecmp.cmp(new_fasta_file,
            cfg.convert_fasta_file, shallow=False))
        utils.remove_temp_file(new_fasta_file)     
        
    def test_fastq_to_fasta(self):
        """
        Test the fastq_to_fasta function with a set of sequences
        which have the @ quality score as the first score
        This tests that the sequence id and sequence are selected correctly
        even though the @ starts a sequence id line and a quality score line
        """
        
        new_fasta_file=utilities.fastq_to_fasta(
            cfg.convert_fastq_at_character_file)
        self.assertTrue(filecmp.cmp(new_fasta_file,
            cfg.convert_fasta_file, shallow=False))
        utils.remove_temp_file(new_fasta_file)  
                       

    def test_double_sort(self):
        """
        Test the double_sort function
        """
        
        pathways={}
        pathways["pathA"]=3
        pathways["pathB"]=3
        pathways["pathC"]=2
        pathways["pathD"]=1
        pathways["pathE"]=1
        
        double_sorted_pathways_keys_check=["pathA","pathB","pathC",
            "pathD","pathE"]
        
        double_sorted_pathways_keys=utilities.double_sort(pathways)
        
        self.assertEqual(double_sorted_pathways_keys,double_sorted_pathways_keys_check)
        
    def test_double_sort_float(self):
        """
        Test the double_sort function with a float value
        """
        
        pathways={}
        pathways["pathA"]=0.1
        pathways["pathB"]=3
        pathways["pathC"]=2
        pathways["pathD"]=1
        pathways["pathE"]=1
        
        double_sorted_pathways_keys=utilities.double_sort(pathways)       
        
        double_sorted_pathways_keys_check=["pathB","pathC",
            "pathD","pathE","pathA"]      
        
        self.assertEqual(double_sorted_pathways_keys, double_sorted_pathways_keys_check)
        
    def test_determine_file_format_fasta(self):
        """
        Test the determine_file_format function with a fasta file
        """
        
        format=utilities.determine_file_format(cfg.small_fasta_file)
        
        self.assertEqual(format,"fasta")
        
    def test_determine_file_format_fastq(self):
        """
        Test the determine_file_format function with a fastq file
        """
        
        format=utilities.determine_file_format(cfg.small_fastq_file)
        
        self.assertEqual(format,"fastq")
        
    def test_determine_file_format_sam_with_header(self):
        """
        Test the determine_file_format function with a sam file with header
        """
        
        format=utilities.determine_file_format(cfg.sam_file_with_header)
        
        self.assertEqual(format,"sam")        
        
    def test_determine_file_format_sam_without_header(self):
        """
        Test the determine_file_format function with a sam file without header
        """
        
        format=utilities.determine_file_format(cfg.sam_file_without_header)
        
        self.assertEqual(format,"sam")    
        
    def test_determine_file_format_bam(self):
        """
        Test the determine_file_format function with a bam file
        """
        
        format=utilities.determine_file_format(cfg.bam_file)
        
        self.assertEqual(format,"bam")   
        
    def test_determine_file_format_sam_without_header_with_tags(self):
        """
        Test the determine_file_format function with a sam file without header
        with tags
        """
        
        format=utilities.determine_file_format(cfg.sam_file_without_header_with_tags)
        
        self.assertEqual(format,"sam") 
        
        
    def test_determine_file_format_rapsearch2_with_header(self):
        """
        Test the determine_file_format function with a rapsearch2 output file
        with the standard header format
        """
        
        format=utilities.determine_file_format(cfg.rapsearch2_output_file_with_header)
        
        self.assertEqual(format,"blastm8") 
        
    def test_determine_file_format_rapsearch2_without_header(self):
        """
        Test the determine_file_format function with a rapsearch2 output file
        without the standard header format
        """
        
        format=utilities.determine_file_format(cfg.rapsearch2_output_file_without_header)
        
        self.assertEqual(format,"blastm8")
        
    def test_determine_file_format_fasta_gzip(self):
        """
        Test the determine_file_format function with a fasta file that is gzipped
        """
        
        # create a small gzipped fasta file
        # read in the small fasta file
        file_handle=open(cfg.small_fasta_file,"rt")
        
        # create a temp file
        file_out, gzip_fasta_file=tempfile.mkstemp(suffix=".gz")
        os.close(file_out)
        
        # write the gzipped file
        file_handle_gzip=gzip.open(gzip_fasta_file,"wt")
        shutil.copyfileobj(file_handle, file_handle_gzip)
        file_handle.close()
        file_handle_gzip.close()
           
        format=utilities.determine_file_format(gzip_fasta_file)
        
        # remove the temp gzipped file
        utils.remove_temp_file(gzip_fasta_file)
        
        self.assertEqual(format,"fasta.gz")

    def test_determine_file_format_fastq_gzip(self):
        """
        Test the determine_file_format function with a fastq file that is gzipped
        """
        
        # create a small gzipped fastq file
        # read in the small fastq file
        file_handle=open(cfg.small_fastq_file,"rt")
        
        # create a temp file
        file_out, gzip_fastq_file=tempfile.mkstemp(suffix=".gz")
        os.close(file_out)
        
        # write the gzipped file
        file_handle_gzip=gzip.open(gzip_fastq_file,"wt")
        shutil.copyfileobj(file_handle, file_handle_gzip)
        file_handle.close()
        file_handle_gzip.close()
           
        format=utilities.determine_file_format(gzip_fastq_file)
        
        # remove the temp gzipped file
        utils.remove_temp_file(gzip_fastq_file)
        
        self.assertEqual(format,"fastq.gz")
        
    def test_determine_file_format_biom(self):
        """
        Test the determine_file_format function with a biom file
        """
        
        format=utilities.determine_file_format(cfg.biom_file)
        
        self.assertEqual(format,"biom")   

    def test_determine_file_format_genetable(self):
        """
        Test the determine_file_format function with a gene table tsv file
        """
        
        format=utilities.determine_file_format(cfg.genetable_file)
        
        self.assertEqual(format,"genetable")   
        
    def test_gunzip_file(self):
        """
        Test the gunzip_file function
        """
        
        # create a small gzipped fastq file
        # read in the small fastq file
        file_handle=open(cfg.small_fastq_file,"rt")
        
        # create a temp file
        file_out, gzip_fastq_file=tempfile.mkstemp(suffix=".gz")
        os.close(file_out)
        
        # write the gzipped file
        file_handle_gzip=gzip.open(gzip_fastq_file,"wt")
        shutil.copyfileobj(file_handle, file_handle_gzip)
        file_handle.close()
        file_handle_gzip.close()
        
        # Redirect stdout
        sys.stdout=open(os.devnull,"w")

        # create the ungzipped file
        new_file=utilities.gunzip_file(gzip_fastq_file)
    
        # Undo stdout redirect
        sys.stdout=sys.__stdout__
    
        # remove the temp gzipped file
        utils.remove_temp_file(gzip_fastq_file)
    
        # check the files are the same
        self.assertTrue(filecmp.cmp(new_file,
            cfg.small_fastq_file, shallow=False))
        
        # remove the temp gunzipped file
        utils.remove_temp_file(new_file)

    def test_add_length_annotation(self):
        """
        Test the add_length_annotation function
        """
        
        id="UniRef50_ABC"
        length=100
        
        expected_result="UniRef50_ABC|100"
        
        self.assertEqual(utilities.add_length_annotation(id, length),expected_result)
        
    def test_add_length_annotation_spaces(self):
        """
        Test the add_length_annotation function with spaces in the id
        """
        
        id="UniRef50_ABC space here"
        length=100
        
        expected_result="UniRef50_ABC|100"
        
        self.assertEqual(utilities.add_length_annotation(id, length),expected_result)
        
    def test_get_length_annotation(self):
        """
        Test the get_length_annotation function
        """
        
        annotated_id="UniRef50_ABC|100"
        expected_result=("UniRef50_ABC",100)
        
        result=utilities.get_length_annotation(annotated_id)
        
        self.assertEqual(result[0],expected_result[0])
        self.assertEqual(result[1],expected_result[1])
        
    def test_get_length_annotation_no_length(self):
        """
        Test the get_length_annotation function from an id without a length
        """
        
        annotated_id="UniRef50_ABC"
        expected_result=("UniRef50_ABC",1)
        
        result=utilities.get_length_annotation(annotated_id)
        
        self.assertEqual(result[0],expected_result[0])
        self.assertEqual(result[1],expected_result[1])
        
    def test_get_length_annotation_non_int_length(self):
        """
        Test the get_length_annotation function from an id without an int length
        """
        
        annotated_id="UniRef50_ABC|A"
        expected_result=("UniRef50_ABC|A",1)
        
        result=utilities.get_length_annotation(annotated_id)
        
        self.assertEqual(result[0],expected_result[0])
        self.assertEqual(result[1],expected_result[1])

    def test_get_length_annotation_multiple_annotations(self):
        """
        Test the get_length_annotation function from an id with multiple annotations
        """
        
        annotated_id="UniRef50_ABC|A|C|100"
        expected_result=("UniRef50_ABC|A|C",100)
        
        result=utilities.get_length_annotation(annotated_id)
        
        self.assertEqual(result[0],expected_result[0])
        self.assertEqual(result[1],expected_result[1])         
    
        
