import unittest
import subprocess
import tempfile
import os
import filecmp
import shutil

import cfg
import utils

class TestFunctionalHumannEndtoEnd(unittest.TestCase):
    """
    Test humann with end to end functional tests
    """

    def test_humann_fastq(self):
        """
        Test the standard humann flow on a fastq input file
        """
        
        # create a temp directory for output
        tempdir = utils.create_temp_folder("fastq")
        
        # run humann test
        command = ["humann","--input",cfg.demo_fastq,"--output",tempdir]
        utils.run_humann(command)
        
        # check the output files are as expected
        for expression, message in utils.check_output(cfg.expected_demo_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    def test_humann_fasta(self):
        """
        Test the standard humann flow on a fasta input file
        """
        
        # create a temp directory for output
        tempdir = utils.create_temp_folder("fasta")
        
        # run humann test
        command = ["humann","--input",cfg.demo_fasta,"--output",tempdir]
        utils.run_humann(command)
        
        # check the output files are as expected
        for expression, message in utils.check_output(cfg.expected_demo_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    def test_humann_sam(self):
        """
        Test the standard humann flow on a sam input file
        """
        
        # create a temp directory for output
        tempdir = utils.create_temp_folder("sam")
        
        # run humann test
        command = ["humann","--input",cfg.demo_sam,"--output",tempdir]
        utils.run_humann(command)
        
        # check the output files are as expected
        for expression, message in utils.check_output(cfg.expected_demo_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    def test_humann_m8(self):
        """
        Test the standard humann flow on a m8 input file
        """
        
        # create a temp directory for output
        tempdir = utils.create_temp_folder("m8")
        
        # run humann test
        command = ["humann","--input",cfg.demo_m8,"--output",tempdir]
        utils.run_humann(command)
        
        # check the output files are as expected
        for expression, message in utils.check_output(cfg.expected_demo_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    def test_humann_gene_families(self):
        """
        Test the standard humann flow on a gene families output file as input
        """
        
        # create a temp directory for output
        tempdir = utils.create_temp_folder("gene_families")
        
        # run humann test
        command = ["humann","--input",cfg.demo_gene_families,"--output",tempdir]
        utils.run_humann(command)
        
        # check the output files are as expected
        # it will include all output files except the gene families output file
        # since this file was used as input
        for expression, message in utils.check_output(cfg.expected_demo_output_files_genefamilies_input, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    def test_humann_fastq_bypass_nucleotide_search(self):
        """
        Test the standard humann flow on a fastq input file
        Test with bypassing nucleotide search
        """
        
        # create a temp directory for output
        tempdir = utils.create_temp_folder("fastq_bypass_nucleotide_search")
        
        # run humann test
        command = ["humann","--input",cfg.demo_fastq,"--output",tempdir,"--bypass-nucleotide-search"]
        utils.run_humann(command)
        
        # check the output files are as expected
        for expression, message in utils.check_output(cfg.expected_demo_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    def test_humann_fasta_bypass_nucleotide_search(self):
        """
        Test the standard humann flow on a fasta input file
        Test with bypassing nucleotide search
        """
        
        # create a temp directory for output
        tempdir = utils.create_temp_folder("fasta_bypass_nucleotide_search")
        
        # run humann test
        command = ["humann","--input",cfg.demo_fasta,"--output",tempdir,"--bypass-nucleotide-search"]
        utils.run_humann(command)
        
        # check the output files are as expected
        for expression, message in utils.check_output(cfg.expected_demo_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    def test_humann_fastq_bypass_translated_search(self):
        """
        Test the standard humann flow on a fastq input file
        Test with bypassing translated search
        """
        
        # create a temp directory for output
        tempdir = utils.create_temp_folder("fastq_bypass_translated_search")
        
        # run humann test
        command = ["humann","--input",cfg.demo_fastq,"--output",tempdir,"--bypass-translated-search"]
        utils.run_humann(command)
        
        # check the output files are as expected
        for expression, message in utils.check_output(cfg.expected_demo_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    def test_humann_fasta_bypass_translated_search(self):
        """
        Test the standard humann flow on a fasta input file
        Test with bypassing translated search
        """
        
        # create a temp directory for output
        tempdir = utils.create_temp_folder("fasta_bypass_translated_search")
        
        # run humann test
        command = ["humann","--input",cfg.demo_fasta,"--output",tempdir,"--bypass-translated-search"]
        utils.run_humann(command)
        
        # check the output files are as expected
        for expression, message in utils.check_output(cfg.expected_demo_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    def test_humann_fastq_bypass_prescreen(self):
        """
        Test the standard humann flow on a fastq input file
        Test with bypassing prescreen
        """
        
        # create a temp directory for output
        tempdir = utils.create_temp_folder("fastq_bypass_prescreen")
        
        # run humann test
        command = ["humann","--input",cfg.demo_fastq,"--output",tempdir,"--bypass-prescreen"]
        utils.run_humann(command)
        
        # check the output files are as expected
        for expression, message in utils.check_output(cfg.expected_demo_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    def test_humann_fasta_bypass_prescreen(self):
        """
        Test the standard humann flow on a fasta input file
        Test with bypassing prescreen
        """
        
        # create a temp directory for output
        tempdir = utils.create_temp_folder("fasta_bypass_prescreen")
        
        # run humann test
        command = ["humann","--input",cfg.demo_fasta,"--output",tempdir,"--bypass-prescreen"]
        utils.run_humann(command)
        
        # check the output files are as expected
        for expression, message in utils.check_output(cfg.expected_demo_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)

    def test_humann_fastq_custom_taxonomic_profile(self):
        """
        Test the standard humann flow on a fastq input file
        Test with a custom taxonomic profile
        """
        
        # create a temp directory for output
        tempdir = utils.create_temp_folder("fastq_custom_taxonomic_profile")
        
        # run humann test
        command = ["humann","--input",cfg.demo_fastq,"--output",tempdir,"--taxonomic-profile",
                   cfg.demo_bugs_list]
        utils.run_humann(command)
        
        # check the output files are as expected
        for expression, message in utils.check_output(cfg.expected_demo_output_files, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
