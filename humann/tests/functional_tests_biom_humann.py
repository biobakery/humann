import unittest
import subprocess
import tempfile
import os
import filecmp
import shutil

import cfg
import utils

class TestFunctionalHumannEndtoEndBiom(unittest.TestCase):
    """
    Test humann with end to end functional tests
    """

    def test_humann_fastq_biom_output(self):
        """
        Test the standard humann flow on a fastq input file
        Test biom output is written
        """
        
        # create a temp directory for output
        tempdir = utils.create_temp_folder("fastq")
        
        # run humann test
        command = ["humann","--input",cfg.demo_fastq,"--output",tempdir,
                   "--output-format", "biom"]
        utils.run_humann(command)
        
        # check the output files are as expected
        for expression, message in utils.check_output(cfg.expected_demo_output_files_biom, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    def test_humann_fastq_biom_output_pathways(self):
        """
        Test the standard humann flow on a fastq input file
        Test biom output is written
        Test the expected pathways are identified
        """
        
        # create a temp directory for output
        tempdir = utils.create_temp_folder("fastq")
        
        # run humann test
        command = ["humann","--input",cfg.demo_fastq,"--output",tempdir,
                   "--output-format", "biom", "--gap-fill", "off"]
        utils.run_humann(command)
        
        # check the output file of pathway abundance has the expected pathways
        pathways_file_tsv=utils.read_biom_table(os.path.join(tempdir,"demo_pathabundance.biom"))
        pathways_found=set([x.split("\t")[0].split(":")[0] for x in filter(lambda x: "PWY" in x, pathways_file_tsv)])
        
        self.assertEqual(pathways_found,cfg.expected_demo_output_files_biom_pathways)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
        
    def test_humann_gene_families_biom_input(self):
        """
        Test the standard humann flow on a gene families output file as input
        Test with the biom format of the gene families file
        """
        
        # create a temp directory for output
        tempdir = utils.create_temp_folder("gene_families")
        
        # run humann test
        command = ["humann","--input",cfg.demo_gene_families_biom,"--output",tempdir]
        utils.run_humann(command)
        
        # check the output files are as expected
        # it will include all output files except the gene families output file
        # since this file was used as input
        for expression, message in utils.check_output(cfg.expected_demo_output_files_genefamilies_input, tempdir):
            self.assertTrue(expression,message)

        # remove the temp directory
        utils.remove_temp_folder(tempdir)
