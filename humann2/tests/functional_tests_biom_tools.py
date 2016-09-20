import unittest
import tempfile
import os

import cfg
import utils

class TestFunctionalHumann2ToolsBiom(unittest.TestCase):
    """
    Test humann2.tools
    """

    def test_humann2_join_tables_biom(self):
        """
        Test joining biom files with humann2_join_tables
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # join the files
        utils.run_command(["humann2_join_tables","--input",
                           cfg.data_folder,"--output",new_file,"--file_name",
                           cfg.multi_sample_genefamilies_split_basename_biom,"--verbose"])
        
        # check the joined file is as expected
        self.assertTrue(utils.check_output(new_file))

        # remove the temp file
        utils.remove_temp_file(new_file)
 
    def test_humann2_split_tables_tsv(self):
        """
        Test splitting a tsv file with humann2_split_tables
        """

        input_file=cfg.multi_sample_genefamilies_biom

        # create a temp directory
        temp_directory=utils.create_temp_folder("split_tables_biom")

        # split the file
        utils.run_command(["humann2_split_table","--input", input_file,
                           "--output",temp_directory,"--verbose"])

        # test the split files are as expected
        for file in cfg.multi_sample_split_files_biom:
            self.assertTrue(utils.check_output(file,temp_directory))

        # remove the temp folder
        utils.remove_temp_folder(temp_directory)
        

    def test_humann2_regroup_table_uniref50_rxn_biom(self):
        """
        Test regrouping the biom file with humann2_regroup_table
        Test with uniref50 to reactions mappings
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_regroup_table","--input",cfg.regroup_input_biom,"--output",
                           new_file,"--groups","uniref50_rxn"])
        
        # check the output is as expected
        self.assertTrue(utils.check_output(new_file))

        # remove the temp file
        utils.remove_temp_file(new_file)
       
    def test_humann2_rename_table_uniref50_biom(self):
        """
        Test renaming the biom file entries with humann2_rename_table
        Test with uniref50 names
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_rename_table","--input",cfg.rename_input_biom,"--output",
                           new_file,"--names","uniref50"])
        
        # check the output is as expected
        self.assertTrue(utils.check_output(new_file))

        # remove the temp file
        utils.remove_temp_file(new_file)     
            
    def test_humann2_renorm_table_cpm_biom(self):
        """
        Test renorm the biom file entries with humann2_renorm_table
        Test with cpm
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_renorm_table","--input",cfg.renorm_input_biom,"--output",
                           new_file,"--units","cpm"])
        
        # check the output is as expected
        self.assertTrue(utils.check_output(new_file))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_renorm_table_cpm_biom_output(self):
        """
        Test renorm the biom file entries with humann2_renorm_table
        Test with cpm
        Test with biom output
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp",suffix=".biom")
        
        # run the command
        utils.run_command(["humann2_renorm_table","--input",cfg.renorm_input_biom,"--output",
                           new_file,"--units","cpm"])
        
        # check the output is as expected
        self.assertTrue(utils.check_output(new_file))

        # remove the temp file
        utils.remove_temp_file(new_file)
