import unittest
import subprocess
import tempfile
import os
import filecmp
import shutil

import cfg
import utils

class TestFunctionalHumann2Tools(unittest.TestCase):
    """
    Test humann2.tools
    """

    def test_humann2_join_tables_tsv(self):
        """
        Test joining tsv files with humann2_join_tables
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # join the files
        utils.run_command(["humann2_join_tables","--input",
                           cfg.data_folder,"--output",new_file,"--file_name",
                           cfg.multi_sample_genefamilies_split_basename,"--verbose"])
        
        # check the joined file is as expected
        self.assertTrue(filecmp.cmp(new_file, cfg.multi_sample_genefamilies, 
            shallow=False))

        # remove the temp file
        utils.remove_temp_file(new_file)
 
    def test_humann2_split_tables_tsv(self):
        """
        Test splitting a tsv file with humann2_split_tables
        """

        input_file=cfg.multi_sample_genefamilies

        # create a temp directory
        temp_directory=utils.create_temp_folder("split_tables_tsv")

        # split the file
        utils.run_command(["humann2_split_table","--input", input_file,
                           "--output",temp_directory,"--verbose"])

        # test the split files are as expected
        output_files=os.listdir(temp_directory)

        # sort the output files
        file_pairs=[]
        for file in output_files:
            filebasename=os.path.basename(file)
            # get the sample number for the file
            file=os.path.join(temp_directory,file)
            if filebasename[-1] == 1:
                file_pairs.append([file,cfg.multi_sample_genefamilies_split1])
            elif filebasename[-1] == 2:
                file_pairs.append([file,cfg.multi_sample_genefamilies_split2])

        for temp_file, file in file_pairs:
            self.assertTrue(filecmp.cmp(temp_file, file, shallow=False))

        # remove the temp folder
        utils.remove_temp_folder(temp_directory)
        

    def test_humann2_regroup_table_uniref50_rxn_tsv(self):
        """
        Test regrouping the tsv file with humann2_regroup_table
        Test with uniref50 to reactions mappings
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_regroup_table","--input",cfg.regroup_input,"--output",
                           new_file,"--groups","uniref50_rxn"])
        
        # check the output is as expected
        self.assertTrue(filecmp.cmp(new_file, cfg.regroup_rxn_output, shallow=False))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_regroup_table_uniref50_rxn_tsv(self):
        """
        Test regrouping the tsv file with humann2_regroup_table
        Test with uniref50 to reactions mappings
        Test with the mean instead of sum output
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_regroup_table","--input",cfg.regroup_input,"--output",
                           new_file,"--groups","uniref50_rxn","--function","mean"])
        
        # check the output is as expected
        self.assertTrue(filecmp.cmp(new_file, cfg.regroup_rxn_mean_output, shallow=False))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_regroup_table_uniref50_ec_tsv(self):
        """
        Test regrouping the tsv file with humann2_regroup_table
        Test with uniref50 to ec mappings
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_regroup_table","--input",cfg.regroup_input,"--output",
                           new_file,"--groups","uniref50_ec"])
        
        # check the output is as expected
        self.assertTrue(filecmp.cmp(new_file, cfg.regroup_ec_output, shallow=False))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_regroup_table_uniref50_go_tsv(self):
        """
        Test regrouping the tsv file with humann2_regroup_table
        Test with uniref50 to go mappings
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_regroup_table","--input",cfg.regroup_input,"--output",
                           new_file,"--groups","uniref50_go"])
        
        # check the output is as expected
        self.assertTrue(filecmp.cmp(new_file, cfg.regroup_go_output, shallow=False))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_regroup_table_uniref50_ko_tsv(self):
        """
        Test regrouping the tsv file with humann2_regroup_table
        Test with uniref50 to ko mappings
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_regroup_table","--input",cfg.regroup_input,"--output",
                           new_file,"--groups","uniref50_ko"])
        
        # check the output is as expected
        self.assertTrue(filecmp.cmp(new_file, cfg.regroup_ko_output, shallow=False))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_regroup_table_custom_grouping_tsv(self):
        """
        Test regrouping the tsv file with humann2_regroup_table
        Test with custom mappings
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_regroup_table","--input",cfg.regroup_custom_input,"--output",
                           new_file,"--custom",cfg.regroup_custom_groups])
        
        # check the output is as expected
        self.assertTrue(filecmp.cmp(new_file, cfg.regroup_custom_groups_output,
                                    shallow=False))

        # remove the temp file
        utils.remove_temp_file(new_file)
       
    def test_humann2_rename_table_uniref50_tsv(self):
        """
        Test renaming the tsv file entries with humann2_rename_table
        Test with uniref50 names
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_rename_table","--input",cfg.rename_input,"--output",
                           new_file,"--names","uniref50"])
        
        # check the output is as expected
        self.assertTrue(filecmp.cmp(new_file, cfg.rename_uniref50_output, shallow=False))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_rename_table_ko_tsv(self):
        """
        Test renaming the tsv file entries with humann2_rename_table
        Test with ko names
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_rename_table","--input",cfg.rename_ko_input,"--output",
                           new_file,"--names","ko"])
        
        # check the output is as expected
        self.assertTrue(filecmp.cmp(new_file, cfg.rename_ko_output, shallow=False))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_rename_table_ec_tsv(self):
        """
        Test renaming the tsv file entries with humann2_rename_table
        Test with ec names
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_rename_table","--input",cfg.rename_ec_input,"--output",
                           new_file,"--names","ec"])
        
        # check the output is as expected
        self.assertTrue(filecmp.cmp(new_file, cfg.rename_ec_output, shallow=False))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_rename_table_rxn_tsv(self):
        """
        Test renaming the tsv file entries with humann2_rename_table
        Test with rxn names
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_rename_table","--input",cfg.rename_rxn_input,"--output",
                           new_file,"--names","metacyc-rxn"])
        
        # check the output is as expected
        self.assertTrue(filecmp.cmp(new_file, cfg.rename_rxn_output, shallow=False))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_rename_table_pathways_tsv(self):
        """
        Test renaming the tsv file entries with humann2_rename_table
        Test with pathways names
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_rename_table","--input",cfg.rename_pathway_input,"--output",
                           new_file,"--names","metacyc-pwy"])
        
        # check the output is as expected
        self.assertTrue(filecmp.cmp(new_file, cfg.rename_pathway_output, shallow=False))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
    def test_humann2_rename_table_custom_tsv(self):
        """
        Test renaming the tsv file entries with humann2_rename_table
        Test with custom names file
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # run the command
        utils.run_command(["humann2_rename_table","--input",cfg.rename_input,"--output",
                           new_file,"--custom",cfg.rename_custom_mapping])
        
        # check the output is as expected
        self.assertTrue(filecmp.cmp(new_file, cfg.rename_custom_output, shallow=False))

        # remove the temp file
        utils.remove_temp_file(new_file)
            