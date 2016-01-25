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
       
            