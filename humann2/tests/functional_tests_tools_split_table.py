import unittest
import subprocess
import tempfile
import os
import filecmp
import shutil

import cfg
import utils

class TestFunctionalHumann2SplitTable(unittest.TestCase):
    """
    Test humann2.tools.split_table
    """

    def test_split_tsv(self):
        """
        Test splitting a tsv file
        """
        
        input_file=cfg.multi_sample_genefamilies
        
        # create a temp directory
        temp_directory=tempfile.mkdtemp(prefix="humann2_test")
        
        # split the file
        try:
            stdout=subprocess.check_output(["humann2_split_table","--input",
                input_file,"--output",temp_directory])
        except EnvironmentError:
            print("Warning: Unable to execute humann2_split_table in test.")
            
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


