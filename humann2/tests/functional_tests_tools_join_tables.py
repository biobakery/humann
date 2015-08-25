import unittest
import subprocess
import tempfile
import os
import filecmp
import shutil

import cfg
import utils

class TestFunctionalHumann2JoinTables(unittest.TestCase):
    """
    Test humann2.tools.join_tables
    """

    def test_join_tsv(self):
        """
        Test joining tsv files
        """
        
        # create a temp file
        file_out, new_file=tempfile.mkstemp(prefix="humann2_temp")
        
        # join the files
        try:
            stdout=subprocess.check_output(["humann2_join_tables","--input",
                cfg.data_folder,"--output",new_file,"--file_name",
                cfg.multi_sample_genefamilies_split_basename])
        except EnvironmentError:
            print("Warning: Unable to execute humann2_join_table in test.")
        
        # check the joined file is as expected
        self.assertTrue(filecmp.cmp(new_file, cfg.multi_sample_genefamilies, 
            shallow=False))

        # remove the temp file
        utils.remove_temp_file(new_file)
        
            