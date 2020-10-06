import unittest
import filecmp
import logging
import tempfile

import cfg
import utils

from humann import utilities
from humann import config

class TestAdvancedHumannUtilitiesFunctions(unittest.TestCase):
    """
    Test the functions found in humann.utilities
    """
    
    def setUp(self):
        config.unnamed_temp_dir=tempfile.gettempdir()
        
        # set up nullhandler for logger
        logging.getLogger('humann.utilities').addHandler(logging.NullHandler())

    def test_fastq_to_fasta_with_pick_frames(self):
        """
        Test the fastq_to_fasta function with pick frames
        """
         
        new_fasta_file=utilities.fastq_to_fasta(
            cfg.convert_fastq_file, apply_pick_frames=True)
        self.assertTrue(filecmp.cmp(new_fasta_file,
            cfg.convert_fasta_pick_frames_file, shallow=False))
        utils.remove_temp_file(new_fasta_file)     
        
    def test_pick_frames_from_fasta(self):
        """
        Test the pick_frames_from_fasta function
        """
         
        new_fasta_file=utilities.pick_frames_from_fasta(
            cfg.convert_fasta_multiline_file)
        self.assertTrue(filecmp.cmp(new_fasta_file,
            cfg.convert_fasta_pick_frames_file, shallow=False))
        utils.remove_temp_file(new_fasta_file)     
                     
