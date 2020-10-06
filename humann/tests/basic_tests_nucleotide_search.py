import unittest
import logging

import cfg
import utils
import tempfile

from humann.search import nucleotide
from humann import config

class TestBasicHumannNucleotideSearchFunctions(unittest.TestCase):
    """
    Test the functions found in humann.search.nucleotide
    """
    
    def setUp(self):
        config.unnamed_temp_dir=tempfile.gettempdir()
        config.temp_dir=tempfile.gettempdir()
        config.file_basename="HUMAnN_test"
        
        # set up nullhandler for logger
        logging.getLogger('humann.search.nucleotide').addHandler(logging.NullHandler())
        
    def test_calculate_percent_identity(self):
        """
        Test the calculate percent identity function
        """
        
        cigar_string="100S84M16S"
        md_field="MD:Z:27A5G2T6G4T1A6T2C0A7A1A0A0C0G1G1A5"
        
        expected_identity= 100.0 * ( 68 / 84.0 )
        
        identity, alignment_length, reference_length =nucleotide.calculate_percent_identity(cigar_string,md_field)
        
        self.assertEqual(identity, expected_identity)

    def test_calculate_percent_identity_insert_delete(self):
        """
        Test the calculate percent identity function
        Test cigar string with insert and deletes
        """
        
        cigar_string="100I84M2D"
        md_field="MD:Z:27A5G2T6G4T1A6T2C0A7A1A0A0C0G1G1A5"
        
        expected_identity= 100.0 * ( 68 / 186.0 )
        
        identity, alignment_length, reference_length =nucleotide.calculate_percent_identity(cigar_string,md_field)
        
        self.assertEqual(identity, expected_identity)

    def test_calculate_percent_identity_two_match_mismatch_identifiers(self):
        """
        Test the calculate percent identity function
        Test cigar string with the two match and mistmatch identifiers (X and =)
        """

        cigar_string="100I84M1X2=2D"
        md_field="MD:Z:27A5G2T6G4T1A6T2C0A7A1A0A0C0G1G1A5"

        expected_identity= 100.0 * ( 68 / 189.0 )

        identity, alignment_length, reference_length =nucleotide.calculate_percent_identity(cigar_string,md_field)

        self.assertEqual(identity, expected_identity)

    def test_calculate_percent_identity_multiple_M_cigar_fields(self):
        """
        Test the calculate percent identity function
        Test with multiple M fields
        """

        cigar_string="100S84M16S10M"
        md_field="MD:Z:27A5G2T6G4T1A6T2C0A7A1A0A0C0G1G1A5"

        expected_identity= 100.0 * ( 68 / 94.0 )

        identity, alignment_length, reference_length =nucleotide.calculate_percent_identity(cigar_string,md_field)

        self.assertEqual(identity, expected_identity)

    def test_calculate_percent_identity_multiple_M_cigar_fields_reference_length(self):
        """
        Test the calculate percent identity function
        Test with multiple M fields to compute reference length
        """

        cigar_string="100S84M16S10M5D2N3S4I3=3X"
        md_field="MD:Z:27A5G2T6G4T1A6T2C0A7A1A0A0C0G1G1A5"

        expected_length= 84+10+5+2+3+3

        identity, alignment_length, reference_length =nucleotide.calculate_percent_identity(cigar_string,md_field)

        self.assertEqual(reference_length, expected_length)

    def test_calculate_percent_identity_simple_md_field(self):
        """
        Test the calculate percent identity function
        Test with a simple md field where all bases match
        """
        
        cigar_string="100S84M16S"
        md_field="MD:Z:84"
        
        expected_identity= 100.0 * ( 84 / 84.0 )
        
        identity, alignment_length, reference_length =nucleotide.calculate_percent_identity(cigar_string,md_field)
        
        self.assertEqual(identity, expected_identity)
        
    def test_calculate_percent_identity_simple_cigar_string(self):
        """
        Test the calculate percent identity function
        Test with a simple cigar string
        """
        
        cigar_string="84M"
        md_field="MD:Z:27A5G2T6G4T1A6T2C0A7A1A0A0C0G1G1A5"
        
        expected_identity= 100.0 * ( 68 / 84.0 )
        
        identity, alignment_length, reference_length=nucleotide.calculate_percent_identity(cigar_string,md_field)
        
        self.assertEqual(identity, expected_identity)

    def test_calculate_percent_identity_missing_cigar_identifier(self):
        """
        Test the calculate percent identity function
        Test with a cigar string that does not have the match/mismatch identifier
        """
        
        cigar_string="*"
        md_field="MD:Z:84"
        
        expected_identity=0.0
        
        identity, alignment_length, reference_length=nucleotide.calculate_percent_identity(cigar_string,md_field)
        
        self.assertEqual(identity, expected_identity)
        
    def test_calculate_percent_identity_missing_md_field(self):
        """
        Test the calculate percent identity function
        Test with an empty md field
        """
        
        cigar_string="100S84M16S"
        md_field=""
        
        expected_identity=0.0
        
        identity, alignment_length, reference_length=nucleotide.calculate_percent_identity(cigar_string,md_field)
        
        self.assertEqual(identity, expected_identity)
        
    def test_find_md_field(self):
        """
        Test the find md field function
        """
        
        info=['r99983|640963016.fna|116533|116684|_from_', '16', 
            'gi|223667726|ref|NZ_DS264546.1|:115458-117362|46503|g__Parabacteroides.s__Parabacteroides_merdae|UniRef90_R6Y4Q3|UniRef50_D6D6J1|1905', 
            '927', '42', '151M', '*', '0', '0', 
            'AAACAATGGAAAAGTACAAGAGCGACTATGGCAAATGTTTTACAACATGGATGCGATTCAACAGGAATGGGTTTGGTTTACCACTCGTGATAAAGATACGGGCTGGTCTGGCGATGTCTTGCCTCCTTCCCAAAATGGTCATGCCCGTCAA',
            '0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000',
            'AS:i:0', 'XN:i:0', 'XM:i:0', 'XO:i:0', 'XG:i:0', 'NM:i:0', 'MD:Z:151', 'YT:Z:UU']
        
        expected_result=info[-2]
        
        result=nucleotide.find_md_field(info)
        
        self.assertEqual(result, expected_result)
        
    def test_find_md_field_missing_column(self):
        """
        Test the find md field function
        Test with an example that does not include the md field column
        """
        
        info=['r99983|640963016.fna|116533|116684|_from_', '16', 
            'gi|223667726|ref|NZ_DS264546.1|:115458-117362|46503|g__Parabacteroides.s__Parabacteroides_merdae|UniRef90_R6Y4Q3|UniRef50_D6D6J1|1905', 
            '927', '42', '151M', '*', '0', '0', 
            'AAACAATGGAAAAGTACAAGAGCGACTATGGCAAATGTTTTACAACATGGATGCGATTCAACAGGAATGGGTTTGGTTTACCACTCGTGATAAAGATACGGGCTGGTCTGGCGATGTCTTGCCTCCTTCCCAAAATGGTCATGCCCGTCAA',
            '0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000',
            'AS:i:0', 'XN:i:0', 'XM:i:0', 'XO:i:0', 'XG:i:0', 'NM:i:0', 'YT:Z:UU']
        
        expected_result=""
        
        result=nucleotide.find_md_field(info)
        
        self.assertEqual(result, expected_result)
        
        
        
