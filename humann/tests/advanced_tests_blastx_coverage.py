import unittest
import logging

import cfg
import utils
import tempfile
import re

from humann2.search import blastx_coverage
from humann2 import store
from humann2 import config

class TestBasicHumann2Blastx_CoverageFunctions(unittest.TestCase):
    """
    Test the functions found in humann2.search.nucleotide
    """
    
    def setUp(self):
        config.unnamed_temp_dir=tempfile.gettempdir()
        config.temp_dir=tempfile.gettempdir()
        config.file_basename="HUMAnN2_test"
        
        # set up nullhandler for logger
        logging.getLogger('humann2.search.blastx_coverage').addHandler(logging.NullHandler())
        
    def test_blastx_coverage_gene_names_default(self):
        """
        Test the blastx_coverage function
        Test the gene names
        Test without filter
        """
        
        # set the coverage threshold to zero so as to not test with filter on
        current_coverage_threshold=config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold=0
        
        # get the set of allowed proteins
        allowed_proteins = blastx_coverage.blastx_coverage(cfg.rapsearch2_output_file_without_header_coverage,
            config.translated_subject_coverage_threshold, log_messages=True)
        
        # load the blastm8-like output
        file_handle=open(cfg.rapsearch2_output_file_without_header_coverage)
        
        all_proteins = set()
        for line in file_handle:
            if not re.search("^#",line):
                data=line.strip().split(config.blast_delimiter)
                protein_name=data[config.blast_reference_index].split(config.chocophlan_delimiter)[0]
                all_proteins.add(protein_name)

        # reset the coverage threshold
        config.translated_subject_coverage_threshold=current_coverage_threshold
        
        # check the expected proteins are found
        self.assertEqual(sorted(all_proteins),sorted(allowed_proteins))
        
    def test_blastx_coverage_gene_names_custom_annotation(self):
        """
        Test the blastx_coverage function
        Test the gene names with custom annotation
        Test without filter
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        
        # set the coverage threshold to zero so as to not test with filter on
        current_coverage_threshold=config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold=0
        
        # get the set of allowed proteins
        allowed_proteins = blastx_coverage.blastx_coverage(cfg.rapsearch2_output_file_without_header_coverage_custom_annotations,
            config.translated_subject_coverage_threshold, alignments, log_messages=True)
        
        # load the blastm8-like output
        file_handle=open(cfg.rapsearch2_output_file_without_header_coverage_custom_annotations)
        
        all_proteins = set()
        for line in file_handle:
            if not re.search("^#",line):
                data=line.strip().split(config.blast_delimiter)
                protein_name=data[config.blast_reference_index].split(config.chocophlan_delimiter)[0]
                all_proteins.add(protein_name)

        # reset the coverage threshold
        config.translated_subject_coverage_threshold=current_coverage_threshold
        
        # check the expected proteins are found
        self.assertEqual(sorted(all_proteins),sorted(allowed_proteins))
        
    def test_blastx_coverage_gene_names_chocophlan_annoation(self):
        """
        Test the blastx_coverage function
        Test the gene names with chocophlan annotations
        Test without filter
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        
        # set the coverage threshold to zero so as to not test with filter on
        current_coverage_threshold=config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold=0
        
        # get the set of allowed proteins
        allowed_proteins = blastx_coverage.blastx_coverage(cfg.rapsearch2_output_file_without_header_coverage_chocophlan_annotations,
            config.translated_subject_coverage_threshold, alignments, log_messages=True)
        
        # load the blastm8-like output
        file_handle=open(cfg.rapsearch2_output_file_without_header_coverage_chocophlan_annotations)
        
        all_proteins = set()
        for line in file_handle:
            if not re.search("^#",line):
                data=line.strip().split(config.blast_delimiter)
                protein_name=data[config.blast_reference_index].split(config.chocophlan_delimiter)[config.chocophlan_gene_indexes[0]]
                all_proteins.add(protein_name)

        # reset the coverage threshold
        config.translated_subject_coverage_threshold=current_coverage_threshold
        
        # check the expected proteins are found
        self.assertEqual(sorted(all_proteins),sorted(allowed_proteins))
        
    def test_blastx_coverage_gene_names_id_mapping(self):
        """
        Test the blastx_coverage function
        Test the gene names with chocophlan annotations
        Test without filter
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        
        # process the id mapping
        alignments.process_id_mapping(cfg.coverage_id_mapping_file)
        
        # set the coverage threshold to zero so as to not test with filter on
        current_coverage_threshold=config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold=0
        
        # get the set of allowed proteins
        allowed_proteins = blastx_coverage.blastx_coverage(cfg.rapsearch2_output_file_without_header_coverage,
            config.translated_subject_coverage_threshold, alignments, log_messages=True)
        
        # load the blastm8-like output
        file_handle=open(cfg.rapsearch2_output_file_without_header_coverage)
        
        all_proteins = set()
        for line in file_handle:
            if not re.search("^#",line):
                data=line.strip().split(config.blast_delimiter)
                # just like the id mapping, remove the UniRef50_
                protein_name=data[config.blast_reference_index].split(config.chocophlan_delimiter)[0]
                protein_name = protein_name.replace("UniRef50_","")
                all_proteins.add(protein_name)

        # reset the coverage threshold
        config.translated_subject_coverage_threshold=current_coverage_threshold
        
        # check the expected proteins are found
        self.assertEqual(sorted(all_proteins),sorted(allowed_proteins))
        
    def test_blastx_coverage(self):
        """
        Test the coverage filter
        Test with one protein with one alignment passing threshold
        Test with one protein with two alignments passing threshold (does not pass with only one alignment)
        Test with other proteins with one more more alignments not passing threshold
        """
        
        # create a set of alignments
        alignments=store.Alignments()
        
        # set the coverage threshold to a small value so as to have some alignments pass
        current_coverage_threshold=config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold=50.0
        
        # get the set of allowed proteins
        allowed_proteins = blastx_coverage.blastx_coverage(cfg.rapsearch2_output_file_without_header_coverage,
            config.translated_subject_coverage_threshold, alignments, True)
        
        # load the blastm8-like output
        file_handle=open(cfg.rapsearch2_output_file_without_header_coverage)
        
        found_proteins=set()
        for line in file_handle:
            if not re.search("^#",line):
                data=line.strip().split(config.blast_delimiter)
                
                referenceid=data[config.blast_reference_index]
                gene, length, bug = alignments.process_reference_annotation(referenceid)
                queryid=data[config.blast_query_index]
                identity=float(data[config.blast_identity_index])
                alignment_length=float(data[config.blast_aligned_length_index])
            
                # the proteins that pass have "_coverage50" as part of their names
                if "_coverage50" in gene:
                    found_proteins.add(gene)   
        file_handle.close()
        
        # reset the coverage threshold
        config.translated_subject_coverage_threshold=current_coverage_threshold
        
        # check the values are unchanged
        self.assertEqual(sorted(allowed_proteins), sorted(found_proteins))
