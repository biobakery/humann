import unittest
import logging

import cfg
import utils
import tempfile
import re

from humann.search import blastx_coverage
from humann import store
from humann import config

class TestBasicHumannBlastx_CoverageFunctions(unittest.TestCase):
    """
    Test the functions found in humann.search.nucleotide
    """
    
    def setUp(self):
        config.unnamed_temp_dir=tempfile.gettempdir()
        config.temp_dir=tempfile.gettempdir()
        config.file_basename="HUMAnN_test"
        
        # set up nullhandler for logger
        logging.getLogger('humann.search.blastx_coverage').addHandler(logging.NullHandler())

    def gene_names(self, allowed_references, alignments=None):
        """
        blastx_coverage returns the full reference annotations, since coverage is
        computed per reference sequence, so map them back to gene family names
        """

        if alignments is None:
            alignments=store.Alignments()

        return set(alignments.process_reference_annotation(reference)[0]
            for reference in allowed_references)

    def test_blastx_coverage_gene_names_default(self):
        """
        Test the blastx_coverage function
        Test the gene names
        Test without filter
        """
        
        # set the coverage threshold to zero so as to not test with filter on
        current_coverage_threshold=config.translated_subject_coverage_threshold
        config.translated_subject_coverage_threshold=0
        
        # get the set of allowed reference sequences
        allowed_references = blastx_coverage.blastx_coverage(cfg.rapsearch2_output_file_without_header_coverage,
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
        self.assertEqual(sorted(all_proteins),sorted(self.gene_names(allowed_references)))
        
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
        
        # get the set of allowed reference sequences
        allowed_references = blastx_coverage.blastx_coverage(cfg.rapsearch2_output_file_without_header_coverage_custom_annotations,
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
        self.assertEqual(sorted(all_proteins),sorted(self.gene_names(allowed_references)))
        
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
        
        # get the set of allowed reference sequences
        allowed_references = blastx_coverage.blastx_coverage(cfg.rapsearch2_output_file_without_header_coverage_chocophlan_annotations,
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
        self.assertEqual(sorted(all_proteins),sorted(self.gene_names(allowed_references)))
        
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
        
        # get the set of allowed reference sequences
        allowed_references = blastx_coverage.blastx_coverage(cfg.rapsearch2_output_file_without_header_coverage,
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
        self.assertEqual(sorted(all_proteins),sorted(self.gene_names(allowed_references,alignments)))
        
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
        
        # get the set of allowed reference sequences
        allowed_references = blastx_coverage.blastx_coverage(cfg.rapsearch2_output_file_without_header_coverage,
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
        self.assertEqual(sorted(self.gene_names(allowed_references,alignments)), sorted(found_proteins))

    def write_temp_alignments(self, alignments_list):
        """
        Write a blastm8-like file of the alignments provided as
        ( reference, subject start, subject stop ) tuples
        Return the name of the file written
        """

        file_handle=tempfile.NamedTemporaryFile(mode="w", suffix=".m8", delete=False)
        for index, (reference, subject_start, subject_stop) in enumerate(alignments_list):
            data=[""] * config.blast_total_columns
            data[config.blast_query_index]="r"+str(index)
            data[config.blast_reference_index]=reference
            data[config.blast_identity_index]="100.0"
            data[config.blast_aligned_length_index]=str(subject_stop - subject_start)
            data[config.blast_query_start_index]="1"
            data[config.blast_query_end_index]=str(subject_stop - subject_start)
            data[config.blast_subject_start_index]=str(subject_start)
            data[config.blast_subject_end_index]=str(subject_stop)
            data[config.blast_evalue_index]="0"
            file_handle.write(config.blast_delimiter.join(data)+"\n")
        file_handle.close()

        return file_handle.name

    def test_blastx_coverage_repeated_gene_family_uses_per_sequence_length(self):
        """
        Test the coverage filter with the same gene family present in two species
        Test that each reference sequence is normalized by its own length
        Test that the result does not depend on the order of the alignments
        """

        # the same UniRef90 in two species, with very different sequence lengths
        short_reference="1__X__gene1|g__Bacteroides.s__Bacteroides_dorei|UniRef90_X|UniRef50_X|100"
        long_reference="2__X__gene1|g__Bacteroides.s__Bacteroides_uniformis|UniRef90_X|UniRef50_X|1000"

        # both sequences are hit over the same 59 positions, which covers 59% of the
        # short sequence but only 5.9% of the long sequence
        alignments_list=[(short_reference, 1, 60), (long_reference, 1, 60)]

        for ordered_alignments in [alignments_list, list(reversed(alignments_list))]:
            alignment_file=self.write_temp_alignments(ordered_alignments)

            allowed_references = blastx_coverage.blastx_coverage(alignment_file,
                50.0, log_messages=True, nucleotide=True)

            utils.remove_temp_file(alignment_file)

            # only the sequence that is actually covered should pass
            self.assertEqual(sorted(allowed_references), [short_reference])

            # the gene family is still resolved from the allowed reference sequence
            self.assertEqual(sorted(self.gene_names(allowed_references)), ["UniRef90_X"])

    def test_blastx_coverage_repeated_gene_family_hits_not_merged(self):
        """
        Test the coverage filter with the same gene family present in two species
        Test that the hit positions of the two sequences are not merged
        """

        first_reference="1__X__gene1|g__Bacteroides.s__Bacteroides_dorei|UniRef90_X|UniRef50_X|100"
        second_reference="2__X__gene1|g__Bacteroides.s__Bacteroides_uniformis|UniRef90_X|UniRef50_X|100"

        # each sequence is covered over a different portion of its length, so neither
        # meets the threshold on its own even though together they would span 89%
        alignment_file=self.write_temp_alignments([(first_reference, 1, 45),
            (second_reference, 50, 95)])

        allowed_references = blastx_coverage.blastx_coverage(alignment_file,
            50.0, log_messages=True, nucleotide=True)

        utils.remove_temp_file(alignment_file)

        self.assertEqual(sorted(allowed_references), [])

    def test_blastx_coverage_repeated_gene_family_both_covered(self):
        """
        Test the coverage filter with the same gene family present in two species
        Test that both sequences are allowed when both are covered
        """

        first_reference="1__X__gene1|g__Bacteroides.s__Bacteroides_dorei|UniRef90_X|UniRef50_X|100"
        second_reference="2__X__gene1|g__Bacteroides.s__Bacteroides_uniformis|UniRef90_X|UniRef50_X|200"

        alignment_file=self.write_temp_alignments([(first_reference, 1, 60),
            (second_reference, 1, 160)])

        allowed_references = blastx_coverage.blastx_coverage(alignment_file,
            50.0, log_messages=True, nucleotide=True)

        utils.remove_temp_file(alignment_file)

        self.assertEqual(sorted(allowed_references), sorted([first_reference, second_reference]))
