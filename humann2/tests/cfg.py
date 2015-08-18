
import os

verbose=False

data_folder=os.path.join(os.path.dirname(os.path.abspath(__file__)),"data")

small_fasta_file=os.path.join(data_folder,"file.fasta")
small_fasta_file_total_sequences=3
small_fastq_file=os.path.join(data_folder,"file.fastq")
small_fastq_file_total_sequences=2

convert_fastq_file=os.path.join(data_folder,"convert_file.fastq")
convert_fastq_at_character_file=os.path.join(data_folder,"convert_file_at_character.fastq")
convert_fasta_file=os.path.join(data_folder,"convert_file.fasta")
convert_fasta_multiline_file=os.path.join(data_folder,"convert_file_multiline.fasta")
convert_fasta_pick_frames_file=os.path.join(data_folder,"convert_file_pick_frames.fasta")

small_fasta_file_single_line_sequences=os.path.join(data_folder,"file_single_line_sequences.fasta")

# The pathways dat file is a test file for the script that creates the structured pathways file
# The result of running this script should be the file named "test_pathways_results.tsv"
# This file is in the structured format expected by the humann2 PathwaysDatabase store class
pathways_file_dat=os.path.join(data_folder, "test_pathways.dat")
pathways_file=os.path.join(data_folder,"test_pathways_results.tsv")
pathways_flat_file=os.path.join(data_folder, "test_pathways_flat.tsv")

reactions_file=os.path.join(data_folder, "reactions.tsv")

usearch_file=os.path.join(data_folder, "usearch_output.tsv")
usearch_uniref50_file=os.path.join(data_folder, "usearch_with_uniref50_output.tsv")
usearch_file_bug_list=["bug1","bug2","bug3","unclassified"]
usearch_file_gene_list=["UniRef50_unknown","gene1","gene2","gene3","gene4","gene5"]

gene_familes_file=os.path.join(data_folder, "gene_families.tsv")
gene_familes_uniref50_with_names_file=os.path.join(data_folder, "gene_families_uniref50_with_names.tsv")
larger_gene_families_uniref50_with_names_file=os.path.join(data_folder, "demo_genefamilies_with_names.tsv")
gene_families_to_names_file=os.path.join(data_folder, "small_map_uniref50_name.txt.gz")

demo_pathabundance_file=os.path.join(data_folder, "demo_pathabundance_with_names.tsv")
demo_pathcoverage_file=os.path.join(data_folder, "demo_pathcoverage_with_names.tsv")

sam_file_with_header=os.path.join(data_folder, "file_with_header.sam")
sam_file_without_header=os.path.join(data_folder, "file_without_header.sam")
sam_file_without_header_with_tags=os.path.join(data_folder, "file_without_header_with_tags.sam")
sam_file_unaligned_reads=os.path.join(data_folder,"2_aligned_3_unaligned.sam")
sam_file_unaligned_reads_total_aligned=2
sam_file_unaligned_reads_total_unaligned=3

bam_file=os.path.join(data_folder,"file.bam")

biom_file=os.path.join(data_folder,"genefamilies.biom")
genetable_file=os.path.join(data_folder,"genefamilies_biom_match.tsv")

genetable_file_bug_scores={
    "s__Bacteroides_xylanisolvens": {"UniRef50_A0A015Q7T6": 6.80272108844 },
    "s__Parabacteroides_merdae" : { "UniRef50_A7ALT1" : 5.74712643678 },
    "s__Bacteroides_uniformis" : { "UniRef50_Q5L8C7": 5.64971751412 },
    "s__Faecalibacterium_prausnitzii" : { "UniRef50_R7EH44" : 5.64971751412 },
    "all" : {
        "UniRef50_A0A015Q7T6" : 6.80272108844,
        "UniRef50_A7ALT1" : 5.74712643678,
        "UniRef50_Q5L8C7" : 5.64971751412,
        "UniRef50_R7EH44" : 5.64971751412,
        "UniRef50_B0NLN0" : 5.26315789474 }
    }

sam_file_annotations=os.path.join(data_folder,"annotations.sam")
rapsearch_file_annotations=os.path.join(data_folder,"annotations.m8")

rapsearch2_output_file_with_header=os.path.join(data_folder, "rapsearch2_output_with_header.m8")
rapsearch2_output_file_without_header=os.path.join(data_folder, "rapsearch2_output_without_header.m8")
rapsearch2_output_file_with_header_no_log=os.path.join(data_folder, "rapsearch2_output_with_header_no_log.m8")

id_mapping_file=os.path.join(data_folder, "id_mapping.tsv")
id_mapping_gene_table_file=os.path.join(data_folder, "id_mapping_gene_table.tsv")

genetable_file_bug_scores_id_mapping={
    "bug1": {"gene1": 6.80272108844 },
    "s__Parabacteroides_merdae" : { "gene2" : 5.74712643678 },
    "s__Bacteroides_uniformis" : { "UniRef50_Q5L8C7": 5.64971751412 },
    "s__Faecalibacterium_prausnitzii" : { "UniRef50_R7EH44" : 5.64971751412 },
    "all" : {
        "gene1" : 6.80272108844,
        "gene2" : 5.74712643678,
        "gene5" : 6.80272108844,
        "UniRef50_A7ALT1" : 5.74712643678,
        "UniRef50_Q5L8C7" : 5.64971751412,
        "UniRef50_R7EH44" : 5.64971751412,
        "UniRef50_B0NLN0" : 5.26315789474 }
    }

