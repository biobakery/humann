
import os

verbose=False

data_folder=os.path.join(os.path.dirname(os.path.abspath(__file__)),"data")

demo_fastq=os.path.join(data_folder,"demo.fastq")
demo_fasta=os.path.join(data_folder, "demo.fasta")
demo_sam=os.path.join(data_folder, "demo.sam")
demo_m8=os.path.join(data_folder,"demo.m8")
demo_gene_families=os.path.join(data_folder,"demo_genefamilies.tsv")
demo_bugs_list=os.path.join(data_folder, "demo_metaphlan_bugs_list.tsv")

expected_demo_output_files=["demo_genefamilies.tsv","demo_pathabundance.tsv","demo_pathcoverage.tsv"]
expected_demo_output_files_genefamilies_input=["demo_genefamilies_pathabundance.tsv","demo_genefamilies_pathcoverage.tsv"]

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
rapsearch2_output_file_without_header_coverage=os.path.join(data_folder,"rapsearch2_output_without_header_coverage_filter.m8")
rapsearch2_output_file_without_header_coverage_custom_annotations=os.path.join(data_folder,"rapsearch2_output_without_header_coverage_filter_custom_annotation.m8")
rapsearch2_output_file_without_header_coverage_chocophlan_annotations=os.path.join(data_folder,"rapsearch2_output_without_header_coverage_filter_chocophlan_annotation.m8")
rapsearch2_output_file_with_header_no_log=os.path.join(data_folder, "rapsearch2_output_with_header_no_log.m8")
coverage_id_mapping_file=os.path.join(data_folder,"id_mapping_coverage_filter.tsv")

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

multi_sample_genefamilies = os.path.join(data_folder, "multi_sample_genefamilies.tsv")
multi_sample_genefamilies_split1=os.path.join(data_folder, "multi_sample_genefamilies_sample1.tsv")
multi_sample_genefamilies_split2=os.path.join(data_folder, "multi_sample_genefamilies_sample2.tsv")
multi_sample_genefamilies_split_basename="multi_sample_genefamilies_"

# Files to test the regroup script
regroup_folder="tooltest-regroup_table"
regroup_input=os.path.join(data_folder, regroup_folder, "regroup_table-builtin_input.txt")
regroup_rxn_output=os.path.join(data_folder, regroup_folder, "regroup_table-builtin_uniref50_rxn_output.txt")
regroup_rxn_mean_output=os.path.join(data_folder, regroup_folder, "regroup_table-builtin_uniref50_rxn_mean_output.txt")
regroup_ec_output=os.path.join(data_folder, regroup_folder, "regroup_table-builtin_uniref50_ec_output.txt")
regroup_go_output=os.path.join(data_folder, regroup_folder, "regroup_table-builtin_uniref50_go_output.txt")
regroup_ko_output=os.path.join(data_folder, regroup_folder, "regroup_table-builtin_uniref50_ko_output.txt")

regroup_custom_input=os.path.join(data_folder, regroup_folder, "regroup_table-custom_input.txt")
regroup_custom_groups=os.path.join(data_folder, regroup_folder, "regroup_table-custom_groups.txt")
regroup_custom_groups_output=os.path.join(data_folder, regroup_folder, "regroup_table-custom_output.txt")

# Files to test the rename script
rename_folder="tooltest-rename_table"
rename_input=os.path.join(data_folder, rename_folder, "rename_table-input.txt")
rename_uniref50_output=os.path.join(data_folder, rename_folder, "rename_table-builtin_output.txt")

rename_pathway_input=os.path.join(data_folder, rename_folder, "pathabundance_without_names.txt")
rename_pathway_output=os.path.join(data_folder, rename_folder, "rename_table-builtin_pathways_output.txt")

rename_ec_input=os.path.join(data_folder, rename_folder, "rename_table-input_ec.txt")
rename_ec_output=os.path.join(data_folder, rename_folder, "rename_table-builtin_ec_output.txt")

rename_ko_input=os.path.join(data_folder, rename_folder, "rename_table-input_ko.txt")
rename_ko_output=os.path.join(data_folder, rename_folder, "rename_table-builtin_ko_output.txt")

rename_rxn_input=os.path.join(data_folder, rename_folder, "rename_table-input_rxn.txt")
rename_rxn_output=os.path.join(data_folder, rename_folder, "rename_table-builtin_rxn_output.txt")

rename_custom_output=os.path.join(data_folder, rename_folder, "rename_table-custom_output.txt")
rename_custom_mapping=os.path.join(data_folder, rename_folder, "rename_table-custom.txt.gz")

# Files to test the renorm script
# Output files of cpm and relab were manually verified
renorm_folder="tooltest-renorm_table"
renorm_input=os.path.join(data_folder, renorm_folder, "renorm_table-input.txt")
renorm_relab_output=os.path.join(data_folder, renorm_folder, "renorm_table-relab_output.txt")
renorm_cpm_output=os.path.join(data_folder, renorm_folder, "renorm_table-cpm_output.txt")

# Files to test the rna/dna script
rna_dna_norm_folder="tooltest-rna_dna_norm"
rna_dna_norm_dna_input=os.path.join(data_folder, rna_dna_norm_folder, "rna_dna_norm-dna.txt")
rna_dna_norm_rna_input=os.path.join(data_folder, rna_dna_norm_folder, "rna_dna_norm-rna.txt")

rna_dna_norm_file_names=["-relative_expression.tsv","-smoothed_dna.tsv","-smoothed_rna.tsv"]

rna_dna_norm_laplace_output_files=[os.path.join(data_folder, rna_dna_norm_folder, "rna_dna_norm-laplace"+file)
                                   for file in rna_dna_norm_file_names]

rna_dna_norm_witten_bell_output_files=[os.path.join(data_folder, rna_dna_norm_folder, "rna_dna_norm-witten-bell"+file)
                                   for file in rna_dna_norm_file_names]

rna_dna_norm_log_output_files=[os.path.join(data_folder, rna_dna_norm_folder, "rna_dna_norm-log"+file)
                                   for file in rna_dna_norm_file_names]

rna_dna_norm_log_10_output_files=[os.path.join(data_folder, rna_dna_norm_folder, "rna_dna_norm-log-10"+file)
                                   for file in rna_dna_norm_file_names]

# Files to test the strain profile script
strain_profile_folder="tooltest-strain_profiler"
strain_profile_input=os.path.join(data_folder, strain_profile_folder, "strain_profiler-input.txt")
strain_profile_file_names=["s1-strain_profile.tsv","s2-strain_profile.tsv"]
strain_profile_m1_n2_output_files=[os.path.join(data_folder, strain_profile_folder, file) for file in strain_profile_file_names]

# Files to test the merge abundance script
merge_abundance_folder="tooltest-merge_abundance_tables"
merge_abundance_genefamilies_input=os.path.join(data_folder, merge_abundance_folder, "demo_genefamilies.tsv")
merge_abundance_pathways_input=os.path.join(data_folder, merge_abundance_folder, "demo_pathabundance.tsv")
merge_abundance_output=os.path.join(data_folder, merge_abundance_folder, "merge_output.tsv")
merge_abundance_remove_taxonomy_output=os.path.join(data_folder, merge_abundance_folder, "merge_output_remove_taxonomy.tsv")

# Get the locations of the demo chocophlan and uniref databases
example_demo_data_folder=os.path.join(os.path.dirname(os.path.abspath(__file__)),os.pardir,"data")
chocophlan_example_demo_folder=os.path.join(example_demo_data_folder, "chocophlan_DEMO")
uniref_example_demo_folder=os.path.join(example_demo_data_folder, "uniref_DEMO")

# Get the locations of the files for the biom tests
expected_demo_output_files_biom=["demo_genefamilies.biom","demo_pathabundance.biom","demo_pathcoverage.biom"]
demo_gene_families_biom=os.path.join(data_folder,"demo_genefamilies.tsv")
expected_demo_output_files_biom_pathways=set(["PWY-6305","PWY490-3"])
regroup_input_biom=os.path.join(data_folder, regroup_folder, "regroup_table-builtin_input.biom")
rename_input_biom=os.path.join(data_folder, rename_folder, "rename_table-input.biom")
renorm_input_biom=os.path.join(data_folder, renorm_folder, "renorm_table-input.biom")
renorm_cpm_output_biom=os.path.join(data_folder, renorm_folder, "renorm_table-cpm_output.biom")

multi_sample_genefamilies_biom = os.path.join(data_folder, "multi_sample_genefamilies.biom")
multi_sample_genefamilies_split_basename_biom="multi_sample_biom_genefamilies_"
multi_sample_split_files_biom=["multi_sample_genefamilies_sample1.biom","multi_sample_genefamilies_sample2.biom"]
