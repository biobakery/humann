"""
Configuration settings
"""

import os

# software run modes
debug = False
verbose = False
bypass_prescreen = False
bypass_nucleotide_index = False

# translated alignment options
translated_alignment_choices = ["usearch","rapsearch"]
translated_alignment_selected = translated_alignment_choices[1]

# file naming
temp_dir=""
file_basename=""

bugs_list_name="_metaphlan_bugs_list.tsv"
metaphlan_bowtie2_name="_metaphlan_bowtie2.txt"

chocophlan_custom_database_name="_custom_chocophlan_database.ffn"
bowtie2_index_name="_bowtie2_index"
chocophlan_alignment_name="_bowtie2_aligned.sam"

nucleotide_unaligned_reads_name_no_ext="_bowtie2_unaligned"
nucleotide_aligned_reads_name_tsv="_bowtie2_aligned.tsv"

translated_alignment_name="_aligned.tsv"
translated_unaligned_reads_name_no_ext="_unaligned"

reactions_name="_reactions.tsv"
pathways_minpath="_pathways.tsv"

pathabundance_file="_pathabundance.tsv"
pathcoverage_file="_pathcoverage.tsv"
genefamilies_file="_genefamilies.tsv"

# metaphlan options
prescreen_threshold=0.01
metaphlan_opts=["-t","rel_ab"]
metaphlan_pkl_file="db_v20/mpa_v20_m200.pkl"
metaphlan_mpa_index="db_v20/mpa_v20_m200"

# bowtie2 options and threshold
bowtie2_large_index_threshold=4000000000
bowtie2_index_ext_list=[".1.bt2",".2.bt2",".3.bt2",".4.bt2",
    ".rev.1.bt2",".rev.2.bt2"]
bowtie2_large_index_ext=".1.bt2l"

bowtie2_build_opts=[]
bowtie2_align_opts=["--very-sensitive"]

# usearch options
usearch_database_extension=".udb"
id_threshold_default=0.4
usearch_max_seqs=10000
usearch_version="usearch v7.0.1001"
usearch_opts=[]

# rapsearch options
rapsearch_database_extension=".info"
rapsearch_opts=["-v",-1]

# humann1 scripts
humann1_scripts="src/humann1/"
humann1_script_bam_to_hits="bam2hits.py"
humann1_script_blast_to_hits="blast2hits.py"
humann1_script_hits_to_reactions="hits2metacyc.py"
humann1_script_minpath="MinPath1.2hmp.py"
humann1_script_enzymes_to_pathways="enzymes2pathways_mp.py"

# data files
data_folder="data/"
metacyc_gene_to_reactions="mcc"
metacyc_reactions_to_pathways="mcpc"

# MinPath
minpath_file="minpath1.2.tar.gz"
minpath_url=("http://omics.informatics.indiana.edu/mg/get.php?" + 
    "justdoit=yes&software=" + minpath_file)
minpath_folder="MinPath"
