"""
Configuration settings
"""

# software run modes
debug = False
verbose = False

# translated alignment options
translated_alignment_choices = ["usearch","rapsearch"]
translated_alignment_selected = translated_alignment_choices[1]

# file naming
bugs_list_name="_bugs_list.tsv"
metaphlan_bowtie2_name="_bowtie2_out.txt"

bowtie2_index_name="_bowtie2_index"

chocophlan_custom_database_name="_custom_database.ffn"
chocophlan_alignment_name="_chocophlan_align.sam"

unaligned_reads_name_no_ext="_unaligned_reads"
aligned_reads_name_tsv="_aligned_reads.tsv"

translated_alignment_name="_translated_aligned"

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
id_threshold_default=0.4
usearch_max_seqs=10000
usearch_version="v7.0.1001"
usearch_opts=[]

# rapsearch options
rapsearch_opts=["-v",-1]

# humann1 scripts
humann1_scripts="src/humann1/"
