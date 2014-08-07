"""
Configuration settings
"""

# software run modes
debug = False
verbose = False

# temp file naming
bugs_list_name="_bugs_list.tsv"
metaphlan_bowtie2_name="_bowtie2_out.txt"

bowtie2_index_name="_bowtie2_index"

chocophlan_custom_database_name="_custom_database.ffn"
chocophlan_alignment_name="_chocophlan_align.sam"

unaligned_reads_name_no_ext="_unaligned_reads"
aligned_reads_name_tsv="_aligned_reads.tsv"

translated_alignment_name="_translated_aligned"

# metaphlan options
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
usearch_max_seqs=10000
usearch_version="v7.0.1001"
usearch_opts=["-target_cov","0.9","-query_cov","0.9","-match",
    "1.0","-mismatch","-1.0","-gapopen","3.0I/0.0E","-gapext",
    "0.5I/0.0E"]

# humann1 scripts
humann1_scripts="src/humann1/"
