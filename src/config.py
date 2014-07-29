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
chocophlan_alignment_name="_chocophlan_align.sam"
unaligned_reads_name_no_ext="_unaligned_reads"
aligned_reads_name_tsv="_aligned_reads.tsv"
translated_alignment_name="_translated_aligned"
usearch_name_temp=".usearch.tmp"


# metaphlan options
metaphlan_opts=["-t","rel_ab"]
metaphlan_pkl_file="db_v20/mpa_v20_m200.pkl"
metaphlan_mpa_index="db_v20/mpa_v20_m200"

chocophlan_custom_database_name="_custom_database.ffn"

# bowtie2 options and threshold
bowtie2_large_index_threshold=4000000000
bowtie2_index_ext=".1.bt2"
bowtie2_large_index_ext=".1.bt2l"

bowtie2_build_opts=[]
bowtie2_align_opts=["--very-sensitive"]

# usearch options
usearch_opts=[]
