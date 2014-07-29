"""
Configuration settings
"""

# software run modes
debug = False
verbose = False

# metaphlan options
metaphlan_opts=["-t","rel_ab"]

# bowtie2 options and threshold
bowtie2_large_index_threshold=4000000000

bowtie2_build_opts=[]
bowtie2_align_opts=["--very-sensitive"]

# usearch options
usearch_opts=[]
