"""
HUMAnN2: config module
Configuration settings

Copyright (c) 2014 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import os
import ConfigParser

# required python version
required_python_version_major = 2
required_python_version_minor = 7

# software run modes
resume = False
verbose = False
bypass_prescreen = False
bypass_nucleotide_index = False
bypass_nucleotide_search = False
bypass_translated_search = False

# number of threads
threads=1

# evalue threshold
evalue_threshold=1.0

# output file decimal places
output_max_decimals=10

# average read length
average_read_length=1

# log options
log_level_choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"]
log_level=log_level_choices[0]

# turn on/off computations
toggle_choices=["on","off"]
xipe_toggle = "off"
minpath_toggle = "on"
pick_frames_toggle = "off"

# file format
output_format_choices=["tsv", "biom"]
output_format=output_format_choices[0]
input_format_choices=["fastq","fastq.gz","fasta","fasta.gz","sam","bam","blastm8","genetable","biom"]

# translated alignment options
translated_alignment_choices = ["usearch","rapsearch","diamond"]
translated_alignment_selected = translated_alignment_choices[1]

# file naming
temp_dir=""
unnamed_temp_dir=""
file_basename=""
fasta_extension=".fa"

bugs_list_name="_metaphlan_bugs_list.tsv"
metaphlan_bowtie2_name="_metaphlan_bowtie2.txt"

chocophlan_custom_database_name="_custom_chocophlan_database.ffn"
bowtie2_index_name="_bowtie2_index"
chocophlan_alignment_name="_bowtie2_aligned.sam"

nucleotide_unaligned_reads_name_no_ext="_bowtie2_unaligned"
nucleotide_unaligned_reads_picked_frames_name_no_ext="_bowtie2_unaligned_picked_frames"
nucleotide_aligned_reads_name_tsv="_bowtie2_aligned.tsv"

translated_alignment_name="_aligned.tsv"
translated_unaligned_reads_name_no_ext="_unaligned"

pathabundance_file="_pathabundance"
pathcoverage_file="_pathcoverage"
genefamilies_file="_genefamilies"

# metaphlan options
prescreen_threshold=0.01
metaphlan_opts=["-t","rel_ab"]
metaphlan_pkl_file="db_v20/mpa_v20_m200.pkl"
metaphlan_mpa_index="db_v20/mpa_v20_m200"

# chocophlan formatting
chocophlan_delimiter="|"
chocophlan_bug_index=-3
chocophlan_gene_indexes=[-1]
chocophlan_location_index=-5
chocophlan_location_delimiter="-"
chocophlan_location_extra_characters="[:|c]"
chocophlan_multiple_location_delimiter=","

# uniref formatting
uniref_delimiter="|"
uniref_gene_index=-2
uniref_length_index=-1
uniref_gene_filters=["UniRef50_unknown","UniRef90_unknown"]

# bowtie2 options and threshold
bowtie2_large_index_threshold=4000000000
bowtie2_index_ext_list=[".1.bt2",".2.bt2",".3.bt2",".4.bt2",
    ".rev.1.bt2",".rev.2.bt2"]
bowtie2_large_index_ext=".1.bt2l"

bowtie2_build_opts=[]
bowtie2_align_opts=["--very-sensitive"]

#set the locations of data in the sam file
sam_read_name_index=0
sam_flag_index=1
sam_reference_index=2
sam_mapq_index=4
sam_read_index=9
sam_read_quality=10
sam_unmapped_flag=0x4
sam_delimiter="\t"

#set the locations of data in a tabulated blast formatted file
# all translated alignment files will be of the tabulated blast format
blast_delimiter="\t"
blast_query_index=0
blast_reference_index=1
blast_identity_index=2
blast_aligned_length_index=3
blast_evalue_index=10
blast_total_columns=12

# output file formats
output_file_column_delimiter="\t"
output_file_category_delimiter="|"
remove_stratified_output=False

# gene table file formats
gene_table_comment_indicator="^#"
gene_table_delimiter="\t"
gene_table_category_delimiter="|"
gene_table_gene_index=0
gene_table_value_index=1
gene_table_total_columns=2

# id mapping formats
id_mapping_comment_indicator="^#"
id_mapping_delimiter="\t"
id_mapping_reference_index=0
id_mapping_gene_index=1
id_mapping_gene_length_index=2
id_mapping_bug_index=3

# usearch options
usearch_database_extension=".udb"
identity_threshold=0.4
usearch_max_seqs=10000
usearch_version="usearch v7.0.1001"
usearch_opts=["-maxhits",20]

# rapsearch options
rapsearch_database_extension=".info"
rapsearch_opts=["-v",20]
rapsearch_output_file_extension=".m8"

# diamond options
diamond_database_extension=".dmnd"
diamond_opts=["--max-target-seqs",20,"--sensitive"]
diamond_cmmd_protein_search="blastp"
diamond_cmmd_nucleotide_search="blastx"

# data files
chocophlan=""
uniref=""

metacyc_gene_to_reactions="databases/pathways/metacyc_reactions.uniref"
metacyc_reactions_to_pathways="databases/pathways/metacyc_pathways"
pathway_identifier="PWY"
pathways_recursion=True

unipathway_database_part1="databases/pathways/unipathway_uniprots.uniref"
unipathway_database_part2="databases/pathways/unipathway_pathways"

reactions_database_delimiter="\t"
pathways_database_delimiter="\t"


# pathways databases
pathways_database_part1=metacyc_gene_to_reactions
pathways_database_part2=metacyc_reactions_to_pathways

# MinPath
minpath_file="minpath1.2.tar.gz"
minpath_url=("http://omics.informatics.indiana.edu/mg/get.php?" + 
    "justdoit=yes&software=" + minpath_file)
minpath_folder="MinPath"
minpath_script="MinPath12hmp.py"
minpath_original_script="MinPath1.2.py"
minpath_reaction_index=0
minpath_pathway_index=7
minpath_pathway_identifier="^path"
minpath_pathway_delimiter=" "
minpath_update_script="MinPath12hmp.py"

# Xipe
xipe_script="xipe.py"
xipe_delimiter="\t"
xipe_percent=str(0.1)
xipe_probability=0.9
xipe_bin=1

# User config file
user_edit_config_file="humann2.cfg"
config_database_section="database_folders"
config_uniref_key="uniref"
config_chochoplan_key="chocophlan"

def get_chochophlan_folder_location():
    """
    Return the location of the chocophlan folder from the user edit config
    """
    
    requested_keys={config_chochoplan_key:""}
    requested_keys=read_user_edit_config_file(config_database_section, requested_keys)

    return requested_keys.get(config_chochoplan_key,"")

def get_uniref_folder_location():
    """
    Return the location of the uniref folder from the user edit config
    """
    
    requested_keys={config_uniref_key:""}
    requested_keys=read_user_edit_config_file(config_database_section, requested_keys)

    return requested_keys.get(config_uniref_key,"")

def write_user_edit_config_file_database_folders(uniref=None,chocophlan=None):
    """
    Write the two database folders to the user edit config file
    """
    
    # if one of the keys is not provided, find the setting from the config file
    if not uniref or not chocophlan:
        requested_keys={config_uniref_key:"",config_chochoplan_key:""}
        requested_keys = read_user_edit_config_file(config_database_section, requested_keys)
    
        # if uniref not provided, set to current config file
        if not uniref:
            uniref=requested_keys[config_uniref_key]
        
        # if chocophlan not provided, set to the current config file
        if not chocophlan:
            chocophlan=requested_keys[config_chochoplan_key]
    
    new_settings={}
    new_settings[config_database_section]={}
    new_settings[config_database_section][config_uniref_key]=uniref
    new_settings[config_database_section][config_chochoplan_key]=chocophlan
    
    write_user_edit_config_file(new_settings)

def write_user_edit_config_file(new_settings):
    """
    Write the settings to the user editable config file
    """
    
    config = ConfigParser.RawConfigParser()
    
    for section in new_settings:
        config.add_section(section)
        for key,value in new_settings[section].items():
            value=str(value)
            if "file" in section or "folder" in section:
                # convert to absolute path if needed
                if not os.path.isabs(value):
                    value=os.path.abspath(value)
            config.set(section,key,str(value))
            
    full_user_edit_config_file=os.path.join(os.path.dirname(os.path.realpath(__file__)),
        user_edit_config_file)
    
    try:
        file_handle=open(full_user_edit_config_file,"wb")
        config.write(file_handle)
        file_handle.close()
    except EnvironmentError:
        sys.exit("Unable to write to the config file: " + full_user_edit_config_file)
    
def read_user_edit_config_file(section,keys):
    """
    Read the settings from the config file
    """
    
    config = ConfigParser.ConfigParser()
    
    full_user_edit_config_file=os.path.join(os.path.dirname(os.path.realpath(__file__)),
        user_edit_config_file)
    
    try:
        config.read(full_user_edit_config_file)
    except EnvironmentError:
        sys.exit("Unable to read from the config file: " + full_user_edit_config_file)
        
    for key in keys:
        value=str(config.get(section,key))
        if "file" in section or "folder" in section:
            # if not absolute path, then return absolute path relative to this folder
            if not os.path.isabs(value):
                value=os.path.abspath(os.path.join(os.path.dirname(full_user_edit_config_file),value))
        keys[key]=value
        
    return keys

def get_humann2_base_directory():
    """ 
    Return the location of the humann2 base directory
    """
    
    config_file_location=os.path.dirname(os.path.realpath(__file__))
    
    # The humann2 base directory is parent directory of the config file location
    humann2_base_directory=os.path.abspath(os.path.join(config_file_location,os.pardir))
    
    return humann2_base_directory
