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
import sys
import ConfigParser

# User config file
user_edit_config_file="humann2.cfg"

full_path_user_edit_config_file=os.path.join(os.path.dirname(os.path.realpath(__file__)),
    user_edit_config_file)

def update_user_edit_config_file_database_folders(uniref=None,chocophlan=None):
    """
    Update the two database folders to the user edit config file
    """
    
    config_items={}
    config_items[config_database_section]={}
    
    if uniref:
        config_items[config_database_section][config_uniref_key]=uniref
        
    if chocophlan:
        config_items[config_database_section][config_chochoplan_key]=chocophlan
    
    update_user_edit_config_file(config_items)

def update_user_edit_config_file(new_config_items):
    """
    Update the settings to the user editable config file
    """
    
    config = ConfigParser.RawConfigParser()
    
    # start with the current config settings
    config_items = read_user_edit_config_file()
    
    # update with the new config items
    for section in new_config_items:
        for name,value in new_config_items[section].items():
            config_items[section][name]=value
    
    for section in config_items:
        config.add_section(section)
        for name,value in config_items[section].items():
            value=str(value)
            if "file" in section or "folder" in section:
                # convert to absolute path if needed
                if not os.path.isabs(value):
                    value=os.path.abspath(value)
            config.set(section,name,value)
    
    try:
        file_handle=open(full_path_user_edit_config_file,"wb")
        config.write(file_handle)
        file_handle.close()
    except EnvironmentError:
        sys.exit("Unable to write to the config file: " + full_path_user_edit_config_file)
    
def read_user_edit_config_file():
    """
    Read the settings from the config file
    """
    
    config = ConfigParser.ConfigParser()
    
    try:
        config.read(full_path_user_edit_config_file)
    except EnvironmentError:
        sys.exit("Unable to read from the config file: " + full_path_user_edit_config_file)
        
    # read through all of the sections
    config_items = {}
    for section in config.sections():
        config_list = config.items(section)
        config_items[section]={}
        for name,value in config_list:
            if "file" in section or "folder" in section:
                # if not absolute path, then return absolute path relative to this folder
                if not os.path.isabs(value):
                    value=os.path.abspath(os.path.join(os.path.dirname(full_path_user_edit_config_file),value))
            config_items[section][name]=value
        
    return config_items

def get_item(config_items, section, name, type=None):
    """
    Get the item from the dictionary of section/names from the user edit config file
    """
    
    # try to obtain the value from the config dictionary
    try:
        value=config_items[section][name]
    except KeyError:
        sys.exit("CRITICAL ERROR: Unable to load value from " + full_path_user_edit_config_file +
            " . \nItem not found. \nItem should be in section (" + section + ") with name (" + name + ").")
        
    # if present, try to change the value type
    if type:
        try:
            if type == "string":
                value=str(value)
            elif type == "int":
                value=int(value)
            elif type == "float":
                value=float(value)
            elif type == "bool":
                if value in ["False","false","F","f"]:
                    value=False
                elif value in ["True","true","T","t"]:
                    value=True
                else:
                    raise ValueError
        except ValueError:
            sys.exit("CRITICAL ERROR: Unable to load value from " + full_path_user_edit_config_file +
                     " . \nItem found in section (" + section + ") with name (" + name + "). " +
                     "\nItem is not of type (" + type + ").")
        
    return value

# get the base settings from the user edit config file
config_items=read_user_edit_config_file()

# set those items that are included in the user edit config file

# database folders
chocophlan=get_item(config_items, "database_folders" , "chocophlan", "string")
uniref=get_item(config_items, "database_folders", "uniref", "string")
    
# pathways files
metacyc_gene_to_reactions=get_item(config_items, "pathways_files", "metacyc_gene_to_reactions", "string")
metacyc_reactions_to_pathways=get_item(config_items, "pathways_files", "metacyc_reactions_to_pathways", "string")
    
unipathway_database_part1=get_item(config_items, "pathways_files", "unipathway_database_part1", "string")
unipathway_database_part2=get_item(config_items, "pathways_files", "unipathway_database_part2", "string")
    
# run modes
resume=get_item(config_items, "run_modes", "resume", "bool")
verbose=get_item(config_items, "run_modes", "verbose", "bool")
bypass_prescreen=get_item(config_items, "run_modes", "bypass_prescreen", "bool")
bypass_nucleotide_index=get_item(config_items, "run_modes", "bypass_nucleotide_index", "bool")
bypass_nucleotide_search=get_item(config_items, "run_modes", "bypass_nucleotide_search", "bool")
bypass_translated_search=get_item(config_items, "run_modes", "bypass_translated_search", "bool")
    
# number of threads
threads=get_item(config_items, "run_modes", "threads", "int")
    
# evalue threshold
evalue_threshold=get_item(config_items, "alignment_settings", "evalue_threshold", "float")
    
# average read length
average_read_length=get_item(config_items, "alignment_settings", "average_read_length", "int")
    
# prescreen threshold
prescreen_threshold=get_item(config_items, "alignment_settings", "prescreen_threshold", "float")

# translated search identity threshold
identity_threshold=get_item(config_items, "alignment_settings", "identity_threshold", "float")
    
# output file decimal places
output_max_decimals=get_item(config_items, "output_format", "output_max_decimals", "int")
    
# stratified output flag
remove_stratified_output=get_item(config_items, "output_format", "remove_stratified_output", "bool")

# pathways database selection
pathways_database_choices=["metacyc","unipathway"]
pathways_database=pathways_database_choices[0]

# selected pathways
pathways_database_part1=metacyc_gene_to_reactions
pathways_database_part2=metacyc_reactions_to_pathways

# pathways settings
reactions_database_delimiter="\t"
pathways_database_delimiter="\t"

pathway_identifier="NA"
pathways_recursion=False

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

