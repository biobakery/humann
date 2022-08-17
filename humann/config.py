"""
HUMAnN: config module
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

# try to import the python2 ConfigParser
# if unable to import, then try to import the python3 configparser
try:
    import ConfigParser as configparser
except ImportError:
    import configparser

import logging

# name global logging instance
logger=logging.getLogger(__name__)

def log_settings():
    """
    Write to the log file the config settings for the run
    """
    
    lines=[]    
    lines.append("DATABASE SETTINGS")
    lines.append("nucleotide database folder = " + nucleotide_database)
    lines.append("protein database folder = " + protein_database)
    if pathways_database_part1:
        lines.append("pathways database file 1 = " + pathways_database_part1)
        lines.append("pathways database file 2 = " + pathways_database_part2)
    else:
        lines.append("pathways database file = " + pathways_database_part2)
    lines.append("utility mapping database folder = " + utility_mapping_database)
    lines.append("")
    
    lines.append("RUN MODES")
    lines.append("resume = " + str(resume))
    lines.append("verbose = " + str(verbose))
    lines.append("bypass prescreen = " + str(bypass_prescreen))
    lines.append("bypass nucleotide index = " + str(bypass_nucleotide_index))
    lines.append("bypass nucleotide search = " + str(bypass_nucleotide_search))
    lines.append("bypass translated search = " + str(bypass_translated_search))
    lines.append("translated search = " + translated_alignment_selected)
    lines.append("threads = " + str(threads))
    lines.append("")
    
    lines.append("SEARCH MODE")
    lines.append("search mode = " + search_mode)
    lines.append("nucleotide identity threshold = " + str(nucleotide_identity_threshold))
    lines.append("translated identity threshold = " + str(identity_threshold))
    lines.append("")
    
    lines.append("ALIGNMENT SETTINGS")
    lines.append("bowtie2 options = " + str(" ".join(map(str,bowtie2_align_opts))))
    lines.append("diamond options = " + str(" ".join(map(str,diamond_opts))))
    lines.append("evalue threshold = " + str(evalue_threshold))
    lines.append("prescreen threshold = " + str(prescreen_threshold))
    lines.append("translated subject coverage threshold = " + str(translated_subject_coverage_threshold))
    lines.append("translated query coverage threshold = " + str(translated_query_coverage_threshold))
    lines.append("nucleotide subject coverage threshold = " + str(nucleotide_subject_coverage_threshold))
    lines.append("nucleotide query coverage threshold = " + str(nucleotide_query_coverage_threshold))
    lines.append("")
    
    lines.append("PATHWAYS SETTINGS")
    lines.append("minpath = " + minpath_toggle)
    lines.append("xipe = " + xipe_toggle)
    lines.append("gap fill = " + gap_fill_toggle)
    lines.append("")    
    
    lines.append("INPUT AND OUTPUT FORMATS")
    lines.append("input file format = " + input_format)
    lines.append("output file format = " + output_format)
    lines.append("output max decimals = " + str(output_max_decimals))
    lines.append("remove stratified output = " + str(remove_stratified_output))
    lines.append("remove column description output = " + str(remove_column_description_output))
    lines.append("log level = " + log_level)
    lines.append("")
    
    logger.info("\nRun config settings: \n\n" + "\n".join(lines))

# User config file
user_edit_config_file="humann.cfg"

full_path_user_edit_config_file=os.path.join(os.path.dirname(os.path.abspath(__file__)),
    user_edit_config_file)

def update_user_edit_config_file_single_item(section,name,value):
    """
    Update the settings to the user editable config file for one item
    """
    
    new_config_items={section:{ name : value }}
    
    update_user_edit_config_file(new_config_items)
    
    print("HUMAnN configuration file updated: "+ section + " : " + name + " = " + str(value))

def update_user_edit_config_file(new_config_items):
    """
    Update the settings to the user editable config file
    """
    
    config = configparser.RawConfigParser()
    
    # start with the current config settings
    config_items = read_user_edit_config_file()
    
    # update with the new config items
    for section in new_config_items:
        for name,value in new_config_items[section].items():
            if section in config_items:
                if name in config_items[section]:
                    config_items[section][name]=value
                else:
                    sys.exit("ERROR: Unable to add new name ( " + name + 
                        " ) to existing section ( " + section + " ) to " +
                        " config file: " + full_path_user_edit_config_file)
            else:
                sys.exit("ERROR: Unable to add new section ( " + section + 
                    " ) to config file: " + full_path_user_edit_config_file)
    
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
        file_handle=open(full_path_user_edit_config_file,"wt")
        config.write(file_handle)
        file_handle.close()
    except EnvironmentError:
        sys.exit("Unable to write to the HUMAnN config file.")
    
def read_user_edit_config_file():
    """
    Read the settings from the config file
    """
    
    config = configparser.ConfigParser()
    
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
nucleotide_database=get_item(config_items, "database_folders" , "nucleotide", "string")
protein_database=get_item(config_items, "database_folders", "protein", "string")
utility_mapping_database=get_item(config_items, "database_folders", "utility_mapping", "string")
        
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
    
# prescreen threshold
prescreen_threshold=get_item(config_items, "alignment_settings", "prescreen_threshold", "float")

# nucletide search identity threshold
nucleotide_identity_threshold = 0.0

# translated search identity threshold
identity_threshold_uniref90_mode = 80.0
identity_threshold_uniref50_mode = 50.0

# translated search coverage thresholds
translated_subject_coverage_threshold=get_item(config_items, "alignment_settings", "translated_subject_coverage_threshold", "float")
translated_query_coverage_threshold=get_item(config_items, "alignment_settings", "translated_query_coverage_threshold", "float")
    
# nucleotide search coverage thresholds
nucleotide_subject_coverage_threshold=get_item(config_items, "alignment_settings", "nucleotide_subject_coverage_threshold", "float")
nucleotide_query_coverage_threshold=get_item(config_items, "alignment_settings", "nucleotide_query_coverage_threshold", "float")

# output file decimal places
output_max_decimals=get_item(config_items, "output_format", "output_max_decimals", "int")
    
# stratified output flag
remove_stratified_output=get_item(config_items, "output_format", "remove_stratified_output", "bool")

# column description flag
remove_column_description_output=get_item(config_items, "output_format", "remove_column_description_output", "bool")

# pathways files
humann_install_directory=os.path.dirname(os.path.abspath(__file__))
metacyc_gene_to_reactions=os.path.abspath(os.path.join(humann_install_directory,"data","pathways","metacyc_reactions_level4ec_only.uniref.bz2"))
metacyc_reactions_to_pathways=os.path.abspath(os.path.join(humann_install_directory,"data","pathways","metacyc_pathways_structured_filtered_v24"))
    
unipathway_database_part1=os.path.abspath(os.path.join(humann_install_directory,"data","pathways","unipathway_uniprots.uniref.bz2"))
unipathway_database_part2=os.path.abspath(os.path.join(humann_install_directory,"data","pathways","unipathway_pathways"))

# pathways and gene families name mapping files
gene_family_name_mapping_file=os.path.abspath(os.path.join(humann_install_directory,"data","misc","map_uniref50_name.txt.bz2"))
pathway_name_mapping_file=os.path.abspath(os.path.join(humann_install_directory,"data","misc","map_metacyc-pwy_name.txt.gz"))
name_mapping_file_delimiter="\t"
name_mapping_join=": "

# pathways database selection
pathways_database_choices=["metacyc","unipathway"]
pathways_database=pathways_database_choices[0]

# selected pathways
pathways_database_part1=metacyc_gene_to_reactions
pathways_database_part2=metacyc_reactions_to_pathways
pathways_ec_column=True

# pathways settings
reactions_database_delimiter="\t"
pathways_database_delimiter="\t"
pathways_database_stucture_delimiter=" "
pathway_AND="+"
pathway_OR=","
pathway_reaction_optional="-"

# sgb association file
sgb_to_species_file=os.path.abspath(os.path.join(humann_install_directory,"data","misc","mpa_vJan21_CHOCOPhlAnSGB_202103.tsv"))

# query annotation
query_length_annotation_delimiter="|"

# memory use
memory_use_options=["minimum","maximum"]
memory_use=memory_use_options[0]

# log options
log_level_choices=["DEBUG","INFO","WARNING","ERROR","CRITICAL"]
log_level=log_level_choices[0]

# turn on/off computations
toggle_choices=["on","off"]
xipe_toggle = "off"
minpath_toggle = "on"
gap_fill_toggle = "on"
pick_frames_toggle = "off"

# file format
output_format_choices=["tsv", "biom"]
output_format=output_format_choices[0]
input_format_choices=["fastq","fastq.gz","fasta","fasta.gz","sam","bam","blastm8","genetable","biom"]
input_format=""

# translated alignment options
translated_alignment_choices = ["usearch","rapsearch","diamond"]
translated_alignment_selected = translated_alignment_choices[2]

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
metaphlan_version={
    "flag" : "--version",
    "major" : 3,
    "minor" : 0,
    "line" : -1,
    "column" : 2}

metaphlan_v3_db_version="v3"
metaphlan_v4_db_version="vJan21"
metaphlan_v3_db_matching_uniref="v201901_v31"

matching_uniref="201901b"

# chocophlan formatting
chocophlan_delimiter="|"
chocophlan_bug_index=1
chocophlan_bug_species_index=-1
chocophlan_bug_genera_index=-2
chocophlan_gene_indexes=[3]
chocophlan_gene_indexes_uniref90_mode=[2]
chocophlan_gene_indexes_uniref50_mode=[3]
chocophlan_location_delimiter="-"
chocophlan_location_extra_characters="[:|c]"
chocophlan_multiple_location_delimiter=","
chocophlan_length_index=4

# uniref formatting
uniref_delimiter="|"
uniref_gene_index=-2
uniref_length_index=-1

# bowtie2 options and threshold
bowtie2_large_index_threshold=4000000000
bowtie2_index_ext_list=[".1.bt2",".2.bt2",".3.bt2",".4.bt2",
    ".rev.1.bt2",".rev.2.bt2"]
bowtie2_large_index_ext=".1.bt2l"
bowtie2_version={
    "flag" : "--version",
    "major" : 2,
    "minor" : 2,
    "line" : 0,
    "column" : 2}

bowtie2_build_opts=[]
bowtie2_align_opts=["--very-sensitive"]

#set the locations of data in the sam file
sam_read_name_index=0
sam_flag_index=1
sam_reference_index=2
sam_pos_index=3
sam_mapq_index=4
sam_cigar_index=5
sam_tlen_index=8
sam_read_index=9
sam_start_optional_index=11

sam_read_quality=10
sam_unmapped_flag=0x4
sam_delimiter="\t"

sam_cigar_match_mismatch_indel_identifiers=["M","=","X","I","D"]
sam_cigar_add_to_reference_identifiers=["M","D","N","=","X"]
sam_md_field_identifier="MD:Z:"

#set the locations of data in a tabulated blast formatted file
# all translated alignment files will be of the tabulated blast format
blast_delimiter="\t"
blast_query_index=0
blast_reference_index=1
blast_identity_index=2
blast_aligned_length_index=3
blast_query_start_index=6
blast_query_end_index=7
blast_subject_start_index=8
blast_subject_end_index=9
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
usearch_version={
    "flag" : "--version",
    "major" : 7,
    "minor" : 0,
    "line" : 0,
    "column" : 1}
usearch_opts=["-maxhits",20]

# rapsearch options
rapsearch_database_extension=".info"
rapsearch_opts=["-v",20]
rapsearch_output_file_extension=".m8"
rapsearch_version={
    "flag" : "-v",
    "major" : 2,
    "minor" : 21,
    "line" : 1,
    "column" : 1}

# diamond options
diamond_database_extension=".dmnd"
diamond_options_custom=False
diamond_opts_uniref50=["--top","1","--sensitive","--outfmt","6"]
diamond_opts_uniref90=["--top","1","--outfmt","6"]
diamond_cmmd_protein_search="blastp"
diamond_cmmd_nucleotide_search="blastx"
diamond_version={
    "flag" : "--version",
    "major" : 0,
    "minor" : 9,
    "second minor" : 36,
    "line" : 0,
    "column" : 2}

# MinPath
minpath_script="MinPath12hmp.py"
minpath_reaction_index=0
minpath_pathway_index=7
minpath_pathway_identifier="^path"
minpath_pathway_delimiter=" "

# Xipe
xipe_script="xipe.py"
xipe_delimiter="\t"
xipe_percent=str(0.1)
xipe_probability=0.9
xipe_bin=1

# Alignment Score defaults
default_reference_length=1000
match_power=2

# Output file tags
unmapped_gene_name = "UNMAPPED"
unmapped_pathway_name = "UNMAPPED"
unintegrated_pathway_name = "UNINTEGRATED"

# Max arguments
max_arguments=500

# set the search mode options
search_mode_uniref90="uniref90"
search_mode_uniref50="uniref50"
search_mode=search_mode_uniref90

