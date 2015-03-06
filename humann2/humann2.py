#!/usr/bin/env python

"""
HUMAnN2 : HMP Unified Metabolic Analysis Network 2

HUMAnN2 is a pipeline for efficiently and accurately determining 
the coverage and abundance of microbial pathways in a community 
from metagenomic data. Sequencing a metagenome typically produces millions 
of short DNA/RNA reads.

Dependencies: MetaPhlAn2, ChocoPhlAn, Bowtie2, Rapsearch2 or Usearch

To Run: humann2 -i <input.fastq> -o <output_dir>

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

import sys
    
# Try to load one of the humann2 src modules to check the installation
try:
    from . import check
except ImportError:
    sys.exit("CRITICAL ERROR: Unable to find the HUMAnN2 python package." +
        " Please check your install.") 

# Check the python version
check.python_version()
    
import argparse
import subprocess
import os
import time
import tempfile
import re
import logging  

from . import config
from . import store
from . import utilities
from .search import prescreen
from .search import nucleotide
from .search import translated
from .quantify import families
from .quantify import modules

# name global logging instance
logger=logging.getLogger(__name__)

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "HUMAnN2 : HMP Unified Metabolic Analysis Network 2\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-v","--verbose", 
        help="additional output is printed\n", 
        action="store_true",
        default=config.verbose)
    parser.add_argument(
        "-r","--resume", 
        help="bypass commands if the output files exist\n", 
        action="store_true",
        default=config.resume)
    parser.add_argument(
        "--bypass-prescreen", 
        help="bypass the prescreen step and run on the full ChocoPhlAn database\n", 
        action="store_true",
        default=config.bypass_prescreen)
    parser.add_argument(
        "--bypass-nucleotide-index", 
        help="bypass the nucleotide index step and run on the indexed ChocoPhlAn database\n", 
        action="store_true",
        default=config.bypass_nucleotide_index)
    parser.add_argument(
        "--bypass-translated-search", 
        help="bypass the translated search step\n", 
        action="store_true",
        default=config.bypass_translated_search)
    parser.add_argument(
        "--bypass-nucleotide-search", 
        help="bypass the nucleotide search steps\n", 
        action="store_true",
        default=config.bypass_nucleotide_search)
    parser.add_argument(
        "-i", "--input", 
        help="input file of type {" +",".join(config.input_format_choices)+ "} \n[REQUIRED]", 
        metavar="<input.fastq>", 
        required=True)
    parser.add_argument(
        "-o", "--output", 
        help="directory to write output files\n[REQUIRED]", 
        metavar="<output>", 
        required=True)
    parser.add_argument(
        "-c", "--chocophlan",
        help="directory containing the ChocoPhlAn database\n[DEFAULT: " 
            + config.chocophlan + "]", 
        metavar="<chocophlan>")
    parser.add_argument(
        "--chocophlan-gene-index",
        help="the index of the gene in the sequence annotation\n[DEFAULT: " 
            + ",".join(str(i) for i in config.chocophlan_gene_indexes) + "]", 
        metavar="<"+",".join(str(i) for i in config.chocophlan_gene_indexes)+">",
        default=",".join(str(i) for i in config.chocophlan_gene_indexes))
    parser.add_argument(
        "-u", "--uniref",
        help="directory containing the UniRef database\n[DEFAULT: " 
            + config.uniref + "]", 
        metavar="<uniref>")
    parser.add_argument(
        "--average-read-length", 
        help="the average length of the reads\n[DEFAULT: " + str(config.average_read_length) + "]", 
        metavar="<" + str(config.average_read_length) + ">", 
        type=int,
        default=config.average_read_length) 
    parser.add_argument(
        "--evalue", 
        help="the evalue threshold\n[DEFAULT: " + str(config.evalue_threshold) + "]", 
        metavar="<" + str(config.evalue_threshold) + ">", 
        type=float,
        default=config.evalue_threshold) 
    parser.add_argument(
        "--metaphlan",
        help="directory containing the MetaPhlAn software\n[DEFAULT: $PATH]", 
        metavar="<metaplhan>")
    parser.add_argument(
        "--o-log", 
        help="log file\n" + 
        "[DEFAULT: temp/sample.log]", 
        metavar="<sample.log>")
    parser.add_argument(
        "--log-level", 
        help="level of messages to display in log\n" + 
        "[DEFAULT: " + config.log_level + "]", 
        default=config.log_level,
        choices=config.log_level_choices)
    parser.add_argument(
        "--remove-temp-output", 
        help="remove temp output files\n" + 
            "[DEFAULT: temp files are not removed]", 
        action="store_true")
    parser.add_argument(
        "--bowtie2",
        help="directory of the bowtie2 executable\n[DEFAULT: $PATH]", 
        metavar="<bowtie2>")
    parser.add_argument(
        "--threads", 
        help="number of threads/processes\n[DEFAULT: " + str(config.threads) + "]", 
        metavar="<" + str(config.threads) + ">", 
        type=int,
        default=config.threads) 
    parser.add_argument(
        "--prescreen-threshold", 
        help="minimum percentage of reads matching a species\n[DEFAULT: "
            + str(config.prescreen_threshold) + "]", 
        metavar="<" + str(config.prescreen_threshold) + ">", 
        type=float,
        default=config.prescreen_threshold) 
    parser.add_argument(
        "--identity-threshold", 
        help="identity threshold to use with the translated search\n[DEFAULT: " 
            + str(config.identity_threshold) + "]", 
        metavar="<" + str(config.identity_threshold) + ">", 
        type=float,
        default=config.identity_threshold) 
    parser.add_argument(
        "--usearch", 
        help="directory containing the usearch executable\n[DEFAULT: $PATH]", 
        metavar="<usearch>")
    parser.add_argument(
        "--rapsearch", 
        help="directory containing the rapsearch executable\n[DEFAULT: $PATH]", 
        metavar="<rapsearch>")
    parser.add_argument(
        "--diamond", 
        help="directory containing the diamond executable\n[DEFAULT: $PATH]", 
        metavar="<diamond>")
    parser.add_argument(
        "--taxonomic-profile", 
        help="a taxonomic profile (the output file created by metaphlan)\n[DEFAULT: file will be created]", 
        metavar="<taxonomic_profile.tsv>")
    parser.add_argument(
        "--id-mapping", 
        help="id mapping file for alignments\n[DEFAULT: alignment reference used]", 
        metavar="<id_mapping.tsv>")
    parser.add_argument(
        "--translated-alignment", 
        help="software to use for translated alignment\n[DEFAULT: " + 
            config.translated_alignment_selected + "]", 
        default=config.translated_alignment_selected,
        choices=config.translated_alignment_choices)
    parser.add_argument(
        "--xipe",
        help="turn on/off the xipe computation\n[DEFAULT: " +
        config.xipe_toggle + "]",
        default=config.xipe_toggle,
        choices=config.toggle_choices)
    parser.add_argument(
        "--minpath",
        help="turn on/off the minpath computation\n[DEFAULT: " + 
        config.minpath_toggle + "]",
        default=config.minpath_toggle,
        choices=config.toggle_choices)
    parser.add_argument(
        "--pick-frames",
        help="turn on/off the pick_frames computation\n[DEFAULT: " + 
        config.pick_frames_toggle + "]",
        default=config.pick_frames_toggle,
        choices=config.toggle_choices)
    parser.add_argument(
        "--output-format",
        help="the format of the output files\n[DEFAULT: " +
        config.output_format + "]",
        default=config.output_format,
        choices=config.output_format_choices)
    parser.add_argument(
        "--output-basename",
        help="the basename for the output files\n[DEFAULT: " +
        "input file basename]",
        default=config.file_basename,
        metavar="<sample_name>")
    parser.add_argument(
        "--remove-stratified-output", 
        help="remove stratification from output\n" + 
            "[DEFAULT: output is stratified]", 
        action="store_true",
        default=config.remove_stratified_output)
    parser.add_argument(
        "--input-format",
        help="the format of the input file\n[DEFAULT: format identified by software]",
        choices=config.input_format_choices)
    parser.add_argument(
        "--pathways-databases",
        help="the two mapping files to use for pathway computations\n[DEFAULT: " +
        config.pathways_database + " databases ]",
        metavar=("<pathways_database_part1.tsv>","<pathways_database_part2.tsv>"),
        nargs=2)
    parser.add_argument(
        "--pathways",
        help="the database to use for pathway computations\n[DEFAULT: " +
        config.pathways_database + "]",
        default=config.pathways_database,
        choices=config.pathways_database_choices)

    return parser.parse_args()
	
	 
def update_configuration(args):
    """
    Update the configuration settings based on the arguments
    """
    
    # Use the full path to the input file
    args.input=os.path.abspath(args.input)

    # If set, append paths executable locations
    if args.metaphlan:
        utilities.add_exe_to_path(os.path.abspath(args.metaphlan))    
    
    if args.bowtie2:
        utilities.add_exe_to_path(os.path.abspath(args.bowtie2))
    
    if args.usearch:
        utilities.add_exe_to_path(os.path.abspath(args.usearch))

    if args.rapsearch:
        utilities.add_exe_to_path(os.path.abspath(args.rapsearch))
        
    if args.diamond:
        utilities.add_exe_to_path(os.path.abspath(args.diamond))
 
    # Set the pathways database selection
    if args.pathways == "metacyc":
        config.pathways_database_part1=config.metacyc_gene_to_reactions
        config.pathways_database_part2=config.metacyc_reactions_to_pathways
        config.pathway_identifier="PWY"
        config.pathways_recursion=True
    elif args.pathways == "unipathway":
        config.pathways_database_part1=config.unipathway_database_part1
        config.pathways_database_part2=config.unipathway_database_part2
        config.pathways_recursion=False
 
    # Set the locations of the pathways databases
    # If provided by the user, this will take precedence over the pathways database selection
    if args.pathways_databases:
        config.pathways_database_part1=os.path.abspath(args.pathways_databases[0])
        config.pathways_database_part2=os.path.abspath(args.pathways_databases[1])
        config.pathways_recursion=False
        
    # Set the locations of the other databases
    if args.chocophlan:
        config.chocophlan=os.path.abspath(args.chocophlan)
        
    if args.uniref:
        config.uniref=os.path.abspath(args.uniref)

    # if set, update the config run mode to resume
    if args.resume:
        config.resume=True  
            
    # if set, update the config run mode to verbose
    if args.verbose:
        config.verbose=True  
   
    # if set, update the config run mode to bypass prescreen step
    if args.bypass_prescreen:
        config.bypass_prescreen=True  
    
    # if set, update the config run mode to bypass nucleotide index step
    if args.bypass_nucleotide_index:
        config.bypass_nucleotide_index=True  
        config.bypass_prescreen=True  
        
    # if set, update the config run mode to bypass translated search step
    # set the pick_frames toggle based on the bypass
    if args.bypass_translated_search:
        config.bypass_translated_search=True
        config.pick_frames_toggle="off"
    else:
        config.pick_frames_toggle=args.pick_frames
        
    # if set, update the config run mode to bypass nucleotide search steps
    if args.bypass_nucleotide_search:
        config.bypass_prescreen=True
        config.bypass_nucleotide_index=True
        config.bypass_nucleotide_search=True
        
    # Update thresholds
    config.prescreen_threshold=args.prescreen_threshold
    config.identity_threshold=args.identity_threshold
    
    # Update threads
    config.threads=args.threads
    
    # Update the average read length
    config.average_read_length=args.average_read_length
    
    # Update the evalue threshold
    config.evalue_threshold=args.evalue
    
    # Update translated alignment software
    config.translated_alignment_selected=args.translated_alignment
        
    # Update the computation toggle choices
    config.xipe_toggle=args.xipe
    config.minpath_toggle=args.minpath
    
    # Update the chocophlan gene indexes
    config.chocophlan_gene_indexes=[]
    for index in args.chocophlan_gene_index.split(","):
        # Look for array range
        if ":" in index:
            split_index=index.split(":")
            start=split_index[0]
            end=split_index[1]
            
            try:
                start=int(start)
                end=int(end)
                config.chocophlan_gene_indexes+=range(start,end)
            except ValueError:
                pass
        else:
            # Convert to int
            try:
                index=int(index)
                config.chocophlan_gene_indexes.append(index)
            except ValueError:
                pass
    
    # Check that the input file exists and is readable
    if not os.path.isfile(args.input):
        sys.exit("CRITICAL ERROR: Can not find input file selected: "+ args.input)
        
    if not os.access(args.input, os.R_OK):
        sys.exit("CRITICAL ERROR: Not able to read input file selected: " + args.input)
        
    # Update the output format
    config.remove_stratified_output=args.remove_stratified_output
     
    # Check that the output directory is writeable
    output_dir = os.path.abspath(args.output)
    
    if not os.path.isdir(output_dir):
        try:
            print("Creating output directory: " + output_dir)
            os.mkdir(output_dir)
        except EnvironmentError:
            sys.exit("CRITICAL ERROR: Unable to create output directory.")
    
    if not os.access(output_dir, os.W_OK):
        sys.exit("CRITICAL ERROR: The output directory is not " + 
            "writeable. This software needs to write files to this directory.\n" +
            "Please select another directory.")
        
    print("Output files will be written to: " + output_dir) 
    
    # Set the basename of the output files if specified as an option
    if args.output_basename:
        config.file_basename=args.output_basename
    else:
        # Determine the basename of the input file to use as output file basename
        input_file_basename=os.path.basename(args.input)
        # Remove gzip extension if present
        if re.search('.gz$',input_file_basename):
            input_file_basename='.'.join(input_file_basename.split('.')[:-1])
        # Remove input file extension if present
        if '.' in input_file_basename:
            input_file_basename='.'.join(input_file_basename.split('.')[:-1])
    
        config.file_basename=input_file_basename
    
    # Set the output format
    config.output_format=args.output_format
    
    # Set final output file names and location
    config.pathabundance_file=os.path.join(output_dir,
            config.file_basename + config.pathabundance_file + "." + 
            config.output_format)
    config.pathcoverage_file=os.path.join(output_dir,
            config.file_basename + config.pathcoverage_file + "." + 
            config.output_format)
    config.genefamilies_file=os.path.join(output_dir,
            config.file_basename + config.genefamilies_file + "." + 
            config.output_format)

    # set the location of the temp directory
    if not args.remove_temp_output:
        config.temp_dir=os.path.join(output_dir,config.file_basename+"_humann2_temp")
        if not os.path.isdir(config.temp_dir):
            try:
                os.mkdir(config.temp_dir)
            except EnvironmentError:
                sys.exit("Unable to create temp directory: " + config.temp_dir)
    else:
        config.temp_dir=tempfile.mkdtemp( 
            prefix=config.file_basename+'_humann2_temp_',dir=output_dir)
        
    # create the unnamed temp directory
    config.unnamed_temp_dir=tempfile.mkdtemp(dir=config.temp_dir)
    
    message="Writing temp files to directory: " + config.temp_dir
    logger.info(message)
    if config.verbose: 
        print("\n"+message+"\n")

    # set the name of the log file 
    log_file=os.path.join(config.temp_dir,config.file_basename+".log")
    
    # change file name if set
    if args.o_log:
        log_file=args.o_log
        
    # configure the logger
    logging.basicConfig(filename=log_file,format='%(asctime)s - %(name)s - %(levelname)s: %(message)s',
        level=getattr(logging,args.log_level), filemode='w', datefmt='%m/%d/%Y %I:%M:%S %p')
    
     
def check_requirements(args):
    """
    Check requirements (file format, dependencies, permissions)
    """

    # Check the pathways database files exist and are readable
    utilities.file_exists_readable(config.pathways_database_part1)
    utilities.file_exists_readable(config.pathways_database_part2)

    # Determine the input file format if not provided
    if not args.input_format:
        args.input_format=utilities.determine_file_format(args.input)
        
        if args.input_format == "unknown":
            sys.exit("CRITICAL ERROR: Unable to determine the input file format." +
                " Please provide the format with the --input_format argument.")
    
    config.input_format = args.input_format
        
    # If the input file is compressed, then decompress
    if args.input_format.endswith(".gz"):
        new_file=utilities.gunzip_file(args.input)
        
        if new_file:
            args.input=new_file
            args.input_format=args.input_format.split(".")[0]
        else:
            sys.exit("CRITICAL ERROR: Unable to use gzipped input file. " + 
                " Please check the format of the input file.")
            
    # If the input format is in binary then convert to sam (tab-delimited text)
    if args.input_format == "bam":
        
        # Check for the samtools software
        if not utilities.find_exe_in_path("samtools"):
            sys.exit("CRITICAL ERROR: The samtools executable can not be found. "
            "Please check the install or select another input format.")
        
        new_file=utilities.bam_to_sam(args.input)
        
        if new_file:
            args.input=new_file
            args.input_format="sam"
        else:
            sys.exit("CRITICAL ERROR: Unable to convert bam input file to sam.")

    # If the input format is in biom then convert to tsv
    if args.input_format == "biom":
        
        # Check for the biom software
        if not utilities.find_exe_in_path("biom"):
            sys.exit("CRITICAL ERROR: The biom executable can not be found. "
            "Please check the install or select another input format.")
        
        new_file=utilities.biom_to_tsv(args.input)
        
        if new_file:
            args.input=new_file
            # determine the format of the file
            args.input_format=utilities.determine_file_format(args.input)
        else:
            sys.exit("CRITICAL ERROR: Unable to convert biom input file to tsv.")
            
    # If the biom output format is selected, check for the biom package
    if config.output_format=="biom":
        if not utilities.find_exe_in_path("biom"):
            sys.exit("CRITICAL ERROR: The biom executable can not be found. "
            "Please check the install or select another output format.")
     
    # If the file is fasta/fastq check for requirements   
    if args.input_format in ["fasta","fastq"]:
        # Check that the chocophlan directory exists
        if not config.bypass_nucleotide_index:
            if not os.path.isdir(config.chocophlan):
                if args.chocophlan:
                    sys.exit("CRITICAL ERROR: The directory provided for the ChocoPhlAn database at " 
                        + args.chocophlan + " does not exist. Please select another directory.")
                else:
                    sys.exit("CRITICAL ERROR: The default ChocoPhlAn database directory of "
                        + config.chocophlan + " does not exist. Please provide the location "
                        + "of the ChocoPhlAn directory using the --chocophlan option.")	
    
        # Check that the files in the chocophlan folder are of the right format
        if not config.bypass_nucleotide_index:
            valid_format_count=0
            for file in os.listdir(config.chocophlan):
                # expect most of the file names to be of the format g__*s__*
                if re.search("^[g__][s__]",file): 
                    valid_format_count+=1
            if valid_format_count == 0:
                sys.exit("CRITICAL ERROR: The directory provided for ChocoPhlAn does not "
                    + "contain files of the expected format (ie \'^[g__][s__]\').")
                
        # Check that the metaphlan2 executable can be found
        if not config.bypass_prescreen and not config.bypass_nucleotide_index:
            if not utilities.find_exe_in_path("metaphlan2.py"): 
                sys.exit("CRITICAL ERROR: The metaphlan2.py executable can not be found. "  
                    "Please check the install.")

        # Check that the bowtie2 executable can be found
        if not config.bypass_nucleotide_search:
            if not utilities.find_exe_in_path("bowtie2"): 
                sys.exit("CRITICAL ERROR: The bowtie2 executable can not be found. "  
                    "Please check the install.")
 
        if not config.bypass_translated_search:
            # Check that the uniref directory exists
            if not os.path.isdir(config.uniref):
                if args.uniref:
                    sys.exit("CRITICAL ERROR: The directory provided for the UniRef database at " 
                        + args.uniref + " does not exist. Please select another directory.")
                else:
                    sys.exit("CRITICAL ERROR: The default UniRef database directory of "
                        + config.uniref + " does not exist. Please provide the location "
                        + "of the UniRef directory using the --uniref option.")            	
    
            # Check that some files in the uniref folder are of the expected extension
            expected_database_extension=""
            if config.translated_alignment_selected == "usearch":
                expected_database_extension=config.usearch_database_extension
            elif config.translated_alignment_selected == "rapsearch":
                expected_database_extension=config.rapsearch_database_extension
            elif config.translated_alignment_selected == "diamond":
                expected_database_extension=config.diamond_database_extension
            
            valid_format_count=0
            database_files=os.listdir(config.uniref)
            for file in database_files:
                if file.endswith(expected_database_extension):
                    # if rapsearch check for the second database file
                    if config.translated_alignment_selected == "rapsearch":
                        database_file=re.sub(config.rapsearch_database_extension+"$","",file)
                        if database_file in database_files:
                            valid_format_count+=1
                    else:
                        valid_format_count+=1
                        
            if valid_format_count == 0:
                sys.exit("CRITICAL ERROR: The UniRef directory provided ( " + config.uniref 
                    + " ) does not contain any files that have been formatted to run with"
                    " the translated alignment software selected ( " +
                    config.translated_alignment_selected + " ). Please format these files so"
                    + " they are of the expected extension ( " + expected_database_extension +" ).")

            # Check for correct usearch version
            if config.translated_alignment_selected == "usearch":
                utilities.check_software_version("usearch","-version",config.usearch_version)
        
            # Check that the translated alignment executable can be found
            if not utilities.find_exe_in_path(config.translated_alignment_selected):
                sys.exit("CRITICAL ERROR: The " +  config.translated_alignment_selected + 
                    " executable can not be found. Please check the install.")
              
def timestamp_message(task, start_time):
    """
    Print and log a message about the task completed and the time
    Log messages are tab delimited for quick task/time access with awk
    Return the new start time
    """
    message="TIMESTAMP: Completed \t" + task + " \t:\t " + \
        str(int(round(time.time() - start_time))) + "\t seconds"
    logger.info(message)
    if config.verbose:
        print("\n"+message.replace("\t","")+"\n")   
        
    return time.time() 
              
def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    # Update the configuration settings based on the arguments
    update_configuration(args)
    
    # Check for required files, software, databases, and also permissions
    check_requirements(args)
    
    # Write the config settings to the log file
    config.log_settings()

    # Initialize alignments and gene scores
    alignments=store.Alignments()
    unaligned_reads_store=store.Reads()
    gene_scores=store.GeneScores()
    
    # If id mapping is provided then process
    if args.id_mapping:
        alignments.process_id_mapping(args.id_mapping)
    
    # Load in the reactions database
    reactions_database=store.ReactionsDatabase(config.pathways_database_part1)

    message="Load pathways database part 1: " + config.pathways_database_part1
    logger.info(message)
    if config.verbose:
        print(message)
    
    # Load in the pathways database
    pathways_database=store.PathwaysDatabase(config.pathways_database_part2, 
        config.pathways_recursion)
    
    message="Load pathways database part 2: " + config.pathways_database_part2
    logger.info(message)
    if config.verbose:
        print(message)

    # Start timer
    start_time=time.time()

    # Process fasta or fastq input files
    if args.input_format in ["fasta","fastq"]:
        # Run prescreen to identify bugs
        bug_file = "Empty"
        if args.taxonomic_profile:
            bug_file = os.path.abspath(args.taxonomic_profile)
        else:
            if not config.bypass_prescreen:
                bug_file = prescreen.alignment(args.input)
                start_time=timestamp_message("prescreen",start_time)
    
        # Create the custom database from the bugs list
        custom_database = ""
        if not config.bypass_nucleotide_index:
            custom_database = prescreen.create_custom_database(config.chocophlan, bug_file)
            start_time=timestamp_message("custom database creation",start_time)
        else:
            custom_database = "Bypass"
    
        # Run nucleotide search on custom database
        if custom_database != "Empty" and not config.bypass_nucleotide_search:
            if not config.bypass_nucleotide_index:
                nucleotide_index_file = nucleotide.index(custom_database)
                start_time=timestamp_message("database index",start_time)
            else:
                nucleotide_index_file = config.chocophlan
                
            nucleotide_alignment_file = nucleotide.alignment(args.input, 
                nucleotide_index_file)
    
            start_time=timestamp_message("nucleotide alignment",start_time)
    
            # Determine which reads are unaligned and reduce aligned reads file
            # Remove the alignment_file as we only need the reduced aligned reads file
            [ unaligned_reads_file_fasta, reduced_aligned_reads_file ] = nucleotide.unaligned_reads(
                nucleotide_alignment_file, alignments, unaligned_reads_store, keep_sam=True)
            
            start_time=timestamp_message("nucleotide alignment post-processing",start_time)
    
            # Print out total alignments per bug
            message="Total bugs from nucleotide alignment: " + str(alignments.count_bugs())
            logger.info(message)
            print(message)
            
            message=alignments.counts_by_bug()
            logger.info("\n"+message)
            print(message)        
    
            message="Total gene families from nucleotide alignment: " + str(alignments.count_genes())
            logger.info(message)
            print("\n"+message)
    
            # Report reads unaligned
            message="Estimate of unaligned reads: " + utilities.estimate_unaligned_reads_stored(
                args.input, unaligned_reads_store) + "%"
            logger.info(message)
            print("\n"+message+"\n")  
        else:
            logger.debug("Custom database is empty")
            reduced_aligned_reads_file = "Empty"
            unaligned_reads_file_fasta=args.input
            unaligned_reads_store=store.Reads(unaligned_reads_file_fasta)
    
        # Do not run if set to bypass translated search in config file
        if not config.bypass_translated_search:
            # Run translated search on UniRef database if unaligned reads exit
            if unaligned_reads_store.count_reads()>0:
                translated_alignment_file = translated.alignment(config.uniref, 
                    unaligned_reads_file_fasta)
        
                start_time=timestamp_message("translated alignment",start_time)
        
                # Determine which reads are unaligned
                translated_unaligned_reads_file_fastq = translated.unaligned_reads(
                    unaligned_reads_store, translated_alignment_file, alignments)
                
                start_time=timestamp_message("translated alignment post-processing",start_time)
        
                # Print out total alignments per bug
                message="Total bugs after translated alignment: " + str(alignments.count_bugs())
                logger.info(message)
                print(message)
            
                message=alignments.counts_by_bug()
                logger.info("\n"+message)
                print(message)
        
                message="Total gene families after translated alignment: " + str(alignments.count_genes())
                logger.info(message)
                print("\n"+message)
        
                # Report reads unaligned
                message="Estimate of unaligned reads: " + utilities.estimate_unaligned_reads_stored(
                    args.input, unaligned_reads_store) + "%"
                logger.info(message)
                print("\n"+message+"\n")  
            else:
                message="All reads are aligned so translated alignment will not be run"
                logger.info(message)
                print(message)
        else:
            message="Bypass translated search"
            logger.info(message)
            print(message)
    
    # Process input files of sam format
    elif args.input_format in ["sam"]:
        # Turn off frame picker if set on
        config.pick_frames_toggle="off"
        
        # Store the sam mapping results
        message="Process the sam mapping results ..."
        logger.info(message)
        print("\n"+message)
            
        [unaligned_reads_file_fasta, reduced_aligned_reads_file] = nucleotide.unaligned_reads(
            args.input, alignments, unaligned_reads_store, keep_sam=True)
        
        start_time=timestamp_message("alignment post-processing",start_time)
            
    # Process input files of tab-delimited blast format
    elif args.input_format in ["blastm8"]:
        
        # Store the blastm8 mapping results
        message="Process the blastm8 mapping results ..."
        logger.info(message)
        print("\n"+message)
        
        translated_unaligned_reads_file_fastq = translated.unaligned_reads(
            unaligned_reads_store, args.input, alignments)
        
        start_time=timestamp_message("alignment post-processing",start_time)
        
    # Clear all of the unaligned reads as they are no longer needed
    unaligned_reads_store.clear()
        
    # Compute or load in gene families
    output_files=[]
    if args.input_format in ["fasta","fastq","sam","blastm8"]:
        # Compute the gene families
        message="Computing gene families ..."
        logger.info(message)
        print("\n"+message)
        
        families_file=families.gene_families(alignments,gene_scores)
        output_files.append(families_file)
    
        start_time=timestamp_message("computing gene families",start_time)

    elif args.input_format in ["genetable"]:
        # Load the gene scores
        message="Process the gene table ..."
        logger.info(message)
        print("\n"+message)
        
        gene_scores.add_from_file(args.input,id_mapping_file=args.id_mapping) 
        
        start_time=timestamp_message("processing gene table",start_time)

    # Handle input files of unknown formats
    else:
        sys.exit("CRITICAL ERROR: Input file of unknown format.")
        
    # Clear all of the alignments data as they are no longer needed
    alignments.clear()
    
    # Identify reactions and then pathways from the alignments
    message="Computing pathways abundance and coverage ..."
    logger.info(message)
    print("\n"+message)
    pathways_and_reactions_store=modules.identify_reactions_and_pathways(
        gene_scores, reactions_database, pathways_database)

    # Compute pathway abundance and coverage
    abundance_file, coverage_file=modules.compute_pathways_abundance_and_coverage(
        pathways_and_reactions_store, pathways_database)
    output_files.append(abundance_file)
    output_files.append(coverage_file)

    start_time=timestamp_message("computing pathways",start_time)

    message="\nOutput files created: \n" + "\n".join(output_files) + "\n"
    logger.info(message)
    print(message)

    # Remove the unnamed temp files
    utilities.remove_directory(config.unnamed_temp_dir)

    # Remove named temp directory
    if args.remove_temp_output:
        utilities.remove_directory(config.temp_dir)
        
