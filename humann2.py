#!/usr/bin/env python

"""
HUMAnN2 : HMP Unified Metabolic Analysis Network 2

HUMAnN2 is a pipeline for efficiently and accurately determining 
the presence/absence and abundance of microbial pathways in a community 
from metagenomic data. Sequencing a metagenome typically produces millions 
of short DNA/RNA reads.

Dependencies: MetaPhlAn, ChocoPhlAn, Bowtie2, Usearch

To Run: ./humann2.py -i <input.fastq> -c <chocophlan/> -u <uniref/>
"""

import argparse, sys, subprocess, os, time, tempfile, re

from src import utilities, prescreen, nucleotide_search
from src import translated_search, config, quantify_families

def parse_arguments (args):
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
        default=False)
    parser.add_argument(
        "-d","--debug", 
        help="bypass commands if the output files exist\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "--bypass_prescreen", 
        help="bypass the prescreen step and run on the full ChocoPhlAn database\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "--bypass_nucleotide_index", 
        help="bypass the nucleotide index step and run on the indexed ChocoPhlAn database\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "-i", "--input", 
        help="fastq/fasta input file\n[REQUIRED]", 
        metavar="<input.fastq>", 
        required=True)
    parser.add_argument(
        "-c", "--chocophlan",
        help="directory containing the ChocoPhlAn database\n[REQUIRED]", 
        metavar="<chocophlan/>",
        required=True)
    parser.add_argument(
        "-u", "--uniref",
        help="directory containing the UniRef database\n[REQUIRED]", 
        metavar="<uniref/>",
        required=True)
    parser.add_argument(
        "--metaphlan",
        help="directory containing the MetaPhlAn software\n[DEFAULT: $PATH]", 
        metavar="<metaplhan/>")
    parser.add_argument(
        "--o_pathabundance", 
        help="output file for pathway abundance\n" + 
        "[DEFAULT: $input_dir/pathabundance.tsv]", 
        metavar="<pathabundance.tsv>")
    parser.add_argument(
        "--o_pathpresence",
        help="output file for pathway presence/absence\n" + 
        "[DEFAULT: $input_dir/pathpresence.tsv]", 
        metavar="<pathpresence.tsv>")
    parser.add_argument(
        "--o_genefamilies", 
        help="output file for gene families\n" + 
        "[DEFAULT: $input_dir/genefamilies.tsv]", 
        metavar="<genefamilies.tsv>")
    parser.add_argument(
        "--temp", 
        help="directory to store temp output files\n" + 
            "[DEFAULT: temp files are removed]", 
        metavar="<temp/>")
    parser.add_argument(
        "--bowtie2",
        help="directory of the bowtie2 executable\n[DEFAULT: $PATH]", 
        metavar="<bowtie2/>")
    parser.add_argument(
        "--threads", 
        help="number of threads to use with bowtie2\n[DEFAULT: 1]", 
        metavar="<1>", 
        type=int,
        default=1) 
    parser.add_argument(
        "--prescreen_threshold", 
        help="minimum percentage of reads matching a species\n[DEFAULT: "
            + str(config.prescreen_threshold) + "]", 
        metavar="<" + str(config.prescreen_threshold) + ">", 
        type=float,
        default=config.prescreen_threshold) 
    parser.add_argument(
        "--identity_threshold", 
        help="identity threshold to use with the translated search\n[DEFAULT: " 
            + str(config.id_threshold_default) + "]", 
        metavar="<" + str(config.id_threshold_default) + ">", 
        type=float,
        default=config.id_threshold_default) 
    parser.add_argument(
        "--usearch", 
        help="directory containing the usearch executable\n[DEFAULT: $PATH]", 
        metavar="<usearch/>")
    parser.add_argument(
        "--rapsearch", 
        help="directory containing the rapsearch executable\n[DEFAULT: $PATH]", 
        metavar="<rapsearch/>")
    parser.add_argument(
        "--metaphlan_output", 
        help="output file created by metaphlan\n[DEFAULT: file will be created]", 
        metavar="<bugs_list.tsv>")
    parser.add_argument(
        "--translated_alignment", 
        help="software to use for translated alignment\n[DEFAULT: " + 
            config.translated_alignment_selected + "]", 
        default=config.translated_alignment_selected,
        choices=config.translated_alignment_choices)

    return parser.parse_args()
	
	 
def check_requirements(args):
    """
    Check requirements (file format, dependencies, permissions)
    """
    # Update translated alignment software if set
    if args.translated_alignment:
        config.translated_alignment_selected=args.translated_alignment

    # Check that the input file exists, is readable, and is fasta/fastq
    if utilities.fasta_or_fastq(args.input) == "error":
        sys.exit("ERROR: The input file is not of a fasta or fastq format.")

    # Check that the chocophlan directory exists
    if not args.bypass_nucleotide_index:
        if not os.path.isdir(args.chocophlan):
            sys.exit("ERROR: The directory provided for ChocoPhlAn at " 
                + args.chocophlan + " does not exist. Please select another directory.")	
    
    # Check that the files in the chocophlan folder are of the right format
    if not args.bypass_nucleotide_index:
        valid_format_count=0
        for file in os.listdir(args.chocophlan):
            # expect most of the file names to be of the format g__*s__*
            if re.search("^[g__][s__]",file): 
                valid_format_count+=1
        if valid_format_count == 0:
            sys.exit("ERROR: The directory provided for ChocoPhlAn does not "
                + "contain files of the expected format.")
 
    # Check that the uniref directory exists
    if not os.path.isdir(args.uniref):
        sys.exit("ERROR: The directory provided for the UniRef database at " 
            + args.uniref + " does not exist. Please select another directory.")	

    # Check that all files in the uniref folder are of *.udb extension or fasta
    # if translated search selected is usearch
    if config.translated_alignment_selected == "usearch":
        for file in os.listdir(args.uniref):
            if not file.endswith(config.usearch_database_extension):
                if utilities.fasta_or_fastq(os.path.join(args.uniref,file)) != "fasta":
                    print os.path.join(args.uniref,file)
                    print utilities.fasta_or_fastq(os.path.join(args.uniref,file))
                    sys.exit("ERROR: The directory provided for the UniRef database "
                        + "contains files of an unexpected format. Only files of the"
                        + " udb or fasta format are allowed.") 

    # Check that some of the database files are of the *.info extension
    valid_format_count=0
    if config.translated_alignment_selected == "rapsearch":
        for file in os.listdir(args.uniref):
            if file.endswith(config.rapsearch_database_extension):
                valid_format_count+=1
        if valid_format_count == 0:
            sys.exit("ERROR: The UniRef directory provided has not been formatted to run with"
                " the rapsearch translated alignment software. Please format these files.")

    # Check for correct usearch version
    if config.translated_alignment_selected == "usearch":
        utilities.check_software_version("usearch","-version",config.usearch_version)
    
    # Check that the translated alignment executable can be found
    if not utilities.find_exe_in_path(config.translated_alignment_selected):
        sys.exit("ERROR: The " +  config.translated_alignment_selected + 
            " executable can not be found. Please check the install.")
    
    # Check that the metaphlan2 executable can be found
    if not args.bypass_prescreen and not args.bypass_nucleotide_index:
        if not utilities.find_exe_in_path("metaphlan2.py"): 
            sys.exit("ERROR: The metaphlan2.py executable can not be found. "  
                "Please check the install.")

    # Check that the bowtie2 executable can be found
    if not utilities.find_exe_in_path("bowtie2"): 
        sys.exit("ERROR: The bowtie2 executable can not be found. "  
            "Please check the install.")

    # Check that the directory that holds the input file is writeable
    input_dir = os.path.dirname(args.input)
    if not input_dir:
        input_dir=os.getcwd()
    if not os.access(input_dir, os.W_OK):
        sys.exit("ERROR: The directory which holds the input file is not " + 
            "writeable. This software needs to write files to this directory.\n" +
            "Please use another directory to hold your input file.") 

    # if set, check that the temp directory location is writeable
    if args.temp:
        if os.path.isdir(args.temp):
            if not os.access(args.temp, os.W_OK):
                sys.exit("ERROR: The directory set to hold the temp files " + 
                    "is not writeable. Please change the permissions or select" +
                    " another directory.")
        else:
            path_to_temp_dir=os.path.dirname(args.temp)
            if not os.path.basename(args.temp):
                path_to_temp_dir=os.path.dirname(os.path.dirname(args.temp))
            if not path_to_temp_dir:
                path_to_temp_dir=os.getcwd()
            if not os.access(path_to_temp_dir, os.W_OK):
                sys.exit("ERROR: The directory set to hold the temp files " + 
                    "is not writeable. Please change the permissions or select" +
                    " another directory.")

    # if set, update the config run mode to debug
    if args.debug:
        config.debug=True  
            
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

    return input_dir	


def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)

    # If set, append paths executable locations
    if args.metaphlan:
        utilities.add_exe_to_path(args.metaphlan)	
	
    if args.bowtie2:
        utilities.add_exe_to_path(args.bowtie2)
	
    if args.usearch:
        utilities.add_exe_to_path(args.usearch)

    if args.rapsearch:
        utilities.add_exe_to_path(args.rapsearch)

    # Add the location of the humann1 scripts to the path
    utilities.add_exe_to_path(
        os.path.join(os.path.dirname(os.path.realpath(__file__)),
            config.humann1_scripts))
	
    # Check for required files, software, databases, and also permissions
    # If all pass, return location of input_dir to write output to
    output_dir=check_requirements(args)

    # if the temp_dir is set by the user then use that directory
    if args.temp:
        config.temp_dir=args.temp
        if not os.path.isdir(config.temp_dir):
            os.mkdir(config.temp_dir)	
    else:
        config.temp_dir=tempfile.mkdtemp( 
            prefix='humann2_temp_', dir=output_dir)

    if config.verbose:
        print "Writing temp files to directory: " + config.temp_dir

    # Set the basename of the temp files to the sample name
    config.file_basename=os.path.splitext(os.path.basename(args.input))[0]

    # Start timer
    start_time=time.time()

    # Run prescreen to identify bugs
    bug_file = "Empty"
    if args.metaphlan_output:
        bug_file = args.metaphlan_output
    else:
        if not config.bypass_prescreen:
            bug_file = prescreen.alignment(args.input, 
                args.threads)

    if config.verbose:
        print str(int(time.time() - start_time)) + " seconds from start"

    # Create the custom database from the bugs list
    custom_database = ""
    if not config.bypass_nucleotide_index:
        custom_database = prescreen.create_custom_database(args.chocophlan, 
            args.prescreen_threshold, bug_file)
    else:
        custom_database = "Bypass"

    if config.verbose:
        print str(int(time.time() - start_time)) + " seconds from start"

    # Run nucleotide search on custom database
    if custom_database != "Empty":
        if not config.bypass_nucleotide_index:
            nucleotide_index_file = nucleotide_search.index(custom_database)
        else:
            nucleotide_index_file = args.chocophlan
        nucleotide_alignment_file = nucleotide_search.alignment(args.input, 
            args.threads, nucleotide_index_file)

        if config.verbose:
            print str(int(time.time() - start_time)) + " seconds from start"

        # Determine which reads are unaligned and reduce aligned reads file
        # Remove the alignment_file as we only need the reduced aligned reads file
        [ unaligned_reads_file_fastq, reduced_aligned_reads_file ] = nucleotide_search.unaligned_reads(
            args.input, nucleotide_alignment_file)

        # Report reads unaligned
        print "\nEstimate of unaligned reads: " + utilities.estimate_unaligned_reads(
            args.input, unaligned_reads_file_fastq) + "%\n"  
    else:
        reduced_aligned_reads_file = "Empty"
        unaligned_reads_file_fastq=args.input

    # Run translated search on UniRef database
    translated_alignment_file = translated_search.alignment(args.uniref, unaligned_reads_file_fastq,
        args.identity_threshold, args.threads)

    if config.verbose:
        print str(int(time.time() - start_time)) + " seconds from start"

    # Determine which reads are unaligned
    translated_unaligned_reads_file_fastq = translated_search.unaligned_reads(
        unaligned_reads_file_fastq, translated_alignment_file)

    # Report reads unaligned
    print "\nEstimate of unaligned reads: " + utilities.estimate_unaligned_reads(
        args.input, translated_unaligned_reads_file_fastq) + "%\n"  

    # Identify the hits from the alignments
    hits_file=quantify_families.find_hits(reduced_aligned_reads_file, translated_alignment_file)

    # Remove temp directory
    if not args.temp:
        utilities.remove_directory(config.temp_dir)

if __name__ == "__main__":
	main()
