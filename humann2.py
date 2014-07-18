#!/usr/bin/env python

"""
HUMAnN2 : HMP Unified Metabolic Analysis Network 2

HUMAnN2 is a pipeline for efficiently and accurately determining 
the presence/absence and abundance of microbial pathways in a community 
from metagenomic data. Sequencing a metagenome typically produces millions 
of short DNA/RNA reads.

Dependencies: MetaPhlAn, ChocoPhlAn, Bowtie2, Usearch, Samtools

To Run: ./humann2.py -i <input.fastq> -m <metaphlan_dir> -c <chocophlan_dir>
"""

import argparse, sys, subprocess, os

from src import utilities, prescreen, nucleotide_search

def parse_arguments (args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "HUMAnN2 : HMP Unified Metabolic Analysis Network 2\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i", "--input", 
        help="fastq/fasta input file.\n[REQUIRED]", 
        metavar="<input.fastq>", 
        required=True)
    parser.add_argument(
        "-m", "--metaphlan",
        help="Location of the MetaPhlAn software.\n[REQUIRED]", 
        metavar="<metaplhan_dir>",
        required=True)
    parser.add_argument(
        "-c", "--chocophlan",
        help="Location of the ChocoPhlAn database.\n[REQUIRED]", 
        metavar="<chocoplhan_dir>",
        required=True)
    parser.add_argument(
        "--o_pathabundance", 
        help="Output file for pathway abundance.\n" + 
        "[DEFAULT: $input_dir/pathabundance.tsv]", 
        metavar="<pathabundance.tsv>")
    parser.add_argument(
        "--o_pathpresence",
        help="Output file for pathway presence/absence.\n" + 
        "[DEFAULT: $input_dir/pathpresence.tsv]", 
        metavar="<pathpresence.tsv>")
    parser.add_argument(
        "--o_genefamilies", 
        help="Output file for gene families.\n" + 
        "[DEFAULT: $input_dir/genefamilies.tsv]", 
        metavar="<genefamilies.tsv>")
    parser.add_argument(
        "--debug", 
        help="Print debug output files to $input_dir/debug/*.\n" + 
            "[DEFAULT: false]", 
        metavar="<false>", 
        type=bool)
    parser.add_argument(
        "--bowtie2",
        help="The directory of the bowtie2 executable.\n[DEFAULT: $PATH]", 
        metavar="<bowtie2/>")
    parser.add_argument(
        "--threads", 
        help="Number of threads to use with bowtie2.\n[DEFAULT: 1]", 
        metavar="<1>", 
        type=int,
        default=1) 
    parser.add_argument(
        "--prescreen_threshold", 
        help="The minimum percentage of reads matching a species.\n[DEFAULT: 0.01]", 
        metavar="<0.01>", 
        type=float,
        default=0.01) 
    parser.add_argument(
        "--usearch", 
        help="The directory of the usearch executable.\n[DEFAULT: $PATH]", 
        metavar="<usearch/>")
    parser.add_argument(
        "--samtools", 
        help="The directory of the samtools executable.\n[DEAFULT: $PATH]", 
        metavar="<samtools/>")
    parser.add_argument(
        "--metaphlan_output", 
        help="The output file created by metaphlan.\n[DEAFULT: Metaphlan will be run to create file.]", 
        metavar="<bugs_list.tsv>")

    return parser.parse_args()
	
	 
def check_requirements(args):
    """
    Check requirements (file format, dependencies, permissions)
    """
    # Check that the input file exists, is readable, and is fasta/fastq
    if utilities.fasta_or_fastq(args.input) == "error":
        sys.exit("ERROR: The input file is not of a fasta or fastq format.")

    # Check that the metphlan directory exists
    if not os.path.isdir(args.metaphlan):
        sys.exit("ERROR: The directory provided for MetaPhlAn at " 
            + args.metaphlan + " does not exist. Please select another directory.")	

    # Check that the chocophlan directory exists
    if not os.path.isdir(args.chocophlan):
        sys.exit("ERROR: The directory provided for ChocoPhlAn at " 
            + args.chocophlan + " does not exist. Please select another directory.")	

    # Check that the bowtie2 executable can be found
    if not utilities.find_exe_in_path("bowtie2"): 
        sys.exit("ERROR: The bowtie2 executable can not be found. "  
            "Please check the install.")

    # Check that the usearch executable can be found
    if not utilities.find_exe_in_path("usearch"):
        sys.exit("ERROR: The usearch executable can not be found. " +  
            "Please check the install.")
	
    # Check that the samtools executable can be found
    if not utilities.find_exe_in_path("samtools"):
        sys.exit("ERROR: The samtools executable can not be found. " +  
            "Please check the install.")
	
    # Check that the directory that holds the input file is writeable
    input_dir = os.path.dirname(args.input)
    if not os.access(input_dir, os.W_OK):
        sys.exit("ERROR: The directory which holds the input file is not " + 
            "writeable. This software needs to write files to this directory.\n" +
            "Please use another directory to hold your input file.") 
	
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
		
    if args.samtools:
        utilities.add_exe_to_path(args.samtools)
				
    # Check for required files, software, databases, and also permissions
    # If all pass, return location of input_dir to write output to
    output_dir=check_requirements(args)

    # Set the location of the debug dir
    debug_dir=os.path.join(output_dir, "debug")

    # Create the debug directory
    if not os.path.isdir(debug_dir):
        os.makedirs(debug_dir)	

    # Run prescreen to identify bugs
    bug_file = ""
    if args.metaphlan_output:
        bug_file = args.metaphlan_output
    else:
        bug_file = prescreen.run_alignment(args.metaphlan, args.input, 
            args.threads, debug_dir)

    # Create the custom database from the bugs list
    custom_database = prescreen.create_custom_database(args.chocophlan, 
        args.prescreen_threshold, bug_file, debug_dir)

    # Run nucleotide search on custom database
    alignment_file = nucleotide_search.alignment(custom_database, args.input, 
        debug_dir, args.threads)

    # Determine which reads are unaligned
    unaligned_reads_fastq = nucleotide_search.unaligned_reads(args.input,
        alignment_file, debug_dir)

if __name__ == "__main__":
	main()
