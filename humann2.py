#!/usr/bin/env python

"""
HUMAnN2 : HMP Unified Metabolic Analysis Network 2

HUMAnN2 is a pipeline for efficiently and accurately determining 
the coverage and abundance of microbial pathways in a community 
from metagenomic data. Sequencing a metagenome typically produces millions 
of short DNA/RNA reads.

Dependencies: MetaPhlAn, ChocoPhlAn, Bowtie2, Rapsearch2 or Usearch

To Run: ./humann2.py -i <input.fastq> -c <chocophlan/> -u <uniref/>
"""

import argparse, sys, subprocess, os, time, tempfile, re

from src import utilities, prescreen, nucleotide_search, store
from src import translated_search, config, quantify_families, quantify_modules

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
        "[DEFAULT: $input_dir/$SAMPLE_pathabundance.tsv]", 
        metavar="<pathabundance.tsv>")
    parser.add_argument(
        "--o_pathcoverage",
        help="output file for pathway coverage\n" + 
        "[DEFAULT: $input_dir/$SAMPLE_pathcoverage.tsv]", 
        metavar="<pathcoverage.tsv>")
    parser.add_argument(
        "--o_genefamilies", 
        help="output file for gene families\n" + 
        "[DEFAULT: $input_dir/$SAMPLE_genefamilies.tsv]", 
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
        help="number of threads/processes\n[DEFAULT: 1]", 
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
    parser.add_argument(
        "--xipe",
        help="turn on/off the xipe computation\n[DEFAULT: " +
        config.xipe_toggle + " ]",
        default=config.xipe_toggle,
        choices=config.toggle_choices)
    parser.add_argument(
        "--minpath",
        help="turn on/off the minpath computation\n[DEFAULT: " + 
        config.minpath_toggle + " ]",
        default=config.minpath_toggle,
        choices=config.toggle_choices)

    return parser.parse_args()
	
	 
def update_configuration(args):
    """
    Update the configuration settings based on the arguments
    """

    # If set, append paths executable locations
    if args.metaphlan:
        utilities.add_exe_to_path(args.metaphlan)    
    
    if args.bowtie2:
        utilities.add_exe_to_path(args.bowtie2)
    
    if args.usearch:
        utilities.add_exe_to_path(args.usearch)

    if args.rapsearch:
        utilities.add_exe_to_path(args.rapsearch)

    # Update data directory to full path
    humann2_fullpath=os.path.dirname(os.path.realpath(__file__))
    config.data_folder=os.path.join(humann2_fullpath,
        config.data_folder)    

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
    
    # Update translated alignment software if set
    if args.translated_alignment:
        config.translated_alignment_selected=args.translated_alignment
        
    # Update the computation toggle choices
    if args.xipe:
        config.xipe_toggle=args.xipe
        
    if args.minpath:
        config.minpath_toggle=args.minpath
        
    # If minpath is set to run, install if not already installed
    if config.minpath_toggle == config.toggle_on:
        utilities.install_minpath()
     
    # Check that the directory that holds the input file is writeable
    # before creating files/directories in that folder
    input_dir = os.path.dirname(args.input)
    if not input_dir:
        input_dir=os.getcwd()
    if not os.access(input_dir, os.W_OK):
        sys.exit("ERROR: The directory which holds the input file is not " + 
            "writeable. This software needs to write files to this directory.\n" +
            "Please use another directory to hold your input file.") 

    # Set the basename of the temp files to the sample name
    config.file_basename=os.path.splitext(os.path.basename(args.input))[0]
    
    # Set final output file names
    if args.o_pathabundance:
        config.pathabundance_file=args.o_pathabundance
    else:
        config.pathabundance_file=os.path.join(input_dir,
            config.file_basename + config.pathabundance_file)

    if args.o_pathcoverage:
        config.pathcoverage_file=args.o_pathcoverage
    else:
        config.pathcoverage_file=os.path.join(input_dir,
            config.file_basename + config.pathcoverage_file)

    if args.o_genefamilies:
        config.genefamilies_file=args.o_genefamilies
    else:
        config.genefamilies_file=os.path.join(input_dir,
            config.file_basename + config.genefamilies_file)

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

    # if the temp_dir is set by the user then use that directory
    if args.temp:
        config.temp_dir=args.temp
        if not os.path.isdir(config.temp_dir):
            os.mkdir(config.temp_dir)    
    else:
        config.temp_dir=tempfile.mkdtemp( 
            prefix='humann2_temp_', dir=input_dir)

    if config.verbose:
        print "\nWriting temp files to directory: " + config.temp_dir + "\n"


     
def check_requirements(args):
    """
    Check requirements (file format, dependencies, permissions)
    """

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

def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    # Update the configuration settings based on the arguments
    update_configuration(args)
	
    # Check for required files, software, databases, and also permissions
    check_requirements(args)

    # Initialize alignments
    alignments=store.alignments()
    
    # Load in the reactions database
    genes_to_reactions=os.path.join(config.data_folder,
        config.metacyc_gene_to_reactions)
    reactions_database=store.reactions_database(genes_to_reactions)

    if config.verbose:
        print "Load reactions from database: " + genes_to_reactions
    
    # Load in the pathways database
    reactions_to_pathways=os.path.join(config.data_folder,
        config.metacyc_reactions_to_pathways)
    pathways_database=store.pathways_database(reactions_to_pathways)
    
    if config.verbose:
        print "Load pathways from database: " + reactions_to_pathways

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
        [ unaligned_reads_file_fasta, unaligned_reads_store, reduced_aligned_reads_file ] = nucleotide_search.unaligned_reads(
            args.input, nucleotide_alignment_file, alignments)

        # Print out total alignments per bug
        print "Total bugs from nucleotide alignment: " + str(alignments.count_bugs())
        alignments.print_bugs()

        print "\nTotal gene families from nucleotide alignment: " + str(alignments.count_genes())

        # Report reads unaligned
        print "\nEstimate of unaligned reads: " + utilities.estimate_unaligned_reads(
            args.input, unaligned_reads_file_fasta) + "%\n"  
    else:
        reduced_aligned_reads_file = "Empty"
        unaligned_reads_file_fasta=args.input
        unaligned_reads_store=store.reads(unaligned_reads_file_fasta)

    # Run translated search on UniRef database
    translated_alignment_file = translated_search.alignment(args.uniref, 
        unaligned_reads_file_fasta, args.identity_threshold, args.threads)

    if config.verbose:
        print str(int(time.time() - start_time)) + " seconds from start"

    # Determine which reads are unaligned
    translated_unaligned_reads_file_fastq = translated_search.unaligned_reads(
        unaligned_reads_store, translated_alignment_file, alignments)

    # Print out total alignments per bug
    print "Total bugs after translated alignment: " + str(alignments.count_bugs())
    alignments.print_bugs()

    print "\nTotal gene families after translated alignment: " + str(alignments.count_genes())

    # Report reads unaligned
    print "\nEstimate of unaligned reads: " + utilities.estimate_unaligned_reads(
        args.input, translated_unaligned_reads_file_fastq) + "%\n"  

    # Compute the gene families
    print "\nComputing gene families ..."
    families_file=quantify_families.gene_families(alignments)

    if config.verbose:
        print str(int(time.time() - start_time)) + " seconds from start"
    
    # Identify reactions and then pathways from the alignments
    print "\nComputing pathways abundance and coverage ..."
    pathways_and_reactions_store=quantify_modules.identify_reactions_and_pathways(
        args.threads, alignments, reactions_database, pathways_database)

    # Compute pathway abundance and coverage
    abundance_file, coverage_file=quantify_modules.compute_pathways_abundance_and_coverage(
        args.threads, pathways_and_reactions_store, pathways_database)

    if config.verbose:
        print str(int(time.time() - start_time)) + " seconds from start"

    output_files=[families_file,abundance_file,coverage_file]
    print "\nOutput files created: \n" + "\n".join(output_files) + "\n"

    # Remove temp directory
    if not args.temp:
        utilities.remove_directory(config.temp_dir)

if __name__ == "__main__":
	main()
