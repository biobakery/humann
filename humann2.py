#!/usr/bin/env python

__title__ = "HUMAnN2 : HMP Unified Metabolic Analysis Network"
__authors__ = "Lauren McIver, George Weingart, Curtis Huttenhower"
__version__ = "0.1"
__date__ = "3 July 2014"
__maintainer__ = "Lauren McIver"
__email__ = "lauren.j.mciver@gmail.com"
__status__ = "WORK_IN_PROGRESS"

import argparse, sys, subprocess, os

## Parse the arguments from the user
def parse_arguments (args):
	parser = argparse.ArgumentParser(description= __title__ + '\n' + "Version: " + __version__ + '\n' 
			 		 + "Authors: " + __authors__ + '\n' 
					 + "\n\nNOTE: This is a " + __status__ + " version with limited function.\n", 
					 formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-i", "--input", help="fastq/fasta input file.\n[REQUIRED]", metavar="<input.fastq>", required=True)
	parser.add_argument("-m", "--metaphlan", help="Location of the MetaPhlAn software.\n[REQUIRED]", 
			    metavar="<metaplhan_dir/>", required=True)
	parser.add_argument("-c", "--chocophlan", help="Location of the ChocoPhlAn database.\n[REQUIRED]", 
			    metavar="<chocoplhan_dir/>", required=True)
	parser.add_argument("--o_pathabundance", help="Output file for pathway abundance.\n[DEFAULT: $input_dir/pathabundance.tsv]", 
			    metavar="<pathabundance.tsv>")
	parser.add_argument("--o_pathpresence", help="Output file for pathway presence/absence.\n[DEFAULT: $input_dir/pathpresence.tsv]", 
			    metavar="<pathpresence.tsv>")
	parser.add_argument("--o_genefamilies", help="Output file for gene families.\n[DEFAULT: $input_dir/genefamilies.tsv]", 
			    metavar="<genefamilies.tsv>")
	parser.add_argument("--debug", help="Print debug output files to $input_dir/debug/*.\n[DEFAULT: false]", 
			    metavar="<false>", type=bool)
	parser.add_argument("--bowtie2", help="Location of the bowtie2 executable.\n[DEFAULT: $PATH/bowtie2]", metavar="<bowtie2.exe>")
	parser.add_argument("--threads", help="Number of threads to use with bowtie2.\n[DEFAULT: 1]", metavar="<1>", type=int) 
	parser.add_argument("--usearch", help="Location of the usearch executable.\n[DEFAULT: $PATH/usearch]", metavar="<userach.exe>")
	parser.add_argument("--samtools", help="Location of the samtools executable.\n[DEAFULT: $PATH/samtools]", metavar="<samtools.exe>")
	return parser.parse_args()

## Check requirements (input file format, required third-party software/databases available, permissions)
def check_requirements(args):

	## Check that the input file exists
	if ( not ( os.path.isfile(args.input) ) ):
		sys.exit("ERROR: The input file provided does not exist at " + args.input + ". Please select another input file.")

	## Check that the metphlan directory exists
	if ( not ( os.path.isdir(args.metaphlan) ) ):
		sys.exit("ERROR: The directory provided for MetaPhlAn at " + args.metaphlan + " does not exist. Please select another directory.")	

	## Check that the chocophlan directory exists
	if ( not ( os.path.isdir(args.chocophlan) ) ):
		sys.exit("ERROR: The directory provided for ChocoPhlAn at " + args.chocophlan + " does not exist. Please select another directory.")	

	## Check that the input file entered is a fastq or fasta file
	if ( not ( args.input.endswith(".fq") or args.input.endswith(".fastq") or args.input.endswith(".fasta") or args.input.endswith(".fa") )):
		sys.exit("ERROR: The input file is not valid. Please provide a fastq or fasta file. Recognized file extensions are " + 
			 "*.fq, *.fastq, *.fasta, and *.fa .")  

	## Check that the bowtie2 executable can be found
	if args.bowtie2:
		try:
			p = subprocess.check_output([args.bowtie2,"-h"])
		except OSError as e:
			sys.exit("ERROR: The bowtie2 executable provided at " + args.bowtie2 + " is not executing properly. Please check the install.")
	else:
		try:
			p = subprocess.check_output(["bowtie2","-h"])
		except OSError as e:
			sys.exit("ERROR: The bowtie2 executable expected to be in your $PATH is not executing properly. \nPlease check the install " +
				 "or provide another path to bowtie2 using the --bowtie2 argument.")

	## Check that the usearch executable can be found
	if args.usearch:
		try:
			p = subprocess.check_output([args.usearch,"-h"])
		except OSError as e:
			sys.exit("ERROR: The usearch executable provided at " + args.usearch + " is not executing properly. Please check the install.")
	else:
		try:
			p = subprocess.check_output(["usearch","-h"])
		except OSError as e:
			sys.exit("ERROR: The usearch executable expected to be in your $PATH is not executing properly. \nPlease check the install " +
				 "or provide another path to usearch using the --usearch argument.")

	## Check that the samtools executable can be found
	if args.samtools:
		try:
			p = subprocess.check_output([args.samtools,"-h"])
		except OSError as e:
			sys.exit("ERROR: The samtools executable provided at " + args.samtools + " is not executing properly. Please check the install.")
	else:
		try:
			p = subprocess.check_output(["samtools","-h"])
		except OSError as e:
			sys.exit("ERROR: The samtools executable expected to be in your $PATH is not executing properly. \nPlease check the install " +
				 "or provide another path to samtools using the --samtools argument.")

	## Check that the directory that holds the input file is writeable
	input_dir = os.path.dirname(args.input)
	if ( not ( os.access(input_dir, os.W_OK) ) ):
		sys.exit("ERROR: The directory which holds the input file is not writeable. This software needs to write files to this directory.\n" +
			 "Please use another directory to hold your input file.") 

def main(argv):

	# Parse arguments from command line
	args=parse_arguments(argv)

	# Check for required files, software, databases, and also permissions
	check_requirements(args)


if __name__ == "__main__":
	main(sys.argv)
	
