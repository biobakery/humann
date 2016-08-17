#!/usr/bin/env python

"""
Create uniref database

This module will create the uniref database used by HUMAnN2. 

To Run: 
$ humann2_create_uniref_database -i <uniref50.fasta> -o <uniref50>

To Run with filtering: 
$ humann2_create_uniref_database -i <uniref50.fasta> -o <uniref50> -f <uncharacterized>

To Run with filtering using a list of ids: 
$ humann2_create_uniref_database -i <uniref50.fasta> -o <uniref50> -l <id_list.tsv>

"""

import argparse
import tempfile
import sys
import re
import shutil
import os
import subprocess
import gzip
import bz2

from humann2 import config

# pathways databases with uniref ids
PATHWAYS_DATABASE1=config.metacyc_gene_to_reactions
PATHWAYS_DATABASE2=config.metacyc_reactions_to_pathways

PATHWAYS_DATABASES=[PATHWAYS_DATABASE1,PATHWAYS_DATABASE2]
PATHWAYS_DELIMITER="\t"
PATHWAYS_UNIREF_IDENTIFIER="uniref50"

# the delimiter between items in the uniref fasta names for each sequence
UNIREF_DELIMITER=" "
# the location of the uniref id in the uniref fasta name
UNIREF_ID_INDEX=0

# the delimiter to use for adding the gene length
ANNOTATION_DELIMITER="|"
# the location to add the gene length annotation
ANNOTATION_INDEX=0

# the delimiter for the file id list
FILTER_ID_DELIMITER="\t"
# the index of the ids in the file
FILTER_ID_INDEX=0

def store_pathways(verbose):
    """
    Store the uniref ids from the pathways file
    """
    
    uniref_pathways={}
    for file in PATHWAYS_DATABASES:
        try:
            if file.endswith(".gz"):
                file_handle = gzip.open(file, "rt")
                readlines = file_handle.readlines()
            elif file.endswith(".bz2"):
                file_handle = bz2.BZ2File(file, "r")
                readlines = [line.decode('utf-8') for line in file_handle.readlines()]
            else:
                file_handle=open(file,"rt")
                readlines = file_handle.readlines()
            if verbose:
                print("Reading data from pathways file: " + file)
        except EnvironmentError:
            sys.exit("Unable to open file: " + file)

        file_handle.close()

        for line in readlines:
            tokens=line.rstrip('\n').split(PATHWAYS_DELIMITER)
            for token in tokens:
                if re.search(PATHWAYS_UNIREF_IDENTIFIER,token,re.IGNORECASE):
                   uniref_pathways[token]=1
                   
        if verbose:
            print("Total ids stored: " + str(len(uniref_pathways.keys())))
        
    return uniref_pathways 

def add_gene_length(sequence_name, sequence):
    """
    Compute the length of the gene and add it to the name
    """
    
    # the length of the sequence is multiplied
    # since the sequence is codons
    gene_length=len(sequence)*3
    
    # add the gene length to the name
    tokens=sequence_name.split(UNIREF_DELIMITER)
    new_annotation=tokens[ANNOTATION_INDEX]+ANNOTATION_DELIMITER+str(gene_length)
    tokens[ANNOTATION_INDEX]=new_annotation
    
    new_sequence_name=UNIREF_DELIMITER.join(tokens)
    
    return new_sequence_name

def check_sequence(sequence_name,uniref_pathways,filter,filter_ids,include_ids):
    """
    Check the name of the sequence to determine if it should
    be removed from the final set
    Also check that the uniref id is not included in a pathways database
    """
    
    filter_value="store"
    #check for the filter sequence
    tokens=sequence_name.split(UNIREF_DELIMITER)
    if filter:
        if re.search(filter, sequence_name, re.IGNORECASE):
            #check that the id is not in a pathway
            if not tokens[UNIREF_ID_INDEX] in uniref_pathways:
                filter_value="remove"
            else:
                filter_value="store_pathways"
    if filter_ids:
        if tokens[UNIREF_ID_INDEX] in filter_ids:
            #check that the id is not in a pathway
            if not tokens[UNIREF_ID_INDEX] in uniref_pathways:
                filter_value="remove"
            else:
                filter_value="store_pathways"
    if include_ids:
        if tokens[UNIREF_ID_INDEX] in include_ids:
            filter_value="store"
        else:
            if not tokens[UNIREF_ID_INDEX] in uniref_pathways:
                filter_value="remove"
            else:
                filter_value="store_pathways"
        
    return filter_value

def process_id_list_file(file):
    """
    Read through the file and store ids
    """
    
    # if set, read in the list of ids from the filter file
    filter_ids={}
    # check file exists
    if not os.path.isfile(file):
        sys.exit("List file can not be found: " + file)
        
    if not os.access(file, os.R_OK):
        sys.exit("List file is not readable: " +  file)
        
    # read through file for ids
    file_handle_read=open(file,"r")
    for line in file_handle_read:
        # bypass comment lines
        if not re.search("^#",line):
            # file format should be a single column of ids
            data=line.rstrip().split(FILTER_ID_DELIMITER)
            filter_ids[data[FILTER_ID_INDEX]]=1
                
    file_handle_read.close()
   
    return filter_ids
    

def process_fasta_file(fasta_file,uniref_pathways,verbose,filter,exclude_list_file, include_list_file):
    """
    Read through the fasta file, storing headers and sequences
    Filter to remove sequences based on string or list
    Add gene lengths to annotations
    """

    if verbose:
        if filter:
            print("Filter set on, using string: " + filter)
        if exclude_list_file:
            print("Using ids to exclude from file: " + exclude_list_file)
        if include_list_file:
            print("Using ids to include from file: " + include_list_file)
            
    # get set of filter ids if file set
    filter_ids={}
    if exclude_list_file:
        filter_ids=process_id_list_file(exclude_list_file)
        if verbose:
            print("Total ids to filter: " + str(len(filter_ids.keys())))
            
    # get set of include ids if file set
    include_ids={}
    if include_list_file:
        include_ids=process_id_list_file(include_list_file)
        if verbose:
            print("Total ids to include: " + str(len(include_ids.keys())))
        
    # check file exists
    if not os.path.isfile(fasta_file):
        sys.exit("Fasta file can not be found: " + fasta_file)
        
    if not os.access(fasta_file, os.R_OK):
        sys.exit("Fasta file is not readable: " +  fasta_file)
    
    file_handle_read=open(fasta_file,"r")

    line=file_handle_read.readline()
    
    # check the file is of fasta format
    if not re.match(">", line):
        file_handle_read.close()
        sys.exit("File provided is not of fasta format: " + fasta_file)
    
    fasta={}
    sequence_name=""
    sequence=""
    
    filtered_count=0
    not_filtered_as_in_pathway=0
    while line:
        if re.match(">",line):
            if sequence_name:
                # store the sequence just read
                # if it is not to be filtered
                filter_type=check_sequence(sequence_name,uniref_pathways,filter,filter_ids,include_ids)
                if "store" in filter_type:
                    new_sequence_name=add_gene_length(sequence_name, sequence)
                    fasta[new_sequence_name]=sequence
                    if "pathways" in filter_type:
                        not_filtered_as_in_pathway+=1
                else:
                    filtered_count+=1
            sequence_name=line.rstrip('\n').replace(">","")
            sequence=""
        else:
            sequence+=line.rstrip('\n')
        line=file_handle_read.readline()

    # store the last sequence if not filtered
    if sequence_name:
        filter_type=check_sequence(sequence_name,uniref_pathways,filter,filter_ids,include_ids)
        if "store" in filter_type:
            new_sequence_name=add_gene_length(sequence_name, sequence)
            fasta[new_sequence_name]=sequence
            if "pathways" in filter_type:
                not_filtered_as_in_pathway+=1
        else:
            filtered_count+=1

    file_handle_read.close()

    if verbose:
        print("Total sequences filtered: " + str(filtered_count))
        print("Total sequences not filtered as in pathway: " + str(not_filtered_as_in_pathway))
        print("Total sequences remaining: " +  str(len(fasta.keys())))

    return fasta

def find_exe_in_path(exe):
    """
    Check that an executable exists in $PATH
    """
    
    paths = os.environ["PATH"].split(os.pathsep)
    for path in paths:
        fullexe = os.path.join(path,exe)
        if os.path.exists(fullexe):
            if os.access(fullexe,os.X_OK):
                return True
    return False


def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Create UniRef database for HUMAnN2\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-v","--verbose", 
        help="additional output is printed\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "-i","--input",
        help="the UniRef fasta file to read\n",
        required=True)
    parser.add_argument(
        "-o","--output",
        help="the UniRef database to write\n",
        required=True)
    parser.add_argument(
        "-f","--filter",
        help="string to use for filtering (example: uncharacterized)\n")
    parser.add_argument(
        "--exclude-list",
        help="file of id list to use for filtering (example: id_list.tsv)\n")
    parser.add_argument(
        "--include-list",
        help="file of id list to use for filtering (example: id_list.tsv)\n")
    parser.add_argument(
        "-d","--format-database",
        choices=["fasta","rapsearch","diamond"],
        default="fasta",
        help="format of output files (default: fasta)\n")

    return parser.parse_args()


def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    # Check for the software to create the database
    if args.format_database == "rapsearch":
        database_software="prerapsearch"
    elif args.format_database == "diamond":
        database_software="diamond"
    else:
        database_software=None
    
    if database_software:
        if not find_exe_in_path(database_software):
            sys.exit("Could not find the location of the software to create the" +
                " database: " + database_software)
    
    # Find the output directory
    args.output=os.path.abspath(args.output)
    output_dir=os.path.dirname(args.output)
    
    if not os.path.isdir(args.output):
        print("Creating output directory: " + args.output)
        try:
            os.mkdir(args.output)
        except EnvironmentError:
            sys.exit("Unable to create output directory.")
    
    # Create temp folder
    temp_dir=tempfile.mkdtemp( 
        prefix='humann2_create_uniref_database_', dir=output_dir)
    if args.verbose:
        print("Temp folder created: " + temp_dir)
        
    if args.verbose:
        print("Storing uniref data from pathways files")
        
    uniref_pathways=store_pathways(args.verbose)
        
    if args.verbose:
        print("Annotating fasta file: " + args.input)
        
    fasta=process_fasta_file(args.input,uniref_pathways,args.verbose,args.filter,args.exclude_list, args.include_list)
    
    # create a temp file of the fasta sequences
    input_file_basename=os.path.splitext(os.path.basename(args.input))[0]
    output_basename=os.path.join(args.output,input_file_basename)
    if not args.format_database == "fasta":
        file_out, temp_fasta_file=tempfile.mkstemp(dir=temp_dir)
    else:
        # if output file is formatted as fasta, then do not write as temp file
        temp_fasta_file=os.path.join(args.output,input_file_basename+".filtered.fasta")
    
    try:
        file_out=open(temp_fasta_file,"w")
    except EnvironmentError:
        sys.exit("ERROR: Unable to open fasta file to write: " + temp_fasta_file)
    
    if args.verbose:
        print("Printing annotated fasta file: " + temp_fasta_file)
    
    for name in fasta.keys():
        file_out.write(">"+name+'\n'+fasta[name]+'\n')
    file_out.close()
        
    # create the database
    cmd=[]
    if database_software == "prerapsearch":
        cmd=[database_software,"-d",temp_fasta_file,"-n",args.output]
    elif database_software == "diamond":
        cmd=[database_software,"makedb","--in",temp_fasta_file,"--db",output_basename]
        
    if cmd:
        print("Formatting database from fasta file")
        if args.verbose:
            print(" ".join(cmd))
        try:
            p_out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        except (EnvironmentError, subprocess.CalledProcessError):
            sys.exit("Error while running the database software command: " + " ".join(cmd))
    
    
    if args.verbose:
        print("Deleting temp files in temp folder: " + temp_dir)
        
    # deleting temp folder with all files
    shutil.rmtree(temp_dir)
    
    print("Database created: " + args.output)

if __name__ == "__main__":
    main()
    
