#!/usr/bin/env python

"""
This script will build a custom database, selecting organisms by genus, and using
id mapping if provided

To Run: 
$ ./humann2_build_custom_database.py --input genes.pep --output custom_database 
  --id-mapping legacy_kegg_idmapping.tsv --format diamond --taxonomic-profile max_taxonomic_profile.tsv

"""

import argparse
import sys
import os
import re
import subprocess

from humann2 import config
from humann2 import store

FASTA_ID_START=">"
FASTA_ID_DELIMITER=" "
WORD_DELIMITER=" "
TAXONOMY_DELIMITER="|"
GENUS_IDENTIFIER="g__"

def format_diamond_database(fasta_file, output_folder):
    """ Format the input file into a diamond database """
    
    database=os.path.join(output_folder,os.path.splitext(os.path.basename(fasta_file))[0])
    
    command=["diamond","makedb","--in",fasta_file,"--db",database]
    
    print("RUNNING: "+" ".join(command))
    try:
        subprocess.check_call(command)
    except (EnvironmentError, subprocess.CalledProcessError):
        sys.exit("ERROR: Unable to execute diamond database build")
        
    print("Created diamond database: " + database)
        
def filter_fasta_file(fasta_file, output_folder, genus, id_mapping):
    """ Read through each sequence in the fasta file filtering by genus """
    
    try:
        file_handle=open(fasta_file, "rt")
    except EnvironmentError:
        sys.exit("ERROR: Unable to read input fasta file: " + fasta_file)
        
    # Create a file in the output folder that is the filtered fasta file
    try:
        input_file_basename=os.path.splitext(os.path.basename(fasta_file))
        new_file=input_file_basename[0]+".filtered"+input_file_basename[1]
    except IndexError:
        new_file=os.path.basename(fasta_file)
        
    new_file=os.path.join(output_folder,new_file)
    
    # Try to open the new file to write
    try:
        file_handle_write=open(new_file,"w")
    except EnvironmentError:
        sys.exit("ERROR: Unable to write output fasta file: " + new_file)
    
    print("Writing filtered fasta file to " + new_file)
    write_sequence=False
    total_sequences_written=0
    for line in file_handle:
        if line[0] == FASTA_ID_START:
            # Get the possible mapping id
            mapping_id=line.split(FASTA_ID_DELIMITER)[0][1:]
            # Apply id mapping if included
            if mapping_id in id_mapping:
                sequence_information=id_mapping[mapping_id][-1].lower().split(WORD_DELIMITER)
            else:
                sequence_information=line[1:].lower().split(WORD_DELIMITER)
                
            # Determine if write sequence based on genus information
            write_sequence=False
            if genus:
                # Check if the sequence taxonomy is an included genus
                # The taxonomy information is a string of words and any word
                # can include the genus information
                if list(filter(lambda x: x in genus, sequence_information)):
                    write_sequence=True
            else:
                write_sequence=True
                
            if write_sequence:
                file_handle_write.write(line)
                total_sequences_written+=1
        else:
            if write_sequence:
                file_handle_write.write(line)
    
    file_handle.close()
    file_handle_write.close()
    print("Finished writing filtered fasta file")

    return new_file, total_sequences_written

def process_taxonomic_profile(taxonomic_profile, abundance_threshold):
    """ Store those genus which pass the threshold """
    
    # The taxonomic profile will be formatted like the metaphlan2 output file
    
    try:
        file_handle=open(taxonomic_profile, "rt")
    except EnvironmentError:
        sys.exit("ERROR: Unable to read taxonomic profile: " + taxonomic_profile)
        
    print("Reading taxonomic profile")
    genus_found=set()
    for line in file_handle:

        # search for the lines that have the genus-level information
        if re.search(GENUS_IDENTIFIER, line):
            # check threshold
            data=line.rstrip().split("\t")
            try:
                read_percent=float(data[1])
            except (IndexError, ValueError):
                read_percent=0
                
            if read_percent >= abundance_threshold:
                genus=list(filter(lambda x: GENUS_IDENTIFIER in x, data[0].split(TAXONOMY_DELIMITER)))[0].replace(GENUS_IDENTIFIER,"").lower()
                if not genus in genus_found:
                    print("Adding genus: " + genus)
                genus_found.add(genus)
        
    file_handle.close()
    
    return genus_found
        
def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Create a custom database file\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i","--input",
        help="the fasta input file\n",
        required=True)
    parser.add_argument(
        "-o","--output",
        help="the output folder\n",
        required=True)
    parser.add_argument(
        "--id-mapping",
        help="the file mapping fasta ids to taxonomy")
    parser.add_argument(
        "--taxonomic-profile",
        help="the file containing the taxonomic profile")
    parser.add_argument(
        "--format",
        help="the final database format",
        default="fasta",
        choices=["fasta","diamond"])
    parser.add_argument(
        "--genus-abundance-threshold",
        help="the minimum abundance for a genus to be included in the database",
        type=float,
        default=config.prescreen_threshold)

    return parser.parse_args()

def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    input_file=os.path.abspath(args.input)
    output_dir=os.path.abspath(args.output)
    
    if not os.path.isdir(output_dir):
        print("Creating output directory: " + output_dir)
        try:
            os.mkdir(output_dir)
        except EnvironmentError:
            sys.exit("Unable to create output directory.")
            
    if not os.access(output_dir, os.W_OK):
        sys.exit("The output directory provided is not writeable. Please update the permissions.")
        
    genus={}
    if args.taxonomic_profile:
        # If taxonomic profile is provided, read in the genus to include in the database
        genus=process_taxonomic_profile(args.taxonomic_profile, args.genus_abundance_threshold)
        if not genus:
            sys.exit("ERROR: Zero genus were identified which pass the abundance threshold.")
        
    id_mapping={}
    if args.id_mapping:
        # If the id mapping file is provided, map the genus to the fasta ids
        print("Reading id mapping file")
        id_mapping=store.store_id_mapping(args.id_mapping)
        
    # Filter the input file to only include those genus selected
    # Apply id mapping if provided
    filtered_fasta_file, total_sequences=filter_fasta_file(input_file, output_dir, genus, id_mapping)
    
    if total_sequences == 0:
        sys.exit("ERROR: When filtering fasta file zero sequences were written.")
    elif genus:
        print("Total sequences included after filtering: " + str(total_sequences))
    
    # If set, format as diamond database
    if args.format == "diamond":
        format_diamond_database(filtered_fasta_file, output_dir)

if __name__ == "__main__":
    main()
