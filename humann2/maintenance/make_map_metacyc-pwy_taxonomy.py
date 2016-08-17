#!/usr/bin/python
"""
Create a file of taxonomy information from MetaCyc flat pathways file (pathways.dat)

Dependencies: MetaCyc flat files and NCBI taxonomy files

To Run: 
$ python make_map_metacyc-pwy_taxonomy.py --input-pathways pathways.dat --input-taxonomy-names names.dmp 
   --input-taxonomy-division division.dmp --input-taxonomy-nodes nodes.dmp 
   --output-names metacyc_pathways_taxonomy_names --output-division metacyc_pathways_taxonomy_division
"""

import sys

try:
    import argparse
except ImportError:
    sys.exit("ERROR: Please upgrade to at least python v2.7")

import os
import re

COMMENT_LINE="#"
METACYC_ID="UNIQUE-ID"
METACYC_PATHWAY_ID="PWY"
METACYC_ID_DELIMITER=" - "
METACYC_TAXONOMIC_RANGE_ID="TAXONOMIC-RANGE"
METACYC_SPECIES_ID="SPECIES"

DELIMITER="\t"

NCBI_DELIMITER="\t|\t"
NCBI_NAMES_ID_COLUMN=0
NCBI_NAMES_NAME_COLUMN=1
NCBI_NAMES_TYPE_COLUMN=3

NCBI_DIVISION_ID_COLUMN=0
NCBI_DIVISION_NAME_COLUMN=2

NCBI_NODES_ID_COLUMN=0
NCBI_NODES_DIVISION_COLUMN=4

def write_output(metacyc_pathways, taxonomy_names, taxonomy_division, 
    taxonomy_nodes, output_names_file, output_division_file):
    """
    Write the output file of pathways with taxonomy 
    """
    
    try:
        file_handle_names=open(output_names_file,"w")
    except EnvironmentError:
        sys.exit("Unable to write to file: " + output_names_file)
    
    try:
        file_handle_division=open(output_division_file,"w")
    except EnvironmentError:
        sys.exit("Unable to write to file: " + output_division_file)
    
    all_division=set()
    for pathway,taxa in metacyc_pathways.items():
        names=set()
        division=set()
        for taxon in taxa:
            if taxon in taxonomy_names:
                names.add(re.sub('[" "|/]',"_",taxonomy_names[taxon]))
            if taxon in taxonomy_nodes:
                division.add(taxonomy_division[taxonomy_nodes[taxon]])
        
        all_division.update(division)
        file_handle_names.write(DELIMITER.join([pathway]+list(names))+"\n")
        file_handle_division.write(DELIMITER.join([pathway]+list(division))+"\n")
        
    file_handle_names.close()
    file_handle_division.close()
    
    print(all_division)
        

def read_ncbi_nodes(nodes_file):
    """
    Read in the nodes for the taxa
    """
    
    taxonomy_nodes={}
    
    try:
        file_handle=open(nodes_file,"r")
        line=file_handle.readline()
    except EnvironmentError:
        sys.exit("Unable to read file: " + nodes_file)
        
    while line:
        data=line.rstrip().split(NCBI_DELIMITER)
        taxonomy_nodes[data[NCBI_NODES_ID_COLUMN]]=data[NCBI_NODES_DIVISION_COLUMN]
        
        line=file_handle.readline()
        
    file_handle.close()
    
    return taxonomy_nodes

def read_ncbi_division(division_file):
    """
    Read in the division for the taxa
    """
    
    taxonomy_division={}
    
    try:
        file_handle=open(division_file,"r")
        line=file_handle.readline()
    except EnvironmentError:
        sys.exit("Unable to read file: " + division_file)
        
    while line:
        data=line.rstrip().split(NCBI_DELIMITER)
        taxonomy_division[data[NCBI_DIVISION_ID_COLUMN]]=data[NCBI_DIVISION_NAME_COLUMN]
        
        line=file_handle.readline()
        
    file_handle.close()
    
    return taxonomy_division

def read_ncbi_names(names_file):
    """
    Read in the names for the taxa
    """
    
    taxonomy_names={}
    
    try:
        file_handle=open(names_file,"r")
        line=file_handle.readline()
    except EnvironmentError:
        sys.exit("Unable to read file: " + names_file)
        
    while line:
        data=line.rstrip().split(NCBI_DELIMITER)
        if "scientific name" in data[NCBI_NAMES_TYPE_COLUMN]:
            taxonomy_names[data[NCBI_NAMES_ID_COLUMN]]=data[NCBI_NAMES_NAME_COLUMN]
        
        line=file_handle.readline()
        
    file_handle.close()
    
    return taxonomy_names

def read_metacyc_pathways(metacyc_pathways_file):
    """
    Process the metacyc pathways file to find the taxonomy information
    """
    
    metacyc_pathways={}
    
    try:
        file_handle=open(metacyc_pathways_file,"r")
        line=file_handle.readline()
    except EnvironmentError:
        sys.exit("Unable to read file: " + metacyc_pathways_file)
        
    pathway=""
    taxonomy_range=set()
    species=set()
    while line:
        if not re.match(COMMENT_LINE, line):
            # find the pathway id
            if re.match(METACYC_ID, line):
                if pathway:
                    # record the last pathway
                    if not taxonomy_range:
                        # if the range is not provided then use the species
                        if species:
                            metacyc_pathways[pathway]=species
                    else:
                        metacyc_pathways[pathway]=taxonomy_range
                        
                pathway=""
                taxonomy_range=set()
                species=set()
                pathway=line.rstrip().split(METACYC_ID_DELIMITER)[-1]
                
            # record the taxonomy range
            elif re.match(METACYC_TAXONOMIC_RANGE_ID, line):
                taxonomy_range.add(re.sub("TAX-","",line.rstrip().split(METACYC_ID_DELIMITER)[-1]))
            elif re.match(METACYC_SPECIES_ID, line):
                species.add(re.sub("TAX-","",line.rstrip().split(METACYC_ID_DELIMITER)[-1]))
                    
        line=file_handle.readline()
        
    # record the last pathway
    if not taxonomy_range:
        # if the range is not provided, then use the species
        if species:
            metacyc_pathways[pathway]=species
    else:
        metacyc_pathways[pathway]=taxonomy_range
        
    file_handle.close()
    
    return metacyc_pathways


def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Find taxonomy for MetaCyc pathways file\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-v","--verbose", 
        help="additional output is printed\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "--input-pathways",
        help="the original MetaCyc pathways file (pathways.dat)\n",
        required=True)
    parser.add_argument(
        "--input-taxonomy-names",
        help="the original NCBI taxonomy names file (names.dmp)\n",
        required=True)
    parser.add_argument(
        "--input-taxonomy-division",
        help="the original NCBI taxonomy division file (division.dmp)\n",
        required=True)
    parser.add_argument(
        "--input-taxonomy-nodes",
        help="the original NCBI taxonomy nodes file (nodes.dmp)\n",
        required=True)
    parser.add_argument(
        "--output-names",
        help="the file to write the pathways with taxonomy names\n",
        required=True)
    parser.add_argument(
        "--output-division",
        help="the file to write the pathways with taxonomy division\n",
        required=True)

    return parser.parse_args()


def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    pathways_file=os.path.abspath(args.input_pathways)
    names_file=os.path.abspath(args.input_taxonomy_names)
    division_file=os.path.abspath(args.input_taxonomy_division)
    nodes_file=os.path.abspath(args.input_taxonomy_nodes)
    output_names_file=os.path.abspath(args.output_names)
    output_divison_file=os.path.abspath(args.output_division)

    # read the metacyc pathways file
    if args.verbose:
        print("Processing metacyc pathways data")
    metacyc_pathways=read_metacyc_pathways(pathways_file)
    
    if args.verbose:
        print("Processing the ncbi names data")
    taxonomy_names=read_ncbi_names(names_file)
    
    if args.verbose:
        print("Processing the ncbi division data")
    taxonomy_division=read_ncbi_division(division_file)
    
    if args.verbose:
        print("Processing the ncbi nodes data")
    taxonomy_nodes=read_ncbi_nodes(nodes_file)
    
    if args.verbose:
        print("Writing pathways with taxonomy file")
    write_output(metacyc_pathways, taxonomy_names, taxonomy_division, taxonomy_nodes, 
        output_names_file, output_divison_file)
         

if __name__ == "__main__":
    main()
