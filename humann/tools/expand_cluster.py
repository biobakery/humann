#! /usr/bin/env python

import argparse
import sys
import os
import gzip

try:
    from humann import config
    from humann.tools import util
except ImportError:
    sys.exit("CRITICAL ERROR: Unable to find the HUMAnN python package." +
        " Please check your install.") 

description = """
HUMAnN utility for expanding clustered table features
=============================================
Given a table of UniRef90 values and a specific UniRef90,
create a table subset that includes all of the UniRef90s
that cluster with the selected UniRef90 in a UniRef50 set.
"""


# path to the mapping file
MAPPING_FILE = os.path.join(config.utility_mapping_database,"map_uniref50_uniref90.txt.gz")
IDENTIFIER = "UniRef90_"

def arg_parse():
    parser = argparse.ArgumentParser(
        description=description, 
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument( 
        "-i", "--input", 
        default=None,
        help="UniRef90 gene families table (tsv format)",
        required=True
        )  
    parser.add_argument( 
        "-g", "--gene", 
        default=None,
        help="Gene family (UniRef90) of interest",
        required=True
        )
    parser.add_argument( 
        "-o", "--output", 
        default=None,
        help="Path for modified output table (tsv format)",
        required=True
        )
    parser.add_argument( 
        "-v", "--verbose", 
        help="Write status information",
        action="store_true"
        )
    args = parser.parse_args()
    return args

def read_mapping(mapping_file):
    # get the uniref90 to uniref50 mappings
    mapping={}
    try:
        for line in gzip.open(mapping_file,"rt"):
            data=line.rstrip().split("\t")
            for uniref in data[1:]:
                if uniref in mapping:
                    mapping[uniref]+=[data[0]]
                else:
                    mapping[uniref]=[data[0]]
    except EnvironmentError:
        sys.exit("Error: Unable to read utility mapping file "+mapping_file)

    return mapping

def match_uniref50s(requested_uniref50s, uniref90, mapping):
    # check if the requested uniref50s are included
    
    match=False
    for uniref50 in mapping.get(uniref90,[]):
        if uniref50 in requested_uniref50s:
            match=True
            break

    return match

def write_genes(uniref50s, mapping, input_file, output_file, verbose):
    # read through the gene family file, writing those uniref90s of interest

    try:
        file_handle=open(output_file,"wt")
    except EnvironmentError:
        sys.exit("Error: Unable to write to output file "+output_file)

    open_function=open
    if input_file.endswith(".gz"):
        open_function=open.gzip

    if verbose:
        print("Writing output file")

    header=""
    for line in open_function(input_file):
         if not header:
             header=line
             file_handle.write(header)
         else:
             uniref90=line.split("\t")[0].split("|")[0]
             if match_uniref50s(uniref50s, uniref90, mapping):
                 file_handle.write(line)

    file_handle.close()

def main( ):
    args = arg_parse()

    # read in the mapping filea
    if args.verbose:
        print("Reading mapping file")
    mapping=read_mapping(MAPPING_FILE)
    if args.verbose:
        print("Done reading mapping file")

    # allow for just the cluster id
    if not args.gene.startswith(IDENTIFIER):
        args.gene=IDENTIFIER+args.gene

    # get uniref50 of interest
    try:
        uniref50s=mapping[args.gene]
    except KeyError:
        sys.exit("Error: No UniRef50s associated with the UniRef90 of interest "+args.gene)

    if args.verbose:
        print("Found UniRef50s associated with selected UniRef90")
        print(uniref50s)

    write_genes(uniref50s, mapping, args.input, args.output, args.verbose)

if __name__ == "__main__":
    main()
