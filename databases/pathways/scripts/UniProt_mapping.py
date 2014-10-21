#!/usr/bin/env python

"""
UniProt_mapping will create gzipped files of UniProt mapping sets for UniRef50 or UniRef90 ids.
It will download the uniprot mapping text file and create the subset map requested.
Mapping output will be a UniProtKB-AC id and the corresponding id_type {"UniRef50","UniRef90"} selected.
"""

import argparse
import urllib
import os
import gzip
import sys

UNIPROT_MAPPING_URL="ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz"
DOWNLOADED_UNIPROT_FILE="idmapping.dat.gz"
ID_INDEX=1
UNIPROT_INDEX=0
OTHER_INDEX=2

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i","--id_type", 
        help="the ids selected for mapping\n",
        choices=["UniRef90","UniRef50"],
        required=True)
    parser.add_argument(
        "--uniprot_mapping",
        help="the location of the full uniprot mapping file [DEFAULT: file downloaded]")
    parser.add_argument(
        "-o","--output",
        help="the file to write the gzipped output",
        required=True)

    return parser.parse_args()

def download_uniprot_mapping():
    """
    Download the file at the url and extract
    """
    
    print("Download from URL:" + UNIPROT_MAPPING_URL)
    print("This file is large (>4G) so this may take some time.")

    try:
        file, headers = urllib.urlretrieve(UNIPROT_MAPPING_URL,DOWNLOADED_UNIPROT_FILE)
    except EnvironmentError:
        sys.exit("Unable to download URL: " + UNIPROT_MAPPING_URL)
   
def main():
    
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
            
    # if the uniprot mapping file is not provided
    # and is not found locally then download
    if not args.uniprot_mapping:
        if not os.path.isfile(DOWNLOADED_UNIPROT_FILE):
            download_uniprot_mapping()        
            
    # open the output file to write in gzipped format
    try:
        file_out=gzip.open(args.output,'wb')
    except EnvironmentError:
        sys.exit("Unable to open the output file for writing: " + args.output)
        
    # open and read the uniprot mapping file 
    try:
        if args.uniprot_mapping:
            file_in=gzip.open(args.uniprot_mapping,'rb')
        else:
            file_in=gzip.open(DOWNLOADED_UNIPROT_FILE,'rb')
    
        for line in file_in:
            data=line.strip().split("\t")
            if data[ID_INDEX] == args.id_type:
                file_out.write(data[UNIPROT_INDEX]+"\t"+data[OTHER_INDEX]+"\n")
            
        file_in.close()
        file_out.close()
    except EnvironmentError:
        sys.exit("Error reading/writing input/output files.")
    
if __name__ == "__main__":
    main()
