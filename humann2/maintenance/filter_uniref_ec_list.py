#!/usr/bin/env python

"""
Filter the uniref ec list.

To Run: 
$ ./filter_uniref_ec_list.py --input map_EC_to_triplet_AC_U50_U90_Swissprot_and_Trembl.txt --uniref 50 --output uniref50_ec4_list.txt --min-ec-level 4 

"""

import argparse
import sys

def read_uniref_ec_input(file, uniref_type, min_ec_level):
    """ Read the input file, storing unirefs based on ec level """
    
    uniref_identifier="UniRef"+uniref_type
    
    filtered_uniref_ids=set()
    with open(file) as file_handle:
        for line in file_handle:
            data=line.rstrip().split("\t")
            if data[0].count(".") >= int(min_ec_level)-1:
                for item in data:
                    if uniref_identifier in item:
                        filtered_uniref_ids.add(item)
    return filtered_uniref_ids

def write_output_file(uniref, file):
    """ Create a file with the list of uniref ids """
    
    with open(file,"w") as file_handle:
        for id in uniref:
            file_handle.write(id+"\n")

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Filter UniRef EC list\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i","--input",
        help="the UniRef EC triplet file [map_EC_to_triplet_AC_U50_U90_Swissprot_and_Trembl.txt]\n",
        required=True)
    parser.add_argument(
        "-o","--output",
        help="the filtered UniRef file to write\n",
        required=True)
    parser.add_argument(
        "-u","--uniref",
        help="the UniRef type\n",
        choices=["50","90"],
        required=True)
    parser.add_argument(
        "-m","--min-ec-level",
        help="the minimum EC level for filtering\n",
        choices=["1","2","3","4"],
        required=True)

    return parser.parse_args()


def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    # Read in the input file, filtering by ec level
    filtered_uniref_ids=read_uniref_ec_input(args.input, args.uniref, args.min_ec_level)
    
    # Write the list of uniref ids
    write_output_file(filtered_uniref_ids, args.output)

if __name__ == "__main__":
    main()
    