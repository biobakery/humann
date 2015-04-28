#!/usr/bin/env python

"""
Join a set of taxonomic profiles into a single profile

Dependencies: None

To Run: 
$ ./join_taxonomic_profiles.py -i <input_dir> -o <taxonomic_profile.tsv>

"""

import argparse
import sys
import os

TABLE_DELIMITER="\t"
TAXON_DELIMITER="|"
TAXON_INDEX=0
VALUE_INDEX=1
        
def join_profile_tables(profile_tables,output,verbose):
    """
    Join the taxonomic profiles into a single profile
    """
    
    taxon_data={}
    taxon_levels={}
    for table in profile_tables:
        # tables are expected to have two columns as follows
        # taxonomy \t percent
        if verbose:
            print("Processing file: " + table)
        with open(table) as file_handle:
            line=file_handle.readline()
            while line:
                data=line.rstrip().split(TABLE_DELIMITER)
                try:
                    taxon=data[TAXON_INDEX]
                    value=float(data[VALUE_INDEX])
                except (ValueError, IndexError):
                    taxon=""
                    
                if taxon:
                    taxon_data[taxon]=taxon_data.get(taxon,0)+value
                    level=taxon.count(TAXON_DELIMITER)
                    if level in taxon_levels:
                        taxon_levels[level].add(taxon)
                    else:
                        taxon_levels[level]=set(taxon)
                    
                line=file_handle.readline()
    
    try:
        file_handle=open(output,"w")
    except EnvironmentError:
        sys.exit("Unable to write file: " + file)  
        
    # write out the taxons by level
    total_profiles=len(profile_tables)
    last_level=max(taxon_levels.keys())
    
    level=0
    while level <= last_level:
        for taxon in taxon_levels.get(level,set()):
            value=taxon_data.get(taxon,0)/float(total_profiles)
            if value:
                file_handle.write(TABLE_DELIMITER.join([taxon,str(value)])+"\n")
        level+=1
    
    file_handle.close()

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Join taxonomic profiles\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-v","--verbose", 
        help="additional output is printed\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "-i","--input",
        help="the directory of taxonomic profiles\n",
        required=True)
    parser.add_argument(
        "-o","--output",
        help="the joined taxonomic profile to write\n",
        required=True)
    parser.add_argument(
        "--file_name",
        help="only join taxonomic profiles with this string included in the file name")

    return parser.parse_args()


def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    # check for format of the gene tables
    input_dir=os.path.abspath(args.input)
    
    # check the directory exists
    if not os.path.isdir(input_dir):
        sys.exit("The input directory provided can not be found." + 
            "  Please enter a new directory.")
    
    profile_tables=[]
    file_list=os.listdir(input_dir)
    
    # filter out files which do not meet the name requirement if set
    if args.file_name:
        reduced_file_list=[]
        for file in file_list:
            if re.search(args.file_name,file):
                reduced_file_list.append(file)
        file_list=reduced_file_list      
    
    for file in file_list:
        profile_tables.append(os.path.join(input_dir,file))
        
    args.output=os.path.abspath(args.output)
        
    if args.verbose:
        print("Joining taxonomic profiles")
        
    join_profile_tables(profile_tables,args.output,args.verbose)
    
    print("Taxonomic profile created: " + args.output)

if __name__ == "__main__":
    main()
