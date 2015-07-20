#!/usr/bin/python
"""
Filter the structured pathways file by reactions

To Run: 
$ ./filter_pathways.py --input-pathways metacyc_structured_pathways --input-reactions metacyc_reactions.uniref --output metacyc_structured_pathways_filtered

"""

import sys

try:
    import argparse
except ImportError:
    sys.exit("ERROR: Please upgrade to at least python v2.7")
    
import re
import os

COLUMN_DELIMITER="\t"
ITEM_DELIMITER=","
OPTIONAL_TAG="-"

def read_reactions(reaction_file):
    """
    Read the reactions file to get EC numbers for reactions
    """
    
    try:
        file_handle=open(reaction_file)
    except EnvironmentError:
        sys.exit("Unable to open reaction file: %s", reaction_file)
        
    ec_numbers={}
    
    line=file_handle.readline()
    reaction=""
    for line in file_handle.readlines():
        data=line.rstrip().split(COLUMN_DELIMITER)
        if len(data) >= 3:
            reaction=data[0]
            ecs=data[1].split(ITEM_DELIMITER)
            ec_numbers[reaction]=ecs
        
    file_handle.close()
        
    return ec_numbers

def filter_pathways(pathways_file, reactions_to_ec, output_file):
    """
    Filter the pathways based on the reactions (using the EC numbers)
    Filter pathways with >10% ECs that are not fully specified (ie 4 levels)
    Filter pathways with less than 4 specific ECs
    """
    
    try:
        file_handle=open(pathways_file)
    except EnvironmentError:
        sys.exit("Unable to open pathways file: %s", pathways_file)
        
    try:
        file_handle_out=open(output_file,"w")
    except EnvironmentError:
        sys.exit("Unable to open output file: %s", output_file)
        
    line=file_handle.readline()
    while line:
        data=line.rstrip().split(COLUMN_DELIMITER)
        total_well_specified_ecs=0
        distinct_ecs=set()
        reactions=set()
        reactions_without_mappings=set()

        for item in data[1].split(" "):
            if not item in ["(","+",",",")",""]:
                
                # get the ec numbers for the reaction
                reactions.add(item)
                ec_numbers=reactions_to_ec.get(item,None)
                if ec_numbers is None:
                    # store those reactions without mappings to uniref50/90
                    reactions_without_mappings.add(item)
                else:
                    # count the number of ECs that are well specified
                    for ec in ec_numbers:
                        ec_level=len(ec.split("."))
                        # check if well specified (ie level >=4)
                        # only add one ec per reaction
                        if ec_level >= 4:
                            total_well_specified_ecs+=1
                            break
                    # count the number of ECs that are distinct
                    for ec in ec_numbers:
                        # add at most 1 EC for reaction
                        if not ec in distinct_ecs:
                            distinct_ecs.add(ec)
                            break
                    
        # check for any ECs that are of the same group with different numbers of levels
        total_distinct_ecs=0
        distinct_ecs=list(distinct_ecs)
        for index, ec in enumerate(distinct_ecs):
            ec_level=len(ec.split("."))
            match=False
            for ec2 in distinct_ecs[(index+1):]:
                ec2_level=len(ec2.split("."))
                # compare ecs at the same level
                if ec2_level < ec_level:
                    reduced_ec=".".join(ec.split(".")[0:ec2_level])
                    if reduced_ec == ec2:
                        match=True
                        break
                elif ec2_level > ec_level:
                    reduced_ec=".".join(ec2.split(".")[0:ec_level])
                    if reduced_ec == ec:
                        match=True
                        break
            if not match:
                total_distinct_ecs+=1   
                    
        # check if this pathway should be filtered
        total_mappable_reactions=len(reactions)-len(reactions_without_mappings)
        if ( total_mappable_reactions >=4 and total_well_specified_ecs/(total_mappable_reactions*1.0) >= 0.10 and
            ( total_distinct_ecs >= 4 or total_distinct_ecs/(total_mappable_reactions*1.0) >= 0.75)):
            # if the pathway is not filtered then make those unmappable reactions optional
            for reaction in reactions_without_mappings:
                line=line.replace(" "+reaction+" "," "+OPTIONAL_TAG+reaction+" ")
                
            file_handle_out.write(line)
            
        line=file_handle.readline()    
    
    file_handle.close()
    file_handle_out.close()
    
def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Filter structured MetaCyc pathways file\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-v","--verbose", 
        help="additional output is printed\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "--input-pathways",
        help="the MetaCyc structured pathways file\n",
        required=True)
    parser.add_argument(
        "--input-reactions",
        help="the MetaCyc reactions file\n",
        required=True)
    parser.add_argument(
        "-o","--output",
        help="the file to write the filtered structured pathways\n",
        required=True)

    return parser.parse_args()


def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
     
    input_file_pathways=os.path.abspath(args.input_pathways)
    input_file_reactions=os.path.abspath(args.input_reactions)
    output_file=os.path.abspath(args.output)
    
    if args.verbose:
        print("Reading MetaCyc reactions file")
    
    reactions_to_ec=read_reactions(input_file_reactions)
    
    if args.verbose:
        print("Total reactions with EC numbers: " + str(len(reactions_to_ec)))
        
    if args.verbose:
        print("Processing MetaCyc structured pathways file")
    
    filter_pathways(input_file_pathways, reactions_to_ec, output_file)
        
if __name__ == "__main__":
    main()
        
     
