#!/usr/bin/python
"""
Filter the structured pathways file by reactions

To Run: 
$ ./filter_pathways.py --input-pathways metacyc_structured_pathways --input-reactions reactions.dat --output metacyc_structured_pathways_filtered

"""

import sys

try:
    import argparse
except ImportError:
    sys.exit("ERROR: Please upgrade to at least python v2.7")
    
import re
import os

DELIMITER="\t"
COMMENT_LINE="#"
METACYC_ID="UNIQUE-ID"
METACYC_ID_DELIMITER=" - "
METACYC_EC_ID="EC-NUMBER"

def read_metacyc_reactions(reaction_file):
    """
    Reaction the MetaCyc reactions file to get EC numbers for reactions
    """
    
    try:
        file_handle=open(reaction_file)
    except EnvironmentError:
        sys.exit("Unable to open reaction file: %s", reaction_file)
        
    ec_numbers={}
    
    line=file_handle.readline()
    reaction=""
    while line:
        if not re.match(COMMENT_LINE, line):
            # store the reaction id
            if re.match(METACYC_ID, line):
                # store the latest set of reaction information
                if reaction and ec:
                    ec_numbers[reaction]=ec
                ec=[]
                reaction=line.rstrip().split(METACYC_ID_DELIMITER)[-1].strip()
            
            # find the ec
            elif re.match(METACYC_EC_ID,line):
                ec.append(line.rstrip().split(METACYC_ID_DELIMITER)[-1].strip())
        line=file_handle.readline()
        
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
        data=line.rstrip().split(DELIMITER)
        total_well_specified_ecs=0
        distinct_ecs=set()
        reactions=set()

        for item in data[1].split(" "):
            if not item in ["(","+",",",")",""]:

                # count the number of ECs that are well specified
                ec_numbers=reactions_to_ec.get(item,[])
                reaction=item
                if not ec_numbers:
                    ec_numbers=reactions_to_ec.get(item[1:],[])
                    reaction=item[1:]
                    
                reactions.add(reaction)
                    
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
                    
        # check if this pathway should be filtered
        if total_well_specified_ecs/(len(reactions)*1.0) >= 0.10 and len(distinct_ecs) >= 4:
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
    
    reactions_to_ec=read_metacyc_reactions(input_file_reactions)
    
    if args.verbose:
        print("Total reactions with EC numbers: " + str(len(reactions_to_ec)))
        
    if args.verbose:
        print("Processing MetaCyc structured pathways file")
    
    filter_pathways(input_file_pathways, reactions_to_ec, output_file)
        
if __name__ == "__main__":
    main()
        
     
