#!/usr/bin/env python

"""
Split a stratified table into two tables: stratified and unstratified

Dependencies: None

To Run: 
$ python split_stratified_table.py -i genefamiles.tsv -o $OUTPUT_FOLDER

"""

import argparse
import sys
import os

STRATIFICATION_DELIMITER="|"
        
def split_table(input,output):
    """
    Split stratified table
    """
    
    # name the two output files
    input_file_name, input_file_extension = os.path.splitext(os.path.basename(input))
    
    output_stratified = os.path.join(output, input_file_name + "_stratified" + input_file_extension)
    output_unstratified = os.path.join(output, input_file_name + "_unstratified" + input_file_extension)
    
    # try to open the output files
    try:
        output_stratified_handle = open(output_stratified, "w")
    except EnvironmentError:
        sys.exit("Unable to open output file for writing: " + output_stratified)

    try:
        output_unstratified_handle = open(output_unstratified, "w")
    except EnvironmentError:
        sys.exit("Unable to open output file for writing: " + output_unstratified)
    
    with open(input, "rt") as file_handle:
        header=file_handle.readline()
        
        # write the header to both output files
        output_stratified_handle.write(header)
        output_unstratified_handle.write(header)
        
        for line in file_handle:
            if STRATIFICATION_DELIMITER in line:
                output_stratified_handle.write(line)
            else:
                output_unstratified_handle.write(line)
                        
    # close the output files
    output_stratified_handle.close()
    output_unstratified_handle.close()
    
    return output_stratified, output_unstratified

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Split stratified table\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i","--input",
        help="the stratified input table\n",
        required=True)
    parser.add_argument(
        "-o","--output",
        help="the output folder\n",
        required=True)

    return parser.parse_args()


def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    input=os.path.abspath(args.input)
    output=os.path.abspath(args.output)
    
    if not os.path.isfile(input):
        sys.exit("The input file provided can not be found." + 
            "  Please enter a new file.")

    # if output folder does not exist, then create
    if not os.path.isdir(output):
        try:
            os.makedirs(output)
        except EnvironmentError:
            sys.exit("ERROR: Unable to create output directory: " + output)
    
    output_stratified, output_unstratified = split_table(args.input,args.output)
    
    print("Split stratified tables created:\n" + "\n".join([output_stratified, output_unstratified]))

if __name__ == "__main__":
    main()
