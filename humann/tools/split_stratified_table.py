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

from humann.tools import util

STRATIFICATION_DELIMITER="|"
COLUMN_DELIMITER="\t"
        
def split_line(line):
    """
    Split the table line into data tokens
    """
    
    return line.split(COLUMN_DELIMITER)
        
def split_table(input,output):
    """
    Split stratified table
    """
    
    # name the two output files
    input_file_split = os.path.basename(input).split(".")
    if input_file_split[-1] in ["bz2", "gz"]:
        input_file_name=".".join(input_file_split[:-2])
        input_file_extension="."+input_file_split[-2]
    else:
        input_file_name=".".join(input_file_split[:-1])
        input_file_extension="."+input_file_split[-1]        
    
    # check for biom output files (based on input file format)
    biom=False
    if input_file_extension == util.BIOM_FILE_EXTENSION:
        biom=True
    
    output_stratified = os.path.join(output, input_file_name + "_stratified" + input_file_extension)
    output_unstratified = os.path.join(output, input_file_name + "_unstratified" + input_file_extension)
    
    # try to open the output files
    if biom: 
        output_stratified_rows = []
        output_unstratified_rows = []
    else:
        try:
            output_stratified_handle = open(output_stratified, "w")
        except EnvironmentError:
            sys.exit("Unable to open output file for writing: " + output_stratified)
    
        try:
            output_unstratified_handle = open(output_unstratified, "w")
        except EnvironmentError:
            sys.exit("Unable to open output file for writing: " + output_unstratified)
    
    # read the input file
    readlines=util.gzip_bzip2_biom_open_readlines(input)    
    header=next(readlines)
    
    # if this is a biom file, there could be a comment line before the header
    if "Constructed from biom file" in header:
        header=next(readlines)
        
    # write the header to both output files
    if biom:
        output_stratified_rows.append(split_line(header))
        output_unstratified_rows.append(split_line(header))
    else:
        output_stratified_handle.write(header+"\n")
        output_unstratified_handle.write(header+"\n")
        
    for line in readlines:
        if STRATIFICATION_DELIMITER in line:
            if biom:
                output_stratified_rows.append(split_line(line))
            else:
                output_stratified_handle.write(line+"\n")
        else:
            if biom:
                output_unstratified_rows.append(split_line(line))
            else:
                output_unstratified_handle.write(line+"\n")
        
    if biom:
        # write the two biom output files
        util.write_biom(output_stratified, iter(output_stratified_rows))
        util.write_biom(output_unstratified, iter(output_unstratified_rows))
    else:                
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
        help="the stratified input table (tsv, tsv.gzip, tsv.bzip2, or biom format)\n",
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
