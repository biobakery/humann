#!/usr/bin/env python

"""
Reduce a table with a specific function

Dependencies: None

To Run: 
$ ./reduce_table.py -i table.tsv -o reduced_table.tsv

"""

import argparse
import sys
import os

from humann.tools import util

TABLE_DELIMITER="\t"
LEVEL_DELIMITER="|"

FUNCTION_OPTIONS={
    "sum":sum,
    "min":min,
    "max":max,
    "mean": lambda scores: sum(scores)/float(len(scores)) if scores else 0
}
DEFAULT_FUNCTION="max"

SORT_OPTIONS={
    # sort the data by the level with the lower levels first
    "level": lambda data: sorted(data, key=lambda k: k.count(LEVEL_DELIMITER)),
    # sort the data by the names
    "name": lambda data: sorted(data),
    # sort data so those with the largest values are printed first
    "value": lambda data: reversed(sorted(data, key=data.get))
}   
        
READ_STATUS_MESSAGE_BLOCK=1000
        
def reduce_table(function,input,output,verbose,sort_by):
    """
    Reduce the table by function
    """
    
    if verbose:
        print("Reading file: " + input)
        print("Using function: " + function)
        
    lines=util.process_gene_table_with_header(input, allow_for_missing_header=True)
    header=next(lines)
    
    if verbose:
        print("Opening output file: " + output)
    
    try:
        file_handle_out=open(output,"w")
    except EnvironmentError:
        sys.exit("Unable to write file: " + output)
        
    # create and write the header
    header_start="# header "
    if header:
        header_start=header.split(TABLE_DELIMITER)[0]
    
    file_handle_out.write(header_start+TABLE_DELIMITER+function+"\n")
    
    store_data={}
    read_line=1
    for line in lines:
        # read in the data and apply the function to each row
        data=line.split(TABLE_DELIMITER)
        item=data.pop(0)
        
        # try to convert the data to floats
        float_data=[]
        for point in data:
            try:
                float_point=float(point)
            except ValueError:
                float_point=0
                
            float_data.append(float_point)
            
        if float_data:
            reduced_data=FUNCTION_OPTIONS[function](float_data)
            
            if sort_by:
                store_data[item]=reduced_data
            else:
                file_handle_out.write(item+TABLE_DELIMITER+str(reduced_data)+"\n")
            
        read_line+=1
        
        if read_line % READ_STATUS_MESSAGE_BLOCK == 0 and verbose:
            print("Finished processing line: " + str(read_line))
    # write the sorted data
    if sort_by:
        for item in SORT_OPTIONS[sort_by](store_data):
            file_handle_out.write(item+TABLE_DELIMITER+str(store_data[item])+"\n")
    
    file_handle_out.close()

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Reduce table\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-v","--verbose", 
        help="additional output is printed\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "-i","--input",
        help="the input table\n",
        required=True)
    parser.add_argument(
        "-o","--output",
        help="the output table\n",
        required=True)
    parser.add_argument(
        "--function",
        choices=list(FUNCTION_OPTIONS),
        default=DEFAULT_FUNCTION,
        help="the function to apply")
    parser.add_argument(
        "--sort-by",
        choices=list(SORT_OPTIONS),
        help="sort the output by the selection")

    return parser.parse_args()


def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    # check for format of the gene tables
    input=os.path.abspath(args.input)
    
    # check the directory exists
    if not os.path.isfile(input):
        sys.exit("The input file provided can not be found." + 
            "  Please enter a new file.")
        
    reduce_table(args.function,args.input,args.output,args.verbose,args.sort_by)
    
    print("Reduced table created: " + args.output)

if __name__ == "__main__":
    main()
