#!/usr/bin/env python

"""
Create a genus level gene families table

Dependencies: None

To Run: 
$ python genefamilies_genus_level.py -i genefamiles.tsv -o genefamilies_genus.tsv

"""

import argparse
import sys
import os

from humann.tools import util

TABLE_DELIMITER="\t"
TAXONOMY_DELIMITER="|"
GENUS_DELIMITER="."
        
def create_table(input,output):
    """
    Create a gene families table with genus level taxonomy
    """
    
    table = util.Table(input)
    values = {}
    for rowhead, row in zip(table.rowheads, table.data):
        if not util.is_special(rowhead):
            feature, name, stratum = util.fsplit(rowhead)
            if stratum is not None:
                stratum = stratum.split(GENUS_DELIMITER)[0]
                rowhead = util.fjoin(feature, name, stratum)
        if rowhead in values:
            values[rowhead] = [float(a) + float(b) for a, b in zip(values[rowhead], row)]
        else:
            values[rowhead] = row[:]
    table.rowheads = util.fsort(values)
    table.data = [[str(value) for value in values[h]] for h in table.rowheads]
    table.write(output)

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Create a genus level gene families file\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i","--input",
        help="the gene families input table\n",
        required=True)
    parser.add_argument(
        "-o","--output",
        help="the output table\n",
        required=True)

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
        
    create_table(args.input,args.output)
    
    print("Genus gene families table created: " + args.output)

if __name__ == "__main__":
    main()
