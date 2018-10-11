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

TABLE_DELIMITER="\t"
TAXONOMY_DELIMITER="|"
GENUS_DELIMITER="."
        
def create_table(input,output):
    """
    Create a gene families table with genus level taxonomy
    """
    
    # read the input gene table
    genus_values={}
    all_values=""
    with open(input, "rt") as file_handle:
        with open(output, "w") as file_handle_write:
            # write the header to the new file
            header=file_handle.readline()
            file_handle_write.write(header)
            
            for line in file_handle:
                if TAXONOMY_DELIMITER in line:
                    # remove the species from the gene/taxonomy
                    data=line.rstrip().split(TABLE_DELIMITER)
                    gene_taxonomy=data[0].split(GENUS_DELIMITER)[0]
                    if gene_taxonomy in genus_values:
                        genus_values[gene_taxonomy]=[float(x) + y for x,y in zip(data[1:], genus_values[gene_taxonomy])]
                    else:
                        genus_values[gene_taxonomy]=[float(x) for x in data[1:]]
                else:
                    # check for empty lines
                    if TABLE_DELIMITER in line:
                        if all_values:
                            file_handle_write.write(all_values)
                            # write genus values for gene, if present
                            for gene_taxonomy, values in genus_values.iteritems():
                                file_handle_write.write(TABLE_DELIMITER.join([gene_taxonomy]+[str(x) for x in values])+"\n")
                            genus_values={}
                        all_values=line

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
