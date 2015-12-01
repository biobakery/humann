#! /usr/bin/env python

"""
This is a HUMAnN2 utility function
* Do a first pass on blastx output
* Identify proteins that were well-covered by reads
* Return them as a dict
* When processing blastx output in HUMAnN2, consider only these proteins
===============================================
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

import sys
import re
import logging
import argparse
from collections import defaultdict

from .. import config

# name global logging instance
logger=logging.getLogger(__name__)

def blastx_coverage( blast6out, min_coverage, alignments=None, log_messages=None ):
    # store protein lengths
    prot_lens = {}
    # store unique positions hit in each protein as sets
    prot_hits = defaultdict( set )
    # track proteins with sufficient coverage
    allowed = set()
    # track alignments unable to compute coverage
    no_coverage=0
    # parse blast6out file
    with open( blast6out ) as fh:
        for line in fh:
            # read if not comment line
            if not re.match("#",line):
                data=line.rstrip().split(config.blast_delimiter)
                reference=data[config.blast_reference_index]
                if alignments:
                    # if the alignment structure is provided use the function with id mapping and
                    # custom annotations to compute annotation information
                    prot_name, gene_len, bug = alignments.process_reference_annotation(reference)
                else:
                    # if alignments are not provided, default to format of "protein | length"
                    reference_information = reference.split("|")
                    prot_name = reference_information[0]
                    try:
                        gene_len = reference_information[1]
                    except IndexError:
                        gene_len = 0
                    
                # the gene length is in nucleotides, so divide by three for protein length
                try:
                    prot_lens[prot_name] = int( gene_len ) / 3
                except ValueError:
                    prot_lens[prot_name] = 0
        
                try:
                    prot_start = int( data[config.blast_protein_start_index] )
                    prot_stop = int( data[config.blast_protein_end_index] )
                    # keep track of unique hit positions in this protein
                    prot_hits[prot_name].update( range( prot_start-1, prot_stop ) )
                except (ValueError,IndexError):
                    no_coverage+=1
    # track proteins without lengths
    no_length=0
    # compute coverage
    for prot_name, hit_positions in prot_hits.items():
        try:
            # compute coverage, with 50 indicating that 50% of the protein is covered
            coverage = len( hit_positions ) / float( prot_lens[prot_name] ) * 100
        except ZeroDivisionError:
            coverage = 0
            no_length+=1
        
        if coverage >= min_coverage:
            allowed.add(prot_name)

    output_messages=["Total alignments without coverage information: "+str(no_coverage)]
    output_messages+=["Total proteins in blastx output: "+str(len( prot_lens ))]
    output_messages+=["Total proteins without lengths: "+str(no_length)]
    output_messages+=["Proteins with coverage greater than threshold ("+str(min_coverage)+"): "+str(len( allowed ))]
    
    # write out informational messages to log or stdout, depending on input parameters
    if log_messages:
        for message in output_messages:
            logger.info(message)
    else:
        print("\n".join(output_messages))
        
    return allowed

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "Compute blastx coverage\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i","--input",
        help="the blastx formatted input file\n",
        required=True)
    parser.add_argument(
        "--coverage-threshold",
        type=float,
        help="the coverage threshold\n[ DEFAULT : "+str(config.coverage_threshold)+" ]", 
        default=config.coverage_threshold)
    parser.add_argument(
        "--print-protein-list",
        action="store_true",
        help="print the list of proteins that meet the coverage threshold")
    
    return parser.parse_args()
    
def main():
    # parse the arguments from the user
    args = parse_arguments(sys.argv)
    
    # run coverage computation
    allowed = blastx_coverage(args.input, args.coverage_threshold)

    if args.print_protein_list:
        print("\n".join(allowed))

if __name__ == "__main__":
    main()
