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

from humann2 import config
from humann2 import utilities
from humann2 import store

# name global logging instance
logger=logging.getLogger(__name__)

def blastx_coverage( blast6out, min_coverage, alignments=None, log_messages=None, apply_filter=None):
    # create alignments instance if none is passed
    if alignments is None:
        alignments=store.Alignments()
    
    # store protein lengths
    protein_lengths = {}
    # store unique positions hit in each protein as sets
    protein_hits = defaultdict( set )
    # track proteins with sufficient coverage
    allowed = set()
    # track alignments unable to compute coverage
    no_coverage=0
    # parse blast6out file, applying filtering as selected
    for alignment_info in utilities.get_filtered_translated_alignments(blast6out, alignments, apply_filter=apply_filter):
        ( protein_name, gene_length, queryid, matches, bug, alignment_length,
          subject_start_index, subject_stop_index) = alignment_info
          
        # divide the gene length by 3 to get protein length from nucleotide length
        gene_length = gene_length / 3
                    
        # store the protein length
        protein_lengths[protein_name] = gene_length
        
        # add the range of the alignment to the protein hits
        protein_range=range(subject_start_index, subject_stop_index)
        if protein_range:
            # keep track of unique hit positions in this protein
            protein_hits[protein_name].update(protein_range)
        else:
            no_coverage+=1
    # track proteins without lengths
    no_length=0
    # compute coverage
    for protein_name, hit_positions in protein_hits.items():
        try:
            # compute coverage, with 50 indicating that 50% of the protein is covered
            coverage = len( hit_positions ) / float( protein_lengths[protein_name] ) * 100
        except ZeroDivisionError:
            coverage = 0
            no_length+=1
        
        if coverage >= min_coverage:
            allowed.add(protein_name)

    output_messages=["Total alignments without coverage information: "+str(no_coverage)]
    output_messages+=["Total proteins in blastx output: "+str(len( protein_lengths ))]
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
        help="the subject coverage threshold\n[ DEFAULT : "+str(config.translated_subject_coverage_threshold)+" ]", 
        default=config.translated_subject_coverage_threshold)
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
