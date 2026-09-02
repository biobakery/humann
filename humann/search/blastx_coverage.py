#! /usr/bin/env python

"""
This is a HUMAnN utility function
* Do a first pass on blastx output
* Identify reference sequences that were well-covered by reads
* Return their annotations as a set
* When processing blastx output in HUMAnN, consider only these reference sequences
===============================================
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

import sys
import re
import logging
import argparse
from collections import defaultdict

from humann import config
from humann import utilities
from humann import store

# name global logging instance
logger=logging.getLogger(__name__)

def blastx_coverage( blast6out, min_coverage, alignments=None, log_messages=None, apply_filter=None, nucleotide = False, query_coverage_threshold=config.translated_query_coverage_threshold, identity_threshold = config.nucleotide_identity_threshold):
    # create alignments instance if none is passed
    if alignments is None:
        alignments=store.Alignments()
    
    # Coverage is tracked per reference sequence rather than per gene family.
    # A gene family (ie a UniRef90 id) appears once per species pangenome in the
    # nucleotide database, and each of those copies has its own sequence length. It can
    # also appear more than once within a pangenome, for the unknown gene family and when
    # running in UniRef50 mode. Keying on the gene family alone would normalize by an
    # arbitrary one of those lengths and would also merge the hit positions of unrelated
    # sequences. The full reference annotation is unique per database sequence (it
    # includes the species and the gene family id) so it is used as the key here.
    # store reference sequence lengths
    reference_lengths = {}
    # store unique positions hit in each reference sequence as sets
    reference_hits = defaultdict( str )
    # track reference sequences with sufficient coverage
    allowed = set()
    # track alignments unable to compute coverage
    no_coverage=0
    # parse blast6out file, applying filtering as selected
    for alignment_info in utilities.get_filtered_translated_alignments(blast6out, alignments, apply_filter=apply_filter, log_filter = log_messages, query_coverage_threshold = query_coverage_threshold, identity_threshold = identity_threshold):
        ( protein_name, gene_length, queryid, matches, bug, alignment_length,
          subject_start_index, subject_stop_index, reference_name) = alignment_info

        # divide the gene length by 3 to get protein length from nucleotide length
        if not nucleotide:
            gene_length = gene_length / 3

        # store the length of this specific reference sequence
        reference_lengths[reference_name] = gene_length

        # add the range of the alignment to the reference sequence hits
        reference_range=range(subject_start_index, subject_stop_index)
        if reference_range:
            # keep track of unique hit positions in this reference sequence
            reference_hits[reference_name]+="{0}-{1};".format(subject_start_index, subject_stop_index)
        else:
            no_coverage+=1
    # track reference sequences without lengths
    no_length=0
    # compute coverage
    for reference_name, hit_positions in reference_hits.items():

        # compile the hit positions
        range_hit_positions = set()
        for alignment_hit in hit_positions.split(";")[:-1]:
             start_index, stop_index = alignment_hit.split("-")
             new_range = range(int(start_index), int(stop_index))
             range_hit_positions.update(new_range)

        try:
            # compute coverage, with 50 indicating that 50% of the sequence is covered
            coverage = len( range_hit_positions ) / float( reference_lengths[reference_name] ) * 100
        except ZeroDivisionError:
            coverage = 0
            no_length+=1

        if coverage >= min_coverage:
            allowed.add(reference_name)

    output_messages=["Total alignments without coverage information: "+str(no_coverage)]
    output_messages+=["Total reference sequences in blastx output: "+str(len( reference_lengths ))]
    output_messages+=["Total reference sequences without lengths: "+str(no_length)]
    output_messages+=["Reference sequences with coverage greater than threshold ("+str(min_coverage)+"): "+str(len( allowed ))]

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
        help="print the list of reference sequences that meet the coverage threshold")
    
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
