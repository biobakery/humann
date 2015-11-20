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
import csv
from collections import defaultdict

def blastx_coverage( blast6out, min_coverage ):
    # store protein lengths
    prot_lens = {}
    # store unique positions hit in each protein as sets
    prot_hits = defaultdict( set )
    # track proteins with sufficient coverage
    allowed = {}
    # parse blast6out file
    with open( blast6out ) as fh:
        for row in csv.reader( fh, dialect="excel-tab" ):
            prot_name, gene_len = row[1].split( "|" )
            prot_lens[prot_name] = int( gene_len ) / 3
            prot_start = int( row[8] )
            prot_stop = int( row[9] )
            # keep track of unique hit positions in this protein
            prot_hits[prot_name].update( range( prot_start-1, prot_stop ) )
    # compute coverage
    for prot_name, hit_positions in prot_hits.items():
        coverage = len( hit_positions ) / float( prot_lens[prot_name] )
        if coverage >= min_coverage:
            allowed[prot_name] = 1
    # done
    print >>sys.stderr, "Total proteins in blastx output:", len( prot_lens )
    print >>sys.stderr, "Proteins with good coverage    :", len( allowed )
    return allowed

if __name__ == "__main__":
    allowed = blastx_coverage( sys.argv[1], 0.5 )
