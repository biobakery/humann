#!/usr/bin/env python

from __future__ import print_function
import sys
import csv
import argparse

def get_args( ):

    parser = argparse.ArgumentParser( )
    parser.add_argument( "dmndout" )
    parser.add_argument( "--protein-coverage", default=0.5, type=float )
    parser.add_argument( "--outfile", default="well_covered.txt" )
    args = parser.parse_args( )
    return args

def interval_overlap_hits( intervals ):

    """
    Given a list of intervals, compute total number of covered sites:
    [[1,5],[3,7],[4,5],[9,11]] => 1:7 9:11 => 7 + 3 => 10
   
    Three types of sorted interval overlaps:

    1. i2 <= leading_edge (no new info):
       ----------L
           1---2

    2. leading_edge < i1:
       ----------L  *****
                    1---2

    3. i1 <= leading_edge < i2:
       ----------L***
                1---2
    """

    hits = 0
    leading_edge = 0
    for i1, i2 in sorted( intervals ):
        if i2 > leading_edge:
            if i1 > leading_edge:
                hits += i2 - i1 + 1
            else:
                hits += i2 - leading_edge
            leading_edge = i2
    return hits

def main( ):
    
    args = get_args( )
    # track the lens of hit proteins
    lens = {}
    # track the start,stop of hits within proteins
    intervals = {}
    
    with open( args.dmndout ) as fh:
        for row in csv.reader( fh, dialect="excel-tab" ):
            # protein name/length
            prot = row[1]
            name, L = prot.split( "|" )
            if name not in lens:
                lens[name] = int( L ) / 3
            # coordinates of the hit within the protein
            i1, i2 = sorted( map( int, row[8:10] ) )
            # uniref len reported as gene equivalent
            intervals.setdefault( name, [] ).append( [i1, i2] )
    
    fh = open( args.outfile, "w" )
    for name in sorted( lens ):
        cover = interval_overlap_hits( intervals[name] ) / float( lens[name] )
        if cover >= args.protein_coverage:
            print( "\t".join( [name, str( cover )] ), file=fh )
    fh.close( )

if __name__ == "__main__":
    main( )
