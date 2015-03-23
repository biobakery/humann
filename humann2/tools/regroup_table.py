#! /usr/bin/env python

"""
HUMAnN2 utility for regrouping TSV files
Run ./regroup_table.py -h for usage help
"""

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import argparse
import sys
import util

def get_args ():
    """ Get args from Argparse """
    parser = argparse.ArgumentParser()
    parser.add_argument( 
        "-i", "--input", 
        default=None,
        help="Original output table (.tsv format); default=[STDIN]",
        )
    parser.add_argument( 
        "-g", "--groups", 
        help="Mapping of features to superfeatures (.tsv or .tsv.gz format)",
        )
    parser.add_argument( 
        "-r", "--reverse",
        action="store_true",
        help="Mapping file is reversed: mapping from superfeatures to features",
    )
    parser.add_argument( 
        "-f", "--function", 
        choices=["sum", "mean"],
        default="sum",
        help="How to combine grouped features; default=[sum]",
        )
    parser.add_argument( 
        "-o", "--output", 
        default=None,
        help="Path for modified output table; default=[STDOUT]",
        )
    args = parser.parse_args()
    return args

def mean( vector ):
    return sum( vector ) / float( len( vector ) )

def regroup( table, groups, function ):
    feature_counts = {}
    function = {"sum":sum, "mean":mean}[function]
    mapping = {}
    for i, rowhead in enumerate( table.rowheads ):
        items = rowhead.split( util.c_strat_delim )
        if items[0] not in feature_counts:
            feature_counts[items[0]] = 0
        # account for previously renamed features
        feature = items[0].split( util.c_name_delim )[0]
        if feature in groups:
            for superfeature in groups[feature]:
                if len( items ) == 1:
                    feature_counts[items[0]] += 1
                # account for stratified feature
                superrowhead = superfeature if len( items ) == 1 \
                                  else util.c_strat_delim.join( [superfeature, items[1]] ) 
                mapping.setdefault( superrowhead, [] ).append( i )
    # rebuild table
    superrowheads = sorted( mapping.keys() )
    superdata = []
    for superrowhead in superrowheads:
        oldrow_index = mapping[superrowhead]
        newrow = [[] for j in range( len( table.colheads ) )]
        for i in oldrow_index:
            for j in range( len( table.colheads ) ):
                newrow[j].append( float( table.data[i][j] ) )
        # collapse groups
        superdata.append( [function( block ) for block in newrow] )
    table.rowheads = superrowheads
    table.data = superdata
    # report
    n = len( feature_counts )
    ungrouped = feature_counts.values().count( 0 )
    grouped_total = n - ungrouped
    grouped_multi = grouped_total - feature_counts.values().count( 1 )
    print( "Features=%d; Grouped 1+ times=%d (%%%.1f); Grouped 2+ times= %d (%%%.1f)" % \
           ( n,
             grouped_total,
             100 * grouped_total / float( n ),
             grouped_multi,
             100 * grouped_multi / float( n ),
         ), file=sys.stderr )

def main ( ):
    args = get_args()
    table = util.Table( args.input )
    groups = util.load_polymap( args.groups )
    if args.reverse:
        groups2 = {}
        for superfeature in groups:
            for feature in groups[superfeature]:
                groups2.setdefault( feature, {} )[superfeature] = 1
        groups = groups2
    regroup( table, groups, args.function )
    fh = open( args.output, "w" ) if args.output is not None else sys.stdout
    table.write( fh )

if __name__ == "__main__":
    main()
