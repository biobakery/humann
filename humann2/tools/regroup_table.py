#! /usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import argparse
import sys
import util

description = """
HUMAnN2 utility for regrouping table features
=============================================
Given a table of feature values and a mapping
of groups to component features, produce a 
new table with group values in place of 
feature values.
"""

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

def get_args ():
    """ Get args from Argparse """
    parser = argparse.ArgumentParser(
        description=description, 
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument( 
        "-i", "--input", 
        default=None,
        help="Original output table (.tsv format); default=[STDIN]",
        )
    parser.add_argument( 
        "-g", "--groups", 
        help="Mapping of groups to features (.tsv or .tsv.gz format)",
        )
    parser.add_argument( 
        "-r", "--reversed",
        action="store_true",
        help="Mapping file is reversed: mapping from features to groups",
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

def regroup( table, map_feature_groups, function ):
    feature_counts = {}
    function = {"sum":sum, "mean":mean}[function]
    mapping = {}
    for i, rowhead in enumerate( table.rowheads ):
        items = rowhead.split( util.c_strat_delim )
        if items[0] not in feature_counts:
            feature_counts[items[0]] = 0
        # account for previously renamed features
        feature = items[0].split( util.c_name_delim )[0]
        if feature in map_feature_groups:
            for group in map_feature_groups[feature]:
                if len( items ) == 1:
                    feature_counts[items[0]] += 1
                # account for stratified feature
                groupname = group if len( items ) == 1 \
                                  else util.c_strat_delim.join( [group, items[1]] ) 
                mapping.setdefault( groupname, [] ).append( i )
    # rebuild table
    groupnames = sorted( mapping.keys() )
    groupdata = []
    for groupname in groupnames:
        oldrow_index = mapping[groupname]
        newrow = [[] for j in range( len( table.colheads ) )]
        for i in oldrow_index:
            for j in range( len( table.colheads ) ):
                newrow[j].append( float( table.data[i][j] ) )
        # collapse groups
        groupdata.append( [function( block ) for block in newrow] )
    table.rowheads = groupnames
    table.data = groupdata
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

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main ( ):
    args = get_args()
    table = util.Table( args.input )
    if not args.reversed:
        map_group_features = util.load_polymap( args.groups )
        map_feature_groups = {}
        for group, fdict in map_group_features.items():
            for feature in fdict:
                map_feature_groups.setdefault( feature, {} )[group] = 1
    else:
        map_feature_groups = util.load_polymap( args.groups )
    regroup( table, map_feature_groups, args.function )
    table.write( args.output )

if __name__ == "__main__":
    main()
