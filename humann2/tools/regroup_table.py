#! /usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
from collections import namedtuple
import argparse
import sys
import os
import util

description = """
HUMAnN2 utility for regrouping table features
=============================================
Given a table of feature values and a mapping
of groups to component features, produce a 
new table with group values in place of 
feature values.

For additional groups files, see https://bitbucket.org/biobakery/humann2/src/tip/humann2/data/misc/
"""

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

p_root = os.path.join( os.path.dirname( os.path.abspath(__file__) ), os.pardir )
Groups = namedtuple( "Groups", ["path", "start", "skip"] )
c_default_groups = {
    "uniref50_rxn": Groups( 
        os.path.join( p_root, "data", "pathways", "metacyc_reactions_level4ec_only.uniref.bz2" ), 0, [1] ),
    "uniref50_ec":  Groups( 
        os.path.join( p_root, "data", "misc", "map_ec_uniref50.txt.gz" ), 0, [] ),
    }
c_protected = [util.c_unmapped, util.c_unintegrated]
c_funcmap = {"sum":sum, "mean":lambda row: sum( row ) / float( len( row ) )}

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
        choices=c_default_groups.keys(),
        default=None,
        help="Built-in grouping options",
        )
    parser.add_argument( 
        "-c", "--custom", 
        default=None,
        help="Custom groups file (.tsv or .tsv.gz format)",
        )
    parser.add_argument( 
        "-r", "--reversed",
        action="store_true",
        help="Custom groups file is reversed: mapping from features to groups",
    )
    parser.add_argument( 
        "-f", "--function", 
        choices=c_funcmap.keys(),
        default="sum",
        help="How to combine grouped features; default=[sum]",
        )
    parser.add_argument( 
        "-u", "--ungrouped",
        default="Y",
        choices=["Y", "N"],
        help="Include an 'UNGROUPED' group to capture features that did not belong to other groups? default=Y",
        )
    parser.add_argument( 
        "-p", "--protected",
        default="Y",
        choices=["Y", "N"],
        help="Carry through protected features, such as 'UNMAPPED'? default=Y",
        )
    parser.add_argument( 
        "-o", "--output", 
        default=None,
        help="Path for modified output table; default=[STDOUT]",
        )
    args = parser.parse_args()
    return args

def regroup( table, map_feature_groups, function, ungrouped=False ):
    feature_counts = {}
    function = c_funcmap[function]
    # index of new group names to old table rows
    mapping = {}
    for i, rowhead in enumerate( table.rowheads ):
        items = rowhead.split( util.c_strat_delim )
        if items[0] not in feature_counts:
            feature_counts[items[0]] = 0
        # account for previously renamed features
        feature = items[0].split( util.c_name_delim )[0]
        if feature in map_feature_groups:
            groups = map_feature_groups[feature]
        elif ungrouped:
            groups = [util.c_ungrouped]
        else:
            groups = []
        for group in groups:
            if len( items ) == 1 and group != util.c_ungrouped:
                feature_counts[items[0]] += 1
            # account for stratified feature
            groupname = group if len( items ) == 1 \
                        else util.c_strat_delim.join( [group, items[1]] ) 
            mapping.setdefault( groupname, [] ).append( i )
    # rebuild table
    # *** sorting changed to force 1|A to come before 11 ***
    groupnames = sorted( mapping.keys(), key=lambda k: k.split( util.c_strat_delim ) )
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
    print( "Original Feature Count: %d; Grouped 1+ times: %d (%%%.1f); Grouped 2+ times: %d (%%%.1f)" % \
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
    # load the table; find all keys in the table
    table = util.Table( args.input )
    allowed_values = {k.split( util.c_name_delim )[0]:1 for k in table.rowheads}
    # decide what grouping file to load and how
    if args.custom is not None:
        print( "Loading custom groups file: {}".format( args.custom ), file=sys.stderr )
        p_groups, start, skip = args.custom, 0, []
    elif args.groups is not None:
        p_groups, start, skip = c_default_groups[args.groups]
    else:
        sys.exit( "Must (i) choose groups option or (ii) provide groups file" )
    # load the grouping file
    map_group_features = util.load_polymap( 
        p_groups, start=start, skip=skip, allowed_values=allowed_values )
    # coerce to features-first format (unless explicitly reversed)
    if not args.reversed:
        map_feature_groups = {}
        for group, fdict in map_group_features.items():
            for feature in fdict:
                map_feature_groups.setdefault( feature, {} )[group] = 1
    else:
        map_feature_groups = map_group_features
    # add protected cases to mapping?
    if args.protected == "Y":
        for feature in c_protected:
            map_feature_groups.setdefault( feature, {} )[feature] = 1
    # perform the table regrouping
    regroup( table, map_feature_groups, args.function, ungrouped=args.ungrouped=="Y" )
    table.write( args.output )

if __name__ == "__main__":
    main()
