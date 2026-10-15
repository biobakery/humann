#! /usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
from collections import namedtuple
import argparse
import sys
import os

try:
    from humann import config
    from humann.tools import util
except ImportError:
    sys.exit("CRITICAL ERROR: Unable to find the HUMAnN python package." +
        " Please check your install.") 

description = """
HUMAnN utility for regrouping table features
=============================================
Given a table of feature values and a mapping
of groups to component features, produce a 
new table with group values in place of 
feature values.
"""

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

p_root = os.path.join( os.path.dirname( os.path.abspath(__file__) ), os.pardir )
Groups = namedtuple( "Groups", ["path", "start", "skip"] )
c_default_groups = {
    "uniref90_rxn": Groups( 
        os.path.join( config.utility_mapping_database, "metacyc_reactions_level4ec_only.uniref.bz2" ), 0, [1] ),
    }

# get a list of all available script mapping files
try:
    all_mapping_files=os.listdir(config.utility_mapping_database)
except EnvironmentError:
    all_mapping_files=[]

# add larger mapping files if they have been downloaded
all_larger_mapping_files=["map_go_uniref90.txt.gz",
                          "map_ko_uniref90.txt.gz",
                          "map_level4ec_uniclust90.txt.gz",
                          "map_pfam_uniref90.txt.gz",
                          "map_eggnog_uniref90.txt.gz",]
larger_mapping_files_found=False
for mapping_file in all_larger_mapping_files:
    if mapping_file in all_mapping_files:
        # get the option name from the file name
        option_to, option_from=mapping_file.split(".")[0].replace("map_","").split("_")
        option="_".join([option_from, option_to])
        c_default_groups[option]=Groups(os.path.join(config.utility_mapping_database, mapping_file), 0, [])
        larger_mapping_files_found=True
    
if not larger_mapping_files_found:
    description+="""
    
For group mapping files, run the following command:
$ humann_databases --download utility_mapping full $DIR
Replacing, $DIR with the directory to download and install the databases."""
    sys.exit("Mapping files not installed.\n"+description)

c_funcmap = {"sum":sum, "mean":lambda row: sum( row ) / float( len( row ) )}

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

def get_args( ):
    """ Get args from Argparse """
    parser = argparse.ArgumentParser(
        description=description, 
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument( 
        "-i", "--input", 
        default=None,
        help="Original output table (tsv or biom format); default=[TSV/STDIN]",
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
        help="How to combine biological groups; residual accounting abundance is always summed; default=sum",
        )
    parser.add_argument( 
        "-e", "--precision",
        default=None,
        type=int,
        help="Decimal places to round to after applying function; default=Don't round",
        )
    parser.add_argument( 
        "-u", "--ungrouped",
        default="Y",
        choices=["Y", "N"],
        help="Collect newly ungroupable features into " + config.ungrouped_reaction_name
             + "? Existing protected accounting rows are retained independently; default=Y",
        )
    parser.add_argument( 
        "-p", "--protected",
        default="Y",
        choices=["Y", "N"],
        help="Preserve accounting rows and bypass custom mappings for them; default=Y",
        )
    parser.add_argument( 
        "-o", "--output", 
        default=None,
        help="Path for modified output table; default=STDOUT",
        )
    args = parser.parse_args()
    return args

# ---------------------------------------------------------------
# regrouping function
# ---------------------------------------------------------------

def regroup( table, map_feature_groups, function, precision, ungrouped=False, protected=True ):
    
    function = c_funcmap[function]
    seen_before = {}
    feature_counts = {}
    # index of new group names to old table rows
    mapping = {}

    for i, rowhead in enumerate( table.rowheads ):
        special = util.is_special(rowhead)
        # Protected input categories must not leak into custom biological groups.
        if protected and special:
            mapping.setdefault(rowhead, []).append(i)
            continue
        feature, name, stratum = util.fsplit( rowhead )
        if not special and feature not in feature_counts:
            feature_counts[feature] = 0
        # decide which groups to use
        if feature in map_feature_groups:
            groups = map_feature_groups[feature]
        elif ungrouped:
            groups = [util.c_ungrouped]
        else:
            groups = []
        # track grouping
        for group in groups:
            if protected and util.is_special(group) and not util.is_ungrouped(group):
                raise ValueError("Custom groups cannot assign biological features to protected category " + group)
            if not special and feature not in seen_before and not util.is_special(group):
                feature_counts[feature] += 1
            # account for stratified feature
            groupname = group 
            if stratum is not None:
                groupname = util.fjoin( groupname, stratum=stratum )
            mapping.setdefault( groupname, [] ).append( i )
        # we have processed an instance of this feature
        seen_before[feature] = 1

    # rebuild table
    groupnames = util.fsort( mapping.keys( ) )
    groupdata = []
    for groupname in groupnames:
        oldrow_index = mapping[groupname]
        newrow = [[] for j in range( len( table.colheads ) )]
        for i in oldrow_index:
            for j in range( len( table.colheads ) ):
                try:
                    newrow[j].append( float( table.data[i][j] ) )
                except ( IndexError, ValueError ):
                    print("WARNING: Unexpected truncated input file")
        # collapse groups
        aggregate = sum if util.is_special(groupname) else function
        newrow = [aggregate( block ) for block in newrow]
        if protected and util.is_special(groupname) and len(oldrow_index) == 1 and util.is_special(table.rowheads[oldrow_index[0]]):
            newrow=table.data[oldrow_index[0]][:]
        if precision is not None and not (protected and util.is_special(groupname)):
            newrow = [round( k, precision ) for k in newrow]
        groupdata.append( newrow )
    table.rowheads = groupnames
    table.data = groupdata

    # report
    n = len( feature_counts )
    ungrouped = list( feature_counts.values( ) ).count( 0 )
    grouped_total = n - ungrouped
    grouped_multi = grouped_total - list(feature_counts.values()).count( 1 )
    print( "Original Feature Count: %d; Grouped 1+ times: %d (%.1f%%); Grouped 2+ times: %d (%.1f%%)" % \
           ( n,
             grouped_total,
             100 * grouped_total / float( n or 1 ),
             grouped_multi,
             100 * grouped_multi / float( n or 1 ),
         ), file=sys.stderr )

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args()
    # load the table; find unique feature ids (no names)
    table = util.Table( args.input )
    features = {k for k in table.rowheads}
    features = {k.split( util.c_strat_delim )[0] for k in features}
    features = {k.split( util.c_name_delim )[0] for k in features}
    # decide what grouping file to load and how
    if args.custom is not None:
        print( "Loading custom groups file: {}".format( args.custom ), file=sys.stderr )
        p_groups, start, skip = args.custom, 0, []
    elif args.groups is not None:
        p_groups, start, skip = c_default_groups[args.groups]
    else:
        sys.exit( "Must specify either 1) built-in groups option [--groups] or 2) custom groups file [--custom]" )
    # load the grouping file
    map_group_features = util.load_polymap( 
        p_groups, start=start, skip=skip, allowed_values=features )
    # coerce to features-first format (unless explicitly reversed)
    if not args.reversed:
        map_feature_groups = {}
        for group, fdict in map_group_features.items():
            for feature in fdict:
                map_feature_groups.setdefault( feature, {} )[group] = 1
    else:
        map_feature_groups = map_group_features
    # perform the table regrouping
    regroup( table, map_feature_groups, args.function, args.precision,
             ungrouped=args.ungrouped=="Y", protected=args.protected=="Y" )
    table.write( args.output )

if __name__ == "__main__":
    main()
