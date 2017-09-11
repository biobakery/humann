#!/usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import os
import sys
import argparse
from collections import namedtuple

try:
    from humann2 import config
    from humann2.tools import util
    #from humann2.tools.better_table import Table
except ImportError:
    sys.exit( "CRITICAL ERROR: Unable to find the HUMAnN2 python package." +
        " Please check your install." )  

import numpy as np
from better_table import Table

description = """
HUMAnN2 utility for regrouping table features
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
        os.path.join( p_root, "data", "pathways", "metacyc_reactions_level4ec_only.uniref.bz2" ), 0, [1] ),
    "uniref50_rxn": Groups( 
        os.path.join( p_root, "data", "pathways", "metacyc_reactions_level4ec_only.uniref.bz2" ), 0, [1] ),
    }

# get a list of all available script mapping files
try:
    all_mapping_files = os.listdir( config.utility_mapping_database )
except EnvironmentError:
    all_mapping_files=[]

"""
Can we autopopulate this by searching the above location?
"""

# add larger mapping files if they have been downloaded
all_larger_mapping_files=["map_go_uniref50.txt.gz","map_go_uniref90.txt.gz",
                          "map_infogo1000_uniref50.txt.gz","map_infogo1000_uniref90.txt.gz",
                          "map_ko_uniref50.txt.gz","map_ko_uniref90.txt.gz",
                          "map_level4ec_uniref50.txt.gz","map_level4ec_uniref90.txt.gz",
                          "map_pfam_uniref50.txt.gz","map_pfam_uniref90.txt.gz",
                          "map_eggnog_uniref50.txt.gz","map_eggnog_uniref90.txt.gz",]

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
    
For additional group mapping files, run the following command:
$ humann2_databases --download utility_mapping full $DIR
Replacing, $DIR with the directory to download and install the databases."""

c_protected = [util.c_unmapped, util.c_unintegrated]

c_functions = {
    "sum":    np.sum,
    "mean":   np.mean,
    "median": np.median,
}

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
        metavar="<path>",
        help="Original output table (tsv or biom format)\nDefault=[STDIN]",
        )  
    parser.add_argument( 
        "-g", "--groups", 
        choices=c_default_groups.keys( ),
        metavar="<choice>",
        default=None,
        help=util.pretty_grid( c_default_groups, desc="Select an available regrouping option:" ),
        )
    parser.add_argument( 
        "-c", "--custom", 
        default=None,
        metavar="<path>",
        help="Custom groups file (.tsv or .tsv.gz format)",
        )
    parser.add_argument( 
        "-r", "--reversed",
        action="store_true",
        help="Custom groups file is reversed: mapping from features to groups",
    )
    parser.add_argument( 
        "-f", "--function", 
        choices=c_functions.keys( ),
        default="sum",
        metavar="<choice>",
        help=util.pretty_grid( c_functions, desc="Select a regrouping function (default=sum):" ),
        )
    """
    parser.add_argument( 
        "-e", "--precision",
        default=None,
        type=int,
        metavar="<int>",
        help="Decimal places to round to after applying function\nDefault=[no rounding]",
        )
        """
    parser.add_argument( 
        "--exclude-ungrouped",
        action="store_true",
        help="Do not include the special 'UNGROUPED' group\nDefault=[include this group]",
        )
    parser.add_argument( 
        "--exclude-special",
        action="store_true",
        help="Do not automatically carry through special features such as 'UNMAPPED'\nDefault=[keep these features]",
        )
    parser.add_argument( 
        "-o", "--output", 
        default=None,
        metavar="<path>",
        help="Path for modified output table\nDefault=[STDOUT]",
        )
    args = parser.parse_args()
    return args

# ---------------------------------------------------------------
# regrouping function
# ---------------------------------------------------------------

def regroup( table, map_feature_groups, function, exclude_ungrouped=False ):
    
    function = c_functions[function]
    # index of which features mapped where
    community_features = {}
    # index of new group names to original features
    mapping = {}
    for f in util.fsort( table.data ):
        fbase, fname, stratum = util.fsplit( f )
        inner = community_features.setdefault( fbase, set( ) )
        # decide which groups to use
        if fbase in map_feature_groups:
            groups = map_feature_groups[fbase]
            inner.update( groups )
        elif not exclude_ungrouped:
            groups = [util.c_ungrouped]
        else:
            groups = []
        # track grouping
        for group in groups:
            groupname = util.fjoin( group, stratum=stratum )
            mapping.setdefault( groupname, set( ) ).add( f )

    # rebuild table
    new_data = {}
    for groupname, features in mapping.items( ):
        # 1. vertically stack rows from the original table to be regrouped
        # 2. collapse the stack to a single row by applying function to columns of the stack
        stack = np.vstack( [table.data[f] for f in features] )
        new_data[groupname] = function( stack, axis=0 )

    # summarize
    fdict = {}
    fdict["GT"] = len( {util.fsplit( group )[0] for group in mapping} )
    fdict["NT"] = len( community_features )
    fdict["N1"] = len( [fbase for fbase, groups in community_features.items( ) \
                            if len( groups ) >= 1] )
    fdict["N2"] = len( [fbase for fbase, groups in community_features.items( ) \
                            if len( groups ) >= 2] )
    fdict["P1"] = 100 * fdict["N1"] / float( fdict["NT"] )
    fdict["P2"] = 100 * fdict["N2"] / float( fdict["NT"] )
    string = ("Regrouping report:\n"
              "  # of original community features: {NT}\n"  
              "  # of features grouped one+ times: {N1} ({P1:.1f}%)\n"
              "  # of features grouped two+ times: {N2} ({P2:.1f}%)\n"
              "  # of final community groups:      {GT}")
    print( string.format( **fdict ), file=sys.stderr )
    
    # new table
    return Table( new_data, headers=table.headers, metadata=table.metadata )

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args( )
    # load the table; find unique feature ids (no names)
    table = Table( args.input )
    features = {util.fsplit( f )[0] for f in util.fsort( table.data )}
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
    # add protected cases to mapping?
    if not args.exclude_special:
        for feature in c_protected:
            map_feature_groups.setdefault( feature, {} )[feature] = 1
    # perform the table regrouping
    new_table = regroup( table, map_feature_groups, args.function, args.exclude_ungrouped )
    new_table.write( args.output, unfloat=True )

if __name__ == "__main__":
    main( )
