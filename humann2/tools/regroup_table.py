#!/usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import os
import sys
import argparse
import re

try:
    from humann2 import config
    from humann2.tools import util
    from humann2.tools.humann2_table import Table
except ImportError:
    sys.exit( "CRITICAL ERROR: Unable to find the HUMAnN2 python package.\n" +
              "Please check your install." )  

try:
    import numpy as np
except ImportError:
    sys.exit( "CRITICAL ERROR: This script requires the python scientific stack (e.g. numpy)" )

description = util.wrap( """
HUMAnN2 utility for regrouping table features

Given a table of feature values and a mapping of groups to component features, produce a 
new table with group values in place of feature values. This is most often used to sum
the abundance of gene families that are annotated to broader functional categories ("groups"),
such as GO terms or Pfam domains.
""" )

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_functions = {
    "sum":    np.sum,
    "mean":   np.mean,
    "median": np.median,
}

c_reactions_map = os.path.join(
    "data",
    "pathways",
    "metacyc_reactions_level4ec_only.uniref.bz2",
    )

# ---------------------------------------------------------------
# handling of regrouping files
# ---------------------------------------------------------------

mapping_options = {}

# metacyc mapping is included with the HUMAnN2 install (same file for uniref90/50)
humann2_install = os.path.join( os.path.dirname( os.path.abspath(__file__) ), os.pardir )
mapping_options["uniref90_rxn"] = os.path.join( humann2_install, c_reactions_map )
mapping_options["uniref50_rxn"] = os.path.join( humann2_install, c_reactions_map )

# additional mapping files
try:
    additional_files = os.listdir( config.utility_mapping_database )
except EnvironmentError:
    additional_files = []

FOUND_SOMETHING = False
for f in additional_files:
    match = re.search( "map_(.*?)_(uniref\d+).txt.gz", f )
    if match:
        FOUND_SOMETHING = True
        key = "{}_{}".format( match.group( 2 ), match.group( 1 ) )
        mapping_options[key] = os.path.join( config.utility_mapping_database, f )

# **** no longer needed ****
for k in ["uniref50_ec", "uniref50_transporter"]:
    if k in mapping_options:
        del mapping_options[k]

# for cli choser
mapping_options_order = sorted( mapping_options,
                                key=lambda x: x.split( "_" )[::-1] )

# tell the user how to get more mapping files
if not FOUND_SOMETHING:
    description += """
    
For additional group mapping files, run the following command:

$ humann2_databases --download utility_mapping full $DIR

Replacing, $DIR with the directory to download and install the databases."""

# ---------------------------------------------------------------
# command-line interface
# ---------------------------------------------------------------

def get_args( ):
    parser = argparse.ArgumentParser(
        description=description, 
        formatter_class=argparse.RawTextHelpFormatter,
    )
    util.attach_common_arguments( parser )
    parser.add_argument( 
        "-g", "--groups", 
        choices=mapping_options.keys( ),
        metavar="<choice>",
        default=None,
        help=util.pretty_grid( mapping_options_order, cols=2, 
                               desc="Select an available regrouping option:" ),
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
        metavar="<" + "/".join( c_functions.keys( ) ) + ">",
        help="Select a regrouping function\n[Default=sum]",
        )
    args = parser.parse_args()
    return args

# ---------------------------------------------------------------
# regrouping function
# ---------------------------------------------------------------

def regroup( table, map_feature_groups, function ):
    
    function = c_functions[function]
    # index of which features mapped where
    community_features = {}
    total_mass = table.zeros( )
    group_mass = table.zeros( )
    
    # process the regrouping of the table's features
    mapping = {}
    for f in util.fsort( table.data ):
        fbase, fname, stratum = util.fsplit( f )
        # decide how to group this row
        if fbase == util.c_unmapped:
            groups = [util.c_unmapped]
        else:
            groups = map_feature_groups.get( fbase, [util.c_ungrouped] )
        # update mapping
        for group in groups:
            groupname = util.fjoin( group, stratum=stratum )
            mapping.setdefault( groupname, set( ) ).add( f )
        # special tracking for community features (for summary)
        if stratum is None and fbase != util.c_unmapped:
            inner = community_features.setdefault( fbase, set( ) )
            total_mass += table.data[f]
            if groups != [util.c_ungrouped]:
                inner.update( groups )
                group_mass += table.data[f]

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
    fdict["MS"] = 100 * np.mean( group_mass / (util.c_eps + total_mass) )
    string = ("Regrouping report:\n"
              "  Original mass regrouped (mean over samples, excludes unclassified): {MS:.1f}%\n"
              "  # of original community features: {NT:,}\n"  
              "  # of features grouped 1+ times: {N1:,} ({P1:.1f}%)\n"
              "  # of features grouped 2+ times: {N2:,} ({P2:.1f}%)\n"
              "  # of final community groups: {GT:,}")
    print( string.format( **fdict ), file=sys.stderr )
    
    # new table
    return Table( new_data, headers=table.headers, metadata=table.metadata )

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args( )
    # decide what grouping file to load and how
    if args.custom is not None:
        print( "Loading custom groups file: {}".format( args.custom ), file=sys.stderr )
        mapping_path = args.custom
    elif args.groups is not None:
        mapping_path = mapping_options[args.groups]
    else:
        sys.exit( "Must specify either 1) built-in groups option [--groups] or 2) custom groups file [--custom]" )
    # load the table; find unique feature ids (no names)
    table = Table( args.input, last_metadata=args.last_metadata )
    features = {util.fsplit( f )[0] for f in util.fsort( table.data )}
    # load the grouping file
    map_group_features = util.load_polymap( 
        mapping_path, start=0, skip=[], allowed_values=features )
    # coerce to features-first format (unless explicitly reversed)
    if not args.reversed:
        map_feature_groups = {}
        for group, fdict in map_group_features.items():
            for feature in fdict:
                map_feature_groups.setdefault( feature, {} )[group] = 1
    else:
        map_feature_groups = map_group_features
    # perform the table regrouping
    new_table = regroup( table, map_feature_groups, args.function )
    new_table.write( args.output, unfloat=True )

if __name__ == "__main__":
    main( )
