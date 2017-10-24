#!/usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import sys
import os
import argparse
import re

try:
    from humann2 import config
    from humann2.tools import util
    from humann2.tools.humann2_table import Table
except ImportError:
    sys.exit( "CRITICAL ERROR: Unable to find the HUMAnN2 python package.\n" +
              "Please check your install." )

description = util.wrap( """
HUMAnN2 utility for renaming table features

Attaches "glosses" (human-readable names) to HUMAnN2 features identified
by codes only (e.g. UniRef IDs). Maintains table stratifications.
""" )

# ---------------------------------------------------------------
# handling of renaming files
# ---------------------------------------------------------------

def name_match( fname ):
    return re.search( "map_(.*?)_name\..*", fname )

names = {}
p_root = os.path.join( os.path.dirname( os.path.abspath(__file__) ), os.pardir )
bundled_name_files = os.listdir( os.path.join( p_root, "data", "misc" ) )
for fname in bundled_name_files:
    match = name_match( fname )
    if match:
        name = match.group( 1 )
        names[name] = os.path.join( p_root, "data", "misc", fname )

# possible additional large mapping files
try:
    extra_name_files=os.listdir( config.utility_mapping_database )
except EnvironmentError:
    extra_name_files=[]
FOUND_EXTRA_FILES = False
for fname in extra_name_files:
    match = name_match( fname )
    if match:
        FOUND_EXTRA_FILES = True
        name = match.group( 1 )
        names[name] = os.path.join( config.utility_mapping_database, fname )
    
if not FOUND_EXTRA_FILES:
    description += """

You are missing the large renaming files. To acquire them, run the following command:

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
        "-n", "--names", 
        choices=names.keys( ),
        default=None,
        metavar="<choice>",
        help=util.pretty_grid( sorted( names.keys( ) ), desc="Select an available renaming option:" ),
        )
    parser.add_argument( 
        "-c", "--custom", 
        default=None,
        metavar="<path>",
        help="Custom mapping of feature IDs to full names (.tsv or .tsv.gz)",
        )
    parser.add_argument( 
        "-s", "--simplify", 
        action="store_true",
        help="Remove non-alphanumeric characters from names",
        )
    args = parser.parse_args()
    return args

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

def rename( table, polymap ):
    seen = set( )
    named = set( )
    new_data = {}
    for f in table.data:
        f_old = f_new = f
        fbase, fname, stratum = util.fsplit( f )
        if fbase not in util.c_topsort:
            seen.add( fbase )
            if fbase not in polymap:
                f_new = util.fjoin( fbase, name=util.c_no_name, stratum=stratum )
            else:
                named.add( fbase )
                choices = sorted( polymap[fbase] )
                if len( choices ) > 1:
                    print( "More than one name for:", fbase, file=sys.stderr )
                f_new = util.fjoin( fbase, name=choices[0], stratum=stratum )
        new_data[f_new] = table.data[f_old]
    # replace
    table.data = new_data
    # report on progress
    print( "Renaming report:", file=sys.stderr )
    a = len( seen )
    print( "  Candidate community features: {}".format( a ), file=sys.stderr )
    a = len( named )
    b = 100 * a / float( len( seen ) )
    print( "  Renamed: {} ({:.1f})%".format( a, b ), file=sys.stderr )
    # in place operation
    return None

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args()
    table = Table( args.input, last_metadata=args.last_metadata )
    allowed_keys = {util.fsplit( f )[0] for f in table.data}
    if args.custom is not None:
        polymap = util.load_polymap( args.custom, allowed_keys=allowed_keys )
    elif args.names is not None:
        polymap = util.load_polymap( names[args.names], allowed_keys=allowed_keys )
    else:
        sys.exit( "Must (i) choose names option or (ii) provide names file" )
    if args.simplify:
        for c, cnames in polymap.items( ):
            cnames = {re.sub( "[^A-Za-z0-9]+", "_", n ) for n in cnames}
            polymap[c] = names
    rename( table, polymap )
    table.write( args.output, unfloat=True )

if __name__ == "__main__":
    main( )
