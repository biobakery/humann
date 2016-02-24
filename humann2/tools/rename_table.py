#! /usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
from collections import namedtuple
import sys
import os
import util
import argparse
import re

description = """
HUMAnN2 utility for renaming table features
===========================================
"""

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

p_root = os.path.join( os.path.dirname( os.path.abspath(__file__) ), os.pardir )
Names = namedtuple( "Names", ["path"] )
c_default_names = {
    "uniref50": Names(
        os.path.join( p_root, "data", "misc", "map_uniref50_name.txt.bz2" ) ),
    "kegg-orthology": Names(
        os.path.join( p_root, "data", "misc", "map_ko_name.txt.gz" ) ),
    "kegg-pathway": Names(
        os.path.join( p_root, "data", "misc", "map_kegg_pathways_name.txt.gz" ) ),
    "kegg-module": Names(
        os.path.join( p_root, "data", "misc", "map_kegg_modules_name.txt.gz" ) ),
    "ec": Names(
        os.path.join( p_root, "data", "misc", "map_ec_name.txt.gz" ) ),
    "metacyc-rxn": Names(
        os.path.join( p_root, "data", "misc", "map_metacyc-rxn_name.txt.gz" ) ),    
    "metacyc-pwy": Names(
        os.path.join( p_root, "data", "misc", "map_metacyc-pwy_name.txt.gz" ) ),
    }

# ---------------------------------------------------------------
# utilities 
# ---------------------------------------------------------------
    
def get_args ():
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
        "-n", "--names", 
        choices=c_default_names.keys(),
        default=None,
        help="Table features that can be renamed with included data files",
        )
    parser.add_argument( 
        "-c", "--custom", 
        default=None,
        help="Custom mapping of feature IDs to full names (.tsv or .tsv.gz)",
        )
    parser.add_argument( 
        "-s", "--simplify", 
        action="store_true",
        help="Remove non-alphanumeric characters from names",
        )
    parser.add_argument( 
        "-o", "--output", 
        default=None,
        help="Path for modified output table; default=[STDOUT]",
        )
    args = parser.parse_args()
    return args

def rename ( table, polymap ):
    seen = {}
    for i, rowhead in enumerate( table.rowheads ):
        items = rowhead.split( util.c_strat_delim )
        # account for previously renamed features
        old_name = items[0].split( util.c_name_delim )[0]
        seen[old_name] = False
        if old_name in polymap:
            new_name = util.c_multiname_delim.join( polymap[old_name].keys() )
            seen[old_name] = True
        else:
            new_name = util.c_str_unknown
        items[0] = util.c_name_delim.join( [old_name, new_name] )
        table.rowheads[i] = util.c_strat_delim.join( items )
    tcount = seen.values().count( True )
    print( "Renamed %d of %d entries (%.2f%%)" \
           % ( tcount, len( seen ), 100 * tcount / float( len( seen ) ) ), 
           file=sys.stderr )

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main ( ):
    args = get_args()
    table = util.Table( args.input )
    allowed_keys = {k.split( util.c_strat_delim )[0]:1 for k in table.rowheads}
    if args.custom is not None:
        polymap = util.load_polymap( args.custom, allowed_keys=allowed_keys )
    elif args.names is not None:
        polymap = util.load_polymap( c_default_names[args.names].path, allowed_keys=allowed_keys )
    else:
        sys.exit( "Must (i) choose names option or (ii) provide names file" )
    if args.simplify:
        for c, ndict in polymap.items():
            ndict = {re.sub( "[^A-Za-z0-9]+", "_", n ):1 for n in ndict}
            polymap[c] = ndict
    rename( table, polymap )
    table.write( args.output )

if __name__ == "__main__":
    main()
