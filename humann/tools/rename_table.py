#! /usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
from collections import namedtuple
import sys
import os
import argparse
import re

try:
    from humann2 import config
    from humann2.tools import util
except ImportError:
    sys.exit("CRITICAL ERROR: Unable to find the HUMAnN2 python package." +
        " Please check your install.") 

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
    "kegg-orthology": Names(
        os.path.join( p_root, "data", "misc", "map_ko_name.txt.gz" ) ),
    "kegg-pathway": Names(
        os.path.join( p_root, "data", "misc", "map_kegg-pwy_name.txt.gz" ) ),
    "kegg-module": Names(
        os.path.join( p_root, "data", "misc", "map_kegg-mdl_name.txt.gz" ) ),
    "ec": Names(
        os.path.join( p_root, "data", "misc", "map_level4ec_name.txt.gz" ) ),
    "metacyc-rxn": Names(
        os.path.join( p_root, "data", "misc", "map_metacyc-rxn_name.txt.gz" ) ),    
    "metacyc-pwy": Names(
        os.path.join( p_root, "data", "misc", "map_metacyc-pwy_name.txt.gz" ) ),
    "pfam": Names(
        os.path.join( p_root, "data", "misc", "map_pfam_name.txt.gz" ) ),
    "eggnog": Names(
        os.path.join( p_root, "data", "misc", "map_eggnog_name.txt.gz" ) ),
    "go": Names(
        os.path.join( p_root, "data", "misc", "map_go_name.txt.gz" ) ),
    # infogo1000 is just a subset of go, but adding a separate entry for consistency
    "infogo1000": Names(
        os.path.join( p_root, "data", "misc", "map_go_name.txt.gz" ) )}

# get a list of all available script mapping files
try:
    all_mapping_files=os.listdir(config.utility_mapping_database)
except EnvironmentError:
    all_mapping_files=[]

# add the options for the larger mapping files if they are present
larger_mapping_files_found=False
if "map_uniref50_name.txt.bz2" in all_mapping_files:
    c_default_names["uniref50"]=Names(os.path.join(config.utility_mapping_database,"map_uniref50_name.txt.bz2"))
    larger_mapping_files_found=True
if "map_uniref90_name.txt.bz2" in all_mapping_files:
    c_default_names["uniref90"]=Names(os.path.join(config.utility_mapping_database,"map_uniref90_name.txt.bz2"))
    larger_mapping_files_found=True
    
if not larger_mapping_files_found:
    description+="""

For additional name mapping files, run the following command:
$ humann2_databases --download utility_mapping full $DIR
Replacing, $DIR with the directory to download and install the databases."""

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
        help="Original output table (tsv or biom format); default=[TSV/STDIN]",
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
        # don't touch special features like UNMAPPED
        if old_name not in util.c_topsort:
            items[0] = util.c_name_delim.join( [old_name, new_name] )
        table.rowheads[i] = util.c_strat_delim.join( items )
    tcount = list(seen.values()).count( True )
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
