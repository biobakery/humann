#!/usr/bin/env python

from __future__ import print_function #PYTHON 2.7+ REQUIRED
import os
import sys
import argparse
import glob
import re

include = """
map_go_name.txt.gz
map_go_uniref50.txt.gz
map_go_uniref90.txt.gz
map_infogo1000_uniref50.txt.gz
map_infogo1000_uniref90.txt.gz
map_ko_name.txt.gz
map_ko_uniref50.txt.gz
map_ko_uniref90.txt.gz
map_level4ec_name.txt.gz
map_level4ec_uniref50.txt.gz
map_level4ec_uniref90.txt.gz
map_pfam_name.txt.gz
map_pfam_uniref50.txt.gz
map_pfam_uniref90.txt.gz
map_eggnog_name.txt.gz
map_eggnog_uniref50.txt.gz
map_eggnog_uniref90.txt.gz
map_uniref50_name.txt.bz2
map_uniref90_name.txt.bz2
map_metacyc-pwy_name.txt.gz
map_metacyc-rxn_name.txt.gz
map_kegg-pwy_name.txt.gz
map_kegg-mdl_name.txt.gz
uniref50-tol-lca.dat.gz
uniref90-tol-lca.dat.gz
"""
include = [k for k in include.split( "\n" ) if k != ""]

# argument parsing (python argparse)
parser = argparse.ArgumentParser()
parser.add_argument( "humann2_path", 
                     help="path to target humann2 repo" )
parser.add_argument( "--tarball", 
                     default="utility_mapping.tar.gz", 
                     help="path to output tarball" )
parser.add_argument( "--ignore-extras", 
                     action="store_true", 
                     help="ignore extra mapping files in data/misc" )
args = parser.parse_args( )

p_data_misc = os.path.join( args.humann2_path, "humann2", "data", "misc" )
p_tools = os.path.join( args.humann2_path, "humann2", "tools" )

# check that all the include files are found
for p in include:
    if not os.path.exists( os.path.join( p_data_misc, p ) ):
        sys.exit( "Expected file <{}> not found in <{}>\n".format( p, p_data_misc ) )
else:
    print( "\nAll expected files found.\n", file=sys.stderr )

# warn if included renames aren't registered with rename
check = {k:False for k in include if "_name." in k}
with open( os.path.join( p_tools, "rename_table.py" ) ) as fh:
    for line in fh:
        for k in check:
            if k in line:
                check[k] = True
for k in check:
    if not check[k]:
        print( "Warning: {} not registered in {}".format(
                k, "rename_table.py" ), file=sys.stderr )

# warn if included regroups aren't registered with rename
check = {k:False for k in include if "map_" in k and "_name." not in k}
with open( os.path.join( p_tools, "regroup_table.py" ) ) as fh:
    for line in fh:
        for k in check:
            if k in line:
                check[k] = True
for k in check:
    if not check[k]:
        print( "Warning: {} not registered in {}".format(
                k, "group_table.py" ), file=sys.stderr )

extras = False
# look for extra naming files
for p in glob.glob( os.path.join( p_data_misc, "map_*_name.*" ) ):
    basename = os.path.split( p )[1]
    if basename not in include:
        print( "Extra naming file:", basename, file=sys.stderr )
        extras = True
# look for extra regroup files
for p in glob.glob( os.path.join( p_data_misc, "map_*_*" ) ):
    basename = os.path.split( p )[1]
    if basename not in include and "_name." not in basename:
        print( "Extra regroup file:", basename, file=sys.stderr )
        extras = True
if extras and not args.ignore_extras:
    sys.exit( """
Extra mapping files found in <{}>.
Edit <{}> to include these or run with "--ignore-extras".
""".format( p_data_misc, sys.argv[0] ) )

# make the tarball
print( "Making tarball: <{}>".format( args.tarball ), file=sys.stderr )
os.system( "tar czfv {} -C {} {}".format( args.tarball, p_data_misc, " ".join( include ) ) )
