#! /usr/bin/env python

"""
Make the mapping from metacyc pathway codes to pathway names
Parses the pathways.day file from metacyc
Works with humann_rename_table
"""

import os, sys, re, glob, argparse

parser = argparse.ArgumentParser()
parser.add_argument( 'pathwaydat', help='' )
args = parser.parse_args()

# constants
p_outfile = "map_metacyc-pwy_name.txt"

# parse out names
names = {}
with open( args.pathwaydat ) as fh:
    for line in fh:
        match = re.search( "^UNIQUE-ID - (.*)", line )
        if match:
            code = match.group( 1 )
        match = re.search( "^COMMON-NAME - (.*)", line )
        if match and code not in names:
            names[code] = match.group( 1 )

# clean up names
for code, name in names.items():
    # remove html
    name = re.sub( "<.*?>", "", name )
    names[code] = name

# output
with open( p_outfile, "w" ) as fh:
    for code in sorted( names ):
        print >>fh, "\t".join( [code, names[code]] )
