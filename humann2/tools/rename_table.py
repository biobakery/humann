#! /usr/bin/env python

"""
HUMAnN2 utility for renaming TSV files
Run ./rename_table.py -h for usage help
"""

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import sys
import util
import argparse
    
def get_args ():
    parser = argparse.ArgumentParser()
    parser.add_argument( 
        "-i", "--input", 
        default=None,
        help="Original output table (.tsv format); default=[STDIN]",
        )
    parser.add_argument( 
        "-n", "--names", 
        help="Mapping of feature IDs to full names (.tsv or .tsv.gz)",
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
        old_name = items[0]
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

def main ( ):
    args = get_args()
    table = util.Table( args.input )
    polymap = util.load_polymap( args.names )
    rename( table, polymap )
    fh = open( args.output, "w" ) if args.output is not None else sys.stdout
    table.write( fh )

if __name__ == "__main__":
    main()
