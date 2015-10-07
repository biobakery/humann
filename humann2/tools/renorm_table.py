#! /usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import argparse
import sys
import util

description = """
HUMAnN2 utility for renormalizing TSV files
===========================================
Each level of a stratified table will be 
normalized using the desired scheme.
"""

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
        "-n", "--norm", 
        choices=["cpm", "relab"],
        default="cpm",
        help="Normalization scheme: copies per million [cpm], relative abundance [relab]; default=[cpm]",
        )
    parser.add_argument( 
        "-o", "--output", 
        default=None,
        help="Path for modified output table; default=[STDOUT]",
        )
    args = parser.parse_args()
    return args

def normalize ( table, cpm=False ):
    divisor = 1e-6 if cpm else 1.0
    # compute totals by delim level (e.g. community vs species)
    totals_by_level = {}
    for i, row in enumerate( table.data ):
        level = len( table.rowheads[i].split( util.c_strat_delim ) )
        if level not in totals_by_level:
            totals_by_level[level] = [0 for k in range( len( table.colheads ) )]
        table.data[i] = [float( k ) for k in row]
        totals_by_level[level] = [k1 + k2 for k1, k2 in zip( totals_by_level[level], table.data[i] )]
    # check for sample / level combinations with zero sum
    for level in sorted( totals_by_level ):
        totals = totals_by_level[level]
        for j, total in enumerate( totals ):
            if total == 0:
                totals[j] = 1
                print( "WARNING: Column {} ({}) has zero sum at level {}".format( \
                        j+1, table.colheads[j], level ), file=sys.stderr )               
    # normalize
    for i, row in enumerate( table.data ):
        level = len( table.rowheads[i].split( util.c_strat_delim ) )
        totals = totals_by_level[level]
        table.data[i] = ["%.6g" % ( row[j] / totals[j] / divisor ) for j in range( len( totals ) )] 

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main ( ):
    args = get_args()
    table = util.Table( args.input )
    normalize( table, cpm=True if args.norm == "cpm" else False )
    table.write( args.output )

if __name__ == "__main__":
    main()
