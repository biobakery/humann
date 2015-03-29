#! /usr/bin/env python

"""
HUMAnN2 utility for renormalizing TSV files
Run ./renorn_table.py -h for usage help
"""

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import argparse
import sys
import util

def get_args ():
    """ Get args from Argparse """
    parser = argparse.ArgumentParser()
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
    divisor = 1 
    if table.is_stratified:
        divisor /= 2.0
    if cpm:
        divisor /= 1e6
    totals = [0 for k in range( len( table.colheads ) )]
    for i, row in enumerate( table.data ):
        table.data[i] = [float( k ) for k in row]
        totals = [k1 + k2 for k1, k2 in zip( totals, table.data[i] )]
    for i, row in enumerate( table.data ):
        table.data[i] = ["%.6g" % ( row[j] / totals[j] / divisor ) for j in range( len( totals ) )] 

def main ( ):
    args = get_args()
    table = util.Table( args.input )
    normalize( table, cpm=True if args.norm == "cpm" else False )
    table.write( args.output )

if __name__ == "__main__":
    main()
