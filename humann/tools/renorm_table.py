#! /usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import argparse
import sys
import re

from humann.tools import util

description = """
HUMAnN utility for renormalizing TSV files
===========================================
Each level of a stratified table will be 
normalized using the desired scheme.
"""

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_default_suffix = "-RPKs"
c_special = [
    util.c_unmapped, 
    util.c_unintegrated, 
    util.c_ungrouped,
    ]

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
        help="Original output table (tsv or biom format); default=[TSV/STDIN]",
        )
    parser.add_argument( 
        "-u", "--units", 
        choices=["cpm", "relab"],
        default="cpm",
        help="Normalization scheme: copies per million [cpm], relative abundance [relab]; default=[cpm]",
        )
    parser.add_argument( 
        "-m", "--mode", 
        choices=["community", "levelwise"],
        default="community",
        help="Normalize all levels by [community] total or [levelwise] totals; default=[community]",
        )
    parser.add_argument( 
        "-s", "--special", 
        choices=["y", "n"],
        default="y",
        help="Include the special features UNMAPPED, UNINTEGRATED, and UNGROUPED; default=[y]",
        )
    parser.add_argument( 
        "-p", "--update-snames", 
        action="store_true",
        help="Update '-RPK' in sample names to appropriate suffix; default=off",
        )
    parser.add_argument( 
        "-o", "--output", 
        default=None,
        help="Path for modified output table; default=[STDOUT]",
        )
    args = parser.parse_args()
    return args

def normalize ( table, cpm=True, levelwise=False, special=True ):
    divisor = 1e-6 if cpm else 1.0
    # remove special features?
    if not special:
        test = [rowhead.split( util.c_strat_delim )[0] not in c_special for rowhead in table.rowheads]
        for flag, rowhead in zip( test, table.rowheads ):
            if not flag:
                print( "Excluding special feature:", rowhead, file=sys.stderr )
        table.rowheads = [rowhead for i, rowhead in enumerate( table.rowheads ) if test[i]]
        table.data = [row for i, row in enumerate( table.data ) if test[i]]
    # compute totals by delim level
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
        # level=1 corresponds to the community level (no strata)
        totals = totals_by_level[level] if levelwise else totals_by_level[1]
        table.data[i] = ["%.6g" % ( row[j] / totals[j] / divisor ) for j in range( len( totals ) )]

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main ( ):
    args = get_args()
    table = util.Table( args.input )
    normalize( 
        table, 
        cpm = args.units=="cpm",
        levelwise = args.mode=="levelwise",
        special = args.special=="y",
        )
    if args.update_snames:
        for i, colhead in enumerate( table.colheads ):
            if re.search( c_default_suffix+"$", colhead ):
                table.colheads[i] = re.sub( c_default_suffix+"$", "-"+args.units.upper(), colhead )
            else:
                table.colheads[i] += "-"+args.units.upper()
    table.write( args.output )

if __name__ == "__main__":
    main()
