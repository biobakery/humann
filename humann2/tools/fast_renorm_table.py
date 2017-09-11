#! /usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import argparse
import sys
import re

#from humann2.tools import util
#from humann2.tools.arraytable import ArrayTable
import util
from better_table import Table

description = """
HUMAnN2 utility for renormalizing TSV files
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

def get_args( ):
    """ Get args from Argparse """
    parser = argparse.ArgumentParser( 
        description=description, 
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument( 
        "-i", "--input",
        metavar="<path>",
        default=None,
        help="Original output table (tsv or biom format); default=[TSV/STDIN]",
        )
    parser.add_argument( 
        "-u", "--units", 
        choices=["cpm", "relab"],
        default="cpm",
        metavar="<choice>",
        help="Normalization scheme: copies per million [cpm], relative abundance [relab]; default=[cpm]",
        )
    parser.add_argument( 
        "--mode", 
        choices=["community", "levelwise"],
        default="community",
        metavar="<choice>",
        help="Normalize all levels by [community] total or [levelwise] totals; default=[community]",
        )
    parser.add_argument( 
        "--exclude-special",
        action="store_true",
        help="Exclude the special features UNMAPPED, UNINTEGRATED, and UNGROUPED; default=[include these features]",
        )
    parser.add_argument( 
        "--update-sample-names", 
        action="store_true",
        help="Update '-RPK' in sample names to appropriate suffix; default=[do not change names]",
        )
    parser.add_argument( 
        "-o", "--output", 
        default=None,
        metavar="<path>",
        help="Path for modified output table; default=[STDOUT]",
        )
    args = parser.parse_args( )
    return args

def normalize( table, cpm=True, levelwise=False, exclude_special=False ):
    new_data = {}
    # compute totals by level
    totals_by_level = {}
    for f in util.fsort( table.data ):
        if exclude_special and f in c_special:
            print( "Excluding special feature:", f, file=sys.stderr )
            continue
        level = len( f.split( util.c_strat_delim ) )
        if level not in totals_by_level:
            totals_by_level[level] = table.zeroes[:]
        totals_by_level[level] += table.data[f]
        new_data[f] = table.data[f]
    # check for zero totals
    levels = [1] if not levelwise else sorted( levels )
    for l in levels:
        totals = totals_by_level[l]
        if min( totals_by_level[l] ) == 0:
            for i, v in enumerate( totals ):
                if v == 0:
                     # avoid divide by zero
                    totals[i] = 1.0
                    print( "WARNING: Sample {} ({}) has zero sum at level {}".format( 
                            i+1, table.headers[i], l ), file=sys.stderr )                   
    # normalize
    for f in util.fsort( new_data ):
        l = 1 if not levelwise else len( f.split( util.c_strat_delim ) )
        new_data[f] /= totals_by_level[l]
        new_data[f] *= 1e6 if cpm else 1
    # new table
    return Table( new_data, metadata=table.metadata, headers=table.headers )

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args()
    table = Table( args.input )
    table = normalize(
        table, 
        cpm = args.units=="cpm",
        levelwise = args.mode=="levelwise",
        exclude_special = args.exclude_special,
        )
    if args.update_sample_names:
        for i, header in enumerate( table.headers ):
            if re.search( c_default_suffix+"$", header ):
                table.headers[i] = re.sub( c_default_suffix+"$", "-"+args.units.upper( ), header )
            else:
                table.headers[i] += "-"+args.units.upper( )
    table.write( args.output, unfloat=True )

if __name__ == "__main__":
    main( )
