#! /usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import argparse
import sys
import re

try:
    from humann2 import config
    from humann2.tools import util
    from humann2.tools.humann2_table import Table
except ImportError:
    sys.exit( "CRITICAL ERROR: Unable to find the HUMAnN2 python package.\n" +
              "Please check your install." )

description = util.wrap( """
HUMAnN2 utility for renormalizing output tables

Each level of a stratified table will be normalized using the desired scheme.
Relative abundance (relab) units causes normalized levels to sum to 1.0.
Copies per million (CPM) units causes normalized levels to sum to 1 million.
CPMs are more convenient for tables with many features (e.g. gene families).
Note that in this case CPM does NOT refer to _counts_ per million.
""" )

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
# command-line interface
# ---------------------------------------------------------------

def get_args( ):
    """ Get args from Argparse """
    parser = argparse.ArgumentParser( 
        description=description, 
        formatter_class=argparse.RawTextHelpFormatter,
    )
    util.attach_common_arguments( parser )
    parser.add_argument( 
        "-u", "--units", 
        choices=["cpm", "relab"],
        default="cpm",
        metavar="<cpm/relab>",
        help="Normalization scheme: copies per million [cpm], relative abundance [relab]\nDefault=[cpm]",
        )
    parser.add_argument( 
        "-m", "--mode", 
        choices=["community", "levelwise"],
        default="community",
        metavar="<community/levelwise>",
        help="Normalize all levels by [community] total or [levelwise] totals\nDefault=[community]",
        )
    parser.add_argument( 
        "-x", "--exclude-special",
        help="Exclude the special features UNMAPPED, UNINTEGRATED, and UNGROUPED\nDefault=[include these features]",
        action="store_true",
        )
    parser.add_argument( 
        "-n", "--update-sample-names", 
        action="store_true",
        help="Update '-RPK' in sample names to appropriate suffix\nDefault=[do not change names]",
        )
    args = parser.parse_args( )
    return args

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

def normalize( table, cpm=True, levelwise=False, exclude_special=False ):
    new_data = {}
    # compute totals by level
    totals_by_level = {}
    for f in util.fsort( table.data ):
        fbase, name, strat = util.fsplit( f )        
        if exclude_special and fbase in c_special:
            print( "Excluding special feature:", f, file=sys.stderr )
            continue
        level = len( f.split( util.c_strat_delim ) )
        if level not in totals_by_level:
            totals_by_level[level] = table.zeros( )
        totals_by_level[level] += table.data[f]
        new_data[f] = table.data[f]
    # check for zero totals
    levels = [1] if not levelwise else sorted( totals_by_level )
    for l in levels:
        totals = totals_by_level[l]
        if min( totals ) == 0:
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
    table = Table( args.input, last_metadata=args.last_metadata )
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
