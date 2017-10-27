#!/usr/bin/env python

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
HUMAnN2 utility for filtering features from a table

Includes options for filtering based on abundance/prevalence,
feature name, and feature stratification status. Does not
perturb table metadata if present. With the exception of stratification
status, all other filters are carried out on community totals, and
stratifications within that total will be kept/discarded if the
corresponding total is kept/discarded.
""" )

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

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
        "-p", "--min-prevalence",
        default=0,
        metavar="<0-1.0>",
        type=float,
        help="Remove features that are detected in less than p (fraction) of samples",
        )
    parser.add_argument( 
        "-n", "--min-samples",
        default=0,
        metavar="<int>",
        type=int,
        help="Remove features that are detected in less than s samples",
        )
    parser.add_argument( 
        "-a", "--abund-detect",
        default=util.c_eps,
        metavar="<float>",
        type=float,
        help="Abundance threshold for detection\n[Default: >0]",
        )
    parser.add_argument( 
        "-x", "--exclude-strata",
        action="store_true",
        help="Remove feature strata (keep only totals)",
        )
    parser.add_argument( 
        "-X", "--exclude-totals",
        action="store_true",
        help="Remove feature totals (keep only strata)",
        )
    parser.add_argument( 
        "-g", "--grep", 
        default=None,
        type=str,
        metavar="<pattern>",
        help="Keep only features whose IDs/names match the specified pattern",
        )
    parser.add_argument( 
        "-G", "--strata-grep", 
        default=None,
        type=str,
        metavar="<pattern>",
        help="Keep only stratified features whose species names match the specified pattern",
        )
    args = parser.parse_args( )
    return args

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args()
    table = Table( args.input, last_metadata=args.last_metadata )
    # process table totals
    allowed = set( )
    for f in util.fsort( table.data ):
        INCLUDE = True
        fcode, fname, stratum = util.fsplit( f )
        ffull = util.fjoin( fcode, fname )
        if stratum is not None:
            continue
        values = table.data[f]
        # check prevalence
        abund_detect = args.abund_detect
        n = sum( [1 for k in values if k >= abund_detect] )
        if n < args.min_samples:
            INCLUDE = False
        if n / float( len( values ) ) < args.min_prevalence:
            INCLUDE = False
        # check pattern
        if args.grep is not None and not re.search( args.grep, ffull ):
            INCLUDE = False
        # decide on this total
        if INCLUDE:
            allowed.add( fcode )
    # process strata
    data2 = {}
    for f in util.fsort( table.data ):
        fcode, fname, stratum = util.fsplit( f )
        if stratum is None:
            if fcode in allowed and not args.exclude_totals:
                data2[f] = table.data[f]
        else:
            if fcode in allowed and not args.exclude_strata:
                if args.strata_grep is None or re.search( args.strata_grep, stratum ):
                    data2[f] = table.data[f]
    # rebuild table
    table = Table( data2, metadata=table.metadata, headers=table.headers )
    table.write( args.output, unfloat=True )

if __name__ == "__main__":
    main( )
