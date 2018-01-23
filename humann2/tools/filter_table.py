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

try:
    import numpy as np
except:
    sys.exit( "This module requires the Python scientific stack: numpy, scipy, and matplotlib." )

# ---------------------------------------------------------------
# description
# ---------------------------------------------------------------

description = util.wrap( """
HUMAnN2 utility for filtering features from a table

Includes options for filtering based on abundance/prevalence,
feature name, feature stratification, and dominance of a feature by
one of its strata. Does not perturb table metadata if present. 
With the exception of stratification status/name, all other filters 
are carried out on community totals, and stratifications within 
that total will be kept/discarded if the corresponding total 
is kept/discarded.
""" )

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# command-line interface
# ---------------------------------------------------------------

def get_args( ):
    parser = argparse.ArgumentParser( 
        description=description, 
        formatter_class=argparse.RawTextHelpFormatter,
    )
    util.attach_common_arguments( parser )
    parser.add_argument( 
        "-p", "--min-prevalence",
        type=float,
        default=0,
        metavar="<0-100>",
        help="Remove features that are detected in < p%% of samples",
        )
    parser.add_argument( 
        "-n", "--min-samples",
        type=int,
        default=0,
        metavar="<int>",
        help="Remove features that are detected in < n samples",
        )
    parser.add_argument( 
        "-a", "--abund-detect",
        type=float,
        default=util.c_eps,
        metavar="<float>",
        help="Abundance threshold for detection\n[Default: non-zero]",
        )
    parser.add_argument( 
        "-d", "--filter-dominated",
        type=float,
        default=None,
        metavar="<0-100>",
        help="REMOVE features that are > p%% explained by a single non-unclassified stratum in > p%% of samples\n[Default: off]",
        )
    parser.add_argument( 
        "-u", "--filter-unclassified",
        type=float,
        default=None,
        metavar="<0-100>",
        help="REMOVE features that are > p%% unclassified in > p%% of samples\n[Default: off]",
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
        type=str,
        default=None,
        metavar="<pattern>",
        help="Keep only stratified features whose species names match the specified pattern",
        )
    parser.add_argument( 
        "-v", "--invert",
        action="store_true",
        help="Return rows that failed to match one or more filters (does not affect metadata)",
        )
    args = parser.parse_args( )
    return args

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

def is_dominated( totals, stratum_values, perc ):
    ret = False
    index = [i for i, v in enumerate( totals ) if v > util.c_eps]
    if len( index ) > 0:
        totals = totals[index]
        stratum_values = stratum_values[index]
        stratum_values /= totals
        if np.percentile( 100 * stratum_values, perc ) >= perc:
            ret = True
    return ret

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args( )
    table = Table( args.input, last_metadata=args.last_metadata )
    # evaluate features
    bad_features = set( )
    for f in util.fsort( table.data ):
        fcode, fname, stratum = util.fsplit( f )
        fbase = util.fjoin( fcode, fname )
        # properties of the total that can kill the total
        if stratum is None:
            values = table.data[f]
            # check prevalence
            abund_detect = args.abund_detect
            n = sum( [1 for k in values if k >= abund_detect] )
            if n < args.min_samples:
                bad_features.add( fcode )
            if 100 * n / float( len( values ) ) < args.min_prevalence:
                bad_features.add( fcode )
            # check pattern
            if args.grep is not None and not re.search( args.grep, fbase ):
                bad_features.add( fcode )
        # properties of a stratum that can kill the total: too much unclassified
        if stratum == util.c_unclassified and args.filter_unclassified is not None:
            perc = args.filter_unclassified
            totals = table.data[fbase]
            stratum_values = table.data[f]
            if is_dominated( totals, stratum_values, perc ):
                bad_features.add( fcode )
        # properties of a stratum that can kill the total: too much something else
        if stratum not in [None, util.c_unclassified] and args.filter_dominated is not None:
            perc = args.filter_dominated
            totals = table.data[fbase]
            stratum_values = table.data[f]
            if is_dominated( totals, stratum_values, perc ):
                bad_features.add( fcode )
        # properties of a stratum that can kill the stratum: pattern
        if stratum is not None and args.strata_grep is not None \
                and not re.search( args.strata_grep, stratum ):
            bad_features.add( f )
    # rebuild table
    data2 = {}
    for f in util.fsort( table.data ):
        fcode, fname, stratum = util.fsplit( f )
        if fcode not in bad_features:
            if stratum is None and not args.exclude_totals:
                data2[f] = table.data[f]
            if stratum is not None and not args.exclude_strata:
                data2[f] = table.data[f]
    if args.invert:
        data2 = {k:row for k, row in table.data.items( ) if k not in data2}
    table = Table( data2, metadata=table.metadata, headers=table.headers )
    table.write( args.output, unfloat=True )

if __name__ == "__main__":
    main( )
