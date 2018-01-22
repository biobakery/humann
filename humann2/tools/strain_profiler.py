#! /usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED

import os
import sys
import argparse
import csv
from math import sqrt

try:
    from scipy.stats.mstats import mquantiles
except:
    sys.exit( "This module requires the Python scientific stack: numpy, scipy, and matplotlib." )

try:
    from humann2.tools import util
    from humann2.tools.humann2_table import Table
except ImportError:
    sys.exit( "CRITICAL ERROR: Unable to find the HUMAnN2 python package.\n" +
              "Please check your install." )

description = util.wrap( """
HUMAnN2 utility for building strain profiles from genefamilies output

The script looks for species in a given sample with a coverage "plateau",
which is expected for well-covered pangenomes with a single, dominant strain.
If inputting unnormalize RPK units, you can additionally enforce a height
for this plateau. For example, 20 RPK ~ 1x coverage assuming 100 nt reads. 
""" )

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_formats_help = """Select an output format for the strain profiles:
values: output the strain profiles in the original abundance units
binary: output the strain profiles in 1/0 (presence/absence) units
hybrid: output binary/original value pairs
[Default=binary]"""
c_binarize_help = """Select a binarization (1=presence; 0=absence) method:
half:    'present' defined as abund > 0.5 * plateau height   :: strict
sqrt:    'present' defined as abund > sqrt( plateau height ) :: lenient
nonzero: 'present' defined as abund > 0                      :: very lenient
[Default=half]"""

# ---------------------------------------------------------------
# utilities 
# ---------------------------------------------------------------

def get_args( ):
    """ Get args from Argparse """
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter
        )
    util.attach_common_arguments( parser, no_output=True )
    parser.add_argument( 
        "-o", "--outdir", 
        metavar="<path>",
        default=".",
        help="Directory to write strain profiles\n[Default=.]",
        )
    parser.add_argument( 
        "-g", "--min-nonzero-genes", 
        metavar="<int>",
        type=int,
        default=500,
        help=("To be considered, a species must recruit reads to this "
              "many gene families in a sample\n[Default=500]"),
        )
    parser.add_argument( 
        "-s", "--plateau-slope", 
        type=float,
        metavar="<float>",
        default=2.0,
        help=("To be considered, Q3/Q1 for non-zero genes cannot exceed this value\n"
              "[Default=2.0]"),
        )
    parser.add_argument( 
        "-a", "--plateau-height", 
        type=float,
        metavar="<float>",
        default=util.c_eps,
        help=("To be considered, non-zero genes must meet this median abundance\n"
              "[Default=No minimum]"),
        )
    parser.add_argument( 
        "-n", "--min-samples",
        type=int,
        metavar="<int>",
        default=1,
        help=("Only write strain profiles for strains detected in at least this many samples\n"
              "[Default=1]"),
        )
    parser.add_argument( 
        "-b", "--binarize",
        metavar="<half/sqrt/nonzero>",
        default="half",
        help=c_binarize_help,
        choices=["half", "sqrt", "nonzero"],
        )
    parser.add_argument( 
        "-z", "--fuzzy",
        action="store_true",
        help="When binarizing, set non-detect:non-zero genes to '0.5' and high outliers to '1.5'",
        )
    parser.add_argument( 
        "-f", "--format",
        metavar="<binary/values/hybrid>",
        default="binary",
        choices=["values", "hybrid", "binary"],
        help=c_formats_help,
        )
    args = parser.parse_args( )
    return args

def strainer( bug_table, args ):
    good_samples = {}
    for i, h in enumerate( bug_table.headers ):
        nonzero = [row[i] for row in bug_table.data.values( )]
        nonzero = [k for k in nonzero if k > 0]
        GOOD = False
        if len( nonzero ) >= args.min_nonzero_genes:
            q1, med, q3 = mquantiles( nonzero )
            if med >= args.plateau_height:
                if q3 / q1 <= args.plateau_slope:
                    good_samples[h] = [q1, med, q3]
    return good_samples

def reformat( bug_table, sample_quartiles, args ):
    for i, h in enumerate( bug_table.headers ):
        q1, med, q3 = sample_quartiles[h]
        upper_inner_fence = q3 + 1.5 * (q3 - q1)
        crit = util.c_eps
        if args.binarize == "sqrt":
            crit = sqrt( med )
        elif args.binarize == "half":
            crit = med / 2.0
        for f in bug_table.data:
            # convert from array to list
            row = bug_table.data[f] = list( bug_table.data[f] )
            binary = 1 if row[i] >= crit else 0
            # adjust for fuzzy binarization
            if args.fuzzy:
                if binary == 0 and row[i] > 0:
                    binary = 0.5
                elif row[i] > upper_inner_fence:
                    binary = 1.5
            # impose desired output format
            if args.format == "binary":
                row[i] = binary
            elif args.format == "hybrid":
                row[i] = "|".join( map( str, [row[i], binary] ) )
            elif args.format == "values":
                row[i] = str( row[i] )
    return None

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args( )
    table = Table( args.input, last_metadata=args.last_metadata )
    # regroup rows by bug
    bug_datas = {}
    for f in util.fsort( table.data ):
        fcode, fname, bug = util.fsplit( f )
        if bug is not None and "s__" in bug and "_unknown" not in fcode:
            inner = bug_datas.setdefault( bug, {} )
            inner[f] = table.data[f]
    # convert to tables
    bug_tables = {}
    for bug, data in bug_datas.items( ):
        if len( data ) >= args.min_nonzero_genes:
            bug_tables[bug] = Table( data, metadata=table.metadata, 
                                     headers=table.headers )
    # process each table
    total = len( bug_tables )
    index = 0
    for bug in sorted( bug_tables ):
        index += 1
        print( "Analyzing {: >3} of {}: {}".format( index, total, bug ), file=sys.stderr )
        filename = ".".join( [bug, "strains", args.format, "tsv"] )
        bt = bug_tables[bug]
        good_samples = strainer( bt, args )
        if len( good_samples ) >= args.min_samples:
            print( "  Identified {} good profiles".format( len( good_samples ) ), file=sys.stderr )
            # restrict to good samples
            bt = bt.resample( good_samples )
            # remove all-zero rows
            bt.data = {k:row for k, row in bt.data.items( ) if sum( row ) > 0}
            # change value encoding
            reformat( bt, good_samples, args )
            # write the table
            bt.write( os.path.join( args.outdir, filename ) )
        else:
            print( "  Not enough good profiles ({})".format( len( good_samples ) ), file=sys.stderr )
    print( "Finished successfully.", file=sys.stderr )

if __name__ == "__main__":
    main( )
