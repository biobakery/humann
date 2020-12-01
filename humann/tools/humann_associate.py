#!/usr/bin/env python

import os
import sys
import csv
import argparse

try:
    import numpy as np
    from scipy.stats import spearmanr
    from scipy.stats.mstats import kruskalwallis
except ImportError:
    sys.exit( "This script requires the Python scientific stack: numpy and scipy." )

try:
    from humann import config
    from humann.tools import util
except ImportError:
    sys.exit( "CRITICAL ERROR: Unable to find the HUMAnN python package." +
              " Please check your install.")

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------
   
c_allowed_impute = 0.2

# ---------------------------------------------------------------
# argument parsing
# ---------------------------------------------------------------

description = """
HUMAnN utility for performing metadata association
===================================================
"""

def get_args( ):
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter,
        )
    parser.add_argument( "-i", "--input", 
                         metavar="<path>",
                         required=True,
                         help="HUMAnN table with metadata rows at the top", )
    parser.add_argument( "-m", "--focal-metadatum",
                         metavar="<str>",
                         required=True,
                         help="Indicate metadatum to test vs. community feature totals", )
    parser.add_argument( "-l", "--last-metadatum",
                         metavar="<str>",
                         required=True,
                         help="Indicate end of metadata rows", )
    parser.add_argument( "-t", "--focal-type",
                         required=True,
                         choices=["continuous", "categorical"],
                         help="Metadatum type", )
    parser.add_argument( "-o", "--output",
                         metavar="<path>",
                         default="associations.tsv",
                         help="Where to save the output", )
    parser.add_argument( "-f", "--fdr",
                         metavar="<float>",
                         type=float,
                         default=0.2,
                         help="FDR threshold (default=0.2)", )
    return parser.parse_args( )

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

def pvalues2qvalues( pvalues ):
    n = len( pvalues )
    # after sorting, index[i] is the original index of the ith-ranked value
    index = range( n )
    index = sorted( index, key=lambda i: pvalues[i] )
    pvalues = sorted( pvalues )
    qvalues = [pvalues[i-1] * n / i for i in range( 1, n+1 )]
    # adjust qvalues to enforce monotonic behavior
    # q( i ) = min( q( i..n ) )
    qvalues.reverse()
    for i in range( 1, n ):
        if qvalues[i] > qvalues[i-1]:
            qvalues[i] = qvalues[i-1]
    qvalues.reverse()
    # rebuild qvalues in the original order
    ordered_qvalues = [None for q in qvalues]
    for i, q in enumerate( qvalues ):
        ordered_qvalues[index[i]] = q
    return ordered_qvalues

def adjust_stats( stats ):
    pvalues = [stat[-1] for stat in stats]
    qvalues = pvalues2qvalues( pvalues )
    for i in range( len( stats ) ):
        stats[i].append( qvalues[i] )
    return sorted( stats, key=lambda stat: stat[-1] )

def spearman_analysis( mvalues, fnames, fvalues ):
    stats = []
    for fname, frow in zip( fnames, fvalues ):
        try:
            rho, p = spearmanr( mvalues, frow )
            stats.append( [fname, "%.4g" % rho, p] )
        except:
            sys.stderr.write( "NOTE: Unable to compute Spearman r with feature: " + fname +"\n" )
    return adjust_stats( stats )

def shatter( cats, values ):
    lists = {}
    for c, v in zip( cats, values ):
        lists.setdefault( c, [] ).append( v )
    return lists

def kruskalwallis_analysis( mvalues, fnames, fvalues ):
    stats = []
    for fname, frow in zip( fnames, fvalues ):
        try:
            lists = shatter( mvalues, frow )
            summary = {k:"%.4g" % ( np.mean( v ) ) for k, v in lists.items( )}
            summary = [":".join( [k, v] ) for k, v in summary.items( )]
            summary = "|".join( summary )
            hstat, p = kruskalwallis( *lists.values( ) )
            stats.append( [fname, summary, p] )
        except:
            sys.stderr.write("NOTE: Unable to compute Kruskal-Wallis with feature: " + fname + "\n")
    return adjust_stats( stats )

def test_float_list( values ):
    good = []
    values2 = []
    for v in values:
        try:
            v = float( v )
            values2.append( v )
            good.append( v )           
        except:
            values2.append( None )
    if len( good ) / float( len( values2 ) ) >= 1 - c_allowed_impute:
        impute = np.mean( good )
        values = [k if k is not None else impute for k in values2]
        return values
    else:
        return None

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args( )
    mname, mvalues = None, []
    fnames, fvalues = [], []
    adding = False
    with open( args.input , "rt" ) as fh:
        for row in csv.reader( fh, csv.excel_tab ):
            header, values = row[0], row[1:]
            if header == args.focal_metadatum:
                mname, mvalues = header, values
            if header == args.last_metadatum:
                adding = True
                continue
            if adding and "|" not in header:
                fnames.append( header )
                fvalues.append( list(map( float, values ) ) )
    # tests
    if not adding:
        sys.exit( "STOPPED: last metadata row <{}> not found.\n".format( args.last_metadatum ) )
    if mname is None:
        sys.exit( "STOPPED: focal metadata row <{}> not found.\n".format( args.focal_metadatum ) )
    if args.focal_type == "categorical" and len( set( mvalues ) ) > np.sqrt( len( mvalues ) ):
        sys.stderr.write( "WARNING: categorical metadata <{}> has many distinct values.\n".format( args.focal_metadatum ) )
    # begin analysis
    fh = open( args.output, "w" ) if args.output is not None else sys.stdout
    if args.focal_type == "continuous":
        mvalues = test_float_list( mvalues )
        if mvalues is None:
            sys.exit( "STOPPED: failed to float many entries in focal row <{}>".format( args.focal_metadatum ) )
        stats = spearman_analysis( mvalues, fnames, fvalues )
        sys.stderr.write( "Performing Spearman analysis vs. metadatum: " + mname + "\n" )
        fh.write( "# Feature\tRho\tP-value\tQ-value\n" )
    elif args.focal_type == "categorical":
        stats = kruskalwallis_analysis( mvalues, fnames, fvalues )
        sys.stderr.write( "Performing Kruskal-Wallis analysis vs. metadatum: " + mname + "\n" )
        fh.write( "# Feature\tLevel means (|ed)\tP-value\tQ-value\n" )
    found_something = False
    for stat in stats:
        if stat[-1] <= args.fdr:
            found_something = True
            stat[-1] = "%.4g" % stat[-1]
            stat[-2] = "%.4g" % stat[-2]
            fh.write( "\t".join( list( map( str, stat ) ) ) + "\n" )
    if not found_something:
        sys.stderr.write( "NO FDR SIGNIFICANT ASSOCIATIONS\n" )
    fh.close( )
    sys.stderr.write( "Finished successfully.\n" )
        
if __name__ == "__main__":
    main( )
