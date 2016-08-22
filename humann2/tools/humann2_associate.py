#!/usr/bin/env python

"""
HUMAnN2 Plotting Tool
===============================================
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

import os
import sys
import csv
import argparse

try:
    import numpy as np
    from scipy.stats import spearmanr
    from scipy.stats.mstats import kruskalwallis
except:
    sys.exit( "This script requires the Python scientific stack: numpy and scipy." )
    
# ---------------------------------------------------------------
# argument parsing
# ---------------------------------------------------------------

def get_args( ):
    parser = argparse.ArgumentParser()
    parser.add_argument( "-i", "--input", 
                         help="HUMAnN2 table with metadata", )
    parser.add_argument( "-l", "--last-metadatum",
                         help="Indicate end of metadata rows", )
    parser.add_argument( "-m", "--focal-metadatum",
                         help="Indicate metadatum to test vs. community feature totals", )
    parser.add_argument( "-t", "--focal-type",
                         choices=["continuous", "categorical"],
                         default="categorical",
                         help="Metadatum type", )
    parser.add_argument( "-o", "--output",
                         default=None,
                         help="Where to save the output", )
    return parser.parse_args()

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
    for fname, fvalue in zip( fnames, fvalues ):
        try:
            rho, p = spearmanr( mvalues, fvalues )
            stats.append( ["%.4g" % rho, p] )
        except:
            sys.stderr.write("Unable to compute spearman r with feature: " + fname +"\n")
    return adjust_stats( stats )

def shatter( cats, values ):
    lists = {}
    for c, v in zip( cats, values ):
        lists.setdefault( c, [] ).append( v )
    return lists

def kruskalwallis_analysis( mvalues, fnames, fvalues ):
    stats = []
    for fname, fvalue in zip( fnames, fvalues ):
        try:
            lists = shatter( mvalues, fvalue )
            summary = {k:"%.4g" % ( np.mean( v ) ) for k, v in lists.items( )}
            summary = [":".join( [k, v] ) for k, v in summary.items( )]
            summary = "|".join( summary )
            hstat, p = kruskalwallis( *lists.values( ) )
            stats.append( [fname, summary, p] )
        except:
            sys.stderr.write("Unable to compute kruskal-wallis with feature: " + fname + "\n")
    return adjust_stats( stats )

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args( )
    fnames, fvalues = [], []
    adding = False
    with open( args.input , "rt") as fh:
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
    fh = open( args.output, "w" ) if args.output is not None else sys.stdout
    if args.focal_type == "continuous":
        mvalues = list( map( float, mvalues ) )
        stats = spearman_analysis( mvalues, fnames, fvalues )
        fh.write("# spearman analysis of metadatum: " + mname + "\n")
        fh.write("# feature\trho\tp-value\tq-value\n")
    elif args.focal_type == "categorical":
        stats = kruskalwallis_analysis( mvalues, fnames, fvalues )
        fh.write("# kruskal-wallis analysis of metadatum: " + mname + "\n")
        fh.write("# feature\tlevel means\tp-value\tq-value\n")
    for stat in stats:
        stat[-1] = "%.4g" % stat[-1]
        stat[-2] = "%.4g" % stat[-2]
        fh.write("\t".join( list(map( str, stat )))+"\n")
            
if __name__ == "__main__":
    main()
