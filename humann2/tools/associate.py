#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import csv
import argparse

try:
    from humann2 import config
    from humann2.tools import util
    from humann2.tools.humann2_table import Table
except ImportError:
    sys.exit( "CRITICAL ERROR: Unable to find the HUMAnN2 python package." +
              " Please check your install.")

try:
    import numpy as np
    from scipy.stats import spearmanr
    from scipy.stats.mstats import kruskalwallis
except ImportError:
    sys.exit( "This script requires the Python scientific stack: numpy and scipy." )

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------
   
c_allowed_impute = 0.2

# ---------------------------------------------------------------
# argument parsing
# ---------------------------------------------------------------

description = util.wrap( """
HUMAnN2 utility for associating functions with metadata

This script is part of the HUMAnN2 demo and is not intended for
serious statistical analyses. Table must have metadata rows at the
top. One continuous or categorical metadata row can be non-parametrically
associated (Spearman/Kruskal-Wallis) with the function totals in 
the input HUMAnN2 table.
""" )

def get_args( ):
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter,
        )
    utils.attach_common_arguments( no_output=True )
    parser.add_argument( "-m", "--focal-metadatum",
                         metavar="<str>",
                         required=True,
                         help="Indicate metadatum to test vs. community feature totals", )
    parser.add_argument( "-t", "--focal-type",
                         required=True,
                         choices=["continuous", "categorical"],
                         metavar="<continuous/categorical>",
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

def metafloat( values ):
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
        sys.exit( "FAILED: Not enough numeric entries in metadata row." )

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

def fdr_correct( stats ):
    H = [h for h in stats]
    P = [stats[h].get( "P-value", 1 ) for h in stats]
    Q = pvalues2qvalues( P )
    for h, q in zip( H, Q ):
        stats[h]["q-value"] = q
    return None

def adjust_stats( stats ):
    pvalues = [stat[-1] for stat in stats]
    qvalues = pvalues2qvalues( pvalues )
    for i in range( len( stats ) ):
        stats[i].append( qvalues[i] )
    return sorted( stats, key=lambda stat: stat[-1] )

def spearman_analysis( T, mvals ):
    print( "Performing Spearman analysis vs. metadatum:" + mname, file=sys.stderr )
    mvals = metafloat( mvals )
    stats = {}
    for f in util.fsort( T.data ):
        fcode, fname, fstrat = util.fsplit( f )
        if fstrat is None:
            fname = util.fjoin( fcode, fname )
            inner = stats.setdefault( fname, {} )
            try:
                rho, p = spearmanr( T.data[f], mvals )
                inner["effect size"] = rho
                inner["p-value"] = p
            except:
                print( "WARNING: Unable to compute Spearman with feature:", fname, file=sys.stderr )
    return stats

def shatter( cats, values ):
    lists = {}
    for c, v in zip( cats, values ):
        lists.setdefault( c, [] ).append( v )
    return lists

def kruskalwallis_analysis( T, mvals ):
    if len( set( mvals ) ) > np.sqrt( len( mvals ) ):
        print( "WARNING: Categorical metadata row has many distinct values.", file=sys.stderr )
    stats = {}
    for f in util.fsort( T.data ):
        fcode, fname, fstrat = util.fsplit( f )
        if fstrat is None:
            fname = util.fjoin( fcode, fname )
            inner = stats.setdefault( fname, {} )
            try:
                lists = shatter( mvals, T.data[f] )
                summary = {k:"%.4g" % ( np.mean( v ) ) for k, v in lists.items( )}
                summary = [":".join( [k, v] ) for k, v in summary.items( )]
                summary = "; ".join( summary )
                hstat, p = kruskalwallis( *lists.values( ) )
                inner["effect size"] = hstat
                inner["p-value"] = p
                inner["level means"] = summary
            except:
                print( "WARNING: Unable to compute Kruskal-Wallis with feature:", fname, file=sys.stderr )
    return stats

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args( )
    T = Table( args.input, last_metadata=args.last_metadata )
    # grab/test metadata
    mname = args.focal_metadatum
    if mname not in T.metadata:
        sys.exit( "FAILED: Requested metadata row <{}> not found.".format( mname ) )
    mvals = T.metadata[mname]
    # headers
    headers = [
        "# feature",
        "effect size",
        "p-value",
        "q-value",
        ]    
    # perform relevant analysis
    if args.focal_type == "continuous":
        stats = spearman_analysis( T, mvals )
    elif args.focal_type == "categorical":
        stats = kruskalwallis_analysis( T, mvals )
        headers.append( "level means" )
    # add qvalues
    fdr_correct( stats )
    # write results
    fh = open( args.output, "w" ) if args.output is not None else sys.stdout
    writer = csv.writer( fh, csv.excel_tab )
    writer.writerow( headers )
    FOUND_SOMETHING = False
    for f in sorted( stats ):
        if stats[f].get( "Q-value" ) <= args.fdr:
            FOUND_SOMETHING = True
            row = [f] + [stats[f].get( h, "#N/A" ) for h in headers[1:]]
            writer.writerow( row )
    if not FOUND_SOMETHING:
        print( "No FDR significant associations found.", file=sys.stderr )
    fh.close( )
    print( "Finished successfully.", file=sys.stderr )
        
if __name__ == "__main__":
    main( )
