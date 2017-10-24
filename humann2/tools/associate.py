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
              " Please check your install." )

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
serious statistical analyses. Given a HUMAnN2 table with metadata rows
at the top, perform repeated, non-parametric comparison between a focal
metadata row and each feature total in the table. A continuous metadata feature
will be compared by the Spearman correlation, and a categorical metadata feature
will be compared by the Kruskal-Wallis test. Resulting p-values are subjected
to Benjamini-Hochberg FDR correction.
""" )

def get_args( ):
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter,
        )
    util.attach_common_arguments( parser, no_output=True )
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
                         default=None,
                         help="Where to write the statistical results\n[Default=STDOUT]", )
    parser.add_argument( "-f", "--fdr",
                         metavar="<float>",
                         type=float,
                         default=0.2,
                         help="Benjamini-Hochberg FDR threshold\n[Default=0.2]", )
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
    qvalues.reverse( )
    for i in range( 1, n ):
        if qvalues[i] > qvalues[i-1]:
            qvalues[i] = qvalues[i-1]
    qvalues.reverse( )
    # rebuild qvalues in the original order
    ordered_qvalues = [None for q in qvalues]
    for i, q in enumerate( qvalues ):
        ordered_qvalues[index[i]] = q
    return ordered_qvalues

def fdr_adjust( stats ):
    H = [h for h in stats]
    P = [stats[h].get( "p-value", 1.0 ) for h in stats]
    Q = pvalues2qvalues( P )
    for h, q in zip( H, Q ):
        stats[h]["q-value"] = q
    return None

def spearman_analysis( T, mvals ):
    mvals = metafloat( mvals )
    stats = {}
    for f in util.fsort( T.data ):
        fcode, fname, fstrat = util.fsplit( f )
        if fstrat is None:
            fname = util.fjoin( fcode, fname )
            try:
                inner = {}
                rho, p = spearmanr( T.data[f], mvals )
                inner["effect_size"] = rho
                inner["p-value"] = p
                stats[fname] = inner
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
            try:
                inner = {}
                lists = shatter( mvals, T.data[f] )
                summary = {k:"%.4g" % ( np.mean( v ) ) for k, v in lists.items( )}
                summary = [":".join( [k, summary[k]] ) for k in sorted( summary )]
                summary = "; ".join( summary )
                hstat, p = kruskalwallis( *lists.values( ) )
                inner["effect_size"] = hstat
                inner["p-value"] = p
                inner["level_means"] = summary
                stats[fname] = inner
            except:
                print( "WARNING: Unable to compute Kruskal-Wallis with feature:", fname, file=sys.stderr )
    return stats

def reformat( stats ):
    for f in stats:
        for k in ["effect_size", "p-value", "q-value"]:
            stats[f][k] = "%.4g" % ( stats[f][k] )
    return None

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
        "effect_size",
        "p-value",
        "q-value",
        ]    
    # perform relevant analysis
    if args.focal_type == "continuous":
        print( "Performing Spearman analysis vs. metadatum:", mname, file=sys.stderr )
        stats = spearman_analysis( T, mvals )
    elif args.focal_type == "categorical":
        print( "Performing Kruskal-Wallis analysis vs. metadatum:", mname, file=sys.stderr )
        stats = kruskalwallis_analysis( T, mvals )
        headers.append( "level_means" )
    # add qvalues / format
    fdr_adjust( stats )
    order = [(stats[f].get( "q-value", 1.0 ), f) for f in stats] 
    order = [(q, f) for q, f in order if q <= args.fdr]
    order.sort( )
    reformat( stats )
    # write results
    fh = open( args.output, "w" ) if args.output is not None else sys.stdout
    writer = csv.writer( fh, csv.excel_tab )
    writer.writerow( [k.upper( ) for k in headers] )
    for q, f in order:
        row = [f] + [stats[f].get( h, "#N/A" ) for h in headers[1:]]
        writer.writerow( row )
    fh.close( )
    # wrap up
    if len( order ) == 0:
        print( "No FDR significant associations found.", file=sys.stderr )
    print( "Finished successfully.", file=sys.stderr )
        
if __name__ == "__main__":
    main( )
