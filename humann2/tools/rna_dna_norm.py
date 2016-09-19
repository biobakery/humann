#! /usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
from math import log
import argparse
import sys

from humann2.tools import util

description = """
HUMAnN2 utility for normalizing combined meta'omic sequencing data
==================================================================
Given a DNA table and a RNA table, produce smoothed RNA and DNA 
values as well as relative expression values. "Smoothing" means
substituting a small value in place of a zero or missing value.
The default method used is "Laplace" (pseudocount) scaling, where
the pseudocount is the sample-specific minimum non-zero value.
(Witten-Bell smoothing is also implemented.)

-- The DNA and RNA columns must be 1:1 and in the same order.

-- If working with stratified data, smoothing is carried out on the
stratified values and then community totals are recomputed.
"""

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_new_dna_extension = "-smoothed_dna.tsv"
c_new_rna_extension = "-smoothed_rna.tsv"
c_norm_rna_extension = "-relative_expression.tsv"

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
        "-d", "--input_dna", 
        help="Original DNA output table (tsv or biom format)",
        )
    parser.add_argument( 
        "-r", "--input_rna", 
        help="Original RNA output table (tsv or biom format)",
        )
    parser.add_argument( 
        "-o", "--output_basename", 
        default = "results",
        help="Path/basename for the three output tables; DEFAULT=results",
        )
    parser.add_argument( 
        "-m", "--method",
        choices=["laplace", "witten_bell"],
        default = "laplace",
        help="Choice of smoothing method; DEFAULT=laplace",
        )
    parser.add_argument( 
        "-l", "--log_transform", 
        action="store_true",
        help="Report log-transformed relative expression values",
        )
    parser.add_argument( 
        "-b", "--log_base", 
        type=float,
        default=2.0,
        help="Base for log transformation (if requested); DEFAULT=2."
        )
    args = parser.parse_args()
    return args

def remove_totals( table ):
    rowheads2, data2 = [], []
    for i, rowhead in enumerate( table.rowheads ):
        if util.c_strat_delim in rowhead:
            rowheads2.append( rowhead )
            data2.append( table.data[i] )
    table.rowheads, table.data = rowheads2, data2

def laplace( table, all_features ):
    alphas = [None for i in table.data[0]]
    for i, row in enumerate( table.data ):
        # float table here
        table.data[i] = list(map( float, row ))
        for j, value in enumerate( table.data[i] ):
            if value > 0:
                if alphas[j] is None or value < alphas[j]:
                    alphas[j] = value
    # compute quick index for table
    rowmap = {rowhead:i for i, rowhead in enumerate( table.rowheads )}
    # rebuild table data
    rowheads2, data2 = [], []
    for feature in all_features:
        rowheads2.append( feature )
        # feature is in the table; still adjust zero values
        if feature in rowmap:
            # augment all feature values by their sample alphas
            i = rowmap[feature]
            table.data[i] = [k1 + k2 for k1, k2 in zip( table.data[i], alphas )]
            data2.append( table.data[i] )
        # feature is absent from the table; use alphas for everyone
        else:
            data2.append( alphas )
    table.rowheads, table.data = rowheads2, data2
    # compute and attach colsums (accounting for new mass)
    colsums = [0 for i in table.data[0]]
    for row in table.data:
        colsums = [k1 + k2 for k1, k2 in zip( colsums, row )]
    table.colsums = colsums

def witten_bell( table, all_features ):
    nonzero = [0 for i in table.data[0]]
    colsums = [0 for i in table.data[0]]
    for i, row in enumerate( table.data ):
        # float table here
        table.data[i] = list(map( float, row ))
        for j, value in enumerate( table.data[i] ):
            nonzero[j] += 1 if value > 0 else 0
            colsums[j] += value
    # save colsums in table object
    table.colsums = colsums
    # compute epsilons/norms
    epsilons, norms = [], []
    for j in range( len( nonzero ) ):
        # total events for column j
        N = colsums[j]
        # total first events
        T = nonzero[j]
        # implied unobserved events
        Z = len( all_features ) - T
        norms.append( N / float( N + T ) )
        epsilons.append( ( norms[-1] * T / float( Z ) ) if Z > 0 else 0 )
    # compute quick index for table
    rowmap = {rowhead:i for i, rowhead in enumerate( table.rowheads )}
    # rebuild table data
    rowheads2, data2 = [], []
    for feature in all_features:
        rowheads2.append( feature )
        # feature is in the table; still adjust zero values
        if feature in rowmap:
            i = rowmap[feature]
            for j, value in enumerate( table.data[i] ):
                if value == 0:
                    table.data[i][j] = epsilons[j]
                else:
                    table.data[i][j] = value * norms[j]
            data2.append( table.data[i] )
        # feature is absent from the table; use epsilon for everyone
        else:
            data2.append( epsilons )
    table.rowheads, table.data = rowheads2, data2

def hsum( table ):
    # look ahead
    groups = {}
    for i, rowhead in enumerate( table.rowheads ):
        groups.setdefault( rowhead.split( util.c_strat_delim )[0], [] ).append( i )
    # rebuild
    rowheads2, data2 = [], []
    for group, poslist in groups.items():
        # attach the sum to the new table
        rowheads2.append( group )
        total = [0 for i in range( len( table.data[0] ) )]
        for i in poslist:
            total = [k1 + k2 for k1, k2 in zip( total, table.data[i] )]
        data2.append( total )
        # attach individual old rows to the new table
        for i in poslist:
            rowheads2.append( table.rowheads[i] )
            data2.append( table.data[i] )
    table.rowheads, table.data = rowheads2, data2

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main ( ):
    args = get_args()
    dna = util.Table( args.input_dna )
    rna = util.Table( args.input_rna )
    method = {"laplace":laplace, "witten_bell":witten_bell}[args.method]
    assert dna.is_stratified == rna.is_stratified, \
        "FAILED: Tables have nonequal stratification status."
    strat_mode = dna.is_stratified
    all_features = sorted( set( dna.rowheads ).__or__( set( rna.rowheads ) ) )
    for t in dna, rna:
        if strat_mode:
            remove_totals( t )
            all_features = [k for k in all_features if util.c_strat_delim in k]
        method( t, all_features )
        if strat_mode:
            hsum( t )
    # write out dna/rna
    dna.write( args.output_basename+c_new_dna_extension, unfloat=True )
    rna.write( args.output_basename+c_new_rna_extension, unfloat=True )
    # normalize rna by dna (account for seq depth [scale]), then write
    scale = [d / r for r, d in zip( rna.colsums, dna.colsums )]
    for i in range( len( dna.data ) ):
        rna.data[i] = [s * r / d for s, r, d in zip( scale, rna.data[i], dna.data[i] )]
        if args.log_transform:
            divisor = log( args.log_base )
            rna.data[i] = list(map( lambda x: log( x ) / divisor, rna.data[i] ))
    rna.write( args.output_basename+c_norm_rna_extension, unfloat=True )

if __name__ == "__main__":
    main()
