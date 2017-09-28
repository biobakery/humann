#!/usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import os
import sys
import argparse
import re

try:
    #from humann2.tools import util
    #from humann2.tools.better_table import Table
    import util
    from better_table import Table
except ImportError:
    sys.exit( "CRITICAL ERROR: Unable to find the HUMAnN2 python package.\n" +
              "Please check your install." )  

try:
    import numpy as np
    import scipy.spatial.distance as spd
except ImportError:
    sys.exit( "CRITICAL ERROR: This script requires the python scientific stack (numpy + scipy)" )

description = """
HUMAnN2 utility for calculating contributional diversity
========================================================

Computes ecological diversity statistics for individual 
functions in a stratified HUMAnN2 profile based on their 
per-species contributions. The analysis is restricted to 
functions that were "well-explained" in the input (i.e. 
a majority of copies were attributed to species in the 
majority of samples).

Diversity is calculated over samples where the feature 
was well-explained and the 'unclassified' stratum is 
excluded. Alpha diversity is calculated with the 
<Gini-Simpson> index. Beta diversity is calculated with 
the <Bray-Curtis> index.
"""

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_eps = 1e-10

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
        "-E", "--min-percent-explained",
        default=75.0,
        type=float,
        metavar="<value in 0-100>",
        help=("Only consider features where species explain >=E%% of copies in >=E%% of samples\n"
              "[Default=75.0]"),
        )
    parser.add_argument( 
        "-A", "--min-abundance",
        default=0.0,
        metavar="<float>",
        help=("Also require total feature abundance >A per sample\n"
              "[Default=0.0, i.e. only exclude samples with zero abundance]")
        )
    args = parser.parse_args( )
    return args

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

def adiv_fast( samples ):
    """ gini-simpson on 2d array with samples as first axis """
    return np.mean( 1 - np.sum( samples**2, axis=1 ) )

def bdiv_fast( samples ):
    """ bray-curtis on 2d array with samples as first axis """
    # pdist returns the non-redundant, non-self distances
    return np.mean( spd.pdist( samples, metric="braycurtis" ) )

def contributional_diversity( table, min_abund=None, min_explained=None ):
    min_explained_frac = min_explained / 100.0
    # total stratified abundance of each function
    total = {}
    # total abundance attributed to species
    explained = {}
    # per-feature abundance vectors
    stacks = {}
    # populate
    for f in table.data:
        fbase, fname, stratum = util.fsplit( f )
        # don't lose track of feature names if present
        fname = util.fjoin( fbase, name=fname )
        if stratum is not None:
            if fname not in total:
                total[fname] = table.zeros( )
                explained[fname] = table.zeros( )
            total[fname] += table.data[f]
            if stratum != util.c_unclassified:
                explained[fname] += table.data[f]
                stacks.setdefault( fname, [] ).append( table.data[f] )
    # determine allowed functions
    allowed = {}
    for f in total:
        my_tot = total[f]
        my_exp = explained[f]
        n = len( my_tot )
        # index of samples we'll use for this function
        index = []
        for i in range( n ):
            if my_tot[i] > min_abund:
                if my_exp[i] / my_tot[i] >= min_explained_frac:
                    index.append( i )
        if len( index ) / float( n ) >= min_explained_frac:
            allowed[f] = index
    # compute stats over passing samples (in index)
    fstats = {}
    for f, index in allowed.items( ):
        inner = fstats.setdefault( f, {} )
        # numpy list slice
        my_tot = total[f][index]
        my_exp = explained[f][index]
        # basic stats (all samples)
        inner["0: Samples used"] = len( index )
        inner["1: Samples used (frac)"] = len( index ) / float( len( total[f] ) )
        inner["2: Mean abundance"] = np.mean( my_tot )
        my_frac = my_exp / my_tot
        inner["3: Mean explained (frac)"] = np.mean( my_frac )
        # stack of per-species contributions normalized to _their_ total
        stack = [row[index] / my_exp for row in stacks[f]]
        # convert to 2d array with samples as first axis
        samples = np.vstack( stack ).transpose( )
        # compute diversity (alpha=gini-simpson / beta=bray-curtis)
        inner["4: Contrib div, alpha"] = adiv_fast( samples )
        inner["5: Contrib div, beta"] = bdiv_fast( samples )
    # report
    print( "Contributional diversity report:", file=sys.stderr )
    tot = len( total )
    print( "  Features considered: {:,}".format( tot ), file=sys.stderr )
    tot = len( fstats )
    frac = 100 * len( fstats ) / float( len( total ) )
    print( "  Required feature abundance >{} and >={}% of copies explained in >={}% of samples".format(
            min_abund, min_explained, min_explained ), file=sys.stderr )
    print( "  Features deemed appropriate for analysis: {:,} ({:.1f}%)".format( tot, frac ), file=sys.stderr )
    # nested dict of per-function stats
    return fstats
        
# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args( )
    table = Table( args.input, last_metadata=args.last_metadata )
    # compute diversity stats
    fstats = contributional_diversity( 
        table, 
        min_abund=args.min_abundance, 
        min_explained=args.min_percent_explained,
        )
    # make output table
    headers = sorted( {v for k in fstats for v in fstats[k]} )
    data = {}
    for f in fstats:
        temp = []
        for h in headers:
            temp.append( fstats[f][h] )
        data[f] = np.array( temp )
    fstats = Table( data, headers=headers )
    fstats.anchor = "# Func \ Stat"
    fstats.write( args.output, unfloat=True )

if __name__ == "__main__":
    main( )
