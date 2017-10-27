#! /usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import sys
import argparse

try:
    from humann2 import config
    from humann2.tools import util
    from humann2.tools.humann2_table import Table
except ImportError:
    sys.exit( "CRITICAL ERROR: Unable to find the HUMAnN2 python package.\n" +
              "Please check your install." )

try:
    import numpy as np
except ImportError:
    sys.exit( "CRITICAL ERROR: This script requires the python scientific stack (e.g. numpy)" )

description = util.wrap( """
HUMAnN2 utility for normalizing combined meta'omic sequencing data

Given HUMAnN2 output for metatranscriptomes (mtx) and metagenomes (mgx)
from the same biosamples, produce a new table of "relative expression"
values by normalizing mtx by their mgx copy number. Normalization can
be by log2-ratio or difference. When using ratios, zero values are 
additively smoothed.
""" )

# ---------------------------------------------------------------
# command-line interface
# ---------------------------------------------------------------

def get_args( ):
    """ Get args from Argparse """
    parser = argparse.ArgumentParser(      
        description=description, 
        formatter_class=argparse.RawTextHelpFormatter,
        )
    parser.add_argument( 
        "-d", "--dna-table", 
        metavar="<path>",
        required=True,
        help="Input metagenomic feature table",
        )
    parser.add_argument( 
        "-r", "--rna-table", 
        metavar="<path>",
        required=True,
        help="Input metatranscriptomic feature table",
        )
    parser.add_argument( 
        "-o", "--output",
        default="relative_expression.tsv",
        metavar="<path>",
        help="Where to write relative expression values\n[Default=relative_expression.tsv]",
        )
    parser.add_argument(
        "-s", "--sample-map",
        metavar="<path>",
        help=("Two columns pairing the DNA samples with RNA samples"
              "If not provided, program will assume that sample N/column N+1"
              "from the two tables should be paired"),
        )
    parser.add_argument( 
        "-m", "--mode",
        choices=["log2ratio", "difference"],
        default= "log2ratio",
        metavar="<log2ratio/difference>",
        help="Method for comparing RNA and DNA values\n[Default=log2ratio]",
        )
    parser.add_argument( 
        "-t", "--threshold",
        default=1e-6,
        metavar="<float>",
        help=("If the relative abundance of a feature is <t in both the RNA and DNA\n"
              "set its log2 ratio to 0.0 (only applies to log2 ratio mode)\n"
              "[Default=1e-6]"),
        )
    args = parser.parse_args( )
    return args

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

def smooth( table ):
    return None

def hsum( table ):
    return None

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args( )
    dna = util.Table( args.dna )
    rna = util.Table( args.rna )
    # normalize rna by dna (account for seq depth [scale]), then write
    scale = [d / r for r, d in zip( rna.colsums, dna.colsums )]
    for i in range( len( dna.data ) ):
        rna.data[i] = [s * r / d for s, r, d in zip( scale, rna.data[i], dna.data[i] )]
        if args.log_transform:
            divisor = log( args.log_base )
            rna.data[i] = list(map( lambda x: log( x ) / divisor, rna.data[i] ))
    rna.write( args.output_basename+c_norm_rna_extension, unfloat=True )

if __name__ == "__main__":
    main( )
