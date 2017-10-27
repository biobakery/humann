#! /usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import argparse
import sys
import csv

try:
    from humann2.tools import util
    from humann2.tools.humann2_table import Table
except ImportError:
    sys.exit( "CRITICAL ERROR: Unable to find the HUMAnN2 python package.\n" +
              "Please check your install." )

description = util.wrap( """
HUMAnN2 utility for building strain profiles

For well-covered species, slice out profiles of gene presence/absence 
as a strain profile across samples. NOTE: Preferred input is a HUMAnN2 gene
families file containing UNNORMALIZED RPK units (for coverage inference).
""" )

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_strain_profile_extension = "-strain_profile.tsv"
c_tax_delim = "."
c_forbidden = ["unclassified"]
c_formats_help = """(You may list more than one.)
abunds: output the strain profiles in the original abundance units
binary: output the strain profiles in 1/0 (presence/absence) units
hybrid: output binary/original value pairs"""

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
        "-o", "--output-directory", 
        metavar="<path>",
        default=".",
        help="Directory to write strain profiles\n[Default=.]",
        )
    parser.add_argument( 
        "-g", "--minimum-nonzero-genes", 
        metavar="<int>",
        type=int,
        default=500,
        help=("To be considered, a species must recruit reads to this "
              "many genes in a sample\n[Default=500]"),
        )
    parser.add_argument( 
        "-c", "--median-coverage", 
        type=float,
        metavar="<float>",
        default=20.0,
        help=("To be considered, non-zero genes must meet this median abundance\n"
              "[Default=20.0; assuming RPK units from 100 nt reads, 20 RPK ~ 2x coverage]"),
        )
    parser.add_argument( 
        "-s", "--minimum-samples",
        type=int,
        metavar="<int>",
        default=2,
        help=("Only write strain profiles for strains detected in at least this many samples\n"
              "[Default=2]"),
        )
    parser.add_argument( 
        "-f", "--output-formats",
        metavar="<abunds/hybrid/binary>",
        #nargs="+",
        default=["binary"],
        choices=["abunds", "hybrid", "binary"],
        help=c_formats_help,
        )
    args = parser.parse_args( )
    return args

class Partition( ):
    def __init__ ( self, name=None ):
        self.name = name
        self.rows = {}
        self.cols = {}
    def add_rows( self, *args ):
        [self.rows.update( [[i, 1]] ) for i in args]
    def add_cols( self, *args ):
        [self.cols.update( [[j, 1]] ) for j in args]
    def del_rows( self, *args ):
        [self.rows.pop( i ) for i in args]    
    def del_cols( self, *args ):
        [self.cols.pop( j ) for j in args]
    def get_rows( self ):
        return sorted( self.rows )
    def get_cols( self ):
        return sorted( self.cols )

def partition_table( table, m, n, pinterval ):
    # first build the partitions
    partitions = {}
    for i, rowhead in enumerate( table.rowheads ):
        table.data[i] = list(map( float, table.data[i] ))
        if util.c_strat_delim in rowhead:
            gene, species = rowhead.split( util.c_strat_delim )
            # disallow the unclassified stratum
            if species not in c_forbidden:
                species = species.split( c_tax_delim )[1]
                if species not in partitions:
                    partitions[species] = Partition( species )
                partitions[species].add_rows( i )
    # limit partitions to subjects with mean non-zero gene abund above threshold
    for name, partition in partitions.items():
        partition.add_cols( *range( len( table.colheads ) ) )
        for j in partition.get_cols():
            values = [table.data[i][j] for i in partition.get_rows()]
            nonzero = [k for k in values if k > 0]
            if len( nonzero ) < n or sum( nonzero ) / float( len( nonzero ) ) < m:
                partition.del_cols( j )
    # limit partitions to interesting features (if user changed f cutoff)
    for name, partition in partitions.items():
        for i in partition.get_rows():
            values = [table.data[i][j] for j in partition.get_cols()]
            nonzero = [k for k in values if k > 0]
            prevalence = len( nonzero ) / ( util.c_eps + float( len( values ) ) )
            # print( name, i, values, nonzero, prevalence, partition.name )
            if not pinterval[0] <= prevalence <= pinterval[1]:
                partition.del_rows( i )
    return partitions

def write_partition ( table, partition, outfile ):
    matrix = [["HEADERS"] + [table.colheads[j] for j in partition.get_cols()]]
    for i in partition.get_rows():
        row = [table.rowheads[i]]
        row += [table.data[i][j] for j in partition.get_cols()]
        matrix.append( row )
    with open( outfile, "w" ) as fh:
        writer = csv.writer( fh, dialect='excel-tab' )
        for row in matrix:
            writer.writerow( row )

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main ( ):
    args = get_args()
    table = util.Table( args.input )
    partitions = partition_table( 
        table,
        args.critical_mean,
        args.critical_count,
        args.pinterval,
        )
    for name, partition in partitions.items():
        if len( partition.get_cols() ) >= args.critical_samples and \
                len( partition.get_rows() ) >= 1 and \
                ( args.limit is None or args.limit in name ):
            write_partition( table, partition, name+c_strain_profile_extension )
            
if __name__ == "__main__":
    main()
