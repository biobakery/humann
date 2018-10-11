#! /usr/bin/env python

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import argparse
import sys
import csv

from humann2.tools import util

description = """
HUMAnN2 utility for making strain profiles
==========================================
Based on the principle of detecting variable 
presence and absence of gene families within a species
that is otherwise well-covered in multiple samples.
"""

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_strain_profile_extension = "-strain_profile.tsv"
c_tax_delim = "."
c_epsilon = 1e-10
c_forbidden = ["unclassified"]

# ---------------------------------------------------------------
# utilities 
# ---------------------------------------------------------------

def get_args ():
    """ Get args from Argparse """
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter
        )
    parser.add_argument( 
        "-i", "--input", 
        default=None,
        help="Original output table (tsv or biom format); default=[TSV/STDIN]",
        )
    parser.add_argument( 
        "-m", "--critical_mean", 
        type=float,
        default=10.0,
        help="Default mean non-zero gene abundance for inclusion; default=10.0",
        )
    parser.add_argument( 
        "-n", "--critical_count", 
        type=int,
        default=500,
        help="Default non-zero number of genes for inclusion; default=500",
        )
    parser.add_argument( 
        "-p", "--pinterval",
        type=float,
        default=[c_epsilon, 1.0],
        nargs=2,
        help="Only genes with prevalence in this interval are allowed; default=[1e-10, 1]",
        )
    parser.add_argument( 
        "-s", "--critical_samples",
        type=int,
        default=2,
        help="Threshold number of samples having strain; default=2",
        )
    parser.add_argument( 
        "-l", "--limit",
        default=None,
        help="Limit output to species matching a particular pattern, e.g. 'Streptococcus'; default=OFF",
        )
    args = parser.parse_args()
    return args

class Partition ( ):
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
            prevalence = len( nonzero ) / ( c_epsilon + float( len( values ) ) )
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
