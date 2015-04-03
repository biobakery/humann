#! /usr/bin/env python

"""
HUMAnN2 utility for making strain profiles
Run ./strain_profiles.py -h for usage help
"""

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import argparse
import sys
import csv
import util

# constants
c_strain_profile_extension = "-strain_profile.tsv"
c_tax_delim = "."
c_epsilon = 1e-10
c_forbidden = ["unclassified"]

def get_args ():
    """ Get args from Argparse """
    parser = argparse.ArgumentParser()
    parser.add_argument( 
        "-i", "--input", 
        default=None,
        help="Original output table (.tsv format); default=[STDIN]",
        )
    parser.add_argument( 
        "-m", "--critical_mean", 
        type=float,
        default=10.0,
        help="Default mean non-zero gene abundance for inclusion",
        )
    parser.add_argument( 
        "-n", "--critical_count", 
        type=int,
        default=100,
        help="Default non-zero number of genes for inclusion",
        )
    parser.add_argument( 
        "-f", "--critical_frequency",
        type=float,
        default=0.0,
        help="Threshold prevalence (or absence) for defining an interesting gene",
        )
    parser.add_argument( 
        "-s", "--critical_samples",
        type=int,
        default=2,
        help="Threshold number of samples with strain for printing",
        )
    args = parser.parse_args()
    return args

class Partition ( ):
    def __init__ ( self, name=None ):
        self.name = name
        self.rows = {}
        self.cols = {}
    def add_row( self, i ):
        self.rows[i] = 1
    def add_col( self, j ):
        self.cols[j] = 1
    def del_row( self, i ):
        del self.rows[i]    
    def del_col( self, j ):
        del self.cols[j]
    def get_rows( self ):
        return sorted( self.rows )
    def get_cols( self ):
        return sorted( self.cols )

def partition_table( table, m, n, p ):
    partitions = {}
    for i, rowhead in enumerate( table.rowheads ):
        table.data[i] = map( float, table.data[i] )
        if util.c_strat_delim in rowhead:
            gene, species = rowhead.split( util.c_strat_delim )
            # disallow the unclassified stratum
            if species not in c_forbidden:
                species = species.split( c_tax_delim )[1]
                if species not in partitions:
                    partitions[species] = Partition( species )
                partitions[species].add_row( i )
    # limit partitions to subjects with mean non-zero above threshold
    for name, partition in partitions.items():
        partition.cols = {k:1 for k in range( len( table.colheads ) )}
        for j in partition.get_cols():
            values = [table.data[i][j] for i in partition.get_rows()]
            nonzero = [k for k in values if k > 0]
            if len( nonzero ) < n or sum( nonzero ) / len( nonzero ) < m:
                partition.del_col( j )
    # limit partitions to interesting features
    for name, partition in partitions.items():
        for i in partition.get_rows():
            values = [table.data[i][j] for j in partition.get_cols()]
            nonzero = [k for k in values if k > 0]
            prevalence = len( nonzero ) / ( c_epsilon + float( len( values ) ) )
            # print( name, i, values, nonzero, prevalence, partition.name )
            if prevalence < p or prevalence > 1 - p:
                partition.del_row( i )
    return partitions

def write_partition ( table, partition, outfile ):
    matrix = [["STRAINS"] + [table.colheads[j] for j in partition.get_cols()]]
    for i in partition.get_rows():
        row = [table.rowheads[i]]
        row += [table.data[i][j] for j in partition.get_cols()]
        matrix.append( row )
    with open( outfile, "w" ) as fh:
        writer = csv.writer( fh, dialect='excel-tab' )
        for row in matrix:
            writer.writerow( row )

def main ( ):
    args = get_args()
    table = util.Table( args.input )
    partitions = partition_table( 
        table,
        args.critical_mean,
        args.critical_count,
        args.critical_frequency,
        )
    for name, partition in partitions.items():
        if len( partition.get_cols() ) >= args.critical_samples and \
                len( partition.get_rows() ) >= 1:
            write_partition( table, partition, name+c_strain_profile_extension )

if __name__ == "__main__":
    main()
