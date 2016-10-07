#!/usr/bin/env python

from __future__ import print_function # Python 2.7+ required
import os
import sys
import csv
import re
import argparse

try:
    from humann2 import config
    from humann2.tools import util
except ImportError:
    sys.exit("CRITICAL ERROR: Unable to find the HUMAnN2 python package." +
        " Please check your install.")

description = """
HUMAnN2 utility for inferring "unclassified" taxonomy
=====================================================
Based on the lowest common ancestor (LCA) annotation
of each UniRef50/90 cluster, infer approximate taxonomy 
for unclassified features at a target level of resolution 
(default=Family). Will modify features of known genus/species 
to match target level."""

# ---------------------------------------------------------------
# utility mapping interface
# ---------------------------------------------------------------

# get a list of all available utility mapping files
try:
    all_mapping_files=os.listdir( config.utility_mapping_database )
except EnvironmentError:
    all_mapping_files=[]

databases = {
    "uniref50": "uniref50-tol-lca.dat.gz",
    "uniref90": "uniref90-tol-lca.dat.gz",
    }

missing = False
for key, value in databases.items( ):
    if value not in all_mapping_files:
        missing = True
    else:
        databases[key] = os.path.join( config.utility_mapping_database, value )
if missing:
    sys.exit( """
This tool requires the HUMAnN2 utility data files.
To add these to your installation, please execute:

$ humann2_databases --download utility_mapping full $DIR

Replacing, $DIR with the directory to install the databases.\n""" )

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_levels = [
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
]
c_tmode = "totals"
c_umode = "unclassified"
c_smode = "stratified"
c_tol_header = "# TOL"
c_lca_header = "# LCA"
c_unclassified = "unclassified"
c_bypass = "AmbiguousBypass"
c_na = "-"

# ---------------------------------------------------------------
# helper objects
# ---------------------------------------------------------------

class Taxon:
    def __init__( self, name, common, rank, pname, status ):
        self.name = name
        self.common = common
        self.rank = rank
        self.pname = pname
        self.status = status

class TreeOfLife:
    def __init__( self ):
        self.nodes = {}
        self.connect = {}
    def attach( self, node ):
        self.nodes[node.name] = node
        if node.status != c_bypass:
            self.connect[node.common] = node.name
    def get_lineage( self, common ):
        lineage = []
        if common in self.connect:
            name = self.connect[common]
            while name in self.nodes and name != c_na:
                node = self.nodes[name]
                lineage.append( [node.rank, node.common] )
                name = node.pname
        return lineage

# ---------------------------------------------------------------
# argument parsing
# ---------------------------------------------------------------

def get_args( ):
    """ Get args from Argparse """
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter
        )
    parser.add_argument( "-i", "--input", 
                         required=True,
                         help="HUMAnN2 output table" )
    parser.add_argument( "-o", "--output", 
                         default=None,
                         help="Destination for modified table; default=STDOUT" )
    parser.add_argument( "-l", "--level",
                         choices=c_levels,
                         default="Family",
                         help="Desired level for taxonomic estimation/summation; default=Family" )
    parser.add_argument( "-r", "--resolution",
                         choices=databases.keys( ),
                         required=True,
                         help="Resolution of the input file." )
    parser.add_argument( "-m", "--mode", 
                         choices=[c_tmode, c_umode, c_smode],
                         default=c_tmode,
                         help="Which rows to include in the estimation/summation; default=totals" )
    parser.add_argument( "-t", "--threshold", 
                         type=float, 
                         default=1e-3, 
                         help="Minimum frequency for a new taxon to be included; default=1e-3" )
    parser.add_argument( "-d", "--dev",
                         help="Manually specify a development database" )
    args = parser.parse_args()
    return args

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

def simplify( name ):
    return re.sub( "[^A-Za-z0-9]+", "_", name )

def build_taxmap( features, target_rank, p_datafile ):
    unirefs = {util.fsplit( k )[0] for k in features}
    unirefs = {k for k in unirefs if "UniRef" in k}
    # load tree of life, subset uniref lca annotation and add to taxmap
    tol = TreeOfLife()
    taxmap = {}
    tol_mode = False
    lca_mode = False
    with util.try_zip_open( p_datafile ) as fh:
        print( "Loading taxonomic data from: "+p_datafile, file=sys.stderr )
        for row in csv.reader( fh, csv.excel_tab ):
            if row[0] == c_tol_header:
                print( "  Loading TOL data", file=sys.stderr )
                tol_mode = True
                continue
            if row[0] == c_lca_header:
                print( "  Loading LCA data", file=sys.stderr )
                tol_mode = False
                lca_mode = True
                continue
            if tol_mode:
                tol.attach( Taxon( *row ) )
            elif lca_mode:
                uni, lca = row
                if uni in unirefs:
                    for rank, common in tol.get_lineage( lca ):
                        if rank == target_rank:
                            taxmap[uni] = rank.lower()[0] + "__" + simplify( common )
                            break
    # augment taxmap with genus-level lineage information for stratified features
    for feature in features:
        feature, name, stratum = util.fsplit( feature )
        if stratum is not None and "g__" in stratum:
            genus = stratum.split( util.c_taxon_delim )[0]
            if target_rank == "Genus":
                taxmap[stratum] = genus
            else:
                genus = genus.replace( "g__", "" )
                for rank, common in tol.get_lineage( genus ):
                    if rank == target_rank:
                        taxmap[stratum] = rank.lower()[0] + "__" + simplify( common )
                        break
    return taxmap

def tax_connect( feature, taxmap ):
    old = feature
    feature, name, stratum = util.fsplit( feature )
    if stratum is None or stratum == c_unclassified:
        stratum2 = taxmap.get( feature, c_unclassified )
    else:
        stratum2 = taxmap.get( stratum, c_unclassified )
    return util.fjoin( feature, name, stratum2 )

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args( )
    tbl = util.Table( args.input )
    # build the taxmap
    print( "Building taxonomic map for input table", file=sys.stderr )
    p_datafile = args.dev if args.dev is not None else databases[args.resolution]
    taxmap = build_taxmap( tbl.rowheads, args.level, p_datafile )
    # refine the taxmap
    counts = {}
    for old, new in taxmap.items():
        counts[new] = counts.get( new, 0 ) + 1
    total = float( sum( counts.values( ) ) )
    count = {k:v/total for k, v in counts.items()}
    taxmap = {old:new for old, new in taxmap.items() if count[new] >= args.threshold}
    # reindex the table
    print( "Reindexing the input table", file=sys.stderr )
    ticker = util.Ticker( tbl.rowheads )
    index = {}
    for i, rowhead in enumerate( tbl.rowheads ):
        ticker.tick()
        feature, name, stratum = util.fsplit( rowhead )
        new_rowhead = tax_connect( rowhead, taxmap )
        # unmapped is never stratified
        if feature == util.c_unmapped:
            index.setdefault( rowhead, [] ).append( i )
        # outside of unclassified mode, keep totals
        elif stratum is None and args.mode != c_umode:
            index.setdefault( rowhead, [] ).append( i )
            # in totals mode, guess at taxonomy from uniref name
            if args.mode == c_tmode:
                index.setdefault( new_rowhead, [] ).append( i )
        elif stratum == c_unclassified and args.mode == c_umode:
            # in unclassified mode, make a new row for the total...
            index.setdefault( util.fjoin( feature, name, None ), [] ).append( i )
            # ...then replace "unclassified" with inferred taxonomy
            index.setdefault( new_rowhead, [] ).append( i )
        elif stratum is not None and args.mode == c_smode:
            index.setdefault( new_rowhead, [] ).append( i )
    # rebuild the table
    print( "Rebuilding the input table", file=sys.stderr )
    rowheads2, data2 = [], []
    ticker = util.Ticker( index )
    for rowhead in util.fsort( index ):
        ticker.tick()
        rowheads2.append( rowhead )
        newrow = [0 for k in tbl.colheads]
        for i in index[rowhead]:
            oldrow = map( float, tbl.data[i] )
            newrow = [a + b for a, b in zip( newrow, oldrow )]
        data2.append( newrow )
    tbl.rowheads = rowheads2
    tbl.data = data2
    # output
    print( "Writing new table", file=sys.stderr )
    tbl.write( args.output, unfloat=True )
    # report on performance
    success, total = 0, 0
    for rowhead in tbl.rowheads:
        feature, name, stratum = util.fsplit( rowhead )
        if stratum is not None:
            total += 1
            if stratum != c_unclassified:
                success += 1
    print( "Summary: Of {TOTAL} stratifications, {SUCCESS} mapped at {TARGET} level ({PERCENT}%)".format( 
            TOTAL=total, 
            SUCCESS=success, 
            TARGET=args.level, 
            PERCENT=round( 100 * success / float( total ), 1 ),
            ), file=sys.stderr,
           )

if __name__ == "__main__":
    main( )
