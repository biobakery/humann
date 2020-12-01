#!/usr/bin/env python

from __future__ import print_function # Python 2.7+ required

import os
import sys
import csv
import re
import argparse
from collections import Counter

try:
    from humann import config
    from humann.tools import util
except ImportError:
    sys.exit("CRITICAL ERROR: Unable to find the HUMAnN python package." +
        " Please check your install.")

description = """
HUMAnN utility for inferring "unclassified" taxonomy
=====================================================
Based on the lowest common ancestor (LCA) annotation
of each UniRef50/90 cluster, infer approximate taxonomy 
for unclassified features at a target level of resolution. 
Will modify features of known genus/species to match 
target level."""

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

c_lca_choices = [
    "source_tax",
    "uniref_lca",
    "humann_lca",
    ]

c_tmode = "totals"
c_umode = "unclassified"
c_smode = "stratified"

c_tol_header = "TAXID"
c_lca_header = "FAMILY"

c_root = "Root"
c_unclassified = util.c_unclassified

# ---------------------------------------------------------------
# utility mapping interface
# ---------------------------------------------------------------

# get a list of all available utility mapping files
try:
    all_mapping_files=os.listdir( config.utility_mapping_database )
except EnvironmentError:
    all_mapping_files=[]

# find the tol-lca files
databases = {}
for f in all_mapping_files:
    if "tol-lca" in f:
        k = f.split( "." )[0]
        v = os.path.join( config.utility_mapping_database, f )
        databases[k] = v

# warn if tol-lca files are missing
if len( databases ) == 0 and "--devdb" not in sys.argv:
    sys.exit( 
"""
You are missing HUMAnN utility mapping files required to use
this script. To add these to your installation, please execute:

$ humann_databases --download utility_mapping full $DIR

Replacing, $DIR with the directory to install the databases.
""" )

# ---------------------------------------------------------------
# helper objects
# ---------------------------------------------------------------

class Taxon:
    
    def __init__( self, rowdict ):
        self.taxid   = rowdict["TAXID"]
        self.sciname = rowdict["NAME"]
        self.rank    = rowdict["RANK"]
        self.parent  = rowdict["PARENT"]

class TreeOfLife:
    
    def __init__( self ):
        self.nodes = {}
        self.genus_map = {}
    
    def add_node( self, node ):
        self.nodes[node.taxid] = node
        if node.rank == "Genus":
            inner = self.genus_map.setdefault( node.sciname, [] )
            inner.append( node.taxid )
    
    def get_lineage( self, taxid ):
        lineage = []
        while taxid in self.nodes and taxid != c_root:
            node = self.nodes[taxid]
            lineage.append( [node.rank, node.sciname] )
            taxid = node.parent
        return lineage
    
    def get_genus_taxid( self, genus ):
        # genus name MIGHT be ambiguous (compensate)
        genera = self.genus_map.get( genus, [] )
        return genera[0] if len( genera ) == 1 else None

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
                         help="HUMAnN genefamilies table" )
    parser.add_argument( "-o", "--output", 
                         default=None,
                         help="Destination for modified table; default=STDOUT" )
    parser.add_argument( "-l", "--level",
                         choices=c_levels,
                         default="Family",
                         help="Desired level for taxonomic estimation/summation; default=Family" )
    parser.add_argument( "-d", "--database",
                         choices=databases.keys( ),
                         help="UniRef-specific taxonomy database" )
    parser.add_argument( "-m", "--mode", 
                         choices=[c_tmode, c_umode, c_smode],
                         default=c_tmode,
                         help="Which rows to include in the estimation/summation; default=totals" )
    parser.add_argument( "-c", "--lca-choice", 
                         choices=c_lca_choices,
                         default="humann_lca",
                         help="Which per-gene taxonomic annotation to consider; default=humann_lca" )
    parser.add_argument( "-t", "--threshold", 
                         type=float, 
                         default=1e-3, 
                         help="Minimum frequency for a new taxon to be included; default=1e-3" )
    parser.add_argument( "--devdb",
                         help="Manually specify a development database" )
    args = parser.parse_args( )
    return args

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

def simplify( sciname ):
    return re.sub( "[^A-Za-z0-9]+", "_", sciname )

def rowdict( headers, row ):
    return {h:r for h, r in zip( headers, row )}

def build_taxmap( features, target_rank, lca_choice, p_tol_lca ):
    # pull out relevant uniref ids (avoids full uniref taxonomy in memory)
    unirefs = {util.fsplit( k )[0] for k in features}
    unirefs = {k for k in unirefs if "UniRef" in k}
    # tree-of-life object
    tol = TreeOfLife()
    # mapping from uniref id to taxon sciname or unclassified
    taxmap = {}
    tol_mode = False
    lca_mode = False
    headers = None
    with util.try_zip_open( p_tol_lca ) as fh:
        print( "  Loading taxonomic data from:", p_tol_lca, file=sys.stderr )
        for row in csv.reader( fh, csv.excel_tab ):
            # determine which parsing mode we're in
            if row[0] == c_tol_header:
                print( "  Loading TOL data", file=sys.stderr )
                tol_mode = True
                headers = row
                continue
            if row[0] == c_lca_header:
                print( "  Loading LCA data", file=sys.stderr )
                tol_mode = False
                lca_mode = True
                headers = row
                continue
            # row is a taxon, add to tree
            if tol_mode:
                R = rowdict( headers, row )
                tol.add_node( Taxon( R ) )
            # row is a uniref lca entry, add to taxmap if relevant
            elif lca_mode:
                R = rowdict( headers, row )
                u = R["FAMILY"]
                if u in unirefs:
                    lca = R[lca_choice.upper( )]
                    for rank, sciname in tol.get_lineage( lca ):
                        if rank == target_rank:
                            taxmap[u] = rank.lower()[0] + "__" + simplify( sciname )
                            break
    # augment taxmap with genus-level lineage information for stratified features
    for feature in features:
        feature, name, stratum = util.fsplit( feature )
        if stratum is not None and "g__" in stratum:
            genus = stratum.split( util.c_taxon_delim )[0]
            if target_rank == "Genus":
                # directly assign stratified genus
                taxmap[stratum] = genus
            else:
                # try to find genus ancestor based on name
                genus = genus.replace( "g__", "" )
                taxid = tol.get_genus_taxid( genus )
                if taxid is not None:
                    for rank, sciname in tol.get_lineage( taxid ):
                        if rank == target_rank:
                            taxmap[stratum] = rank.lower()[0] + "__" + simplify( sciname )
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
    T = util.Table( args.input )

    # build the taxmap (uniref -> taxon mapping)
    print( "Building taxonomic map for input table", file=sys.stderr )
    if args.devdb is not None:
        p_datafile = args.devdb
    elif args.database in databases:
        p_datafile = databases[args.database]
    else:
        sys.exit( "Must specify a valid database (from utility mapping or --devdb)" )
    taxmap = build_taxmap( T.rowheads, args.level, args.lca_choice, p_datafile )

    # refine the taxmap (remove rare taxa)
    counts = Counter( taxmap.values( ) )
    total  = float( sum( counts.values( ) ) )
    counts = {k:v/total for k, v in counts.items()}
    taxmap = {old:new for old, new in taxmap.items() if counts[new] >= args.threshold}

    # reindex the table (which rows to keep, which rows to pseudo-stratify)
    print( "Reindexing the input table", file=sys.stderr )
    ticker = util.Ticker( T.rowheads )
    index = {}
    for i, rowhead in enumerate( T.rowheads ):
        ticker.tick( )
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
        # in unclassified mode, make a new row for the total...
        elif stratum == c_unclassified and args.mode == c_umode:
            index.setdefault( util.fjoin( feature, name, None ), [] ).append( i )
            # ...then replace "unclassified" with inferred taxonomy
            index.setdefault( new_rowhead, [] ).append( i )
        # update strata in stratified mode
        elif stratum is not None and args.mode == c_smode:
            index.setdefault( new_rowhead, [] ).append( i )

    # rebuild the table
    print( "Rebuilding the input table", file=sys.stderr )
    rowheads2, data2 = [], []
    ticker = util.Ticker( index )
    for rowhead in util.fsort( index ):
        ticker.tick( )
        rowheads2.append( rowhead )
        newrow = [0 for k in T.colheads]
        for i in index[rowhead]:
            oldrow = map( float, T.data[i] )
            newrow = [a + b for a, b in zip( newrow, oldrow )]
        data2.append( newrow )
    T.rowheads = rowheads2
    T.data = data2
    print( "Writing new table", file=sys.stderr )
    T.write( args.output, unfloat=True )

    # report on performance
    success, total = 0, 0
    for rowhead in T.rowheads:
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

# run in script mode
if __name__ == "__main__":
    main( )
