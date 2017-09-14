#!/usr/bin/env python

from __future__ import print_function # Python 2.7+ required
import os
import sys
import csv
import re
import argparse

try:
    from humann2 import config
    #from humann2.tools import util
    #from humann2.tools.better_table import Table
    import util
    from better_table import Table
except ImportError:
    sys.exit( "CRITICAL ERROR: Unable to find the HUMAnN2 python package.\n" +
              "Please check your install." )

description = """
HUMAnN2 utility for inferring "unclassified" taxonomy
=====================================================
Based on the lowest common ancestor (LCA) annotation
of each UniRef50/90 cluster, infer approximate taxonomy 
for unclassified features at a target level of resolution 
(default=Family). Will modify features of known genus/species 
to match target level."""

# ---------------------------------------------------------------
# check that user has the required database file
# ---------------------------------------------------------------

try:
    all_mapping_files = os.listdir( config.utility_mapping_database )
except EnvironmentError:
    all_mapping_files = []

databases = {
    "uniref50": "uniref50-tol-lca.dat.gz",
    "uniref90": "uniref90-tol-lca.dat.gz",
    }

SOMETHING_MISSING = False
for key, value in databases.items( ):
    if value not in all_mapping_files:
        SOMETHING_MISSING = True
    else:
        databases[key] = os.path.join( config.utility_mapping_database, value )

if SOMETHING_MISSING:
    sys.exit( """
This tool requires the HUMAnN2 utility data files.
To add these to your installation, please execute:

$ humann2_databases --download utility_mapping full $DIR

Replacing, $DIR with the directory to install the databases.
""" )

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
# command-line interface
# ---------------------------------------------------------------

def get_args( ):
    """ Get args from Argparse """
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter
        )
    parser.add_argument( "-i", "--input", 
                         required=True,
                         metavar="<path>", 
                         help="HUMAnN2 genefamilies.tsv file",
                         )
    parser.add_argument( "-o", "--output", 
                         default=None,
                         metavar="<path>",
                         help="Destination for modified table\n[Default=STDOUT]",
                         )
    parser.add_argument( "-l", "--level",
                         choices=c_levels,
                         metavar="<choice>",
                         default="Family",
                         help=util.pretty_grid( c_levels, cols=7,
                                                desc="Level for taxonomic summarization [default=Family]:" ), 
                         )
    parser.add_argument( "-r", "--resolution",
                         metavar="<choice>",
                         choices=databases.keys( ),
                         required=True,
                         help=util.pretty_grid( databases.keys( ), "Input UniRef type:" ),
                         )
    parser.add_argument( "-m", "--mode", 
                         choices=[c_tmode, c_umode, c_smode],
                         default=c_tmode,
                         metavar="<choice>",
                         help=util.pretty_grid( [c_tmode, c_umode, c_smode],
                                           cols=3, desc="Which rows to operate on [Default={}]:".format( c_tmode ) ),
                         )
    parser.add_argument( "-t", "--threshold", 
                         type=float, 
                         default=1e-3,
                         metavar="<float>",
                         help="Minimum frequency for a new taxon to be included\n[Default=0.001]",
                         )
    parser.add_argument( "-d", "--dev",
                         metavar="<path>",
                         help="Manually specify a development database",
                         )
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
    tol = TreeOfLife( )
    taxmap = {}
    tol_mode = False
    lca_mode = False
    with util.try_zip_open( p_datafile ) as fh:
        print( "Loading taxonomic data from: " + p_datafile, file=sys.stderr )
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
    table = Table( args.input )
    # build the taxmap
    print( "Building taxonomic map for input table", file=sys.stderr )
    p_datafile = args.dev if args.dev is not None else databases[args.resolution]
    taxmap = build_taxmap( table.data.keys( ), args.level, p_datafile )
    # refine the taxmap
    counts = {}
    for old, new in taxmap.items():
        counts[new] = counts.get( new, 0 ) + 1
    total = float( sum( counts.values( ) ) )
    count = {k:v/total for k, v in counts.items()}
    taxmap = {old:new for old, new in taxmap.items() if count[new] >= args.threshold}
    # reindex the table
    print( "Reindexing the table", file=sys.stderr )
    mapping = {}
    for f in table.data:
        fbase, fname, stratum = util.fsplit( f )
        new_f = tax_connect( f, taxmap )
        # unmapped is never stratified
        if f == util.c_unmapped:
            mapping.setdefault( rowhead, set( ) ).add( f )
        # outside of unclassified mode, keep original totals
        elif stratum is None and args.mode != c_umode:
            mapping.setdefault( f, set( ) ).add( f )
            # in totals mode, guess at taxonomy from uniref name
            if args.mode == c_tmode:
                mapping.setdefault( new_f, set( ) ).add( f )
        elif stratum == c_unclassified and args.mode == c_umode:
            # in unclassified mode, make a new row for the total...
            new_tot = util.fjoin( fbase, fname, None )
            mapping.setdefault( new_tot, set( ) ).add( f )
            # ...then replace "unclassified" with inferred taxonomy
            mapping.setdefault( new_f, set( ) ).add( f )
        elif stratum is not None and args.mode == c_smode:
            # work with taxonomy of known stratifications
            mapping.setdefault( new_f, set( ) ).add( f )
    # rebuild the table
    print( "Rebuilding the input table", file=sys.stderr )
    new_data = {}
    for f in mapping:
        new_data[f] = table.zeros( )
        for f2 in mapping[f]:
            new_data[f] += table.data[f2]
    new_table = Table( new_data, metadata=table.metadata, headers=table.headers )
    # report on performance
    success, total = 0, 0
    for f in new_table.data:
        fbase, fname, stratum = util.fsplit( f )
        if stratum is not None:
            total += 1
            if stratum != c_unclassified:
                success += 1
    success_rate = 100 * success / float( total )
    print( "Reclassification summary:", file=sys.stderr )
    print( "  Level: {}".format( args.level ), file=sys.stderr )
    print( "  Features: {:,}".format( total ), file=sys.stderr )
    print( "  Mapped at target level: {:,} ({:.1f})%".format( 
            success, success_rate ), file=sys.stderr )
    # output
    print( "Writing new table", file=sys.stderr )
    new_table.write( args.output, unfloat=True )

if __name__ == "__main__":
    main( )
