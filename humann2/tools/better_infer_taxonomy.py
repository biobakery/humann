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

try:
    import numpy as np
except ImportError:
    sys.exit( "CRITICAL ERROR: This script requires the python scientific stack (e.g. numpy)" )

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

c_modes = [
    "totals:       operate on original totals, discard original strata",
    "unclassified: operate on unclassified, discard original total + species",
    "stratified:   operate on all original strata, including unclassified",
    "species:      operate on original species only, pass-through unclassified",
]
c_mode_names = [k.split( ":" )[0] for k in c_modes]

c_tol_header = "# TOL"
c_lca_header = "# LCA"
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
    util.attach_common_arguments( parser )   
    parser.add_argument( "-l", "--taxonomic-level",
                         choices=c_levels,
                         metavar="<choice>",
                         default="Family",
                         help=util.pretty_grid( c_levels, cols=7,
                                                desc="Level for taxonomic summarization [Default=Family]:" ), 
                         )
    parser.add_argument( "-r", "--resolution",
                         metavar="<choice>",
                         choices=databases.keys( ),
                         help=util.pretty_grid( databases.keys( ), "Required outside of 'species' mode:" ),
                         )
    parser.add_argument( "-m", "--mode",
                         choices=c_mode_names,
                         default=["stratified"],
                         metavar="<choice>",
                         help=util.pretty_grid( c_modes, cols=1, desc="Operating mode [Default=stratified]" ),
                         )
    parser.add_argument( "-t", "--threshold", 
                         type=float, 
                         default=None,
                         metavar="<float>",
                         help="Minimum frequency for a new taxon to be included\n[Default=0/include everything]",
                         )
    parser.add_argument( "-w", "--taxonomy-report",
                         metavar="<path>",
                         default=None,
                         help="Write a taxonomy report at the specified path.\n[Default=no report]",
                         )
    parser.add_argument( "-d", "--dev",
                         metavar="<path>",
                         help="Manually specify a development database",
                         )
    args = parser.parse_args( )
    return args

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

def simplify( name ):
    return re.sub( "[^A-Za-z0-9]+", "_", name )

def genus_taxmap( features ):
    taxmap = {}
    for feature in features:
        fbase, fname, stratum = util.fsplit( feature )
        if stratum is not None and stratum != util.c_unclassified:
            genus = stratum.split( util.c_taxon_delim )[0]
            taxmap[stratum] = genus
    return taxmap

def complete_taxmap( features, target_rank, p_datafile ):
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
                        taxmap[stratum] = rank.lower( )[0] + "__" + simplify( common )
                        break
    return taxmap

def tax_connect( feature, taxmap ):
    old = feature
    feature, name, stratum = util.fsplit( feature )
    if stratum is None or stratum == util.c_unclassified:
        stratum2 = taxmap.get( feature, util.c_unclassified )
    else:
        stratum2 = taxmap.get( stratum, util.c_unclassified )
    return util.fjoin( feature, name, stratum2 )

def generate_mapping( table, taxmap, mode ):
    mapping = {}
    for f in table.data:
        fbase, fname, stratum = util.fsplit( f )
        new_f = tax_connect( f, taxmap )
        # community total
        if stratum is None:
            # always keep UNMAPPED
            if f == util.c_unmapped:
                mapping.setdefault( f, set( ) ).add( f )
            # keep original totals unless in "unclassified" mode
            elif mode != "unclassified":
                mapping.setdefault( f, set( ) ).add( f )
                # total becomes a new stratum in "totals" mode
                if mode == "totals":
                    mapping.setdefault( new_f, set( ) ).add( f )
        # unclassified stratum
        elif stratum == util.c_unclassified:
            # create a new total in unclassified mode and infer
            if mode == "unclassified":
                new_tot = util.fjoin( fbase, fname, None )
                mapping.setdefault( new_tot, set( ) ).add( f )
                mapping.setdefault( new_f, set( ) ).add( f )
            # infer in stratified mode
            elif mode == "stratified":
                mapping.setdefault( new_f, set( ) ).add( f )
            # just pass through in species mode
            elif mode == "species":
                mapping.setdefault( f, set( ) ).add( f )
        # this must be a known-species stratum
        elif "s__" in stratum:
            if mode in ["stratified", "species"]:
                mapping.setdefault( new_f, set( ) ).add( f )
    return mapping

def tax_report( table, path ):
    stacks = {}
    for f in table.data:
        fbase, fname, stratum = util.fsplit( f )
        if fbase == util.c_unmapped:
            stacks[fbase] = [table.data[f]]
        elif stratum is not None:
            stacks.setdefault( stratum, [] ).append( table.data[f] )
    means = {}
    for stratum, stack in stacks.items( ):
        means[stratum] = np.mean( np.vstack( stack ) )
    total = sum( means.values( ) )
    with util.try_zip_open( path, "w" ) as fh:
        for s in sorted( means, key=lambda x: -means[x] ):
            print( "{}\t{:.6f}".format( s, means[s] / float( total ) ), file=fh )

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args( )
    table = Table( args.input, last_metadata=args.last_metadata )
    # make a taxmap
    print( "Building taxonomic map for input table", file=sys.stderr )
    # use default genus.species annotation to regroup at genus level
    if args.mode == "species" and args.taxonomic_level == "Genus":
        taxmap = genus_taxmap( table.data.keys( ) )
    # load a full taxmap
    elif args.dev is not None:
        taxmap = complete_taxmap( table.data.keys( ), args.taxonomic_level, args.dev )
    elif args.resolution not in databases:
        print( ("EXITING: Outside of 'species' mode you must specify your UniRef resolution\n"
                "using the -r/--resolution flag: ").format( databases.keys( ) ), file=sys.stderr )
        sys.exit( )
    else:
        taxmap = complete_taxmap( table.data.keys( ), args.taxonomic_level, databases[args.resolution] )
    # optionally forget very rare taxa in the taxmap
    if args.threshold > 0:
        counts = {}
        for old, new in taxmap.items( ):
            counts[new] = counts.get( new, 0 ) + 1
        total = float( sum( counts.values( ) ) )
        counts = {k:v/total for k, v in counts.items( )}
        for old, new in taxmap.items( ):
            if counts[new] < args.threshold:
                taxmap[old] = util.c_unclassified
    # reindex the table
    print( "Reindexing the table", file=sys.stderr )
    mapping = generate_mapping( table, taxmap, args.mode )
    # rebuild the table
    print( "Rebuilding the input table", file=sys.stderr )
    new_data = {}
    for f in mapping:
        new_data[f] = table.zeros( )
        for f2 in mapping[f]:
            new_data[f] += table.data[f2]
    new_table = Table( new_data, metadata=table.metadata, headers=table.headers )
    # report on performance
    total = set( )
    unclass = set( )
    for f in mapping:
        fbase, fname, stratum = util.fsplit( f )
        if stratum is not None:
            total.add( fbase )
            if stratum == util.c_unclassified:
                unclass.add( fbase )
    # success -> feature no longer has unclassified level
    success = total - unclass
    success_rate = 100 * len( success ) / float( len( total ) )
    print( "Reclassification summary:", file=sys.stderr )
    print( "  Level: {}".format( args.taxonomic_level ), file=sys.stderr )
    print( "  Features considered: {:,}".format( len( total ) ), file=sys.stderr )
    print( "  Fully classified at target level: {:,} ({:.1f})%".format( 
            len( success ), success_rate ), file=sys.stderr )
    # output
    print( "Writing new table", file=sys.stderr )
    new_table.write( args.output, unfloat=True )
    # write tax report?
    if args.taxonomy_report is not None:
        print( "Writing taxonomy report at:", args.taxonomy_report, file=sys.stderr )
        tax_report( new_table, args.taxonomy_report )

if __name__ == "__main__":
    main( )
