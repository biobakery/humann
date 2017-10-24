#!/usr/bin/env python

from __future__ import print_function # Python 2.7+ required
import os
import sys
import csv
import re
import argparse
from collections import defaultdict

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
HUMAnN2 utility for inferring taxonomy

Given a gene families file, this script can infer the taxonomy of
unclassified features based on known lowest command ancestor (LCA)
annotations of UniRef50/90 clusters. This script can also be applied 
to _any_ HUMAnN2 table to regroup species-level stratifications to 
broader taxonomic levels ("species" mode).
""" )

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
This script requires the HUMAnN2 utility data files.
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
    "stratified:   adjust all strata (UniRef50/90 only)",
    "species:      adjust species strata only",
    "unclassified: adjust unclassified, discarding totals + species (UniRef50/90 only)",
    "totals:       adjust totals, discarding strata (UniRef50/90 only)",
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
        formatter_class=argparse.RawTextHelpFormatter,
        )
    util.attach_common_arguments( parser )   
    parser.add_argument( "-l", "--taxonomic-level",
                         choices=c_levels,
                         metavar="<choice>",
                         default="Genus",
                         help=util.pretty_grid( c_levels, cols=7,
                                                desc="Level for taxonomic summarization [Default=Genus]:" ), 
                         )
    parser.add_argument( "-r", "--resolution",
                         metavar="<choice>",
                         choices=databases.keys( ),
                         help=util.pretty_grid( databases.keys( ), "Required outside of 'species' mode:" ),
                         )
    parser.add_argument( "-m", "--mode",
                         choices=c_mode_names,
                         default="stratified",
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
    """Get species->genus map from HUMAnN2 stratifications"""
    taxmap = {}
    for feature in features:
        fbase, fname, stratum = util.fsplit( feature )
        if stratum is not None and stratum != util.c_unclassified:
            genus = stratum.split( util.c_taxon_delim )[0]
            taxmap[stratum] = genus
    return taxmap

def complete_taxmap( features, target_rank, p_datafile ):
    """Load full taxonomy from the TOL file"""
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
                            taxmap[uni] = rank.lower( )[0] + "__" + simplify( common )
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
    """Adjust a feature's taxonomy based on the taxmap"""
    old = feature
    feature, name, stratum = util.fsplit( feature )
    # get taxonomy based on UniRef annotation (UniRefs only)
    if stratum is None or stratum == util.c_unclassified:
        stratum2 = taxmap.get( feature, util.c_unclassified )
    # get taxonomy based on HUMAnN2 taxonony (any type of feature)
    else:
        stratum2 = taxmap.get( stratum, util.c_unclassified )
    return util.fjoin( feature, name, stratum2 )

def generate_mapping( table, taxmap, mode ):
    """Determine which original-table rows to keep and how to recombine them;
    varies substantially with the --mode option""" 
    mapping = defaultdict( set )
    for f in table.data:
        fbase, fname, stratum = util.fsplit( f )
        # new_f is f with new taxonomy (if found) else "unclassified"
        new_f = tax_connect( f, taxmap )
        # community total
        if stratum is None:
            # always keep UNMAPPED
            if f == util.c_unmapped:
                mapping[f].add( f )
            # keep original totals unless in "unclassified" mode
            elif mode != "unclassified":
                mapping[f].add( f )
                # total becomes a new stratum in "totals" mode
                if mode == "totals":
                    mapping[new_f].add( f )
        # unclassified stratum
        elif stratum == util.c_unclassified:
            # create a new total in unclassified mode and infer
            if mode == "unclassified":
                new_tot = util.fjoin( fbase, fname, None )
                mapping[new_tot].add( f )
                mapping[new_f].add( f )
            # infer in stratified mode
            elif mode == "stratified":
                mapping[new_f].add( f )
            # just pass through in species mode
            elif mode == "species":
                mapping[f].add( f )
        # this must be a known-species stratum
        elif "s__" in stratum:
            if mode in ["stratified", "species"]:
                mapping[new_f].add( f )
    return mapping

def tax_report( table, path ):
    """Write a summary of the new taxa in the output table"""
    # s --> stratum throughout this function
    stacks = {}
    for f in table.data:
        fbase, fname, s = util.fsplit( f )
        if s is not None:
            stacks.setdefault( s, [] ).append( table.data[f] )
    totals = table.zeros( )
    # sum within-sample, normalize to sample totals
    masses = {}
    for s, stack in stacks.items( ):
        masses[s] = np.sum( np.vstack( stack ), axis=0 )
        totals += masses[s]
    masses = {s:np.mean( row/totals ) for s, row in masses.items( )}
    # report
    with util.try_zip_open( path, "w" ) as fh:
        print( "Taxon\tMean % of new abundance", file=fh )
        for s in sorted( masses, key=lambda x: -masses[x] ):
            if masses[s] > 0:
                print( "{}\t{:.1f}".format( s, 100 * masses[s] ), file=fh )

def load_appropriate_taxmap( table, args ):
    # infer taxmap from humann2 stratifications
    if args.mode == "species" and args.taxonomic_level == "Genus":
        taxmap = genus_taxmap( table.data.keys( ) )
    # get from the uniref50 tol file, although no uniref data needed
    elif args.mode == "species":
        taxmap = complete_taxmap( table.data.keys( ), args.taxonomic_level, databases["uniref50"] )
    # load a dev taxmap
    elif args.dev is not None:
        taxmap = complete_taxmap( table.data.keys( ), args.taxonomic_level, args.dev )
    # fail if the specific taxmap doesn't exist
    elif args.resolution not in databases:
        sys.exit( ("CRITICAL ERROR: Outside of 'species' mode you must specify your UniRef resolution\n"
                   "using the -r/--resolution flag") )
    # typical case
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
    # done
    return taxmap

def summarize_success( mapping, args ):
    """print a summary of features that are no longer unclassified after mapping"""
    total = set( )
    unclass = set( )
    for f in mapping:
        fbase, fname, stratum = util.fsplit( f )
        if stratum is not None:
            total.add( fbase )
            if stratum == util.c_unclassified:
                unclass.add( fbase )
    success = total - unclass
    success_rate = 0
    if len( total ) > 0:
        success_rate = 100 * len( success ) / float( len( total ) )
    print( "Reclassification summary:", file=sys.stderr )
    print( "  Level: {}".format( args.taxonomic_level ), file=sys.stderr )
    print( "  Features considered: {:,}".format( len( total ) ), file=sys.stderr )
    print( "  Fully classified at target level: {:,} ({:.1f})%".format( 
            len( success ), success_rate ), file=sys.stderr )

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args( )
    table = Table( args.input, last_metadata=args.last_metadata )
    # make a taxmap
    print( "Building taxonomic map for input table", file=sys.stderr )
    taxmap = load_appropriate_taxmap( table, args )
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
    summarize_success( mapping, args )
    # output
    new_table.write( args.output, unfloat=True )
    # write tax report?
    if args.taxonomy_report is not None:
        print( "Writing taxonomy report to <{}>".format( args.taxonomy_report ), file=sys.stderr )
        tax_report( new_table, args.taxonomy_report )

if __name__ == "__main__":
    main( )
