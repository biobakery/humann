#!/usr/bin/env python

from __future__ import print_function
import sys
import csv
import re
import argparse

from zopy.utils import try_open, qw, die
from zopy.dag import Node, DAG

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_na = "-"

# ---------------------------------------------------------------
# classes
# ---------------------------------------------------------------

class Clade( Node ):
    def __init__( self, name ):
        Node.__init__( self, name )
        self.weight = 0
        self.common = name
        self.rank = c_na
        self.lineage = c_na
        self.status = c_na
    def get_parent( self ):
        # Node class allows multiple parents
        p = None
        if len( self.parents ) > 1:
            die( "Clade with more than one parent" )
        elif len( self.parents ) == 1:
            p = list( self.parents )[0]
        return p

class TOL( DAG ):
    def get( self, name ):
        if name not in self.node_dict:
            self.add( Clade( name ) )
        return self.node_dict[name]
        
# ---------------------------------------------------------------
# cli
# ---------------------------------------------------------------

# argument parsing (python argparse)
parser = argparse.ArgumentParser()
parser.add_argument( "ncbi_taxonomy", help="" )
parser.add_argument( "uniref_headers", help="" )
args = parser.parse_args()

# ---------------------------------------------------------------
# load and curate full taxonomy as dag
# ---------------------------------------------------------------

print( "Loading tree of life", file=sys.stderr )

"""
 0 Taxon
 1 Mnemonic
 2 Scientific name
 3 Common name
 4 Synonym
 5 Other Names
 6 Reviewed
 7 Rank
 8 Lineage
 9 Parent
10 Virus hosts
"""

def check( string, default=c_na ):
    if string == "":
        string = default
    return string

tol = TOL( )
common_tax = {}
headers = False
with try_open( args.ncbi_taxonomy ) as fh:
    for row in csv.reader( fh, csv.excel_tab ):
        if not headers:
            headers = True
        else:
            tax       = row[0]
            common    = check( row[2], tax )
            rank      = check( row[7] )
            lineage   = check( row[8] )
            parent    = check( row[9] )
            n         = tol.get( tax )
            n.common  = common
            n.rank    = rank
            n.lineage = lineage
            n.weight  = 0
            p         = tol.get( parent )
            n.add_parent( p )
            common_tax.setdefault( common, [] ).append( tax ) 
tol.update( )

# ---------------------------------------------------------------
# disambiguate common names
# ---------------------------------------------------------------

for clade in tol.nodes( ):
    if clade.is_leaf:
        for lin in clade.get_lineages( ):
            for ancestor in lin:
                ancestor.weight += 1

for common, taxes in common_tax.items( ):
    if len( taxes ) > 1:
        dominant = None
        # determine if one use of common name dominates
        # ex. kingdom "Bacteria" versus stick insect genus "Bacteria"
        if dominant is None:
            choices = []
            for tax in taxes:
                choices.append( [tol.get( tax ).weight, tax] )
            choices.sort( )
            if choices[-1][0] / (1 + choices[-2][0]) > 10:
                dominant = choices[-1][1]
        # if no one is dominant by weight, take most specific
        # this resolves bacterial phylogeny where e.g. p__, c__, o__ have same name
        if dominant is None:
            choices = []
            for tax in taxes:
                lin = tol.get( tax ).get_lineages( )[0]
                choices.append( [len( lin ), tax] )
            choices.sort( )
            if choices[-1][0] > choices[-2][0]:
                dominant = choices[-1][1]
        # enact
        for tax in taxes:
            clade = tol.get( tax )
            if tax == dominant:
                clade.status = "AmbiguousConnect"
            else:
                clade.status = "AmbiguousBypass"

# ---------------------------------------------------------------
# clean up ranks, output TOL
# ---------------------------------------------------------------

for node in tol.nodes( ):
    # what we call kingdom ncbi calls superkingdom
    if node.rank == "Kingdom":
        node.rank = c_na
    elif node.rank == "Superkingdom":
        node.rank = "Kingdom"
    # helpful to label major virus subdivisions as phyla
    if not node.is_root and node.get_parent( ).common == "Viruses":
        node.rank = "Phylum"

print( "# TOL" )
for node in tol.nodes( ):
    if node.name == c_na:
        continue
    outline = [
        node.name,
        node.common,
        node.rank,
        node.get_parent( ).name if node.get_parent( ) is not None else c_na,
        node.status,
        ]
    print( "\t".join( [str( k ) for k in outline] ) )

# ---------------------------------------------------------------
# load uniref taxonomy
# ---------------------------------------------------------------

print( "Loading uniref taxonomy", file=sys.stderr ) 

print( "# LCA" )
uni_tax = {}
with try_open( args.uniref_headers ) as fh:
    for line in fh:
        # >UniRef50_Q6GZX4 Putative transcription factor 001R n=13 Tax=Ranavirus RepID=001R_FRG3G
        if line[0] == ">":
            uni = line[1:].split( )[0]
            tax = re.search( "Tax=(.*) RepID", line ).group( 1 )
            print( "\t".join( [uni, tax] ) )
