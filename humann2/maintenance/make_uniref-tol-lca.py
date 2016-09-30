#!/usr/bin/env python

import sys
import csv
import re
import argparse
from zopy.utils import try_open

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

levels = """
Kingdom
Phylum
Class
Order
Family
Genus
Species
"""
ranks = [k for k in levels.split( "\n" ) if k != ""]

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

print("Loading tree of life")

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
#
map_name_tax = {}
map_tax_parent = {}
map_parent_tax = {}
map_tax_rank = {}
headers = False
with try_open( args.ncbi_taxonomy ) as fh:
    for row in csv.reader( fh, csv.excel_tab ):
        if headers:
            map_tax_parent[row[0]] = row[9]
            map_parent_tax[row[9]] = row[0]
            map_name_tax.setdefault( row[2], [] ).append( row[0] )
            map_tax_rank[row[0]] = row[7]
        else:
            headers = True

# clean taxonomy
for name, taxes in map_name_tax.items():
    # **** sloppy: if the name isn't unique, assign to lowest taxid ****
    map_name_tax[name] = sorted( taxes, key=lambda k: int( k ) )[0]
map_tax_name = {v:k for k, v in map_name_tax.items()}

# update parent
map_name_parent = {}
for c, p in map_tax_parent.items():
    if c in map_tax_name and p in map_tax_name:
        map_name_parent[map_tax_name[c]] = map_tax_name[p]

# update rank
map_name_rank = {}
for c, r in map_tax_rank.items():
    if r == "Kingdom":
        r = "Subkingdom"
    elif r == "Superkingdom":
        r = "Kingdom"
    if c in map_tax_name and r in ranks:
        map_name_rank[map_tax_name[c]] = r

print("# TOL")
for name in sorted( map_name_tax ):
    print("\t".join( [name, map_name_rank.get( name, "n/a" ), map_name_parent.get( name, "root" )] ))

# ---------------------------------------------------------------
# load uniref taxonomy
# ---------------------------------------------------------------

print("Loading uniref taxonomy")

print("# LCA")
map_uni_tax = {}
with try_open( args.uniref_headers ) as fh:
    for line in fh:
        # >UniRef50_Q6GZX4 Putative transcription factor 001R n=13 Tax=Ranavirus RepID=001R_FRG3G
        if line[0] == ">":
            uni = line[1:].split( )[0]
            tax = re.search( "Tax=(.*) RepID", line ).group( 1 )
            print("\t".join( [uni, tax] ))
