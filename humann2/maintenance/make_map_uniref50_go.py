#!/usr/bin/env python

"""
Create a uniref50 to go mapping.
Isolates only those uniref50 entires with go annotation.
Downloads the uniprot GOA mapping if you don't have it (5-10GB).
====
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

import os, sys, re, argparse, subprocess

# constants
url_uniprotgoa = "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz"

# argument parsing (python argparse)
parser = argparse.ArgumentParser()
parser.add_argument( "--uniref50_fasta", \
                     help="path to the uniref50 fasta file" )
parser.add_argument( "--uniprot_goa", default=None, \
                     help="path to the uniprotgoa file [can download if missing" )
parser.add_argument( "--output", default="map_uniref50_go.txt", \
                     help="output map file" )
args = parser.parse_args()

# derived constants
path_uniref50_fasta = args.uniref50_fasta
path_uniprot_goa = args.uniprot_goa
path_output = args.output

# check/get uniprot_goa
if ( args.uniprot_goa is None ) or \
   ( not os.path.exists( args.uniprot_goa ) ):
    response = raw_input( "uniprot_goa file doesn't exist; download? [Y/n]" )
    if response[0] == "Y":
        os.system( "wget %s" % ( url_uniprot_goa ) )
        path_uniprot_goa = "gene_association.goa_uniprot.gz"
    else:
        sys.exit( "cannot continue" )

# derive the list of uniref50 headers
print >>sys.stderr, "Loading uniref50 ids from the fasta file...",
cmd = subprocess.Popen( 
    "grep '>' %s | perl -pe 's/>UniRef50_(.*?) .*/\\1/'" % ( path_uniref50_fasta ), 
    shell=True, 
    stdout=subprocess.PIPE, )
map_uniref50_go = {}
for line in cmd.stdout:
    map_uniref50_go[line.strip()] = {}
print >>sys.stderr, "complete.", "Found", len( map_uniref50_go ), "UniRef50 IDs."

# process goa file, save lines where uniprotid is uniref50
print >>sys.stderr, "Parsing GOA file for uniref50 entries...",

"""
Notes:
The GOA file is tab-delimited. 
Col2 is the uniprot id (a superset of uniref50).
Col5 is the Gene Ontology annotation.
Col4 is a logical modifier of the uniprot->go mapping.
Must exclude the cases where this is "NOT".
"""

cmd = subprocess.Popen( 
    "zcat %s | cut -f2,4,5" % ( path_uniprot_goa ),
    shell=True, 
    stdout=subprocess.PIPE, )
for line in cmd.stdout:
    if line[0] == "!":
        continue
    uniprot_id, modifier, go_id = line.strip().split( "\t" )
    if uniprot_id in map_uniref50_go:
        if "NOT" not in modifier:
            map_uniref50_go[uniprot_id][go_id] = 1
print >>sys.stderr, "complete."

# write the output file
print >>sys.stderr, "Writing the UniRef50 to GO mapping...",
number_annotated = 0
with open( path_output, "w" ) as fh:
    for uniprot_id, go_ids in map_uniref50_go.items():
        if len( go_ids ) > 0:
            number_annotated += 1
            print >>fh, "\t".join( ["UniRef50_"+uniprot_id] + go_ids.keys() )
print >>sys.stderr, "complete.", "There were", number_annotated, "annotated entries."
