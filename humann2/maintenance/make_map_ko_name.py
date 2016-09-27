#!/usr/bin/env python

"""
Make script for building KO name map
Adds legacy KO names from HUMAnN1
Does not replace new names with old names
=======================================
Eric Franzosa (eric.franzosa@gmail.com)
"""

import os
import sys
import re
import gzip
import argparse
import csv

# argument parsing (python argparse)
parser = argparse.ArgumentParser()
parser.add_argument( "humann1_map_kegg", help="legacy KO names" )
args = parser.parse_args()

# constants
c_kegg_ko_url = "http://www.genome.jp/kegg-bin/download_htext?htext=ko00001.keg&format=htext&filedir="
p_tempfile    = "temp-ko.txt"
p_outfile     = "map_ko_name.txt.gz"

# check if KEGG file exists
if not os.path.exists( p_tempfile ):
    print >>sys.stderr, "wget-ing raw data to", p_tempfile
    os.system( "wget {} --output-document {}".format( c_kegg_ko_url, p_tempfile ) )
else:
    print >>sys.stderr, p_tempfile, "exists, so using it; delete to force download"

"""
Format of the KEGG file.
Note that KOs are repeated under different subtrees.
==== 

+D      KO
#<h2><a href="/kegg/kegg2.html"><img src="/Fig/bget/kegg3.gif" align="middle" border=0></a> &nbsp; KEGG Orthology (KO)</h2>
#<!---
#ENTRY       ko00001
#NAME        KO
#DEFINITION  KEGG Orthology (KO)
#--->
!
A<b>Metabolism</b>
B
B  <b>Overview</b>
C    01200 Carbon metabolism [PATH:ko01200]
D      K00844  HK; hexokinase [EC:2.7.1.1]
D      K12407  GCK; glucokinase [EC:2.7.1.2]
D      K00845  glk; glucokinase [EC:2.7.1.2]
"""

# pairs from the web file
ko2name = {}
with open( p_tempfile ) as fh:
    for line in fh:
        match = re.search( "D      (K[0-9]{5})  (.*)", line )
        if match:
            ko, name = match.groups( )
            # strip gene ids from front of name
            name = name.split( "; " )[-1]
            ko2name[ko] = name

# pairs from the legacy file
with open( args.humann1_map_kegg ) as fh:
    for code, name in csv.reader( fh, csv.excel_tab ):
        if re.search( "K[0-9]{5}", code ):
            ko = code
            # do NOT override new names
            if ko not in ko2name:
                ko2name[ko] = name

# dump to file
with gzip.GzipFile( p_outfile, "w" ) as fh:
    for ko in sorted( ko2name ):
        name = ko2name[ko]
        print >>fh, "\t".join( [ko, name] )
print >>sys.stderr, "output written to", p_outfile
