#!/usr/bin/env python

from __future__ import print_function # python 2.7+ required
import subprocess
import argparse
import csv
import gzip

c_halflen = 30
c_maxlen = 2 * c_halflen

parser = argparse.ArgumentParser()
parser.add_argument( 'version', help='eggnog version, e.g. 4.5' )
args = parser.parse_args()

url = "http://eggnogdb.embl.de/download/eggnog_{}/data/NOG/NOG.annotations.tsv.gz".format( args.version )

command = "curl {} | zcat | cut -f2,6".format( url )
ph = subprocess.Popen( command, shell=True, stdout=subprocess.PIPE )

data = {}
for row in csv.reader( ph.stdout, csv.excel_tab ):
    ident, name = row
    if name != "NA":
        if len( name ) > c_maxlen:
            name = name[0:c_halflen] + "[...]" + name[-c_halflen:]
        data[ident] = name

with gzip.GzipFile( "map_eggnog_name.txt.gz", "w" ) as fh:
    for ident in sorted( data ):
        print( ident+"\t"+data[ident], file=fh )
