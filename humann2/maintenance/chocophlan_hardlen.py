#!/usr/bin/env python

import os
import sys
import gzip
from Bio import SeqIO

modified = []
with gzip.GzipFile( sys.argv[1] ) as fh:
    for record in SeqIO.parse( fh, "fasta" ):
        record.id += "|" + str( len( record ) )
        modified.append( record )

for record in modified:
    record.name, record.description = "", ""

newname = os.path.split( sys.argv[1] )[1].replace( ".ffn.m8.gz", ".v0.1.1.ffn.gz" )
with gzip.GzipFile( newname, "w" ) as fh:
    SeqIO.write( modified, fh, "fasta" )
