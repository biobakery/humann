#!/usr/bin/env python

"""
I/O pipe for converting DNA reads to peptides.
Ignores any reading frame from the read with a
non-amino acid character (notably "*" = STOP)

Usage: cat reads.fasta | python <this>
-----------------------------------------------
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

import sys, re
from Bio import SeqIO

# constants
bad_aa_pattern = "[^ACDEFGHIKLMNPQRSTVWY]"

# read fasta dna from stdin / write fasta protein to stdout
for record in SeqIO.parse( sys.stdin, "fasta" ):
    header = str( record.name )
    frames = []
    # 3 translations from the forward direction
    for i in range( 3 ):
        frames.append( ["+%d" % ( i ), str( record.seq[i:].translate() )] )
    # 3 translations from the reverse direction
    for i in range( 3 ):
        frames.append( ["-%d" % ( i ), str( record.seq.reverse_complement()[i:].translate() )] )
    # check for valid peptides
    for index, seq in frames:
        if not re.search( bad_aa_pattern, seq ):
            """
            # this option would append the reading frame to the read name
            e.g. ">read123" -> ">read123_+0", ">read123_+1", etc.
            print ">%s_%s" % ( header, index )
            """
            print ">%s" % ( header )
            print seq
