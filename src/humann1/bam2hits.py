#!/usr/bin/env python

"""
Description: Takes BAM files (produced from Bowtie2 or bwa alignments) and produces a pickled file of the hits.
Point in the pipeline: Raw BAM input --> 00
Program called before: None (Bowtie2 or bwa outside of HUMAnN)
Program called after: hits2enzymes.py, hits2metacyc.py, and/or hits2metarep.py
"""

import hits # Part of HUMAnN, located in src dir
import sys
import math

pHits = hits.CHits( ) # Instantiate a CHits object (see src/hits.py)
for strLine in sys.stdin: # Loop through the lines of the .bam file (bowtie) provided as input
	astrLine = strLine.rstrip( ).split( "\t" ) # Split each line into an array of strings, spliting along the tabs
	try:
		strTo, strFrom, strE = ( astrLine[2], astrLine[0], astrLine[4] ) # load strTo, strFrom and strE from the third, first, and fifth columns of the input file
	except IndexError: 
		sys.stderr.write( "%s\n" % astrLine ) # If there aren't at least 5 columns in the line, print the line's text to console and skip it 
		continue
	try:
		dE = math.pow( 10.0, ( float( strE ) )/( -10.0 ) ) # Local vs. Global
	except ValueError:
		continue
	pHits.add( strTo, strFrom, dE, 1, 1 ) # Adds strTo, strFrom, and dE to the appropriate member variable arrays and hashtables of the pHits object
pHits.save( sys.stdout ) # Pickle all of the gathered data (see hits.py) to sys.stdout
