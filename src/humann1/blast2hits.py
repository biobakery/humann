#!/usr/bin/env python

"""
Description: Takes .txt files of blastx, mblastx, or mapx alignments, and produces a pickled file of the hits.
Point in the pipeline: Raw TXT input --> 00
Program called before: None (blastx, mblastx, or mapx outside of HUMAnN)
Program called after: hits2enzymes.py, hits2metacyc.py, and/or hits2metarep.py
"""

import hits # Part of HUMAnN, located in src dir
import sys

c_strID		= "%identical"

strType = "blastx" if ( len( sys.argv ) <= 1 ) else sys.argv[1] # strType is the first arguement given in the function call.
if strType in ("-h", "-help", "--help"): # If user input flags for help, return a usage message
	raise Exception( "Usage: blast2hits.py [type={blastx,mblastx,mapx}] [filter] < <blast.txt>" )
dFilter = 0 if ( len( sys.argv ) <= 2 ) else float(sys.argv[2]) # dFilter is set equal to the float of the [filter] arguement inputted to the program  if it exists, else 0.

pHits = hits.CHits( ) # Instantiate a CHits object (see src/hits.py)
iID = 2 if ( strType == "blastx" ) else ( 4 if ( strType == "mblastx" ) else None ) # iID records the type of blast in a number: 2 = "blastx", 4 = "mblastx", or none if no strType is given.
for strLine in sys.stdin:  # Open the file passed to stdin, loop through the lines
	astrLine = strLine.rstrip( ).split( "\t" ) # For each line, split it into columns along tabs and remove whitespace.
	if not astrLine[0]:
		continue
	if astrLine[0].startswith( "#" ):
		if iID == None: # If there is no blastx vs mblastx identifier given,
			for i in range( len( astrLine ) ): # Loop through all the columns in this line.
				if astrLine[i] == c_strID: # If any column contains the c_strID defined above in this file,
					iID = i # Then the blastx vs mblastx identifier is equal to the index of that particular column.
					break # Thus, if c_strID appears in the third column, then it is a blastx file. If it appears in the fifth column, it is an mblastx file.
		continue
	if iID == None: # If the line does not start with a pound sign and no blastx vs mblastx identifier is given, then skip this line.
		continue
	try:
		if strType == "mblastx":
			strTo, strFrom, strID, strE, strCov = (astrLine[1], astrLine[0], astrLine[iID],
				astrLine[2], astrLine[5]) # Read the strings representing To, From, ID, E value, and Coverage from the appropriate columns for an mblastx file.
		elif strType == "mapx":
			strTo, strFrom, strID, strE, strCov = (astrLine[0], astrLine[2], astrLine[iID],
				astrLine[-1], astrLine[iID + 1]) # Read the strings representing To, From, ID, E value, and Coverage from the appropriate columns for an mapx file.
		else:
			strTo, strFrom, strID, strE, strCov = (astrLine[1], astrLine[0], astrLine[iID],
				astrLine[-2], astrLine[3]) # Read the strings representing To, From, ID, E value, and Coverage from the appropriate columns for an blastx file.
	except IndexError:
		sys.stderr.write( "%s\n" % astrLine )
		continue
	try:
		dE, dID, dCov = (float(s) for s in (strE, strID, strCov)) # Cast the E-value, ID, and Coverage strings to floats.
	except ValueError:
		continue # If there is a problem casting these values to floats, skip this line.
	if dFilter > 0:
		if ( dID / 100 ) >= dFilter:
			continue
	pHits.add( strTo, strFrom, dE, dID, dCov ) # Add To, From, E-value, ID, and Coverage to the pHits object.
pHits.save( sys.stdout ) # Pickle the member variables associated with pHits, and pass them to stdout.
