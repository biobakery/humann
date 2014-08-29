#!/usr/bin/env python

import os
import subprocess
import sys
import tempfile

if len( sys.argv ) != 3:
	raise Exception( "Usage: enzymes2pathways_mp.py <minpath.py> <mapfile> < <enzymes.txt>" )
strMP, strMap = sys.argv[1:]

hashhashAbs = {}
astrComments = []

for strLine in sys.stdin:
	if strLine and ( strLine[0] == "#" ):
		astrComments.append( strLine )
		continue
	astrLine = strLine.strip( ).split( "\t" )
	if ( astrLine[0] == "GID" ):
		fOrg = len(astrLine) > 2
		continue
	strOrg = astrLine[1] if fOrg else None
	dScore = float( astrLine[2] ) if fOrg else float( astrLine[1] )
	hashhashAbs.setdefault( strOrg, {} )[astrLine[0]] = d = float( astrLine[2] ) if fOrg else float( astrLine[1] )

hashhashPaths = {}
for strOrg, hashAbs in hashhashAbs.items( ):
	iIn, strIn = tempfile.mkstemp( )
	for strID, dAb in hashAbs.items( ):
		os.write( iIn, "%s	%g\n" % (strID, dAb) )
	os.close( iIn )
	iOut, strOut = tempfile.mkstemp( )
	iTmp, strTmp = tempfile.mkstemp( )
	os.close( iOut )
	subprocess.call( [strMP, "-any", strIn, "-map", strMap,
		"-report", "/dev/null", "-details", strOut, "-mps", strTmp],
		env = {"MinPath": os.path.dirname( strMP )}, stdout = sys.stderr )
	os.unlink( strIn )

	strPath = None
	for strLine in open( strOut ):
		astrLine = strLine.strip( ).split( " " )
		if strLine[0:4] == "path":
			strPath = astrLine[7]
		else:
			hashhashPaths.setdefault( strOrg, {} ).setdefault( astrLine[0], [] ).append( strPath )
	os.unlink( strOut )
	os.close( iTmp )

sys.stdout.write( "GID	" )
if fOrg:
	sys.stdout.write( "Organism	" )
print( "Pathway	Abundance" )
sys.stdout.write( "".join( astrComments ) )
for strOrg, hashAbs in hashhashAbs.items( ):
	hashPaths = hashhashPaths.get( strOrg ) or {}
	for strID, dAb in hashAbs.items( ):
		astrPaths = hashPaths.get( strID ) or [""]
		strAb = str(dAb) # / len( astrPaths ))
		if fOrg:
			for strPath in ( astrPaths or [""] ):
				print( "\t".join( [strID, strOrg, strPath, strAb] ) )
		else:
			for strPath in ( astrPaths or [""] ):
				print( "\t".join( [strID, strPath, strAb] ) )
