#!/usr/bin/env python

import re
import sys

if len( sys.argv ) != 1:
	raise Exception( "Usage: metacyc2mcpc.py < <pathways.dat>" )

strID = astrGenes = None
for strLine in sys.stdin:
	mtch = re.search( '^UNIQUE-ID\s+-\s+(\S+)', strLine )
	if mtch:
		if astrGenes:
			print( "\t".join( [strID] + astrGenes ) )
		strID = mtch.group( 1 )
		astrGenes = []
		continue
	mtch = re.search( '^REACTION-LIST\s+-\s+(\S+)', strLine )
	if mtch:
		astrGenes.append( mtch.group( 1 ) )
if astrGenes:
	print( "\t".join( [strID] + astrGenes ) )
