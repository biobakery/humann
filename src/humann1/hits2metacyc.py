#!/usr/bin/env python

import hits
import math
import sys

if len( sys.argv ) < 2:
	raise Exception( "Usage: hits2metacyc.py <mcc> [topn] < <hits.bin>" )
strMCC = sys.argv[1]
iTopN = -1 if ( len( sys.argv ) <= 2 ) else int(sys.argv[2])

hashCCM = {}
hashMCC = {}
for strLine in open( strMCC ):
	astrLine = strLine.rstrip( ).split( "\t" )
	strID, astrGenes = astrLine[0], astrLine[1:]
	hashMCC[strID] = astrGenes
	for strGene in astrGenes:
		hashCCM.setdefault( strGene, set() ).add( strID )

pHits = hits.CHits( )
pHits.open( sys.stdin )
hashGenes = {}
for iFrom in range( pHits.get_froms( ) ):
	aiScores = pHits.get_scores( iFrom )
	if iTopN > 0:
		aiScores = sorted( aiScores, lambda iOne, iTwo: cmp( pHits.get_score( iOne ), pHits.get_score( iTwo ) ) )
		aiScores = aiScores[:iTopN]
# Keep only hits that correspond to at least one reaction
	aiScores = filter( lambda i: hashCCM.get( pHits.get_to( pHits.get_scoreto( i ) ) ), aiScores )
	astrTos = [pHits.get_to( pHits.get_scoreto( i ) ) for i in aiScores]
	adScores = [math.exp( -pHits.get_dic( i )[0] ) for i in aiScores]
	dSum = sum( adScores )
	for i in range( len( astrTos ) ):
		strTo, dCur = (a[i] for a in (astrTos, adScores))
		hashGenes[strGene] = ( dCur / dSum ) + hashGenes.get( strTo, 0 )

hashScores = {}
dSum = 0
for strMC, astrGenes in hashMCC.items( ):
	dScore = 0
	for strGene in astrGenes:
		dScore += hashGenes.get( strGene, 0 )
	if dScore > 0:
		hashScores[strMC] = dScore
		dSum += dScore

print( "GID	Abundance" )
for strMC, dScore in hashScores.items( ):
	print( "\t".join( (strMC, str(dScore)) ) ) # / dSum)) ) )
