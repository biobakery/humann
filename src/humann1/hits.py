#!/usr/bin/env python

import array
import cPickle as pickle # use the cPickle object, which is faster than standard pickle
import sys

class CHits:

	def __init__( self ):
		""" Clear and instantiate all member variables on init. """

		self._clear( )

	def _clear( self ):
		""" Clears any member variables and instantiates empty member variables of the appropriate format """
				
		self.m_hashFroms = {}
		self.m_hashTos = {}
		self.m_astrTos = []
		self.m_astrFroms = []
		self.m_apScores = []
		self.m_pTos = array.array( "L" ) # Array of long integers
		self.m_pEs = array.array( "f" ) # Array of floats
		self.m_pIDs = array.array( "f" )
		self.m_pCovs = array.array( "f" )
		
	def _enhash( self, strID, hashIDs, astrIDs, apScores = None ):
		"""
		Returns value of the index integer coresponding to strID in the hashIDs hashtable.
		
		Adds strID to the hashIDs hashtable if not present. Asigns it the next highest available integer ID.
		Also appends strID to astrIDs.
		"""
		iID = hashIDs.get( strID )
		if iID == None:
			hashIDs[strID] = iID = len( hashIDs ) # Set value of hashIDS[strID] to its index number in hashIDs
			astrIDs.append( strID ) # Append strID to astrIDs
			if apScores != None:
				apScores.append( array.array( "L" ) ) # If apScores is provided, appends empty int array 
			
		return iID # Return the index number of strID in the hashIDs hashtable

	def _repopulate( self, astrIDs, hashIDs ):
		"""
			Creates a dict mapping each ID to an integer index of its position in astrIDs.
			Input: an array of IDs (astrIDs), and a dict (hashIDs).
			Output: a dict mapping each ID to its index in astrIDs.
		"""
		
		for i in range( len( astrIDs ) ):
			hashIDs[astrIDs[i]] = i

	def add( self, strTo, strFrom, dE, dID, dCov ):
		""" Adds strTo, strFrom, dE, dID and dCov to the coresponding member variables """
		
		for astrCur, hashCur in ((self.m_astrTos, self.m_hashTos), (self.m_astrFroms, self.m_hashFroms)):
			if astrCur and ( not hashCur ):
				self._repopulate( astr, hash )

		iTo = self._enhash( strTo, self.m_hashTos, self.m_astrTos )
		iFrom = self._enhash( strFrom, self.m_hashFroms, self.m_astrFroms, self.m_apScores )
		iScore = len( self.m_pTos )
		self.m_apScores[iFrom].append( iScore )
		self.m_pTos.append( iTo )
		self.m_pEs.append( dE )
		self.m_pIDs.append( dID )
		self.m_pCovs.append( dCov )

	def get_froms( self ):
		""" Returns the size of the m_astrFroms member variable """
		
		return len( self.m_astrFroms )
	
	def get_from( self, iFrom ):
		""" Returns the string in the m_astrFroms member variable coresponding to the index iFrom. """
		
		return self.m_astrFroms[iFrom]
	
	def get_to( self, iTo ):
		""" Returns the string in the m_astrTos member variable coresponding to the index iTo. """
		
		return self.m_astrTos[iTo]
	
	def get_tos( self ):
		""" Returns the size of the m_astrTos member variable """
		
		return len( self.m_astrTos )

	def get_scoreto( self, iScore ):
		""" Returns the value associated with the index iScore from the member variable m_pTos. """

		return self.m_pTos[iScore]

	def get_scores( self, iFrom ):
		""" Returns the value associated with the index iFrom from the member variable m_apScores. """

		return self.m_apScores[iFrom]
	
	def get_dic( self, iScore ):
		"""
		Returns an array of member variable arrays, all at the score (iScore) passed in.
		Input: a score integer (iScore)
		Output: an array of tuples, [(m_pEs[iScore], m_pIDs[iScore], m_pCovs[iScore])]
		"""
		
		return [a[iScore] for a in (self.m_pEs, self.m_pIDs, self.m_pCovs)] # Returns an array of tuples
	
	def save( self, fileOut ):
		""" Pickles all of the member variables to fileOut """
	
		pickle.dump( self.m_astrFroms, fileOut, pickle.HIGHEST_PROTOCOL )
		pickle.dump( self.m_astrTos, fileOut, pickle.HIGHEST_PROTOCOL )
		pickle.dump( self.m_apScores, fileOut, pickle.HIGHEST_PROTOCOL )
		pickle.dump( self.m_pTos, fileOut, pickle.HIGHEST_PROTOCOL )
		pickle.dump( self.m_pEs, fileOut, pickle.HIGHEST_PROTOCOL )
		pickle.dump( self.m_pIDs, fileOut, pickle.HIGHEST_PROTOCOL )
		pickle.dump( self.m_pCovs, fileOut, pickle.HIGHEST_PROTOCOL )
		
	def open( self, fileIn ):
		""" Depickles any member variables stored in fileIn, asigns values to coresponding member variables of CHits object """
		
		self._clear( )
		self.m_astrFroms = pickle.load( fileIn )
		self.m_astrTos = pickle.load( fileIn )
		self.m_apScores = pickle.load( fileIn )
		self.m_pTos = pickle.load( fileIn )
		self.m_pEs = pickle.load( fileIn )
		self.m_pIDs = pickle.load( fileIn )
		self.m_pCovs = pickle.load( fileIn )

if __name__ == "__main__":
	pHits = CHits( )
	pHits.open( sys.stdin )
#	pHits.save( sys.stdout )
	for i in pHits.m_pTos: # Print all values in the pTos array (after depickling member vars from the stdin stream)
		print( i )