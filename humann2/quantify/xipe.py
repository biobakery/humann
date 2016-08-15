#!/usr/bin/env python

# Author: Beltran Rodriguez-Mueller
# Early-Late-2010 Version (Oct-2)
# License: GPL Version 3
# See the text of the license at:
# http://www.gnu.org/licenses/agpl-3.0.html

# Use notes:
# the quick and easy way to use it is:
# xipe --file1 aaa.txt --file2 bbb.txt -outfile ccc.txt --nrepetitions 1000 --samplesize 1000
# additionally the --paranoid 1 
# flag can be used to "be paranoid" and do everything in triplicate and take
# the worst confidence level in the outcome.

# input file format
# category1 (tab) count (enter)
# category2 (tab) count (enter)
# ...

# output file format
# category (confidence level, 1 or 2 )
# 1 -> file 1
# 2 -> file 2
# say, 
# category1 (98,1)
# category1 is overrepresented in file 1 with .98

# Additional note:
# dealing with a set of categories that has a low abundance
# element, it might be easier to just duplicate the file as is,
# and set the category to zero to test the confidence level at
# which it can be told from zero.
# A couple of functions in this code are involved in that process:
# writeSamples
# and
# doManyCases
# doManyCases right now has quite a few elements hardcoded, so 
# it might need some tweaking before becoming a generally useful
# function, but its functionality can easily be replicated with
# a wrapper on doOneCase.
# doManyCases was mainly employed to generate contour plots.
import random
import sys
import time
import copy
from optparse import OptionParser
random.seed()
isDebug=False
# functions
# Xipe			 # Xipe3test
# readSample	   #
# homogenizeDicts  # testHomogenizeDicts
def writeSamples( outfile1, outfile2, ncats, catsize, weird, weirdsize ):
	"""writes two files, differing only on the value of one category"""
	try:
		f = open( outfile1, "w" )
	except IOError:
		print("For some reason I cannot open "+outfile1)
		print("Check directory permissions and/or file name")
		exit(-1)
	try:
		g = open( outfile2, "w" )
	except:
		print("For some reason I cannot open "+outfile1)
		print("Check directory permissions and/or file name.")
		exit(-1)
	for i in range(ncats):
		f.write("cat_"+str(i)+"\t"+str(catsize)+"\n" )
		g.write("cat_"+str(i)+"\t"+str(catsize)+"\n" )
	f.write( weird+"\t0\n" )
	g.write( weird+"\t"+str(weirdsize)+"\n" )
	f.close()
	g.close()
	
def readSample( infile, dictOne = None ):
	"""reads an input file in the: category \t count \n format"""
	resDict = {}
	if not infile:
		return resDict
	if dictOne != None:
		try:
			dPerc = float(infile)
			aaOne = sorted( dictOne.items( ), key = lambda a: a[1] )
			i = int(round( dPerc * len( aaOne ) ))
			for a in aaOne[:i]:
				sys.stderr.write( "Removing:	%s\n" % a[0] )
			for a in aaOne[i:]:
				resDict[a[0]] = a[1]
			return resDict
		except ValueError:
	   		pass
	try:
		f = infile if ( type( infile ) == file ) else open( infile )
	except IOError:
		print("For some reason I cannot open "+infile)
		print("Check directory permissions and/or file name.")
		exit(-1)
	for rawline in f:
		line = rawline.strip()
		if line[0] != "#":
			pieces = line.split("\t")
			if len(pieces) < 2:
				continue
			theKey = pieces[0].upper()
			try:
				theVal = float( pieces[1] )
			except ValueError:
				continue
			if theVal > 0:
				resDict[theKey] = theVal
	f.close()
	return resDict
def homogenizeDicts( dict1, dict2, theKeys=0 ):
	new1 = {}
	new2 = {}
	if theKeys == 0:
		allKeys = list( set ( dict1.keys() + dict2.keys() ) )
	else:
		allKeys = theKeys
	for key in allKeys:
		if key in dict1:
			new1[key] = dict1[key]
		else:
			new1[key] = 0
		if key in dict2:
			new2[key] = dict2[key]
		else:
			new2[key] = 0
	return new1, new2
def homogenizeDictsN( allDicts, allKeys=0 ):
	"""Homogenizes more than 3 dictionaries
	
	Not in use right now, for polyxipe"""
	ndicts = len(allDicts)
	if allKeys == 0:
		allKeys = []
		for dict in allDicts:
			allKeys = allKeys + dict.keys()
		allKeys = list( set( allKeys ) )
	newDicts = []
	for dict in allDicts:
		thisDict = {}
		for key in allKeys:
			if key in dict:
				thisDict[key] = dict[key]
			else:
				thisDict[key] = 0
		newDicts.append( thisDict )
	return newDicts
		
def createMix( dict1, dict2 ):
	"""Takes two dictionaries, and makes sure that they share all the keys"""
	resDict = {}
	for key in dict1:
		if key in resDict:
			resDict[key] += dict1[key]
		else:
			resDict[key] = dict1[key]
	for key in dict2:
		if key in resDict:
			resDict[key] += dict2[key]
		else:
			resDict[key] = dict2[key]
	return resDict
def sample( inDict, ss ):
	"""Takes a list of elements and samples ss elements"""

	theSample = []
	if not inDict:
		return theSample
	dSum = reduce( lambda dSum, d: dSum + d, inDict.values( ) )
	while len(theSample) < ss:
		dTarget = random.random( ) * dSum
		d = 0
		for strCur, dCur in inDict.items( ):
			d += dCur
			if d >= dTarget:
				break
		theSample.append( strCur )
	return theSample
def sampleMix( inDict1, inDict2, sampleSize ):
	outDict = {}
	for curDict in (inDict1, inDict2):
		for strKey, dValue in curDict.items( ):
			d = outDict.get( strKey, 0 )
			outDict[strKey] = d + dValue
	return sample( outDict, sampleSize )
def getCountDict( inlist ):
	"""list -> count dict"""
	outCD = {}
	for element in inlist:
		if element in outCD:
			outCD[element] += 1
		else:
			outCD[element] = 1
	return outCD

def mainStep1( sample1, sample2, sampleSize, nreps ):
	sample1h,sample2h = homogenizeDicts( sample1, sample2)
	allKeys = sample1h.keys()

	listDeltas = []
	baseDeltas = []
	for i in range(nreps):
		list1x = sample( sample1, sampleSize )
		list2x = sample( sample2, sampleSize )
		l1xcd = getCountDict( list1x )
		l2xcd = getCountDict( list2x )
		l1xcd2, l2xcd2 = homogenizeDicts( l1xcd, l2xcd, allKeys )
		# the loop below initializes the delta arrays
		# it should NOT be necessary, but I'm out of ideas
		# as for why the code was not working so :/
		deltas1 = {}
		deltas2 = {}
		for key in allKeys:	
			deltas1[key] = 0
			deltas2[key] = 0
		for key in deltas1:
			deltas1[ key ] = l1xcd2[key] - l2xcd2[key]
		listDeltas.append( deltas1 )
		# ^ appends each delta to the main list of dicts
		# the print below is for debugging purposes.
		# print i,listDeltas
		
		b1x = sampleMix( sample1, sample2, sampleSize )
		b2x = sampleMix( sample2, sample1, sampleSize )
		b1xcd = getCountDict( b1x )
		b2xcd = getCountDict( b2x )
		b1xcd2, b2xcd2 = homogenizeDicts( b1xcd, b2xcd, allKeys )
		for key in deltas2:
			deltas2[ key ] = b1xcd2[key] - b2xcd2[key]
		baseDeltas.append( deltas2 )
		# ^ appends the space delta to the main list of space delta dicts
		if isDebug:
			print("This iteration ")
			for key in allKeys:
				outStr = key + "\t" + str(l1xcd2[key]) + " -" + str( l2xcd2[key])
				outStr+= " = " + str(deltas1[key]) + "\t||\t"
				outStr+= str(b1xcd2[key]) + " - " + str(b2xcd2[key])
				outStr+= " = " + str(deltas2[key])
				print(outStr)
		# the print below is for debugging purposes.
		# print i,"-->",listDeltas
		
	#print "End Results ====== "
	returnDict = {}
	for key in allKeys:
		returnList = [ (0,0), (0,0) ]
		thisSpace = []
		thisMedian = []
		for i in range( nreps ):
			thisSpace.append( baseDeltas[i][key ] )
			thisMedian.append( listDeltas[i][key] )
		thisSpace.sort(reverse=True)
		thisMedian.sort(reverse=True)
		#print "key = ",key," === "
		#print thisSpace
		#print " ^^^ "
		#print thisMedian
		if nreps % 2 == 1:
			theMedian = thisMedian[int(round(nreps/2))]
		else:
			theMedian = thisMedian[int(round(nreps/2))]*.5 +  thisMedian[int(round(nreps/2))+1]*.5
		the50a = thisSpace[ int(round(nreps/100)*25 ) ]
		the50b = thisSpace[ int(round(nreps/100)*75 ) ]
		the60a = thisSpace[ int(round(nreps/100)*20 ) ]
		the60b = thisSpace[ int(round(nreps/100)*80 ) ]
		the70a = thisSpace[ int(round(nreps/100)*15 ) ]
		the70b = thisSpace[ int(round(nreps/100)*85 ) ]
		the80a = thisSpace[ int(round(nreps/100)*10 ) ]
		the80b = thisSpace[ int(round(nreps/100)*90 ) ]
		the90a = thisSpace[ int(round(nreps/100)*5 ) ]
		the90b = thisSpace[ int(round(nreps/100)*95 ) ]
		the91a = thisSpace[ int(round(nreps/100)*4.5 ) ]
		the91b = thisSpace[ int(round(nreps/100)*95.5 ) ]
		the92a = thisSpace[ int(round(nreps/100)*4 ) ]
		the92b = thisSpace[ int(round(nreps/100)*96 ) ]
		the93a = thisSpace[ int(round(nreps/100)*3.5 ) ]
		the93b = thisSpace[ int(round(nreps/100)*96.5 ) ]
		the94a = thisSpace[ int(round(nreps/100)*3 ) ]
		the94b = thisSpace[ int(round(nreps/100)*97 ) ]
		the95a = thisSpace[ int(round(nreps/100)*2.5 ) ]
		the95b = thisSpace[ int(round(nreps/100)*97.5 ) ]
		the96a = thisSpace[ int(round(nreps/100)*2 ) ]
		the96b = thisSpace[ int(round(nreps/100)*98 ) ]
		the97a = thisSpace[ int(round(nreps/100)*1.5 ) ]
		the97b = thisSpace[ int(round(nreps/100)*98.5 ) ]	
		the98a = thisSpace[ int(round(nreps/100)*1 ) ]
		the98b = thisSpace[ int(round(nreps/100)*99 ) ]
		the99a = thisSpace[ int(round(nreps/100)*.5 ) ]
		the99b = thisSpace[ int(round(nreps/100)*99.5 ) ]
		#print "99|98|97|96|95|94|93|92|91|90|80|70|60|50|o|50|60|70|80|90|91|92|92|94|95|96|97|98|99"
		strOut = ""
		if theMedian > the99a:
			strOut+= "+ |"
			returnList.append( (99,1) )
		else:
			strOut+= "  |"
		if theMedian > the98a:
			strOut+= "+ |"
			returnList.append( (98,1) )			
		else:
			strOut+= "  |"
		if theMedian > the97a:
			strOut+= "+ |"
			returnList.append( (97,1) )
		else:
			strOut+= "  |"
		if theMedian > the96a:
			returnList.append( (96,1) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian > the95a:
			returnList.append( (95,1) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian > the94a:
			returnList.append( (94,1) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian > the93a:
			returnList.append( (93,1) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian > the92a:
			returnList.append( (92,1) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian > the91a:
			returnList.append( (91,1) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian > the90a:
			returnList.append( (90,1) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian > the80a:
			returnList.append( (80,1) )
			strOut+= "+ |"
		else:
			strOut+= "  |"			
		if theMedian > the70a:
			returnList.append( (70,1) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian > the60a:
			returnList.append( (60,1) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian > the50a:
			returnList.append( (50,1) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		strOut += "o|"
		if theMedian < the50b:
			returnList.append( (50,2) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian < the60b:
			returnList.append( (60,2) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian < the70b:
			returnList.append( (70,2) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian < the80b:
			returnList.append( (80,2) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian < the90b:
			returnList.append( (90,2) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian < the91b:
			returnList.append( (91,2) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian < the92b:
			returnList.append( (92,2) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian < the93b:
			returnList.append( (93,2) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian < the94b:
			returnList.append( (94,2) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian < the95b:
			strOut+= "+ |"
			returnList.append( (95,2) )
		else:
			strOut+= "  |"
		if theMedian < the96b:
			returnList.append( (96,2) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian < the97b:
			returnList.append( (97,2) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian < the98b:
			returnList.append( (98,2) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		if theMedian < the99b:
			returnList.append( (99,2) )
			strOut+= "+ |"
		else:
			strOut+= "  |"
		#print strOut + key
		returnList.sort(reverse=True)
		# element 0 of return list is the one we actually care for
		returnDict[key] = returnList
	return returnDict
def doOneCase(file1, file2, outfile, sampleSize=1000,nreps=2000, paranoid=1):
	"""Performs analysis for a pair of files"""
	sample1 = readSample( file1 or sys.stdin )
	sample2 = readSample( file2, sample1 )
	isDebug = False
	allResults = {}
	# I am being extremely paranoid here and using the worst
	# confidence value out of three tries.
	# If speed is an issue or anybody is feeling lucky or less
	# paranoid than me, turn paranoid to False :)
	if paranoid==1:
		ntries = 3
	else:
		ntries = 1
	for i in range(ntries):
		result = mainStep1( sample1, sample2, sampleSize, nreps )
		for key in result:
			if key in allResults:
				r1 = result[key][0]
				r2 = allResults[key][0]
				if r1 < r2:
					allResults[key] = result[key][0]
			else:
				allResults[key] = result[key][0]
	allkeys = allResults.keys()
	allkeys.sort()
	f = open(outfile, "w" ) if outfile else sys.stdout
	for key in allkeys:
		f.write( key+"\t"+str( allResults[key] )+"\n" )
	f.close()
def doManyCases():
	isDebug = False
	catsize = 70
	file1 = "uniform_c"+str(catsize)+"_1.txt"
	file2 = "uniform_c"+str(catsize)+"_2.txt"
	sampleSize = 1000
	nreps = 1500
	ofi = open( "uniform_c"+str(catsize)+".txt", "w")
	for ncats in [5,7,9,11,15,22,30,40]:
		for weirdsize in range(5,31,5):
			writeSamples( file1, file2, ncats, catsize, 'XXX', weirdsize )
			resultado = mainStep1( file1, file2, sampleSize, nreps )
			r1 = resultado['XXX'][0][0]
			resultado = mainStep1( file1, file2, sampleSize, nreps )
			r2 = resultado['XXX'][0][0]
			resultado = mainStep1( file1, file2, sampleSize, nreps )
			r3 = resultado['XXX'][0][0]
			r123 = r1/3.0 + r2/3.0 + r3/3.0
			ofi.write( str(ncats*catsize)+" "+str(weirdsize)+" "+str(r123)+"\n")
			print("".join([str(i) for i in [ncats, weirdsize, resultado['XXX'][0]]]))
		ofi.write("\n")
	timeEnd = time.time()
	print("took " + str(timeEnd-timeStart) + "seconds")
	ofi.close()
if __name__=="__main__":
	parser = OptionParser()
	parser.add_option("--file1", dest="filename1",
		   help="read first data set from FILE", metavar="FILE")
	parser.add_option("--file2", dest="filename2",
		   help="read second data set from FILE", metavar="FILE")
	parser.add_option("--samplesize", dest="sampleSize",
		   default=100,
		   help="sample size, NUMBER", metavar="NUMBER" )
	parser.add_option("--nrepetitions", dest="nreps",
		   default=100,
		   help="number of repetitions, NUMBER", metavar="NUMBER")
	parser.add_option("--outfile", dest="outfile",
		   help="write output in FILE", metavar="FILE")
	parser.add_option("--paranoid", dest="paranoid", default=1,
		   help="paranoid? 1 = yes, otherwise no. Will do everything in triplicate and choose the worst confidence level for each category (defaults to yes).")
	parser.add_option("-q", "--quiet",
		   action="store_false", dest="verbose", default=True,
		   help="don't print status messages to stdout")
	options, args = parser.parse_args()
	#print "options = ", options
	#doManyCases()
	doOneCase( options.filename1, options.filename2, options.outfile,
			   int(options.sampleSize), int(options.nreps), int(options.paranoid) )
