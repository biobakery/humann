#!/usr/bin/env python
 
import sys,string
import os
from pprint import pprint
import sys, os
from pprint import *
import math
import pdb
import tempfile 
import argparse



#********************************************************************************************
#    Map Pathways to Uniprot IDs  and Uniref50 and 90                                       *
#                                                                                           *
#    The objective of this program is to map  Swissprot Pathways to Uniprot ACs             *
#    and Uniref50 and 90                                                                    *
#                                                                                           *
# This program reads the Swissprot pathways file                                            *
# /n/huttenhower_lab_nobackup/downloads/uniprot_pathways/2014_10/pathway.txt                *
# that looks as follows:                                                                    *
#****                                                                                       *
# Alkaloid biosynthesis; 3alpha(S)-strictosidine biosynthesis; 3alpha(S)-strictosidine      *
#     STS1_ARATH  (P94111)    , STS3_ARATH  (P92976)    , STSY_CATRO  (P18417)    ,         *
#     STSY_RAUMA  (P68174)    , STSY_RAUSE  (P68175)                                        *
#Alkaloid biosynthesis; ajmaline biosynthesis                                               *
#     PNAE_RAUSE  (Q9SE93)                                                                  *
#****                                                                                       *
#And builds the relations: AC --> Reaction and Reaction --> AC                              *
#It also builds an extract file controlled by the parameter --o ValidACs                    *
#which contains a list of the ACs that were output                                          *
#*This means that if all files need to be generated, the first step must be run first*      *
#                                                                                           *
#At this point,  it generates the unipathway_pathways file and it can complete here.        *
#However, it has the option to generate also a file with relations: Reaction-->Uniref50and90*
#If so,  it proceeds to read the Uniref50 and Uniref90 files in the same fashion as         * 
#ReadSwisport.py does (Unzip the 50, 90 files, glue them) and treats, like in the case of   *
#ReadSwissprot.py,  the AC table, as a transaction file and runs a Transaction vs.          *
#Master process (AC Table vs. U5090 file of 80 million recs)  and this way updates the U50  *
#and U90 for the particular AC and generates                                                *
#the extract:  Reaction, AC{s}, U50{s},U90{s}                                               *
#                                                                                           *
#                                                                                           *
#  -----------------------------------------------------------------------------------------*
#  Invoking the program:                                                                    *
#  ---------------------                                                                    *
#  python Build_mapping_Pathways_Uniprot.py 
# --i /n/huttenhower_lab_nobackup/downloads/uniprot_pathways/2014_10/pathway.txt \
# --uniref50gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz\
# --uniref90gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef90.dat.gz\
# --oPathwaysACs  unipathway_pathways \
# --oValidACs  ../list_of_ACs \
# --oPathwaysUniref5090 PathwaysUniref5090                                            
#                                                                                           *
#   Where:                                                                                  *
#    --i_reactions, is the pathways  file, which is currently located at                    *
#    /n/huttenhower_lab_nobackup/downloads/uniprot_pathways/2014_10/pathway.txt             *
#  and it was downloaded from  the site:  http://www.uniprot.org/help/pathway               *
#                                                                                           *
#   --uniref50gz and --uniref90gz are the Uniref50 and 90 mappings,                         *
#      Currently located at /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz and 90
#                                                                                           *
#   --oPathwaysACs  is the Output file containing the Pathways --> ACs relations            *
#                                                                                           *
#  --oPathwaysUniref5090 is the output file containing the Pathways --> Uniref50/90 relations
#     ****NOTE****  If this parameter is not supplied,  this file is not created            *
#  
#   Written by George Weingart  Oct. 20, 2014   george.weingart@gmail.com                   *
#********************************************************************************************



#*************************************************************************************
#* Parse Input parms                                                                 *
#*************************************************************************************
def read_params(x):
	CommonArea = dict()	
	parser = argparse.ArgumentParser(description='Build mapping Pathways to Uniprot Ids')
	parser.add_argument('--oPathwaysACs', action="store", dest='oPathwaysACs',nargs='?')
	parser.add_argument('--oPathwaysUniref5090', action="store", dest='oPathwaysUniref5090',nargs='?')
	parser.add_argument('--i', action="store", dest='i',nargs='?')
	parser.add_argument('--uniref50gz', action="store", dest='Uniref50gz',nargs='?')
	parser.add_argument('--uniref90gz', action="store", dest='Uniref90gz',nargs='?')
	parser.add_argument('--oValidACs', action="store", dest='oValidACs',nargs='?')

	
	CommonArea['parser'] = parser
	return  CommonArea
 
 
#**************************************************************
#*   Build the second database                                *
#**************************************************************
def BuildPathwaysToUniref(CommonArea):
   	#***************************************
	# Processing Uniref50 90 files         *
	#***************************************
	dInputFiles =  InitializeProcess(CommonArea['Uniref50gz'] ,  CommonArea['Uniref90gz'])  # Invoke initialization
	CommonArea['strInput5090'] = dInputFiles["File5090"]		#Name of the Uniref5090 file
	print("Starting the load of the 5090 table\n")
	
	CommonArea = TxnVsMaster(CommonArea)   #Process the 5090 file against the ACs
 
	cmd_remove_tempdir = "rm -r /" + dInputFiles["TempDirName"]		# Remove the temporary directory
	os.system(cmd_remove_tempdir)
	 
	OutputFile2 = open(CommonArea['oPathwaysUniref5090'],'w') 

	lSortedPathwasyACs = sorted(CommonArea['dPathwaysACs'].iteritems())
	for entPathwayACs in  lSortedPathwasyACs:  
		Pathway = entPathwayACs[0]
		lACs = entPathwayACs[1]
		Pathway = Pathway.replace("\t"," ")
		lOutputRecord = [Pathway]
		lU50 = list()
		lU90 = list()
		for AC in lACs:
			try:
				U50 = "Uniref50_" + CommonArea['dUniprotUniref'][AC][0]
				lU50.append(U50)
				U90 = "Uniref90_" + CommonArea['dUniprotUniref'][AC][1]
				lU90.append(U90)
			except:
				continue
		

		if len(lU50) > 0 and len(lU90) > 0:
			sU50 = set(lU50)
			sU90 = set(lU90)
			lU50 = list(sU50)
			lU90 = list(sU90)
			for U50 in sorted(lU50):
				lOutputRecord.append(U50)
			for U90 in sorted(lU90):
				lOutputRecord.append(U90)			


			strBuiltRecord = "\t".join(lOutputRecord) + 	"\n"
			OutputFile2.write(strBuiltRecord )
			
	OutputFile2.close()	
	return CommonArea

	
	
	
#**************************************************************
#*           Match Txn vs Master                              *
#**************************************************************
def TxnVsMaster(CommonArea):

	dUniprotUniref = dict()
	FlagEnd = False
	iTxnIndex = 0	
	iTotalACsProcessed = 0	
	iPrintAfterACReads = 1000   							# Print status after processing this number of ACs	
	iTotalUniref5090RecsRead = 0							# Counter of Uniref5090 recs Read
	iPrintAfterReads = 1000000 								# Print status after read of these number of records
	
	
	lACs = list()
	for AC in CommonArea['sTableACs']:
		lACs.append(AC) 
	CommonArea['lACsSorted'] = sorted(lACs)
	CommonArea['File5090'] = open(CommonArea['strInput5090'])							# Open the file
	MasterLine = CommonArea['File5090'].readline()						# Read the Line
	MKey = MasterLine.split()[0]							# This is the Master Key
	TxnKey = CommonArea['lACsSorted'][0]
	

	while  FlagEnd == False:
		if MKey == TxnKey:
			iTxnIndex = iTxnIndex + 1
			if iTxnIndex >=  len(CommonArea['lACsSorted']) - 1:
				FlagEnd = True
			iTotalACsProcessed+=1							# Count the reads
			if  iTotalACsProcessed %  iPrintAfterACReads == 0:	# If we need to print status
				print("Total of " + str(iTotalACsProcessed) + " ACs Processed against the Uniref5090 file")	
			TxnKey = CommonArea['lACsSorted'][iTxnIndex]
			lInputLineSplit = MasterLine.split() 				# Split the line using space as delimiter
			lEnt5090 = list()									# Initialize list
			lEnt5090.append(lInputLineSplit[1].split("_")[1])	# Entry is of the form UniRef50_Q2FWP1 - We need only the Q2FWP1
			lEnt5090.append(lInputLineSplit[3].split("_")[1])	# Entry is of the form UniRef50_Q2FWP1 - We need only the Q2FWP1
			dUniprotUniref[lInputLineSplit[0]] = lEnt5090		# Post it in the dictionary
			continue
		elif  TxnKey > MKey:
			MasterLine = CommonArea['File5090'].readline() 
			if not MasterLine: 
				FlagEnd = True
			iTotalUniref5090RecsRead+=1							# Count the reads
			if  iTotalUniref5090RecsRead %  iPrintAfterReads == 0:	# If we need to print status
				print("Total of " + str(iTotalUniref5090RecsRead) + " Uniref5090 records read")
			MKey = MasterLine.split()[0]
			continue
		elif  TxnKey < MKey:
			iTxnIndex = iTxnIndex + 1
			if iTxnIndex >=  len(CommonArea['lACsSorted']) -1:
				FlagEnd = True
			TxnKey = CommonArea['lACsSorted'][iTxnIndex]
			iTotalACsProcessed+=1							# Count the reads
			if  iTotalACsProcessed %  iPrintAfterACReads == 0:	# If we need to print status
				print("Total of " + str(iTotalACsProcessed) + " ACs Processed against the Uniref5090 file")	
			continue
			
			
	CommonArea['dUniprotUniref'] = dUniprotUniref
	return CommonArea



 

#*************************************************************************************
#* Map Pathways to UniprotIDs                                                        *
#*************************************************************************************
def Map_Pathways_to_UniprotIDs(CommonArea):
	dPathwaysACs = dict()
	lTableACs = list()
	OutputFile = open(CommonArea['oACsFile'],'w')
	InputFile = open(CommonArea['iFile'])
	dPathwayUniprotAc = dict()
	lOutputLine = list()
	sCurrentPathway = None
	for iLine in InputFile: 
		iLine = iLine.rstrip('\n')
		RC = FilterLine(iLine)
		if RC < 0:
			continue
		if iLine[1] != " ":
			if len(lOutputLine) > 0:
				strBuiltRecord = "\t".join(lOutputLine) + 	"\n"
				strBuiltRecord = strBuiltRecord.replace (" ", "_")   #Modified 20141216
				strBuiltRecord = strBuiltRecord.replace (";_", ";")  #Modified 20141216
				OutputFile.write(strBuiltRecord )
				iLine = iLine.replace("\t"," ",3)
				sCurrentPathway = iLine
			lOutputLine = [iLine]
			if  sCurrentPathway is not None:
				dPathwaysACs[sCurrentPathway] = list()
		else:
			lEntries = iLine.split()
			for sEntry in lEntries:
				if sEntry.startswith("("):
					AC = sEntry.replace("(","").replace(")","")
					lOutputLine.append(AC)
					if sCurrentPathway is not None:
						dPathwaysACs[sCurrentPathway].append(AC)
						lTableACs.append(AC)
	else:
		if len(lOutputLine) > 0:
			strBuiltRecord = "\t".join(lOutputLine) + 	"\n"
			strBuiltRecord = strBuiltRecord.replace (" ", "_")   #Modified 20141216
			strBuiltRecord = strBuiltRecord.replace (";_", ";")  #Modified 20141216
			OutputFile.write(strBuiltRecord )
			
	CommonArea['dPathwaysACs'] = dPathwaysACs
	sTableACs = set(lTableACs) 	#*** Remove dups of ACs
	CommonArea['sTableACs'] = sTableACs
	lListOfACs = list(sTableACs)
	print("There are " + str(len(lListOfACs)) + " unique AC entries in the set of ACs\n")
	

	
	InputFile.close()
	OutputFile.close()
	OutputValidACsFile = open(CommonArea['oValidACs'],'w')    #Create the file of valid ACs to be used in ReadSwissprot.py
	for ValidAC in CommonArea['sTableACs']:
	    OutRecValidAC = str(ValidAC) + "\n"
	    OutputValidACsFile.write(OutRecValidAC)
	OutputValidACsFile.close()
	return CommonArea 


#********************************************************************************************
#*   Initialize the process                                                                 *
#********************************************************************************************
def InitializeProcess(strUniref50gz,  strUniref90gz):

	dInputFiles = dict()									# Initialize the dictionary
	dInputFiles["Uniref50gz"] = strUniref50gz				# Store 1st file name in dictionary
	dInputFiles["Uniref90gz"] = strUniref90gz				# Store 2nd file name in dictionary

	strTempDir = tempfile.mkdtemp()							# Make temporary folder to work in
	dInputFiles["TempDirName"] = strTempDir					# Store the name of the temp dir for future use
	cmd_chmod = "chmod 755 /" + strTempDir					# Change permissions to make usable 
	os.system(cmd_chmod)									# Invoke os
	strUniref50gzFileName = os.path.split(strUniref50gz)[1]
	strUniref90gzFileName = os.path.split(strUniref90gz)[1]
	print("Unzipping uniref50 file")
	cmd_gunzip = "gunzip -c " + strUniref50gz + ">" + strTempDir + "/" + strUniref50gzFileName[:-3] # Build the gunzip command
	os.system(cmd_gunzip)									# Invoke os
	print("Unzipping uniref90 file")
	cmd_gunzip = "gunzip -c " + strUniref90gz + ">" + strTempDir + "/" + strUniref90gzFileName[:-3] # Build the gunzip command
	os.system(cmd_gunzip)									# Invoke os
	print("Pasting Uniref50 to Uniref90")

	cmd_paste =  "paste " +  strTempDir + "/" + strUniref50gzFileName[:-3] + " " +\
						strTempDir + "/" + strUniref90gzFileName[:-3] + ">" +\
						strTempDir + "/" + strUniref50gzFileName[:-7] +  "90"    # Paste the two files together
	os.system(cmd_paste )									# Invoke os
	dInputFiles["File5090"] = strTempDir + "/" + strUniref50gzFileName[:-7] +  "90"  #Post the file created into the Common Area
	return dInputFiles


#*************************************************************************************
#* Filter the records                                                                *
#*************************************************************************************
def FilterLine(iLine):
		if len(iLine) == 0\
		or iLine.startswith("----")\
		or iLine.find("UniProt - Swiss-Prot Protein Knowledgebase") > -1 \
		or iLine.find("SIB Swiss Institute of Bioinformatics; Geneva, Switzerland") > -1\
		or  iLine.find("European Bioinformatics Institute (EBI); Hinxton, United Kingdom") > -1\
		or iLine.find("Protein Information Resource (PIR); Washington DC, USA") > -1\
		or iLine.find("Description: PATHWAY comments: index") > -1\
		or iLine.find("Name:        pathway.txt") > -1\
		or  iLine.find("Release:") > -1\
		or iLine.find("Copyrighted by the UniProt Consortium") > -1\
		or iLine.find("Distributed under the Creative Commons Attribution-NoDerivs License") > -1:
			return -1
		return 0





#*************************************************************************************
#*  Main Program                                                                     *
#*************************************************************************************
print("Program started")

CommonArea = read_params( sys.argv )  # Parse command  
parser = CommonArea['parser'] 
results = parser.parse_args()
CommonArea['oACsFile'] = results.oPathwaysACs
CommonArea['oPathwaysUniref5090'] = results.oPathwaysUniref5090
CommonArea['iFile'] = results.i 
CommonArea['Uniref50gz'] = results.Uniref50gz
CommonArea['Uniref90gz'] = results.Uniref90gz
CommonArea['oValidACs'] = results.oValidACs

 
CommonArea = Map_Pathways_to_UniprotIDs(CommonArea)
 
#**************************************************************
#*   Build the second database                                *
#**************************************************************
if CommonArea['oPathwaysUniref5090'] is not None:
	CommonArea = BuildPathwaysToUniref(CommonArea)


print("Program ended Successfully")
exit(0)
