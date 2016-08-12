#!/usr/bin/env python
import sys,string
import os
import argparse
import tempfile 
from pprint import *
import math
from collections import defaultdict
import operator

#********************************************************************************************
#    Map ECs to Uniref50 and Uniref90                                                       *
#    Logic:                                                                                 *
#                                                                                           *
#  -----------------------------------------------------------------------------------------*
#  Invoking the program:                                                                    *
#  ---------------------                                                                    *
#   python map_ECs_to_Uniref50_90.py                                                        *
#                                                                                           *                                                                                        *
#                                                                                           *           
#                                                                                           *
#                                                                                           *
#                      P R O G R A M      F L O W                                           *
#                     ---------------------------                                           *
#    1. The program reads the Swissprot.dat and builds the relation AC->ECs                 *
#       and does the same for Trembl                                                        *
#    2. It aggregates all the ACs into a sorted list where every AC appears only once       *
#    3. It reads the u50 & u90 zipped files and pastes them together(They are sorted by AC) *
#    4. It reads in parallel the Master file (U5090) vs the TXN file (ACs)  and where there *
#       is a match for the AC in the U5090 file,  it takes the corresponding EC for that AC *
#       and generates a record: EC\t\U50,U90\tU50,U90\t....etc...\n                         *                       
#********************************************************************************************
#   Written by George Weingart - george.weingart@gmail.com   12/20/2015                     *  
#********************************************************************************************



#*************************************************************************************
#* Parse Input parms                                                                 *
#*************************************************************************************
def read_params(x):
	CommonArea = dict()	
	parser = argparse.ArgumentParser(description='Build relation: Reaction --> EC --> UniprotAC --> Uniref5090')
	parser.add_argument('--uniref50gz', action="store", dest='Uniref50gz',nargs='?', default="/n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz")
	parser.add_argument('--uniref90gz', action="store", dest='Uniref90gz',nargs='?' , default = "/n/huttenhower_lab/data/idmapping/map_uniprot_UniRef90.dat.gz")
	parser.add_argument('--o_map1', action="store", dest='o_map1',nargs='?', help='Map EC_ACs_U50U90s', default = "map_EC_to_triplet_AC_U50_U90" )
	parser.add_argument('--o_map2', action="store", dest='o_map2',nargs='?', help='Map EC_U50s' , default = "map_EC_to_U50s" )
	parser.add_argument('--iSwssiprot', action="store", dest='iSwissprot',nargs='?', default="/n/huttenhower_lab/data/uniprot/2015-06/uniprot_sprot.dat")
	parser.add_argument('--iTrembl', action="store", dest='iTrembl',nargs='?', default="/n/hutlab11_nobackup/users/gweingart/uniprotkb/uniprot_trembl/uniprot_trembl.dat")
	CommonArea['parser'] = parser
	return  CommonArea



#*************************************************************************************
#* Read Swissprot                                                                    *
#* Build the AC -->[ECs]  relationship                                               *
#*************************************************************************************
def ReadSwissprot(iFile):
	InputFile = open(iFile)
	LineCntr = 0
	iInterval = 10000000
	dACtoECs = dict()
	lStripChars = ["\n",";",".-",".-",".-"]
	for iLine in InputFile: 
			LineCntr = LineCntr +1
			if LineCntr % iInterval == 0:
				print("Total input records read = " + str(LineCntr)) 

			if iLine.startswith("AC   "):
				lTemp = iLine.rstrip().split().pop().split(";")
				lACs = [var for var in lTemp if var]
				
							
			if iLine.startswith("DE   ") and "EC=" in iLine:
				for AC in lACs:
					if AC not in dACtoECs:
						dACtoECs[AC] = list()
				EC = iLine.split("EC=")[1] 
				EC = EC.split(" ")[0]
				for sStripChar in lStripChars:
					EC = EC.rstrip(sStripChar)

				for AC in lACs:
					dACtoECs[AC].append(EC)
					
	InputFile.close()
	return  dACtoECs

	
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

#********************************************************************************************
#*   Build the cross reference: EC--> [U50,U90]                                            *
#********************************************************************************************
def Build_EC_to_Uniref5090_File(CommonArea):
    OutputFile = open(CommonArea['OutputFile'],'w') 
    OutputFile_Map2 = open(CommonArea['OutputFile_Map2'],'w')
    for EC, lACs in sorted(CommonArea['dEC_to_AC'].iteritems()):
        lECsU50 = list()
        strOutputRec = EC + "\t"
        strOutputRec_Map2 = EC + "\t"
        bFlagAcFound = False        # Set up a flag - do not print if no AC found for EC
        lU50sU90s = list()  # List of U50 and U90s for this EC
        for AC in lACs:
            if AC in CommonArea['dUniprotUniref']:
                bFlagAcFound  = True     # There is an AC for this EC
                lTriplet5090 =  [AC,CommonArea['dUniprotUniref'][AC][0],CommonArea['dUniprotUniref'][AC][1]]  #Triplet: AC +  Uniref5090 for this AC
                lU50sU90s.append(lTriplet5090)
                lECsU50.append(CommonArea['dUniprotUniref'][AC][0])   # Add the UniRef50
        lU50sU90s = sorted( lU50sU90s)   # Sort it
        for eU50U90 in  lU50sU90s: 
            strOutputRec = strOutputRec + eU50U90[0] + "\tUniRef50_" + eU50U90[1] +  "\tUniRef90_" + eU50U90[2] +"\t" 

        strOutputRec = strOutputRec + "\n"
        
 
        lECsU50 = sorted(list(set(lECsU50)))  # Get rid of dups and sort it
        for eU50 in lECsU50:
            strOutputRec_Map2 = strOutputRec_Map2 + "UniRef50_" + eU50 + "\t"
        strOutputRec_Map2 = strOutputRec_Map2 + "\n"
        
        
        
        if bFlagAcFound == True:   # If there was an AC for this EC (And therefore a U50 and 90) 
            OutputFile.write(strOutputRec )
            OutputFile_Map2.write(strOutputRec_Map2)
    
    OutputFile.close()
    OutputFile_Map2.close()
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

#****************************************************************************************************************
#*   Merge the two dictionaries: d1=AC(Swissprot)-->[ECs Swissprot] and d2=AC(Trembl)-->[ECs Trembl]            *
#****************************************************************************************************************
def Merge_Dictionaries(d1,d2):
	d = {}

	for key in set(d1.keys() + d2.keys()):
		try:
			d.setdefault(key,[]).append(d1[key])
		except KeyError:
			pass

		try:
			d.setdefault(key,[]).append(d2[key])
		except KeyError:
			pass

	
	for key, value in d.iteritems():   # Remove EC duplicates - first - flatten the list and set it...
		d[key] = list(set(reduce(operator.add, d[key])) )

	return d







	
#********************************************************************************************
#* Main                                                                                     *
#********************************************************************************************
print("Program started")

CommonArea = read_params( sys.argv )  # Parse command  
parser = CommonArea['parser'] 
results = parser.parse_args()
CommonArea['iFileSwissprot'] = results.iSwissprot
CommonArea['iFileTrembl'] = results.iTrembl

print("Reading Swissprot")
dACtoECsSwissprot = ReadSwissprot(CommonArea['iFileSwissprot']) #Read Swissprot 

print("Reading Trembl")
dACtoECsTrembl = ReadSwissprot(CommonArea['iFileTrembl'])  #Read Trembl

print("Merging the Swissprot and Trembl dictionairies")
dMerged = Merge_Dictionaries(dACtoECsSwissprot,dACtoECsTrembl) # Merge the two dictionaries


CommonArea['dEC_to_AC']  = dict()     #This is the dictionary of EC --> ACs
CommonArea['lACsSorted'] = list()  
for AC, lECs  in sorted(dMerged.iteritems()):
	CommonArea['lACsSorted'].append(AC)  #Build the Selection ACs
	for EC in lECs:
		if EC not in CommonArea['dEC_to_AC']: # And at the same time also invert the relation AC->ECs to ECs->ACs
			CommonArea['dEC_to_AC'][EC] = list()
		CommonArea['dEC_to_AC'][EC].append(AC)



CommonArea['OutputFile'] = results.o_map1
CommonArea['OutputFile_Map2'] = results.o_map2	


#***************************************
# Processing Uniref50 90 files         *
#***************************************

CommonArea['Uniref50gz'] = results.Uniref50gz
CommonArea['Uniref90gz'] = results.Uniref90gz
dInputFiles =  InitializeProcess(CommonArea['Uniref50gz'] ,  CommonArea['Uniref90gz'])  # Invoke initialization
CommonArea['strInput5090'] = dInputFiles["File5090"]		#Name of the Uniref5090 file


CommonArea = TxnVsMaster(CommonArea)   #Process the 5090 file against the ACs
cmd_remove_tempdir = "rm -r /" + dInputFiles["TempDirName"]		# Remove the temporary directory
os.system(cmd_remove_tempdir)
 

CommonArea = Build_EC_to_Uniref5090_File(CommonArea)    #  Build the cross reference

 
print("Program ended Successfully")
exit(0)
