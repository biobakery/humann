#!/usr/bin/env python
from cStringIO import StringIO
import sys,string
import os
import argparse
import tempfile 
from pprint import *
import math
 

#********************************************************************************************
#    Map ECs to Uniref50 and Uniref90                                                       *
#    Logic:                                                                                 *
#                                                                                           *
#  -----------------------------------------------------------------------------------------*
#  Invoking the program:                                                                    *
#  ---------------------                                                                    *
#   python map_ECs_to_Uniref50_90.py    --EnzymesFile  enzyme.dat --uniref50gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz --uniref90gz   /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef90.dat.gz  --o_map1 Mapping_EC_AC_U5090s --o_map2 Map_ECs_U50
#                                                                                           *                                                                                        *
#   Where:                                                                                  *
#  -EnzymesFile  is reactions.dat                                                           *
#  That file can be downloaded from:                                                        *
#   ftp://ftp.expasy.org/databases/enzyme/                                                  * 
#  and  more particularly from                                                              *                         
#  ftp://ftp.expasy.org/databases/enzyme/enzyme.dat                                         * 

#    --i_reactions, is the reactions file, which is currently located at                    *
#    /n/huttenhower_lab_nobackup/downloads/metacyc/18.1/reactions.dat                       *
#                                                                                           *           
#                                                                                           *
#    uniref50gz and uniref90gz are the uniref mappings (Uniref50 --> Uniprot AC)            * 
#     currently located at                                                                  *
#     /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz                         *
#                                                                                           *
#                      P R O G R A M      F L O W                                           *
#                     ---------------------------                                           *
#    1. The program reads the enzyme.dat file and loads a dictionary of EC --> [ACs]        *
#    2. It aggregates all the ACs into a sorted list where every AC appears only once       *
#    3. It reads the u50 & u90 zipped files and pastes them together(They are sorted by AC) *
#    4. It reads in parallel the Master file (U5090) vs the TXN file (ACs)  and where there *
#       is a match for the AC in the U5090 file,  it takes the corresponding EC for that AC *
#       and generates a record: EC\t\U50,U90\tU50,U90\t....etc...\n                         *                       
#********************************************************************************************
#   Written by George Weingart - george.weingart@gmail.com   12/02/2015                     *  
#********************************************************************************************



#*************************************************************************************
#* Parse Input parms                                                                 *
#*************************************************************************************
def read_params(x):
	CommonArea = dict()	
	parser = argparse.ArgumentParser(description='Build relation: Reaction --> EC --> UniprotAC --> Uniref5090')
	parser.add_argument('--uniref50gz', action="store", dest='Uniref50gz',nargs='?')
	parser.add_argument('--uniref90gz', action="store", dest='Uniref90gz',nargs='?')
	parser.add_argument('--o_map1', action="store", dest='o_map1',nargs='?', help='Map EC_ACs_U50U90s' )
        parser.add_argument('--o_map2', action="store", dest='o_map2',nargs='?', help='Map EC_U50s' )
	parser.add_argument('--EnzymesFile', action="store", dest='EnzymesFile',nargs='?')
	CommonArea['parser'] = parser
	return  CommonArea

	
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
				print "Total of ", iTotalACsProcessed, " ACs Processed against the Uniref5090 file"	
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
				print "Total of ", iTotalUniref5090RecsRead, " Uniref5090 records read"
			MKey = MasterLine.split()[0]
			continue
		elif  TxnKey < MKey:
			iTxnIndex = iTxnIndex + 1
			if iTxnIndex >=  len(CommonArea['lACsSorted']) -1:
				FlagEnd = True
			TxnKey = CommonArea['lACsSorted'][iTxnIndex]
			iTotalACsProcessed+=1							# Count the reads
			if  iTotalACsProcessed %  iPrintAfterACReads == 0:	# If we need to print status
				print "Total of ", iTotalACsProcessed, " ACs Processed against the Uniref5090 file"	
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
	print "Unzipping uniref50 file"
	cmd_gunzip = "gunzip -c " + strUniref50gz + ">" + strTempDir + "/" + strUniref50gzFileName[:-3] # Build the gunzip command
	os.system(cmd_gunzip)									# Invoke os
	print "Unzipping uniref90 file"
 	cmd_gunzip = "gunzip -c " + strUniref90gz + ">" + strTempDir + "/" + strUniref90gzFileName[:-3] # Build the gunzip command
	os.system(cmd_gunzip)									# Invoke os
	print "Pasting Uniref50 to Uniref90"

	cmd_paste =  "paste " +  strTempDir + "/" + strUniref50gzFileName[:-3] + " " +\
						strTempDir + "/" + strUniref90gzFileName[:-3] + ">" +\
						strTempDir + "/" + strUniref50gzFileName[:-7] +  "90"    # Paste the two files together
	os.system(cmd_paste )									# Invoke os
	dInputFiles["File5090"] = strTempDir + "/" + strUniref50gzFileName[:-7] +  "90"  #Post the file created into the Common Area
	return dInputFiles










	
#********************************************************************************************
#*   Read the Dictionary of EC--->[ACs]  and aggregate the ACs                              *
#********************************************************************************************	
def AggregateACs(CommonArea):

    lAggregatedACs = list()
    for EC, lACs  in CommonArea['dEC_to_AC'].iteritems():
        for AC in lACs:
            lAggregatedACs.append(AC)
    CommonArea['lACsSorted'] = sorted(list(set(lAggregatedACs)))
    return CommonArea 

#********************************************************************************************
#*   Read the enzymes.dat file and create a dictionary that for each EC                     *
#********************************************************************************************
def 	ReadEnzymesFile(CommonArea):	
        CommonArea['dEC_to_AC'] = dict()
	
	InputFile = open(CommonArea['EnzymesFile'])
 
	for iLine in InputFile: 
		iLine = iLine.rstrip('\n')
		if iLine.startswith("ID   "):
			EC = iLine.split()[1].rstrip()
			if EC not in CommonArea['dEC_to_AC']:
			    CommonArea['dEC_to_AC'][EC] = list()
 
		if iLine.startswith("DR   "):
			lAC_Entries =  iLine[2:].replace(" ","").split(";")
 
                        for AC_Entry in lAC_Entries:
                            if len(AC_Entry) > 0:
                                AC = AC_Entry.split(",")[0]
                                if AC not in CommonArea['dEC_to_AC'][EC]:
                                    CommonArea['dEC_to_AC'][EC].append(AC)			

	return CommonArea
		
	
#********************************************************************************************
#* Main                                                                                     *
#********************************************************************************************
print "Program started"

CommonArea = read_params( sys.argv )  # Parse command  
parser = CommonArea['parser'] 
results = parser.parse_args()
CommonArea['EnzymesFile'] = results.EnzymesFile	
CommonArea['OutputFile'] = results.o_map1
CommonArea['OutputFile_Map2'] = results.o_map2	
CommonArea = ReadEnzymesFile(CommonArea)   #  Read the Enzymes file and get EC and all the ACs for each EC
CommonArea = AggregateACs(CommonArea)   # Aggregate the ACs    


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

 
print "Program ended Successfully"
exit(0)