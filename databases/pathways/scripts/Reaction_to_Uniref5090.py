#!/usr/bin/env python
from cStringIO import StringIO
import sys,string
import sys, os
import argparse
import tempfile 



#********************************************************************************************
#    Map Reactions to Uniref5090                                                            *
#                                                                                           *
#    The objective of this program is to map reactions from the metacyc file to             *
#       uniref50/90                                                                         *
#                                                                                           *
#    Logic:                                                                                 *
#    1. Read the reactions file (See location below) and build relation: REACTION--> EC     *
#    2. Read Swissprot file (See location below) and build relations EC--> Swissprot AC     *
#    3. Build the relations REACTIONs --> UniprotKb ACs                                     *
#    4. Build and print the relations REACTIONS --> UniRef50, 90                            *
#                                                                                           *
#  -----------------------------------------------------------------------------------------*
#  Invoking the program:                                                                    *
#  ---------------------                                                                    *
#  python Reaction_to_Uniref5090.py --i_reactions /n/huttenhower_lab_nobackup/downloads/metacyc/18.1/data/reactions.dat  --i_sprot /n/huttenhower_lab/data/uniprot/2014-09/uniprot_sprot.dat  --uniref50gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz --uniref90gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef90.dat.gz  --o mapping_reactions_to_uniref5090
#                                                                                           *
#   Where:                                                                                  *
#    --i_reactions, is the reactions file, which is currently located at                    *
#    /n/huttenhower_lab_nobackup/downloads/metacyc/18.1/reactions.dat                       *
#                                                                                           *           
#   --i_sprot input_file is the UniprotKB Swissprot text file, which can be downloaded from *
#    ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz *
#   The current downloaded i_sprot file, which serves as input,  resides on hutlab3 in      *
#    /n/huttenhower_lab/data/uniprot/2014-09/uniprot_sprot.dat                              *
#                                                                                           *
#    uniref50gz and uniref90gz are the uniref mappings (Uniref50 --> Uniprot AC)            * 
#     currently located at                                                                  *
#     /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz                         *

#   Written by George Weingart - george.weingart@gmail.com   10/08/2014                     *  
#********************************************************************************************








#*************************************************************************************
#* Read Reactions file and build relation: REACTION --> ECs                          *
#*************************************************************************************
def ReadReactions(iFileReactions):
	bFlagEC = False
	lECs = list()
	InputFile = open(iFileReactions)
	LineReactionsCntr = 0
	ReactionsToECCntr = 0
	dReactionsToECs = dict()
	for iLine in InputFile: 
			LineReactionsCntr = LineReactionsCntr +1
			if iLine.startswith("UNIQUE-ID"):
				ReactionName = iLine.split(" ")[2].replace("\n","")
  
			if iLine.startswith("EC-NUMBER"):
				EC = iLine.split("-")[3].replace("\n","")
				bFlagEC = True
				lECs.append(EC)
				
 
			if  bFlagEC == True and iLine.startswith("//"):
				if ReactionName not in dReactionsToECs:
						ReactionsToECCntr = ReactionsToECCntr + 1
						dReactionsToECs[ReactionName] = list()
				for EC in  lECs:
					dReactionsToECs[ReactionName].append(EC)
					bFlagEC = False
					lECs = list()
	print "Read " + str(LineReactionsCntr) + " Input Lines from the reactions file"
	print "The table Reactions --> ECs contains " + str(ReactionsToECCntr) + " Reactions"
	InputFile.close()
	CommonArea["dReactionsToECs"] =  dReactionsToECs
	return CommonArea







#*************************************************************************************
#* Read Swissprot  and build relation: EC --> UniprotACs                             *
#*************************************************************************************
def ReadSwissprot(i):
	iFile = i
	bFlagEC = False
	lECs = list()
	InputFile = open(iFile)
	LineCntr = 0
	ECsToACsCntr = 0
	dECsToUniprotACs = dict()
	for iLine in InputFile: 
			LineCntr = LineCntr +1
			if iLine.startswith("AC   "):
				lTemp = iLine.rstrip().split().pop().split(";")
				lACs = [var for var in lTemp if var]
  
			if iLine.startswith("DE   ") and "EC=" in iLine:
				bFlagEC = True
				iLine1 = iLine.replace(";"," ")
				iECLoc = iLine1.find("EC=") + 3
				EC = iLine1[iECLoc:].split(" ")[0]
				EC  =  EC.replace("EC=","")
				EC  =  EC.replace(".-","")
				lECs.append(EC)
 
			if  bFlagEC == True and iLine.startswith("//"):
				   for ProtAC in lACs:
						OutputLine = ProtAC 
						for EC in  lECs:
							if EC not in dECsToUniprotACs:
								ECsToACsCntr = ECsToACsCntr + 1
								dECsToUniprotACs[EC] = list()
							dECsToUniprotACs[EC].append(ProtAC)
							
				   bFlagEC = False
				   lECs = list()
				   lACs = list() 
	print "Read " + str(LineCntr) + " Input Lines from uniprot_sprot.dat"
	print "The table ECs --> UniprotKb ACs contains "  + str(ECsToACsCntr) + " ECs"
	InputFile.close()
	CommonArea["dECsToUniprotACs"] = dECsToUniprotACs
	
	return CommonArea








#*************************************************************************************
#* Parse Input parms                                                                 *
#*************************************************************************************
def read_params(x):
	CommonArea = dict()	
	parser = argparse.ArgumentParser(description='Build relation: Reaction --> EC --> UniprotAC --> Uniref5090')
	parser.add_argument('--i_reactions', action="store", dest='i_reactions',nargs='?')
	parser.add_argument('--i_sprot', action="store", dest='i_sprot',nargs='?')
	parser.add_argument('--uniref50gz', action="store", dest='Uniref50gz',nargs='?')
	parser.add_argument('--uniref90gz', action="store", dest='Uniref90gz',nargs='?')
	parser.add_argument('--o', action="store", dest='o',nargs='?' )
	CommonArea['parser'] = parser
	return  CommonArea


#*************************************************************************************
#* Build relation: REACTION --> ACs                                                  *
#*************************************************************************************
def ResolveReactionsToACs(CommonArea):
    ECsToUniprotACsCntr = 0 
    dReactionsToACs = dict()
    for Reaction, lECs  in CommonArea["dReactionsToECs"].iteritems():
		for EC in lECs:
			 for EC in lECs:
				if EC in   CommonArea["dECsToUniprotACs"]:
					try:
						if Reaction not in dReactionsToACs:
							ECsToUniprotACsCntr = ECsToUniprotACsCntr + 1
							dReactionsToACs[Reaction] = list()
						dReactionsToACs[Reaction] = CommonArea["dECsToUniprotACs"][EC]
					except:
						pass
    print "The table Reactions --> UniprotKb ACs contains "  + str(ECsToUniprotACsCntr) + " Reactions"
    CommonArea['dReactionsToACs'] = dReactionsToACs
    del CommonArea["dECsToUniprotACs"]
	#*****************************************************************
	#*   Select only the ACs that are needed for next step           *
	#*****************************************************************
    lAccumulatedACs = list()
    for Reaction, lACs  in CommonArea["dReactionsToACs"].iteritems():
		for AC in lACs:
			lAccumulatedACs.append(AC)
    CommonArea["sAccumulatedACs"] = set(lAccumulatedACs)		#Get rid of dups
    print "The total number of AC entries is: ",  str(len(CommonArea["sAccumulatedACs"]))
    return CommonArea

	

#********************************************************************************************
#*   Read Uniref 5090 file  and build translation table                                     *
#********************************************************************************************
def ReadUniref5090File(strInput5090,sUniprotIds):
	
	dUniprotUniref = dict()
	iTotalUniref5090RecsLoaded = 0							# Counter of Uniref5090 recs loaded
	iTotalUniref5090RecsRead = 0							# Counter of Uniref5090 recs Read
	iPrintAfter = 1000000									# Print status after multiple of these number of records
	iPrintAfterReads = 100000  								# Print status after read of these number of records
	File5090 = open(strInput5090)							# Open the file
	for strInputLine in File5090:   						# Read Input file
		iTotalUniref5090RecsRead+=1							# Count the reads
		if  iTotalUniref5090RecsRead %  iPrintAfterReads == 0:	# If we need to print status
			print "Total of ", iTotalUniref5090RecsRead, " Uniref5090 records read"
		lInputLineSplit = strInputLine.split() 				# Split the line using space as delimiter

		if lInputLineSplit[0]  not in sUniprotIds:			# If we don't need this ID 
			continue										# skip it 
		lEnt5090 = list()									# Initialize list
		lEnt5090.append(lInputLineSplit[1].split("_")[1])	# Entry is of the form UniRef50_Q2FWP1 - We need only the Q2FWP1
		lEnt5090.append(lInputLineSplit[3].split("_")[1])	# Entry is of the form UniRef50_Q2FWP1 - We need only the Q2FWP1
		dUniprotUniref[lInputLineSplit[0]] = lEnt5090		# Post it in the dictionary
		iTotalUniref5090RecsLoaded+=1						# Increase counter
		if  iTotalUniref5090RecsLoaded %  iPrintAfter == 0:	# If we need to print status
			print "Total of ", iTotalUniref5090RecsLoaded, " Uniref5090 records loaded into the table"
	print "Load Complete - Total of ", iTotalUniref5090RecsLoaded, " Uniref5090 records loaded into the table"
	File5090.close()										# CLose the file
	return dUniprotUniref




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
#* Generate Extract:  Reaction --> Uniref50, 90                                             *
#********************************************************************************************
def GenerateExtract(CommonArea, OutputFileName):
	strNewLine = "\n"
	strTab = "\t"
	ReactionToUnirefCntr = 0
	sReactionsWritten = set()
	
	OutputFile = open(OutputFileName,'w')		#Open the Output file
	for Reaction, lACs  in CommonArea["dReactionsToACs"].iteritems():
		bFlagUnirefFound = False
		lBuiltRecord = [Reaction]
		lBuiltRecord.append(CommonArea["dReactionsToECs"][Reaction][0])  #Post the first EC
		lU50 = list()
		lU90 = list()	
		for AC in lACs:
			if AC in CommonArea['dUniprotUniref']:
				Uniref50 = "UniRef50_" + CommonArea['dUniprotUniref'][AC][0] 
				lU50.append(Uniref50)
				Uniref90 = "UniRef90_" + CommonArea['dUniprotUniref'][AC][1] 
				lU90.append(Uniref90)
 
		lU50Sorted = sorted(list(set(lU50)))
		lU90Sorted = sorted(list(set(lU90)))
 
		for U50 in lU50Sorted:
			bFlagUnirefFound = True
			lBuiltRecord.append(U50)
		for U90 in lU90Sorted:
			bFlagUnirefFound = True
			lBuiltRecord.append(U90)
 
		strBuiltRecord = "\t".join(lBuiltRecord) + 	strNewLine
		lBuiltRecord = list()
			
		if Reaction not in sReactionsWritten  and bFlagUnirefFound == True:
			sReactionsWritten = sReactionsWritten | {Reaction}
			OutputFile.write(strBuiltRecord )
			ReactionToUnirefCntr = ReactionToUnirefCntr + 1
			
	OutputFile.close()	
	print "Total REACTION -> Uniref(50) Uniref(90) relations generated in the file " + OutputFileName + " = " + str(ReactionToUnirefCntr)
	return CommonArea  


#********************************************************************************************
#* Main                                                                                     *
#********************************************************************************************
print "Program started"
CommonArea = read_params( sys.argv )  # Parse command  
parser = CommonArea['parser'] 
results = parser.parse_args()


iFileReactions = results.i_reactions
iFileSProt = results.i_sprot
strUniref50gz = results.Uniref50gz					# The first file is the zipped version of the Uniref50 Translation file
strUniref90gz = results.Uniref90gz					# The 2nd file is the zipped version of the Uniref90 Translation file
OutputFileName = results.o								# The output file
 
CommonArea  = ReadReactions(iFileReactions)			#Read the reactions file to build the relations REACTION -> EC
CommonArea  = ReadSwissprot(iFileSProt)				#Read the Swissprot file to build the relations ECs --> UniprotKb ACs
CommonArea  = ResolveReactionsToACs(CommonArea)		#Build the relations Reactions --> UniprotKb ACs

 



#***************************************
# Processing Uniref50 90 files         *
#***************************************

dInputFiles =  InitializeProcess(strUniref50gz,  strUniref90gz)  # Invoke initialization
strInput5090 =  dInputFiles["File5090"]		#Name of the Uniref5090 file


print "Starting the load of the 5090 table\n"
CommonArea['dUniprotUniref'] = ReadUniref5090File(strInput5090,CommonArea["sAccumulatedACs"])	#Invoke reading of the file

cmd_remove_tempdir = "rm -r /" + dInputFiles["TempDirName"]		# Remove the temporary directory
os.system(cmd_remove_tempdir)

#******************************************
# Generate extract: Reaction: Uniref50 90 *
#******************************************
CommonArea = GenerateExtract(CommonArea,OutputFileName)			#Generate extract



print "Program ended Successfully"
exit(0)
