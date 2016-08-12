#!/usr/bin/env python
import sys,string
import sys, os
import argparse
import tempfile 
import pdb



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
#  python Reaction_to_Uniref5090.py --i_reactions /n/huttenhower_lab_nobackup/downloads/metacyc/19.1/data/reactions.dat 
##        --i_sprot /n/huttenhower_lab/data/uniprot/2015-06/uniprot_sprot.dat 
##        --uniref50gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz 
##        --uniref90gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef90.dat.gz
##        --EnzymesFile   /n/huttenhower_lab_nobackup/downloads/enzymes/enzyme.dat
##         --o mapping_reactions_to_uniref5090
#                                                                                           *
#   Where:                                                                                  *
#    --i_reactions, is the reactions file, which is currently located at                    *
#    /n/huttenhower_lab_nobackup/downloads/metacyc/18.1/reactions.dat                       *
#                                                                                           *           
#   --i_sprot input_file is the UniprotKB Swissprot text file, which can be downloaded from *
#    ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz *
#   The current downloaded i_sprot file, which serves as input,  resides on hutlab3 in      *
#    /n/huttenhower_lab/data/uniprot/2015-06/uniprot_sprot.dat                              *
#                                                                                           *
#    uniref50gz and uniref90gz are the uniref mappings (Uniref50 --> Uniprot AC)            * 
#     currently located at                                                                  *
#     /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz                         *
#                                                                                           *
#   Written by George Weingart - george.weingart@gmail.com   10/08/2014                     *  
#********************************************************************************************


#********************************************************************************************
#                                                                                           *
#  Modification by George Weingart - george.weingart@gmail.com   29/08/2014                 *  
#                                                                                           *
#  Added an option to expand EC Level 3's to all its subclasses, for example                *
#  if EC=6.1.1 is detected,  we will post 6.1.1.1, 6.1.1.2,.....,6.1.1.n                    *
#  This will based on a parm passed by the User --ECLevel3Expand                            *
#  If this parm is received,  we will read the enzymes.dat file, build the table that       *
#  contains the expansion ECLevel3 -> [ECLevel4(1), ECLevel4(2)....] and post it            *
#  The location of the enzymes.dat is specified via the parm --EnzymesFile filename         *
#  << Please note >> That file can be downloaded from:                                      *
#   ftp://ftp.expasy.org/databases/enzyme/                                                  * 
#  and  more particularly from                                                              *                         
#  ftp://ftp.expasy.org/databases/enzyme/enzyme.dat                                         *                     
#********************************************************************************************




#*************************************************************************************
#* Read Reactions file and build relation: REACTION --> ECs                          *
#*                                                                                   *
#*  Modification LOG                                                                 *
#*  George Weingart 2015/07/10                                                       *
#*  When reading the reactions.dat file,  also capture the direct relation           *
#*  reaction --> AC                                                                  *
#* Load it into a dictionary and use it after the relation Reaction -> EC -> AC      *
#* as an additional set of relations                                                 *
#* These new relations are identified via the DBLINKS - (UNIPROT prefix in the       *
#* reactions.dat file                                                                *
#*************************************************************************************
def ReadReactions(iFileReactions):
	bFlagEC = False
	lECs = list()
	InputFile = open(iFileReactions)
	LineReactionsCntr = 0
	ReactionsToECCntr = 0
	dReactionsToECs = dict()
	dDirectReactionToACs = dict()  #*  Modification: GW 20150710 - Get the direct reaction--> AC relations*
	for iLine in InputFile: 
			LineReactionsCntr = LineReactionsCntr +1
			
				
			if iLine.startswith("UNIQUE-ID"):
				ReactionName = iLine.split(" ")[2].replace("\n","")
				
			#***********************************************************************
			#*  Modification: GW 20150710 - Get the direct reaction--> AC relations*
			#***********************************************************************

			if iLine.startswith("DBLINKS - (UNIPROT"):	
				AC = iLine.split('"')[1]					# Get the AC
				if ReactionName not in dDirectReactionToACs: #  If new entry - build the dictionary entry made of a new list
					dDirectReactionToACs[ReactionName] = list() 
				dDirectReactionToACs[ReactionName].append(AC)		# And add the AC entry
  			#***************************************************************************
			#*  End Modification: GW 20150710 - Get the direct reaction--> AC relations*
			#***************************************************************************
  
			try:                       #Modification GW 20150825 - Allow for secondary syntax in the dat file for the EC
				if iLine.startswith("EC-NUMBER"):
					EC = iLine.split("-")[3].replace("\n","").lstrip()
					bFlagEC = True  	# GW 20150709
					lECs.append(EC)   	# GW 20150709
			except:
				EC = iLine.split("-")[2].replace("\n","").lstrip() #Modification GW 20150825 - Allow for secondary syntax in the dat file for the EC
				bFlagEC = True  	#Modification GW 20150825 - Allow for secondary syntax in the dat file for the EC
				lECs.append(EC)   	#Modification GW 20150825 - Allow for secondary syntax in the dat file for the EC

				
 
			if  bFlagEC == True and iLine.startswith("//"):
				if ReactionName not in dReactionsToECs:
						ReactionsToECCntr = ReactionsToECCntr + 1
						dReactionsToECs[ReactionName] = list()
				for EC in  lECs:
					dReactionsToECs[ReactionName].append(EC)
					bFlagEC = False
					lECs = list()
	print("Read " + str(LineReactionsCntr) + " Input Lines from the reactions file")
	print("The table Reactions --> ECs contains " + str(ReactionsToECCntr) + " Reactions")
	InputFile.close()
	CommonArea["dReactionsToECs"] =  dReactionsToECs
	CommonArea["dDirectReactionToACs"] = dDirectReactionToACs #   Modification: GW 20150710
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
				lTemp = iLine.rstrip().split("AC   ").pop().split(";")
				lACs = [var.replace(" ","") for var in lTemp if var]

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
	print("Read " + str(LineCntr) + " Input Lines from uniprot_sprot.dat")
	print("The table ECs --> UniprotKb ACs contains "  + str(ECsToACsCntr) + " ECs")
	InputFile.close()
	CommonArea["dECsToUniprotACs"] = dECsToUniprotACs
	
	return CommonArea








#*************************************************************************************
#* Parse Input parms                                                                 *
#* Midification log:                                                                 *
#* George Weingart: 2016/06/08                                                       *
#* Added a parm to select only centroids of the U50 and U90 where a centroid is      *
#*   the U50 or U90 where U50 or U90 = AC                                            *
#*   --UseU5090CentroidsOnly                                                         *
#* And added defaults for the parameters                                             *
#*************************************************************************************
def read_params(x):
	CommonArea = dict()	
	parser = argparse.ArgumentParser(description='Build relation: Reaction --> EC --> UniprotAC --> Uniref5090')
	parser.add_argument('--i_reactions',
	   action="store",
	   dest='i_reactions',
	   nargs='?',
	   default="/n/huttenhower_lab_nobackup/downloads/metacyc/19.1/data/reactions.dat")
	   
	parser.add_argument('--i_sprot',
   	    action="store",
   	    dest='i_sprot',
   	    nargs='?',
   	    default="/n/huttenhower_lab/data/uniprot/2015-06/uniprot_sprot.dat")
   	    
	parser.add_argument('--uniref50gz',
	   action="store",
	   dest='Uniref50gz',
	   nargs='?',
	   default="/n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz")
	   
	parser.add_argument('--uniref90gz',
	   action="store",
	   dest='Uniref90gz',
	   nargs='?',
	   default="/n/huttenhower_lab/data/idmapping/map_uniprot_UniRef90.dat.gz")

	parser.add_argument('--ECLevel3Expand',
	   action="store", 
	   dest='ECLevel3Expand',
	   nargs='?')
	   
	parser.add_argument('--EnzymesFile',
	   action="store",
	   dest='EnzymesFile',
	   nargs='?',
	   default="/n/huttenhower_lab_nobackup/downloads/enzymes/enzyme.dat")
	   
	parser.add_argument('--o',
	   action="store",
	   dest='o',
	   nargs='?',
	   default="mapping_reactions_to_uniref5090.tab")

	parser.add_argument('--UseU5090CentroidsOnly',
	   action="store",
	   dest='UseU5090CentroidsOnly',
	   nargs='?',
	   default="Y")	   
	   	   	   
	CommonArea['parser'] = parser
	return  CommonArea

#********************************************************************************************
#*  Modification: GW 20150710 - Get the direct reaction--> AC relations                     *
#*                                                                                           *
#********************************************************************************************
#*   Incorporate the direct relations  Reactions--> AC  to the ones retrieved from          *
#*       Reactions --> EC --> AC                                                            *
#*   We will be adding the entries in dDirectReactionsToACs to dReactionsToACs              *
#*   Note:  Many of the reactions don't have a direct reaction --> AC                       *
#*   but we want to take advantage of these relations                                       *
#********************************************************************************************	
def  IncorporateDirectReactionsToACEntries(CommonArea):
	for Reaction, lACs  in CommonArea["dDirectReactionToACs"].iteritems():
		if Reaction  not in CommonArea['dReactionsToACs']:   #If there is not an entry for this reaction - build it
			CommonArea['dReactionsToACs'][Reaction] = list()  # As an empty list
		CommonArea['dReactionsToACs'][Reaction] = list(set(CommonArea['dReactionsToACs'][Reaction])|set(CommonArea["dDirectReactionToACs"][Reaction])) # And merge the lists
	return CommonArea
	
	
	

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
 	
	CommonArea['dReactionsToACs'] = dReactionsToACs
	CommonArea =  IncorporateDirectReactionsToACEntries(CommonArea)  #   Modification: GW 20150710
	print("The table Reactions --> UniprotKb ACs contains "  + str(ECsToUniprotACsCntr) + " Reactions")

 
	del CommonArea["dECsToUniprotACs"]
	#*****************************************************************
	#*   Select only the ACs that are needed for next step           *
	#*****************************************************************
	lAccumulatedACs = list()
	for Reaction, lACs  in CommonArea["dReactionsToACs"].iteritems():
		for AC in lACs:
			lAccumulatedACs.append(AC)
	CommonArea["sAccumulatedACs"] = set(lAccumulatedACs)		#Get rid of dups
	print("The total number of AC entries is: " + str(len(CommonArea["sAccumulatedACs"])))
	

	return CommonArea

	
#********************************************************************************************
#*   Read Uniref 5090 file  and build translation table                                     *
#********************************************************************************************
def ReadUniref5090File(strInput5090,sUniprotIds):
	
	dUniprotUniref = dict()
	iTotalUniref5090RecsLoaded = 0							# Counter of Uniref5090 recs loaded
	iTotalUniref5090RecsRead = 0							# Counter of Uniref5090 recs Read
	iPrintAfter = 10000000									# Print status after multiple of these number of records
	iPrintAfterReads = 1000000  								# Print status after read of these number of records
	File5090 = open(strInput5090)							# Open the file
	for strInputLine in File5090:   						# Read Input file
		iTotalUniref5090RecsRead+=1							# Count the reads
		if  iTotalUniref5090RecsRead %  iPrintAfterReads == 0:	# If we need to print status
			print("Total of " + str(iTotalUniref5090RecsRead) + " Uniref5090 records read")
		lInputLineSplit = strInputLine.split() 				# Split the line using space as delimiter

		if lInputLineSplit[0]  not in sUniprotIds:			# If we don't need this ID 
			continue										# skip it 
		lEnt5090 = list()									# Initialize list
		U50 = lInputLineSplit[1].split("_")[1]     
		U90 = lInputLineSplit[3].split("_")[1]       

		#***************************************************
		#  Modification GW 20160609 Add to centroids       *
		#***************************************************			
		CommonArea['U50Centroids'].append(U50)             # Add it to the centroids table
		CommonArea['U90Centroids'].append(U90)             # Add it to the centroids table
		
		
		lEnt5090.append(U50)	          # Entry is of the form UniRef50_Q2FWP1 - We need only the Q2FWP1
		lEnt5090.append(U90)	          # Entry is of the form UniRef50_Q2FWP1 - We need only the Q2FWP1
		dUniprotUniref[lInputLineSplit[0]] = lEnt5090		# Post it in the dictionary
		

		
		
		
		iTotalUniref5090RecsLoaded+=1						# Increase counter
		if  iTotalUniref5090RecsLoaded %  iPrintAfter == 0:	# If we need to print status
			print("Total of " + str(iTotalUniref5090RecsLoaded) + " Uniref5090 records loaded into the table")
	print("Load Complete - Total of " + str(iTotalUniref5090RecsLoaded) + " Uniref5090 records loaded into the table")
	File5090.close()										# CLose the file
	
	CommonArea['U50Centroids'] = set(CommonArea['U50Centroids'])             # Elininate duplicates
	CommonArea['U90Centroids'] = set(CommonArea['U90Centroids'])             # Elininate duplicates
	
	print("Total number of U50 Centroids loaded: " +  str(len(CommonArea['U50Centroids'])))
	print("Total number of U90 Centroids loaded: " +  str(len(CommonArea['U90Centroids'])))	
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

#********************************************************************************************
#* Generate Extract:  Reaction --> Uniref50, 90                                             *
#********************************************************************************************
def GenerateExtract(CommonArea, OutputFileName):
	strNewLine = "\n"
	strTab = "\t"
	strComma = ","
	ReactionToUnirefCntr = 0
	CntrReactionsWithNoEc = 0  # Counting how many there are 

	sReactionsWritten = set()
	
	OutputFile = open(OutputFileName,'w')		#Open the Output file
	for Reaction, lACs  in CommonArea["dReactionsToACs"].iteritems():
		bFlagUnirefFound = False
		bFlagECFound = True    #Modification by GW  DT20150721: If no EC - don't select the record
		lBuiltRecord = [Reaction]

		try:
			lListOfECsToCheckECLevel = CommonArea["dReactionsToECs"][Reaction]  # Modification by George Weingart 20150710  We are checking the EC Level of all ECs
			lPostedECList = list()	# Modification by George Weingart 20150710
			lListOfECsToCheckECLevel = list(set(lListOfECsToCheckECLevel)) #Modification by George Weingart 20150825 - Select only unique EC
			for ECToCheck in lListOfECsToCheckECLevel:   	# Modification by George Weingart 20150710 
				ECToCheck = ECToCheck.replace(".-","")  	# Modification by George Weingart 20150710 
				iECLevel = ECToCheck.count(".") + 1			# Modification by George Weingart 20150710 
				if iECLevel == 4 or (CommonArea['ECLevel3Expand']  == False and iECLevel == 3):  # Modification by George Weingart 20150829
					lPostedECList.append(ECToCheck)  # Modification by George Weingart 20150710
				if iECLevel == 3 and  CommonArea['ECLevel3Expand']  == True:  # Modification by George Weingart 20150829
					for ECLevel4Entry in CommonArea['dECLevel3ToLevel4'][ECToCheck]:    # Review the ECLevel4 entry and if needed, post it
						if ECLevel4Entry not in lPostedECList:
							lPostedECList.append(ECLevel4Entry)  # Modification by George Weingart 20150710
			if len(lPostedECList) == 0:		# Modification by George Weingart 20150710 - If no ECs because we eliminated levels 1 and 2 - Set empty list
				lPostedECList.append(" ")
				bFlagECFound = False    #Modification by GW  DT20150721: If no EC - don't select the record
			lPostedECList = sorted(list(set(lPostedECList)))     #Modification by GW  DT20150728 - Remove duplicate entries from the EC list 			
			lBuiltRecord.append(strComma.join(lPostedECList))  #Post the list of ECs
		except:
			bFlagECFound = False    #Modification by GW  DT20150721: If no EC - don't select the record
			lBuiltRecord.append(" ")    # Note that there is no EC for this reaction - so we are posing an empty list
			CntrReactionsWithNoEc =  CntrReactionsWithNoEc  # No EC for this Reaction
			
		lU50 = list()
		lU90 = list()	
		
		#********************************************************************************
		#  Modification Log:                                                            *
		#  Modified the code below so that if, requested only centroids,  use only      *
		#  centroids (Defined as AC=U90 or AC=U50                                       *
		#  Otherwise, select all ACs,  like we were doing till now                      *
		#  George Weingart  2016/06/08                                                  *
		#********************************************************************************
		for AC in lACs:
			if AC in CommonArea['dUniprotUniref']:
				if  CommonArea['UseU5090CentroidsOnly'] == "Y": 
					if AC == CommonArea['dUniprotUniref'][AC][0]:  # Is it a U50 Centroid? 
                                    		Uniref50 = "UniRef50_" + CommonArea['dUniprotUniref'][AC][0]
                                    		lU50.append(Uniref50)
					if AC == CommonArea['dUniprotUniref'][AC][1]:  # Is it a U90 Centroid? 
				    		Uniref90 = "UniRef90_" + CommonArea['dUniprotUniref'][AC][1]
				    		lU90.append(Uniref90)    #GW  20160702  Was U50 by mistake - fixed it 
				else:
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
			
		if Reaction not in sReactionsWritten  and bFlagUnirefFound == True  and bFlagECFound == True:    #Modification by GW  DT20150721: If no EC - don't select the record 
			sReactionsWritten = sReactionsWritten | {Reaction}
			OutputFile.write(strBuiltRecord )
			ReactionToUnirefCntr = ReactionToUnirefCntr + 1
			
	OutputFile.close()	
	print("Total REACTION -> Uniref(50) Uniref(90) relations generated in the file " + OutputFileName + " = " + str(ReactionToUnirefCntr))
	return CommonArea  

#********************************************************************************************
#*   Read the enzymes.dat file and create a dictionary that for each EClevel3 contains      *
#*   all the EClevel4 subclases, for example:                                               *
#*   6.1.1 --> [6.1.1.1, 6.1.1.2 ....... 6.1.1.x]                                           *
#********************************************************************************************
def 	ReadEnzymesFile(strEnzymesFile):	
	CommonArea['EnzymesFile'] = strEnzymesFile
	CommonArea['ECLevel3Expand'] = True
	InputFile = open(CommonArea['EnzymesFile'])
	dECLevel3ToLevel4 = dict()
 
	for iLine in InputFile: 
		iLine = iLine.rstrip('\n')
		if iLine.startswith("ID   "):
			EC = iLine.split()[1]
			ECLevel = EC.count(".") +1
			ECNodes = EC.split('.')
			if ECLevel == 4:
			    ECLevel3 = ECNodes[0] + '.' + ECNodes[1] + '.' + ECNodes[2]
			    if ECLevel3 not in dECLevel3ToLevel4:
			       dECLevel3ToLevel4[ECLevel3] = list()
			    dECLevel3ToLevel4[ECLevel3].append(EC)
				
	CommonArea['dECLevel3ToLevel4']	 =  dECLevel3ToLevel4
	return CommonArea
	

#********************************************************************************************
#* Main                                                                                     *
#********************************************************************************************
print("Program started")
CommonArea = read_params( sys.argv )  # Parse command  
parser = CommonArea['parser'] 
results = parser.parse_args()

################################################################################
# Modification by George Weingart 20150829                                     *
# If User requested Expansion of EC level 3 to 4, load the table               *
################################################################################
CommonArea['ECLevel3Expand']  = False   # Set default flag to False
if results.ECLevel3Expand  == "YES":    # If user requested 
	CommonArea['ECLevel3Expand']  = True   # Set default flag to True
	strEnzymesFile = results.EnzymesFile
	CommonArea = ReadEnzymesFile(strEnzymesFile)
################################################################################
# End Modification by George Weingart 20150829                                 *
# If User requested Expansion of EC level 3 to 4, load the table               *
################################################################################

iFileReactions = results.i_reactions
iFileSProt = results.i_sprot
strUniref50gz = results.Uniref50gz					# The first file is the zipped version of the Uniref50 Translation file
strUniref90gz = results.Uniref90gz					# The 2nd file is the zipped version of the Uniref90 Translation file
OutputFileName = results.o								# The output file
CommonArea['UseU5090CentroidsOnly'] = results.UseU5090CentroidsOnly      #Modification GW 20160608: Default - Use only U50 and u90 Centroids (Where AC=U50 or AC=U90)
 
CommonArea  = ReadReactions(iFileReactions)			#Read the reactions file to build the relations REACTION -> EC
CommonArea  = ReadSwissprot(iFileSProt)				#Read the Swissprot file to build the relations ECs --> UniprotKb ACs
CommonArea  = ResolveReactionsToACs(CommonArea)		#Build the relations Reactions --> UniprotKb ACs

CommonArea["U50Centroids"] = list()           # We are going to capture all U50s where U50=AC
CommonArea["U90Centroids"] = list()           # We are going to capture all U90s where U90=AC

 



#***************************************
# Processing Uniref50 90 files         *
#***************************************

dInputFiles =  InitializeProcess(strUniref50gz,  strUniref90gz)  # Invoke initialization
strInput5090 =  dInputFiles["File5090"]		#Name of the Uniref5090 file


print("Starting the load of the 5090 table\n")
CommonArea['dUniprotUniref'] = ReadUniref5090File(strInput5090,CommonArea["sAccumulatedACs"])	#Invoke reading of the file

cmd_remove_tempdir = "rm -r /" + dInputFiles["TempDirName"]		# Remove the temporary directory
os.system(cmd_remove_tempdir)

#******************************************
# Generate extract: Reaction: Uniref50 90 *
#******************************************
CommonArea = GenerateExtract(CommonArea,OutputFileName)			#Generate extract



print("Program ended Successfully")
exit(0)
