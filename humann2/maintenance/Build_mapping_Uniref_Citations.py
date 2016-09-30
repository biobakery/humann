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
import csv
import re
import uuid
import subprocess
from subprocess import call


#********************************************************************************************
#    Build mapping of Uniref90 / 50 to the DR citations in Swissprot                        *
#                                                                                           *
#                                                                                           *
#  -----------------------------------------------------------------------------------------*
#  Invoking the program:                                                                    *
#  ---------------------                                                                    *
#  python Build_mapping_Uniref_Citations.py  
# --uniref50gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz\
# --uniref90gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef90.dat.gz\
# --i_sprot /n/huttenhower_lab/data/uniprot/2015-06/uniprot_sprot.dat\
# --o_extract   OutputExtractMapping
#                                                       
#  Please note there is an option to keep the temporary directory,  but in general,         *
#    it is not needed, unless there is a problem                                            *
#                                                                                           *
#   The program uses as input the Uniref50 and Uniref90 files,  glues them together         *
#   generating approximately 80 million records and then selects only those records         *
#   where U90 = AC  -  that is approximately 20 million records                             *
#   and loads a dictionary with uniref90-->uniref50 for those records                       *
#                                                                                           *
#   It then proceeds to read Swissprot (The  text file) and for each protein where the      *
#     AC number matches one of the selected U90s in the previous step,  it generates a      *
#     citations extract from the DE and DR records (EC and KO, GO and Pfam                  *
#                                                                                           *
#                                                                                           *
#   Written by George Weingart  May 26, 2016    george.weingart@gmail.com                   *
#  -----------------------------------------------------------------------------------------*
#   Modification Log                                                                        *
#   George Weingart 8/21/16 - Added extraction of DR eggnogg entries  such as :             *
#   DR   eggNOG; NOG87981; -.                                                               *
#   DR   eggNOG; NOG87986; -.                                                               *
#********************************************************************************************



#*************************************************************************************
#* Parse Input parms                                                                 *
#*************************************************************************************
def read_params(x):
	CommonArea = dict()
	parser = argparse.ArgumentParser(description='Build mapping of Uniref50 90 to DR citations in Swissprot')
	parser.add_argument('--i_sprot',
            action="store",
            dest='i_sprot',
            help='Name of the input Swissprot or Uniprot datasetin text format',
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
	parser.add_argument('--o_extract',
	   action="store",
	   dest='o_extract',
	   help='Name of the output dataset that contais the extract: Uniref90 \ Uniref50 \ All the GO,   Pfam, EC and KO mappings',
	   nargs='?',
	   default="o_extract")
	parser.add_argument('--KeepTemporaryDatasets',
            action="store_true",  
            dest='KeepTempDatasets',
            default=False,
            help='If set,  program will not clean generated temp directory  Example:  --KeepTemporaryDatasets')
 
	CommonArea['parser'] = parser
	return  CommonArea
	
	
	
#********************************************************************************************
#*   Read the uniref glued file and extract only U90s that the U90 ID = AC                  *
#*   Layout example: A0A008APQ8      UniRef50_O32037 A0A008APQ8      UniRef90_A5ITF3        *
#********************************************************************************************
def Process_Unirefs_File(CommonArea):
    iU5090RecCntr = 0
    iU90SelectedCntr = 0
    iPrintAfterReads = 1000000 								# Print status after read of these number of records
    with open(CommonArea['strInput5090'] , 'rb') as U5090InputFile:
        U5090Recs = csv.reader(U5090InputFile , delimiter='\t')
        for U5090Row in U5090Recs:  # Read the file
            iU5090RecCntr+=1  # Count the recs
            if  iU5090RecCntr %  iPrintAfterReads == 0:	# If we need to print status
                print("Total of " + str(iU5090RecCntr) + " Uniref5090 records read")
            U50ID = U5090Row[1].split("_")[1] 
            U90ID = U5090Row[3].split("_")[1] 
            AC = U5090Row[2]
            
            if  AC == U50ID:    #If this U50 is equal to an AC, then we are going to use it 
                CommonArea['U50_Centroids'].append(U50ID)   
                
            if U90ID != AC:  # We will select only those ones where the AC = U90ID
                continue
            iU90SelectedCntr+=1
            CommonArea['U90s_Centroids'][U90ID] = dict()  # Note that U90ID = AC
            CommonArea['U90s_Centroids'][U90ID] = U50ID 
            
    CommonArea['U50_Centroids']  = set(CommonArea['U50_Centroids'])  # For future use  
    print("".join("*" * 100))
    print("* Read:                             " + "{:,.0f}".format(iU5090RecCntr) + " U5090 records")
    print("* A total of:                       " + "{:,.0f}".format(iU90SelectedCntr) + " records had U90 = AC")
    print("* A total of:                       " + "{:,.0f}".format(len(CommonArea['U50_Centroids'])) + " records had U50 = AC")

    return CommonArea
 
#********************************************************************************************
#*   Process Swissprot                                                                      *
#******************************************************************************************** 
def  ProcessSwissprotRecords(CommonArea):
    lProteinInfo = list()   # We are going to store the records for a single AC in this list and process when find a break ("//")
    InputFile = open(CommonArea['i_sprot'])
    for iLine in InputFile: 
        iLine = iLine.rstrip('\n')
        if iLine.startswith("//"):
            CommonArea = ProcessProtein(lProteinInfo, CommonArea)  ## Decode the rows of Swissprot Entries we have accumulated
            lProteinInfo = list()  # Reinitialize the table for the next AC
        else:
            if iLine[:2] in CommonArea['sValidSissprotRecordTypes']:  # We will look only at these (AC and DR)
                lProteinInfo.append(iLine)
    return CommonArea 
 
 
 
#********************************************************************************************
#*  Post the citation type to generate the cross reference of CitationType U501,U502,...... *
#********************************************************************************************
def PostUnirefCrossReference(CitationType,lCitationEntries,CommonArea,UnirefType, UnirefEntry):
    for CitationSingleEntry in lCitationEntries :
        if CitationSingleEntry not in CommonArea[CitationType  + UnirefType]["CrossReference"]:
            CommonArea[CitationType  + UnirefType]["CrossReference"][CitationSingleEntry] = list() # Initialize the cross reference
        CommonArea[CitationType  + UnirefType]["CrossReference"][CitationSingleEntry].append("UniRef" + UnirefType + "_" + UnirefEntry)  # Add the entry to the cross reference
    return CommonArea 
 
#********************************************************************************************
#*  Decode the Swissprot Entries                                                            *
#********************************************************************************************
def  DecodeSwissprotLines(lProteinInfo):
    dProteinInfo = dict()
    lSelectedACsMatchingU90s = list()  # List of ACs that have a matching U90
    lSelectedACsMatchingU50s = list()  # List of ACs that have a matching U50
    lECs = list()  # Initialize the ECs list
    lGOs = list()  # Initialize the GOs list
    lKOs = list()  # Initialize KO list
    lPfams = list() # Initialize Pfams list
    leggNOGS = list()  #  Initialize this list GW 20160821
    for strSwissprotLine in lProteinInfo:
        if strSwissprotLine.startswith("AC"):  
            lACs = strSwissprotLine[5:].split(";")  #Check the ACs - we will only process ACs that AC=U90
            lACs = list(filter(None, lACs))  # Remove empty string
            for AC in lACs:
                CommonArea['ACsProcessed']+=1  # This is the number of ACs that were in the Swissprot file
                if  AC in CommonArea['U90s_Centroids']:
                    lSelectedACsMatchingU90s.append(AC)   # Post to the table that contains AC=U90
                if  AC in CommonArea['U50_Centroids']:
                    lSelectedACsMatchingU50s.append(AC)   # Post to the table that contains AC=U90
                    
                    
                    
        if strSwissprotLine.startswith("DE            EC="):
            #*************************************************
            #*  DE with an EC                                *
            #************************************************* 
            strEC1 = strSwissprotLine.split("DE            EC=")[1] 
            EC  = re.split('[ ;]', strEC1)[0] 
            ECLevel  =  EC.replace(".-","").count(".")+1   # Check the EC level of this EC
            if ECLevel == 4:                            # We are accumulating only ECs at level 4
                lECs.append(EC)   # And store it in the ECs table for this AC
            
        if strSwissprotLine.startswith("DR"):  
            if   strSwissprotLine.startswith("DR   GO; "):   # GO Entries
                GOEntry = strSwissprotLine.split("DR   GO; ")[1].split(";")[0]  # This is the GO entry
                lGOs.append(GOEntry)  # Add it to the list 

            elif   strSwissprotLine.startswith("DR   KO; "):   # KO Entries
                KOEntry = strSwissprotLine.split("DR   KO; ")[1].split(";")[0]  # This is the KO entry
                lKOs.append(KOEntry)  # Add it to the list 
                               
            elif   strSwissprotLine.startswith("DR   Pfam; "):   # Pfam Entries
                PfamEntry = strSwissprotLine.split("DR   Pfam; ")[1].split(";")[0]  # This is the Pfam entry
                lPfams.append(PfamEntry)  # Add it to the list  
                
            elif   strSwissprotLine.startswith("DR   eggNOG;"):   # eggNog Entries
                eggNOGEntry = strSwissprotLine.split("DR   eggNOG; ")[1].split(";")[0]  # This is the eggNOG entry GW 20160821
                leggNOGS.append(eggNOGEntry)  # Add it to the list      
                
                  
    if lSelectedACsMatchingU90s  == 0  and lSelectedACsMatchingU50s == 0:  # If this protein is no U90 or U50 centroid - we will not process it 
        return     dProteinInfo    
    
    #********************************************************
    #* Build the decode information                         *
    #********************************************************
    dProteinInfo['lSelectedACsMatchingU90s'] = lSelectedACsMatchingU90s   
    dProteinInfo['lSelectedACsMatchingU50s'] = lSelectedACsMatchingU50s  
    dProteinInfo['lECs'] = lECs
    dProteinInfo['lKOs'] = lKOs   
    dProteinInfo['lGOs'] = lGOs
    dProteinInfo['lPfams'] = lPfams 
    dProteinInfo['leggNOGS'] = leggNOGS 
    if len(lECs) + len(lKOs) + len(lKOs) + len(lGOs)  + len(leggNOGS)== 0:  #If there were no citations   #GW 20160821  Addded eggNOGS
        dProteinInfo = dict()
     
    return     dProteinInfo   

 
#********************************************************************************************
#*  Process the Swissprot entries for an AC                                                 *
#     Glue the two uniref files                                                             *
#********************************************************************************************
def ProcessProtein(lProteinInfo, CommonArea):
    strTab = "\t"
    strEndOfLine = "\n"
    dProteinInfo = DecodeSwissprotLines(lProteinInfo)  # It will contain the decoded Swissprot Entry for that protein
    if len(dProteinInfo)  == 0:   #If there was nothin interesting in this protein (No AC=U90 or AC=U50 or there were no citations)
        return CommonArea

    for U90 in dProteinInfo['lSelectedACsMatchingU90s']:   # print record in main extract and store in the U90 cross references
        strBuiltRecord  = "UniRef90_" + U90 + strTab + "UniRef50_" +  CommonArea['U90s_Centroids'][U90]  # Added the literals "UniRef90 and 50" for downstream compatibility
        if len(dProteinInfo['lECs'])  > 0:
            CommonArea['iCntrRecsWithCitationsEC']+=1   # Count the records that had EC citations
            strECs = strTab.join(dProteinInfo['lECs'])
            strBuiltRecord  = strBuiltRecord  + strTab + "EC=" + strECs
            CommonArea = PostUnirefCrossReference("EC",dProteinInfo['lECs'],CommonArea,"90", U90)  # We will post it for cross refeencing the entries          
            
        if len(dProteinInfo['lGOs'])  > 0:
            CommonArea['iCntrRecsWithCitationsGO']+=1   # Count the records that had GO citations
            strGOs = strTab.join(dProteinInfo['lGOs'])
            strBuiltRecord  = strBuiltRecord  + strTab + "GO=" + strGOs
            CommonArea = PostUnirefCrossReference("GO",dProteinInfo['lGOs'],CommonArea,"90", U90)  # We will post it for cross refeencing the entries
  
        if len(dProteinInfo['lKOs'])  > 0:
            strKOs = strTab.join(dProteinInfo['lKOs'])
            CommonArea['iCntrRecsWithCitationsKO']+=1 # Count those
            strBuiltRecord  = strBuiltRecord  + strTab + "KO=" + strKOs
            CommonArea = PostUnirefCrossReference("KO",dProteinInfo['lKOs'],CommonArea,"90", U90)  # We will post it for cross refeencing the entries 
          
            
        if len(dProteinInfo['lPfams'])  > 0:
            CommonArea['iCntrRecsWithCitationsPfam']+=1 # COunt these
            strPfams = strTab.join(dProteinInfo['lPfams'])
            strBuiltRecord  = strBuiltRecord  + strTab + "Pfam=" + strPfams
            CommonArea = PostUnirefCrossReference("Pfam",dProteinInfo['lPfams'],CommonArea,"90", U90)  # We will post it for cross refeencing the entries 
   
        if len(dProteinInfo['leggNOGS'])  > 0:  #  GW20160821  added eggNOGS
            CommonArea['iCntrRecsWithCitationseggNOGS']+=1 # Count these #  GW20160821  added eggNOGS
            streggNOGS = strTab.join(dProteinInfo['leggNOGS']) #  GW20160821  added eggNOGS
            strBuiltRecord  = strBuiltRecord  + strTab + "eggNOG=" + streggNOGS #  GW20160821  added eggNOGS
            CommonArea = PostUnirefCrossReference("eggNOG",dProteinInfo['leggNOGS'],CommonArea,"90", U90)  # We will post it for cross refeencing the entries #  GW20160821  added eggNOGS
                                            
                                                                                  
                                                                                                                                                              
        strBuiltRecord  = strBuiltRecord  + strEndOfLine 
        CommonArea['TemporaryOutputFile'].write(strBuiltRecord )

    #***************************************************************
    #   Post cross reference entries for U50s                      *
    #***************************************************************  
    for U50 in dProteinInfo['lSelectedACsMatchingU50s']: 
        if len(dProteinInfo['lECs'])  > 0:
            CommonArea = PostUnirefCrossReference("EC",dProteinInfo['lECs'],CommonArea,"50", U50)  # We will post it for cross refeencing the entries          
            
        if len(dProteinInfo['lGOs'])  > 0:
            CommonArea = PostUnirefCrossReference("GO",dProteinInfo['lGOs'],CommonArea,"50", U50)  # We will post it for cross refeencing the entries
  
        if len(dProteinInfo['lKOs'])  > 0:
            CommonArea = PostUnirefCrossReference("KO",dProteinInfo['lKOs'],CommonArea,"50", U50)  # We will post it for cross refeencing the entries 
          
            
        if len(dProteinInfo['lPfams'])  > 0:
            CommonArea = PostUnirefCrossReference("Pfam",dProteinInfo['lPfams'],CommonArea,"50", U50)  # We will post it for cross refeencing the entries 
 
        if len(dProteinInfo['leggNOGS'])  > 0:
            CommonArea = PostUnirefCrossReference("eggNOG",dProteinInfo['leggNOGS'],CommonArea,"50", U50)    #  GW20160821  added eggNOGS
                                                                   
                                                                                                                                                                                                      
    return  CommonArea  
 


 
#********************************************************************************************
#*   Initialize the process                                                                 *
#     Glue the two uniref files                                                             *
#     and open a temporary file in the temp directory to store the output before it is      *
#     sorted                                                                                *
#********************************************************************************************
def Glue_Uniref50_Uniref90_Files(strUniref50gz,  strUniref90gz):
	dInputFiles = dict()									# Initialize the dictionary
	dInputFiles["Uniref50gz"] = strUniref50gz				# Store 1st file name in dictionary
	dInputFiles["Uniref90gz"] = strUniref90gz				# Store 2nd file name in dictionary

	strTempDir = tempfile.mkdtemp()							# Make temporary folder to work in
	CommonArea['strTempDir'] =  strTempDir                         # Post the name of the temporary directory in the CommonArea
	
	#*************************************************************************
	#  We will open a temporary file where we post the extract               *
	#     because we are sorting it before delivering it into the real       *
	#     output file                                                        *
	#*************************************************************************

	TempFileNameRandom = uuid.uuid4()
	CommonArea['TemporaryOutputFileName'] = os.path.join(CommonArea['strTempDir'], str(TempFileNameRandom) )  # We are going to post the output here and then sort it by U90 
	CommonArea['TemporaryOutputFile'] = open(CommonArea['TemporaryOutputFileName'],'w')
        
        
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


	cmd_paste =  "paste " +  strTempDir + "/" + strUniref50gzFileName[:-3] + " " +\
						strTempDir + "/" + strUniref90gzFileName[:-3] + ">" +\
						strTempDir + "/" + strUniref50gzFileName[:-7] +  "90"    # Paste the two files together
	os.system(cmd_paste )									# Invoke os
	dInputFiles["File5090"] = strTempDir + "/" + strUniref50gzFileName[:-7] +  "90"  #Post the file created into the Common Area
	return dInputFiles


#*************************************************************************************
#*  Create the dictionaries and filenames for the mappings cross references          *
#*************************************************************************************
def GenerateTablesAndCrossReferenceFilenames(CommonArea):
    for CitationType in CommonArea["Citations"]:
        for UnirefType in CommonArea["Unirefs"]:
            CommonArea[CitationType+UnirefType] = dict()
            CommonArea[CitationType+UnirefType]["OutputFileName"] = "map_" + CitationType + "_uniref" + UnirefType   + ".txt"  # This is the output file name 
            CommonArea[CitationType+UnirefType]["OutputFileName"] = CommonArea[CitationType+UnirefType]["OutputFileName"].lower()  # Convert it to lower case for output
            CommonArea[CitationType+UnirefType]["CrossReference"] = dict()  # This will contain the cross references that will be printed later 
 
    return CommonArea
    
    
#*************************************************************************************
#*  Generate the cross reference files                                               *
#*************************************************************************************
def GenerateCrossReferenceFiles(CommonArea):
    strTab = "\t"
    strEndOfLine = "\n"
    for CitationType in CommonArea["Citations"]:
        for UnirefType in CommonArea["Unirefs"]:
            OutputFile = open(CommonArea[CitationType + UnirefType ]['OutputFileName'],'w')
            for Citation in sorted(CommonArea[CitationType  + UnirefType]["CrossReference"]): 
                strOutputRec = Citation + strTab + strTab.join(sorted(list(set(CommonArea[CitationType  + UnirefType]["CrossReference"][Citation])))) + strEndOfLine  #Eliminate dups by list(set())
                OutputFile.write(strOutputRec )
            OutputFile.close()
            #************************************************
            #*  Gzip the file                               *
            #************************************************
            CurrentStepOutputFileName =  os.path.join(CommonArea['strTempDir'], "gzip" + CitationType + UnirefType + "results.txt") 
            GzipCommand ="gzip " + CommonArea[CitationType + UnirefType ]['OutputFileName']    
            lCommands = GzipCommand.split(" ")
            with open(CurrentStepOutputFileName, 'w') as f:
                        subprocess.call(lCommands, stdout=f)
            #************************************************
            #*                                              *
            #************************************************
    return CommonArea

#*************************************************************************************
#*  Main Program                                                                     *
#*************************************************************************************
print("Program started")

CommonArea = read_params( sys.argv )  # Parse command  
parser = CommonArea['parser'] 
results = parser.parse_args()
CommonArea['i_sprot'] = results.i_sprot
CommonArea['Uniref50gz'] = results.Uniref50gz
CommonArea['Uniref90gz'] = results.Uniref90gz
CommonArea['o_extract'] = results.o_extract
 

#***********************************************************************
#*  Used to generetate the systems Corss reference files               *
#***********************************************************************
CommonArea["Citations"] = ["KO","Pfam","GO","EC","eggNOG"]
CommonArea["Unirefs"] = ["50","90"] 
CommonArea = GenerateTablesAndCrossReferenceFilenames(CommonArea)  # We will generate here the dixtionaries and file names
#***********************************************************************

CommonArea['U50_Centroids'] =  list()     # List of U50s where U50 = AC in Uniprot


CommonArea['U90s_Centroids'] = dict()   # dict of ACs where AC = U90.  Contains U90=AC and U50
CommonArea['sValidSissprotRecordTypes'] = ("AC","DR", "DE")  #These are the Swissprot record types we are accepting

CommonArea['iCntrRecsWithCitationsPfam'] = 0  # Counter for those Records
CommonArea['iCntrRecsWithCitationsKO'] = 0  # Counter for those Records
CommonArea['iCntrRecsWithCitationsGO'] = 0  # Counter for those Records
CommonArea['iCntrRecsWithCitationsEC'] = 0  # Counter for those Records
CommonArea['iCntrRecsWithCitationseggNOGS'] = 0 # GW20160821  counter for eggNOGS
CommonArea['ACsProcessed'] = 0

dInputFiles =  Glue_Uniref50_Uniref90_Files(CommonArea['Uniref50gz'] ,  CommonArea['Uniref90gz'])  # Invoke initialization
CommonArea['strInput5090'] = dInputFiles["File5090"]		#Name of the Uniref5090 file
print("Uniref50 and 90 glued")


CommonArea = Process_Unirefs_File(CommonArea)   #  Process the unirefs file
CommonArea['OutputFile'] = open(CommonArea['o_extract'],'w')
CommonArea = ProcessSwissprotRecords(CommonArea)  # Process Swissprot - Read the records and select the ones where U90 = AC
CommonArea['TemporaryOutputFile'].close()         # Close the temporary output file




#**************************************************************
#* Sort the temporary output file to the output               *
#**************************************************************
 
CurrentStepOutputFileName  = os.path.join(CommonArea['strTempDir'], "SortOutput.txt" ) 
SortCommand ="sort -t\t -k1 " + CommonArea['TemporaryOutputFileName']   
lCommands = SortCommand.split(" ")
with open(CurrentStepOutputFileName, 'w') as f:
     subprocess.call(lCommands, stdout=f)
 
CopyCommand ="cp " + os.path.join(CommonArea['strTempDir'], "SortOutput.txt" )  + " " + CommonArea['o_extract']  
lCommands = CopyCommand.split(" ")        
CurrentStepCopyOutputFileName = os.path.join(CommonArea['strTempDir'], "CopyOutput.txt" )      
with open(CurrentStepCopyOutputFileName, 'w') as f:
     subprocess.call(lCommands, stdout=f)





print("* Number of proteins (ACs) in Swissprot file read:" + "{:,.0f}".format(CommonArea['ACsProcessed']))
print("* Summary of generated records")
print("* ----------------------------")
print("*")
print("* Proteins with ECs from DE records:" + "{:,.0f}".format(CommonArea['iCntrRecsWithCitationsEC']))
print("* Proteins with KO citations       :" + "{:,.0f}".format(CommonArea['iCntrRecsWithCitationsKO']))
print("* Proteins with Pfam citations     :" + "{:,.0f}".format(CommonArea['iCntrRecsWithCitationsPfam']))
print("* Proteins with GO citations       :" + "{:,.0f}".format(CommonArea['iCntrRecsWithCitationsGO']))
print("* Proteins with eggNOG citations   :" + "{:,.0f}".format(CommonArea['iCntrRecsWithCitationseggNOGS']))
print("".join("*" * 100))


CommonArea = GenerateCrossReferenceFiles(CommonArea)  # Generate the citations files


#**************************************************************
#* Remove the temporary directory                             *
#**************************************************************
if results.KeepTempDatasets == False:
    DeleteCommand ="rm -r " + CommonArea['strTempDir'] 
    CurrentStepDelOutputFileName = os.path.join(CommonArea['strTempDir'], "DelOutput.txt" )  
    lCommands = DeleteCommand.split(" ")  
    with open(CurrentStepDelOutputFileName, 'w') as f:
        subprocess.call(lCommands, stdout=f)
else:
    print("Temporary directory will be kept")
    print("Temporary directory " + CommonArea['strTempDir'] + " ****  Will be kept  ****")



print("Program ended Successfully")
exit(0)
