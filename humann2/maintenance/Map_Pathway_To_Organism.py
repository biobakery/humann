#!/usr/bin/env python
import sys,string
import sys, os
import argparse
import tempfile



#********************************************************************************************
#   Map pathway to Organism                                                                 *
#                                                                                           *
#  -----------------------------------------------------------------------------------------*
#  Invoking the program:                                                                    *
#  ---------------------                                                                    *
#   python ReadSwissprot.py  --ic  /n/huttenhower_lab_nobackup/downloads/metacyc/18.1/data/classes.dat \
#    --ip /n/huttenhower_lab_nobackup/downloads/metacyc/18.1/data//pathways.dat \ 
#   --o mapping_patwhays_to_organisms                          
#                                                                                           *
#                                                                                           *
#                                                                                           *
#   Written by George Weingart - george.weingart@gmail.com   8/2/2015                       *  
#                                                                                           *
#********************************************************************************************

#*************************************************************************************
#* Read classes.dat                                                                  *
#* We are interested in the TAXonomy entries that look as follows                    *
#*  UNIQUE-ID - TAX-7038                                                             *
#*  and their corresponding COMMON-NAME  that look as follows:                       *
#*  COMMON-NAME - Locusta migratoria                                                 *
#*  The important item:  Unique ID has to have TAX                                   *
#*************************************************************************************
def ReadClasses(CommonArea):
	dTaxOrganism = dict()
	iFile = CommonArea['iFile']
	InputFile = open(iFile)
	for iLine in InputFile: 
			if iLine.startswith("UNIQUE-ID"):
				bFlagTaxonomy = False			#Set up flag - Only interested in Taxonomy entries
			if iLine.startswith("UNIQUE-ID"):
				bFlagTaxonomy = True			#Set up flag - Only interested in Taxonomy entries
				UId = iLine.split("UNIQUE-ID - ")[1].rstrip()
			if iLine.startswith("COMMON-NAME - ") and  bFlagTaxonomy == True:
				CommonName = iLine.split("COMMON-NAME - ")[1].rstrip()
				dTaxOrganism[UId] = CommonName 
	CommonArea['dTaxOrganism'] = dTaxOrganism
   
	InputFile.close()
	return CommonArea



#*************************************************************************************
#* Read pathwaysdat                                                                  *
#* We are interested in the Unique Id and the Taxonomy entries related to it         *
#* An example of a Unique ID is:                                                     *
#*  UNIQUE-ID - PWY-6261                                                             *
#*  and their corresponding SPECIES    that look as follows:                         *
#*  SPECIES - TAX-2214                                                               *
#*  The important item:  Unique ID has to have TAX                                   *
#*************************************************************************************
def ReadPathways(CommonArea):
	dPathwaysOrganism = dict()
	iPathwaysFile = CommonArea['iFilePathways']
	InputFile = open(iPathwaysFile)
	for iLine in InputFile: 
			if iLine.startswith("UNIQUE-ID"):
				bFlagSpecies = False			#Set up flag - Only interested in Species entries
				lSpecies = list()				#Set up an empty list to store the species
				try:                            #Try to store the previous pathway --> Organisms relations
					dPathwaysOrganism[UPathwayId] = lSpecies
				except:
					pass
				UPathwayId = iLine.split("UNIQUE-ID - ")[1].rstrip()   #Get the patway id 
			if iLine.startswith("SPECIES - "):
				bFlagSpecies = True			#Set up flag - Only interested in Taxonomy entries
			if iLine.startswith("SPECIES - ") and  bFlagSpecies == True:
				Species = iLine.split("SPECIES - ")[1].rstrip()
				SpeciesName = CommonArea["dTaxOrganism"][Species]  #Get the species name
				lSpecies.append(SpeciesName)			#Add the current species
	CommonArea['dPathwaysOrganism'] = dPathwaysOrganism
   
	InputFile.close()
	return CommonArea


#*************************************************************************************
#* Parse Input parms                                                                 *
#*************************************************************************************
def read_params(x):
	CommonArea = dict()	
	parser = argparse.ArgumentParser(description='Read Swissprot and generate an exract with AC and all ECs related to it')
	parser.add_argument('--ic', action="store", dest='ic',nargs='?')
	parser.add_argument('--ip', action="store", dest='ip',nargs='?')
	parser.add_argument('--o', action="store", dest='o',nargs='?',  default='output_analysis')
	CommonArea['parser'] = parser
	return  CommonArea


 

	
 

#**************************************************************
#*      Build the Output file                                 *
#**************************************************************	
def  GenerateOutputFile(CommonArea):
	strTab = "\t"
	strNewLine = "\n"
	OutputLineCntr = 0
	OutputFile = open(CommonArea['oFile'],'w')
	for Pathway in sorted(CommonArea['dPathwaysOrganism'].keys()):
		lOutputRecord = [Pathway]
		strOrganisms = ""
		for Organism in CommonArea['dPathwaysOrganism'][Pathway]:
			strOrganisms = strOrganisms + Organism + "," #  Add the organism
		lOutputRecord.append(strOrganisms)  #Add the Organisms
		strBuiltRecord = strTab.join(lOutputRecord)[:-1] +  strNewLine
		OutputFile.write(strBuiltRecord)
	OutputFile.close()
	return CommonArea
	
#*************************************************************************************
#*  Main Program                                                                     *
#*************************************************************************************
print("Program started")
CommonArea = read_params( sys.argv )  # Parse command  
parser = CommonArea['parser'] 
results = parser.parse_args()

CommonArea['iFile'] = results.ic   #Classes.dat 
CommonArea['iFilePathways'] = results.ip   #Pathways file pathways.dat
CommonArea['oFile'] = results.o
Commonrea = ReadClasses(CommonArea)    #Load the mapping Unique_id --> Organism
Commonrea = ReadPathways(CommonArea)    #Load the mapping Unique_id --> Organism
CommonArea = GenerateOutputFile(CommonArea)  #Print the relationships

print("Program ended Successfully")
exit(0)
