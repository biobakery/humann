#!/usr/bin/env python
import sys,string
import sys, os
import argparse
import tempfile



#********************************************************************************************
#   Map pathway to Organism                                                                 *
#                                                                                           *
#   This program creates extracts of the form: Pathway-->[Organism,Organism]                *
#   OR                                                                                      *
#   Reaction --> [Organism, Organism]                                                       *
#                                                                                           *
#   Depending whether the input extract is from patwhays.dat or reactions.dat               *
#   (See details below) -  but the important item is that the program is the SAME, and      *
#   no need to change a parm,  just point it to the reactions or pathways extract           *
#                                                                                           *
#  -----------------------------------------------------------------------------------------*
#  Invoking the program:                                                                    *
#  ---------------------                                                                    *
#   python Map_Pathway_Or_Reaction_To_Organism_from_Extract.py                              *
#    --i  organisms_pathways_extract \                                                      *                                                    
#   --o mapping_patwhays_to_organisms                                                       *
#                                                                                           *
#   OR,  for the reactions file                                                             *
#                                                                                           *
# python Map_Pathway_Or_Reaction_To_Organism_from_Extract.py                                *
#    --i  organisms_reactions_extract \                                                     *                                                    
#   --o mapping_reactions_to_organisms                                                      *
#  NOTES and EXPLANATIONS                                                                   *
#  ----------------------                                                                   *
#  The program runs against the entire download of the metacyc database which is 58Gb in    *
#      its zipped version and 402Gb unzipped                                                *
#  That directory contains 5,000 directories                                                *
#  We look in each directory and extract from the pathways.dat member in each directory     *
#      the organism name and the pathways associated for that organism                      *
#                                                                                           *
#  Therefore, before this program is run,  one must extract these entries and that is done  *
#      via the command (Ran in the root directory of the metacyc directories)               *
#   egrep -r --include "pathways.dat" '# Species:|UNIQUE-ID'  .>pathways_output_extract     *
#  The file output_extract  is going to be the input of this program                        *
#                                                                                           *
#  Sample entries of that input file look as follows                                        *
#./bbro1208658cyc/pathways.dat:# Species: Bordetella bronchiseptica MO149                   *
#./bbro1208658cyc/pathways.dat:UNIQUE-ID - PWY-7411                                         *
#./bbro1208658cyc/pathways.dat:UNIQUE-ID - PWY-6118                                         *
#                                                                                           *
# <<<<<<<<<<<<<<   NOTE   >>>>>>>>>>>>>>>>>>>>>>>                                           *
#  For the REACTIONS to Organisms file,  one must create the extract as follows             *
#                                                                                           *
#   egrep -r --include "reactions.dat" '# Species:|UNIQUE-ID'  .>reactions_output_extract   *
#                                                                                           *
#  and use the output of that command (reactions_output_extract) as input to this program   *
#  as follows:                                                                              *
#                                                                                           *
# python Map_Pathway_Or_Reaction_To_Organism_from_Extract.py                                *
#    --i  organisms_reactions_extract \                                                     *                                                    
#   --o mapping_reactions_to_organisms                                                      *                   
#                                                                                           *
#                                                                                           *
# <<<<<<<<<<<<<<   NOTE   >>>>>>>>>>>>>>>>>>>>>>>                                           *
#                                                                                           *
#   Written by George Weingart - george.weingart@gmail.com   8/5/2015                       * 
#   Modified on 8/24/15 to use the same program to creat the reactions ->Organisms          *
#   (The program did not change - just the documentation                                    * 
#                                                                                           *
#********************************************************************************************

#*************************************************************************************
#* Read extract                                                                      *
#*************************************************************************************
def ReadExtract(CommonArea):
	dPathwayOrganism = dict()
	iFile = CommonArea['iFile']
	InputFile = open(iFile)

	for iLine in InputFile: 
			iData = iLine.split(":")[1]
			if iData.startswith("# Species"):
				Species = iLine.split("Species: ")[1].rstrip()	# Get the Species
			if iData.startswith("UNIQUE-ID"):
				Pathway = iData.split(" - ")[1].rstrip()	# Get the pathway name
				if  Pathway not in dPathwayOrganism:	# If not there - create it
					dPathwayOrganism[Pathway] = list()	# Create an empty list
				if Species != "MetaCyc":
					dPathwayOrganism[Pathway].append(Species)	# Append the Species
 
	CommonArea['dPathwayOrganism'] = dPathwayOrganism
   
	InputFile.close()
	return CommonArea

#*************************************************************************************
#* Parse Input parms                                                                 *
#*************************************************************************************
def read_params(x):
	CommonArea = dict()	
	parser = argparse.ArgumentParser(description='Read Swissprot and generate an exract with AC and all ECs related to it')
	parser.add_argument('--i', action="store", dest='i',nargs='?')
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
	for Pathway in sorted(CommonArea['dPathwayOrganism'] .keys()):
		lOutputRecord = [Pathway]
		strOrganisms = ""
 
		for Organism in CommonArea['dPathwayOrganism'][Pathway]:
			strOrganisms = strOrganisms + Organism + strTab #  Add the organism
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

CommonArea['iFile'] = results.i   #  The input file
CommonArea['oFile'] = results.o
Commonrea = ReadExtract(CommonArea)    # Read the file  
CommonArea = GenerateOutputFile(CommonArea)  #Print the relationships

print("Program ended Successfully")
exit(0)
