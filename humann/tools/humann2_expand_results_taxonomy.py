#!/usr/bin/env python
import sys,string
import sys, os
import argparse
try:
	from ete2 import NCBITaxa
	ncbi = NCBITaxa()
except:
	print("ete2 is required to run this program")
	print("See http://pythonhosted.org/ete2/tutorial/tutorial_ncbitaxonomy.html")
	print("Please install ete2 and rerun")
	exit()



#********************************************************************************************
#    Expand taxonomy of bugs in humann results                                             *
#                                                                                           *
#    The objective of this program is  to expand the taxonomy for bugs in the               *
#    three humann results files:                                                           *
#    1. pathcoverage.tsv                                                                    *
#    2. pathabundance.tsv                                                                   *
#    3. pathabundance.tsv                                                                   *
#                                                                                           *
#                                                                                           *
#                                                                                           *
#  The program is invoked as follows:                                                       *
#  python humann_expand_results_taxonomy.py --input   input1 --ootput o1                   *
#                                                                                           *
#  Where input1 can be any of the three files mentioned above                               *
#  and o1 will be the output                                                                *
#                                                                                           *
#  It is a standalone program                                                               *  
#                                                                                           *
#  Sample of an input record:                                                               *     
#  UniRef50_Q8A776|g__Bacteroides.s__Bacteroides_thetaiotaomicron	5.89970501479           * 
#  The corresponding output record:                                                         *   
#  UniRef50_Q8A776|Bacteria_Bacteroidetes/Chlorobi group_Bacteroidetes_Bacteroidia_Bacteroidales_Bacteroidaceae_Bacteroides_Bacteroides thetaiotaomicron	5.8997050147
#                                                                                            *
#  Only records that have bugs are modified - the others are written as is                   *
#                                                                                            *
#                                                                                            *
#  The taxonomy of a bug is calculated as follows:                                           *
#  Take the taxonomy name,  for example "Picrophilus torridus"  and use the ete2 ncbi        *
#  taxonimy functions to discern the lineage                                                 *
#  IMPORTANT note:                                                                           *
#  The program   requires ete2  to be installed                                              *
#  Because it invokes                                                                        *
#    from ete2 import NCBITaxa                                                               *
#    ncbi = NCBITaxa()                                                                       *
#                                                                                            *
#   Please see:                                                                              *
#     http://pythonhosted.org/ete2/tutorial/tutorial_ncbitaxonomy.html                       *
#     and                                                                                    *
#     http://pythonhosted.org/ete2/reference/reference_ncbi.html#ete2.NCBITaxa.get_name_translator
#                                                                                            *
# Written by George Weingart  June 8th, 2015   george.weingart@gmail.com                     *
#********************************************************************************************





#*************************************************************************************
#* Parse Input parms                                                                 *
#*************************************************************************************
def read_params(x):
	CommonArea = dict()	
	parser = argparse.ArgumentParser(description='Humann1_kegg extract generator')
	parser.add_argument('--input', action="store", dest='i',nargs='?')
	parser.add_argument('--output', action="store", dest='o',nargs='?')
	CommonArea['parser'] = parser
	return  CommonArea

   
		
		

#************************************************************************
#*  Get the Taxonomy of an Organism                                     *
#************************************************************************
def  GetTaxonomy(OrgIdName):
 
	OrgIDExpandedTaxonomy = OrgIdName  #  Set the default
	try:
		lTaxonomyExlusionList = [1,131567]   #These taxonomies are obvious, so we  are not going to include them
		dTaxId = ncbi.get_name_translator([OrgIdName])	#Get the Tax Id number for the organism
		lLineage =  ncbi.get_lineage( dTaxId[OrgIdName])   #Get the Lineage for the bug
		dRanksOfLineage = ncbi.get_rank(lLineage)   #Get the ranks of the lineages
		lLineageNames = ncbi.get_taxid_translator(lLineage)  # Get the lineage names
		lExpandedOrganism = list()  
		for  LineageEntry  in lLineage:   #Get the Lineage entry in the right order
			if LineageEntry not in lTaxonomyExlusionList:
 
					LineageEntryName = lLineageNames[LineageEntry]    #  Build the list of the lineage names
					if  dRanksOfLineage[LineageEntry] !=  "no rank":  # If there is a rank for this lineage level
						LineageEntryName  =  dRanksOfLineage[LineageEntry]  + "__" + LineageEntryName # Append the lineage level 
						lExpandedOrganism.append(LineageEntryName)   #  add the lineage
		OrgIDExpandedTaxonomy = "|".join(lExpandedOrganism)   # Post to the new org Name
	except:
		pass		
	return  OrgIDExpandedTaxonomy



#*************************************************************************************
#*  Main Program                                                                     *
#*************************************************************************************
print("Program started")


CommonArea = read_params( sys.argv )  # Parse command  
parser = CommonArea['parser'] 
results = parser.parse_args()

CommonArea['ifile'] = results.i 
CommonArea['oFile'] = results.o
OutputFile = open(CommonArea['oFile'],'w')
CommonArea['InputFile'] = open(CommonArea['ifile'], "rt")
 
for sLine in CommonArea['InputFile']: 

	sLineBrokenByPipe = sLine.split('|')
	if len(sLineBrokenByPipe) == 1: 	 #If  it does not contain pipe - Just write the record and continue
		OutputFile.write(sLine)
		continue
	if  not sLineBrokenByPipe[1].startswith('g__'):		#If it does not contain a bug - continue 
		OutputFile.write(sLine)
		continue
 
	sB1 = sLine.split('.s__')[1]    #Split to get the bug name
	sB2 = sB1.split(' ') 
	sB3 = sB2[0].split('\t') 
	sBugName = sB3[0].replace("_"," ")
	BugAndTaxonomy = GetTaxonomy(sBugName)
	RebuiltOutputLine = sLineBrokenByPipe[0] + "\t" + BugAndTaxonomy + "\t" + sB3[1]  
	OutputFile.write(RebuiltOutputLine)
 
	
		
OutputFile.close()
print("Program ended Successfully")
exit(0)
