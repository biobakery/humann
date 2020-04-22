#!/usr/bin/env python
import sys,string
import sys, os
import argparse



#********************************************************************************************
#    Build kegg extract                                                                     *
#                                                                                           *
#    The objective of this program is to build a Kegg extract (Kegg ID Mapping)             *
#    kegg_idmapping                                                                         *
# Format of the Output:                                                                     *
# reference \t gene \t length (nucleotides) \t bug                                          *
#  hsa:399693 \t K0# \t 953*3 \t bug                                                        *
#  takes 3 input files:                                                                     *
#  humann1/data/genels                                                                      *
#  A. Format (first 3 lines bug:gene \t length(aa)):                                        *
#     hsa:399693	953                                                                     *
#     hsa:401052	118                                                                     *
#     hsa:285033	121                                                                     *
#  B. humann1/data/koc (bug:gene to K0)                                                     *
#     Format:                                                                               *
#     K0# \t hsa:# \t hsa:#2                                                                *
#  C. kegg file with bug names: http://www.genome.jp/kegg/catalog/org_list.html             *
#                                                                                           *
#  The program is invoked as follows:                                                       *
#  python humann_humann1_kegg.py --igenels ./data/genels   --ikoc  ./data/koc   --o output_extract *
#                                                                                            *
#
#  The program tries to read the Kegg translation table OrgId -->OrgName that is located in  *
#  ../data/misc/KeggOrgId2OrgNameTable.txt                                                   *
#  And if that fails,  it tries to get that info from Kegg itself ( "http://www.genome.jp/kegg/catalog/org_list.html" )
#                                                                                            *
# Written by George Weingart  March 31, 2015   george.weingart@gmail.com                     *
#                                                                                            *
#                                                                                            *
#                                                                                            *
#   --------------------------------------------------------------------------------------   *
#  PROGRAM MODIFICATION: george.weingart@gmail.com  May 23, 2015                             *
#  ********************                                                                      *
#   --------------------------------------------------------------------------------------   *
#  Modified the program so that it includes the whole taxonomy for each of the bugs          *
#  The taxonomy is calculated as follows:                                                    *
#  Take the taxonomy name,  for example "Picrophilus torridus"  and use the ete2 ncbi        *
#  taxonimy functions to discern the lineage                                                 *
#  IMPORTANT note:                                                                           *
#  The program now requires ete2  to be installed                                            *
#  Because it invokes                                                                        *
#    from ete2 import NCBITaxa                                                               *
#    ncbi = NCBITaxa()                                                                       *
#                                                                                            *
#   Please see:                                                                              *
#     http://pythonhosted.org/ete2/tutorial/tutorial_ncbitaxonomy.html                       *
#     and                                                                                    *
#     http://pythonhosted.org/ete2/reference/reference_ncbi.html#ete2.NCBITaxa.get_name_translator
#                                                                                            *
#********************************************************************************************





#*************************************************************************************
#* Parse Input parms                                                                 *
#*************************************************************************************
def read_params(x):
	CommonArea = dict()	
	parser = argparse.ArgumentParser(description='Humann1_kegg extract generator')
	parser.add_argument('--igenels', action="store", dest='igenels',nargs='?')
	parser.add_argument('--ikeggtrans', action="store", dest='ikeggtrans',nargs='?', default="http://www.genome.jp/kegg/catalog/org_list.html")
	parser.add_argument('--ikoc', action="store", dest='ikoc',nargs='?')
	parser.add_argument('--ikeggOrgId2OrgName', action="store", dest='ikeggOrgId2OrgName',nargs='?')
	parser.add_argument('--o', action="store", dest='o',nargs='?')
	CommonArea['parser'] = parser
	return  CommonArea

	

#************************************************************************
#*   Load the dictionary of OrgId --> OrgName                           *
#* The routine takes its data from the following web-site               *
#*  "http://www.genome.jp/kegg/catalog/org_list.html"                   *
#************************************************************************
def  ReadHtmlKegTable(CommonArea):
		print("Trying to read the online kegg html table to translate OrgId to OrgName : http://www.genome.jp/kegg/catalog/org_list.html")
		from bs4 import BeautifulSoup
		import urllib2
		
		dOrgIdOrgName = dict()
		f1  =  urllib2.urlopen(CommonArea['ikeggtrans'])
		soup1 = BeautifulSoup(f1)

		for child in soup1.recursiveChildGenerator():
				name = getattr(child, "name", None)
				try:
					xchild = str(child)
					xname = str(name)
					if xname.startswith("a")  and xchild.startswith('<a href="/dbget-bin/www_'):
						ychild = str(xchild)
						OrgName = ychild.split(">")[1].split('<')[0]
						OrgNameOrgId = OrgId + "\t" +  OrgName 
						dOrgIdOrgName[OrgId] = OrgName
						
					if xname.startswith("a")     and xchild.startswith('<a href="/kegg-bin/show_organism?org='):
						ychild = str(xchild)
						OrgId = ychild.split(">")[1].split('<')[0]
				except:
					pass
		CommonArea['dOrgIdOrgName'] = dOrgIdOrgName 
		return CommonArea 



#************************************************************************
#*  Try to read the Sequential File                                     *
#************************************************************************
def  ReadSequentialTranslationFile(CommonArea):
		if  CommonArea['ikeggOrgId2OrgName'] is not None:
			CommonArea['ikeggOrgId2OrgName_file'] = open(CommonArea['ikeggOrgId2OrgName'], "rt")
		else:
			FileLocation = os.path.join(os.path.dirname(os.path.abspath(__file__)),os.pardir, "data/misc/KeggOrgId2OrgNameTable.txt")
			CommonArea['ikeggOrgId2OrgName_file'] = open(FileLocation, "rt")  
			
			
		for iKeggOrgLine in CommonArea['ikeggOrgId2OrgName_file']: 
			OrgId = iKeggOrgLine.split("\t")[0]
			OrgName = iKeggOrgLine.split("\t")[1].rstrip('\n')
			CommonArea['dOrgIdOrgName'][OrgId]  = OrgName 
		return CommonArea 
		
		
		
		

#************************************************************************
#*  Get the Taxonomy of an Organism                                     *
#************************************************************************
def  GetTaxonomy(OrgIdName):
	OrgIDExpandedTaxonomy = OrgIdName  #  Set the default
	try:
		lTaxonomyExlusionList = [1,131567]   #These taxonomies are obvious, so we  are not going to include them
		dTaxId = ncbi.get_name_translator([OrgIdName])	#Get the Tax Id number for the organism
		lLineage =  ncbi.get_lineage( dTaxId[OrgIdName])   #Get the Lineage for the bug
		lLineageNames = ncbi.get_taxid_translator(lLineage)  # Get the lineage names
		lExpandedOrganism = list()  
		for  LineageEntry  in lLineage:   #Get the Lineage entry in the right order
			if LineageEntry not in lTaxonomyExlusionList:
					LineageEntryName = lLineageNames[LineageEntry]    #  Build the list of the lineage names
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

CommonArea['igenels'] = results.igenels
CommonArea['igenels_file'] = open(CommonArea['igenels'], "rt")	# Open the igenels file

CommonArea['ikoc'] = results.ikoc
CommonArea['ikoc_file'] = open(CommonArea['ikoc'], "rt")	# Open the ikoc file

CommonArea['ikeggtrans'] = results.ikeggtrans

CommonArea['oFile'] = results.o
OutputFile = open(CommonArea['oFile'],'w')

CommonArea['ikeggOrgId2OrgName'] = results.ikeggOrgId2OrgName




CommonArea['dOrgIdOrgName'] = dict()
try:
	CommonArea  =  ReadSequentialTranslationFile(CommonArea)   #First,  try to read the sequential file in our directory
except:
	print("Reading the Table OrgId--> Name in ../data/misc/KeggOrgId2OrgNameTable.txt Failed - will try the kegg website")
	CommonArea  =  ReadHtmlKegTable(CommonArea)   #But if that fails - try directly from Kegg
 
Flagete2Installed = True
try:
	from ete2 import NCBITaxa
	ncbi = NCBITaxa()
except:
	print("ete2 is not installed and thus we cannot provide the full taxonomy for the bug, providing only the name")
	Flagete2Installed = False


	 
for OrgId, OrgIdName in CommonArea['dOrgIdOrgName'].iteritems():
	if  Flagete2Installed == True:   #If ete2 is installed, try to get the taxonomy
		OrgIDExpandedTaxonomy = GetTaxonomy(OrgIdName)
		if  OrgIDExpandedTaxonomy  !=  OrgIdName:  #If there was a taxonomy list
			CommonArea['dOrgIdOrgName'][OrgId] = OrgIDExpandedTaxonomy  # Update the entry in the dictionary
 







dBugGeneToKO = dict()
#************************************************************************
#*   Load the koc dictionary                                            *
#************************************************************************
for ikocLine in CommonArea['ikoc_file']: 
	kocEntries  = ikocLine.split('\t')
	KO = kocEntries[0]
	kocEntries.pop(0)

	for kocEntry in kocEntries:
		Bug = kocEntry.split("#")[0]
		Gene = kocEntry.split("#")[1].rstrip('\n')
		BugGene = Bug + ":" + Gene
		dBugGeneToKO[BugGene] = KO
		
		
 
for iLine in CommonArea['igenels_file']: 
	iLine = iLine.rstrip('\n')
	iBug = iLine.split(":")[0]
	iGene = iLine.split(":")[1].split("\t")[0]
	iGeneLengthAA = int(iLine.split("\t")[1])
	iGeneLengthNucleotide = 3 * iGeneLengthAA 
	iBugName = "n/a"
	try:
		iBugName = CommonArea['dOrgIdOrgName'][iBug]
	except:
		iBugName = iBug
	
	BugGene = iBug + ":" + iGene
	BugGene  = BugGene.upper()
	try:
		KOEntry = dBugGeneToKO[BugGene]
	except:
		KOEntry = "NoKONum"
		
	
	try:
		OrgName = CommonArea['dOrgIdOrgName'][iBug] 
	except: 
		OrgName = iBug	
		
	if  KOEntry != "NoKONum":
		OutputLine = iBug + ":" + iGene + "\t" + KOEntry + "\t" + str(iGeneLengthNucleotide) + "\t" + OrgName + "\n"
		OutputFile.write(OutputLine)
	
		
OutputFile.close()
print("Program ended Successfully")
exit(0)
