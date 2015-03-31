#!/usr/bin/env python
from cStringIO import StringIO
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
#  python humann1_kegg.py --igenels ./data/genels   --ikoc  ./data/koc --ikeggtrans  "http://www.genome.jp/kegg/catalog/org_list.html"  
#              --ikeggOrgId2OrgName    ../data/pathways/KeggOrgId2OrgNameTable.txt   --o output_extract
#                                                                                            *
#
#  The program tries to read the Kegg translation table OrgId -->OrgName that is located in  *
#  ../data/pathways/KeggOrgId2OrgNameTable.txt                                               *
#  And if that fails,  it tries to get that info from Kegg itself ( "http://www.genome.jp/kegg/catalog/org_list.html" )
#                                                                                            *
# Written by George Weingart  March 31, 2015   george.weingart@gmail.com                     *
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
		print "Trying to read the online kegg html table to translate OrgId to OrgName : http://www.genome.jp/kegg/catalog/org_list.html"
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
		CommonArea['ikeggOrgId2OrgName_file'] = open(CommonArea['ikeggOrgId2OrgName'])
		for iKeggOrgLine in CommonArea['ikeggOrgId2OrgName_file']: 
			OrgId = iKeggOrgLine.split("\t")[0]
			OrgName = iKeggOrgLine.split("\t")[1].rstrip('\n')
			CommonArea['dOrgIdOrgName'][OrgId]  = OrgName 
		return CommonArea 

#*************************************************************************************
#*  Main Program                                                                     *
#*************************************************************************************
print "Program started"

CommonArea = read_params( sys.argv )  # Parse command  
parser = CommonArea['parser'] 
results = parser.parse_args()

CommonArea['igenels'] = results.igenels
CommonArea['igenels_file'] = open(CommonArea['igenels'])	# Open the igenels file

CommonArea['ikoc'] = results.ikoc
CommonArea['ikoc_file'] = open(CommonArea['ikoc'])	# Open the ikoc file

CommonArea['ikeggtrans'] = results.ikeggtrans

CommonArea['oFile'] = results.o
OutputFile = open(CommonArea['oFile'],'w')

CommonArea['ikeggOrgId2OrgName'] = results.ikeggOrgId2OrgName




CommonArea['dOrgIdOrgName'] = dict()
try:
	CommonArea  =  ReadSequentialTranslationFile(CommonArea)   #First,  try to read the sequential file in our directory
except:
	print "Reading the Table OrgId--> Name in ../data/pathways/KeggOrgId2OrgNameTable.txt Failed"
	CommonArea  =  ReadHtmlKegTable(CommonArea)   #But if that fails - try directly from Kegg
 

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
		BugGene = Bug.lower() + ":" + Gene
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
	try:
		KOEntry = dBugGeneToKO[BugGene]
	except:
		KOEntry = "NoKONum"
		
	
	try:
		OrgName = CommonArea['dOrgIdOrgName'][iBug] 
	except: 
		OrgName = iBug	
	OutputLine = iBug + ":" + iGene + "\t" + KOEntry + \
	"\t" + str(iGeneLengthNucleotide) + "\t" + OrgName + "\n"
	OutputFile.write(OutputLine)
	
		
OutputFile.close()
print "Program ended Successfully"
exit(0)
