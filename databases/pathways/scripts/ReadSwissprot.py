#!/usr/bin/env python
from cStringIO import StringIO
import sys,string
import sys, os
import argparse



#********************************************************************************************
#    Read Swissprot program                                                                 *
#    This program reads the unprot.dat file and creates an                                  *
#    extract containing in each line                                                        *
#    The Protein AC and all the ECs related to it                                           *

#  -----------------------------------------------------------------------------------------*
#  Invoking the program:                                                                    *
#  ---------------------                                                                    *
#   python ReadSwissprot.py  --i  input_file  --o output_file                               *
#   Where:                                                                                  *
#   --i input_file is the UniprotKB Swissprot text file, which can be downloaded from       *
#    ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz *
#                                                                                           *
#   The current downloaded file, which serves as input,  resides on hutlab3 in              *
#    /n/huttenhower_lab/data/uniprot/2014-09/uniprot_sprot.dat                              *
#                                                                                           *
#   The current output of thei program is stored in humann2_test directory under the name   *
#   swissprot_mapping_AC_to_EC                                                              * 
#                                                                                           *
#   Written by George Weingart - george.weingart@gmail.com   10/06/2014                     *  
#********************************************************************************************








#*************************************************************************************
#* Read Swissprot                                                                    *
#*************************************************************************************
def ReadSwissprot(i,o):
	iFile = i
	oFile = o
	strTab = "\t"
	strNewLine = "\n"
	bFlagEC = False
	lECs = list()
	InputFile = open(iFile)
	OutputFile = open(oFile,'w')
	LineCntr = 0
	OutputLineCntr = 0
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
							OutputLine =  OutputLine + strTab + EC  
						OutputLine = OutputLine + strNewLine
						OutputFile.write(OutputLine )
						OutputLineCntr = OutputLineCntr + 1
				   bFlagEC = False
				   lECs = list()
				   lACs = list() 
	print "Read " + str(LineCntr) + " Input Lines"
	print "Wrote " + str(OutputLineCntr) + " Output Lines"
	InputFile.close()
	OutputFile.close()
	return 0











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









#*************************************************************************************
#*  Main Program                                                                     *
#*************************************************************************************
print "Program started"
CommonArea = read_params( sys.argv )  # Parse command  
parser = CommonArea['parser'] 
results = parser.parse_args()

iFile = results.i
oFile = results.o
RC  = ReadSwissprot(iFile,oFile)



print "Program ended Successfully"
exit(0)
