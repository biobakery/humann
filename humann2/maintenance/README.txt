HUMAnN2 Maintenance Scripts

These scripts are used to build/maintain HUMAnN2 data files.
The typical user will not need to work with these.

Steps to create a set of HUMAnN2 pathways database files:

1. Create a set of HUMANn2 MetaCyc pathways database files

1.1. Download the meta.tar.gz of flat-files from MetaCyc.
A description of the files along with download instructions 
can be found at http://bioinformatics.ai.sri.com/ptools/flatfile-format.html

1.2. Decompress the download.
$ tar zxvf meta.tar.gz

NOTE: Instructions that follow refer to the metacyc directory as $METACYC.
If you are running on hutlab3, $METACYC=/n/huttenhower_lab_nobackup/downloads/metacyc/

1.3. Create the humann2/data/metacyc_reactions.uniref data file using Uniprot EC mapping.

1.3.1 Download and decompress the UniProtKB SwissProt text file.
$ wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
$ gunzip uniprot_sprot.dat.gz

NOTE: If you are running on hutlab3, this file can be found at /n/huttenhower_lab/data/uniprot/2014-09/ .

1.3.2 Download and decompress the TrEMBL database and then concat with SwissProt for steps that follow.
$ wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz
$ gunzip uniprot_trembl.dat.gz
$ cat uniprot_sprot.dat uniprot_trembl.dat > uniprot_sprot_trembl.dat

NOTE: This is a large file approx 130GB.

1.3.3 Download and create the UniProt to UniRef50/90 mappings from the full set.
$ ./UniProt_mapping.py -i UniRef50 -o map_uniprot_UniRef50.dat.gz
$ ./UniProt_mapping.py -i UniRef90 -o map_uniprot_UniRef90.dat.gz

NOTE: If you are running on hutlab3, these files can be found at /n/huttenhower_lab/data/idmapping/ .

1.3.4 Create the metacyc_reations_level4ec_only.uniref file (default pathways database file1).
$ ./map_reactions_to_uniprot.py --input-reactions $METACYC/18.1/data/reactions.dat --input-enzrxn $METACYC/18.1/data/enzrxns.dat --input-proteins $METACYC/18.1/data/proteins.dat --input-gene-links $METACYC/18.1/data/gene-links.dat --output reactions_level4ec_only.dat
$ python Reaction_to_Uniref5090.py --i_reactions reactions_level4ec_only.dat  --i_sprot uniprot_sprot_trembl.dat  --uniref50gz map_uniprot_UniRef50.dat.gz --uniref90gz map_uniprot_UniRef90.dat.gz  --o metacyc_reations_level4ec_only.uniref

1.4 Create the structured, filtered humann2/data/metacyc_pathways_structured_filtered file (default pathways database file2)
$ ./create_metacyc_structured_pathways_database.py --input $METACYC/18.1/data/pathways.dat --output metacyc_pathways_structured
$ ./filter_pathways.py --input-pathways metacyc_pathways_structured --input-reactions metacyc_reactions_level4ec_only.uniref --output metacyc_structured_pathways_filtered

1.5. (Optional) Create the unstructured humann2/data/metacyc_pathways data file (optional pathways database file2).
$ ./metacyc2mcpc.py < $METACYC/18.1/data/pathways.dat > metacyc_pathways

2. Create a set of HUMANn2 UniPathways pathways database files

2.1 Create the first file (unipathway_pathways) by running Build_mapping_Pathways_Uniprot.py

The objective of this program is to map Swissprot Pathways to Uniprot ACs    

This program reads the Swissprot pathways file /n/huttenhower_lab_nobackup/downloads/uniprot_pathways/2014_10/pathway.txt
that looks as follows:
****
Alkaloid biosynthesis; 3alpha(S)-strictosidine biosynthesis; 3alpha(S)-strictosidine from secologanin and tryptamine: step 1/1
     STS1_ARATH  (P94111)    , STS3_ARATH  (P92976)    , STSY_CATRO  (P18417)    ,
     STSY_RAUMA  (P68174)    , STSY_RAUSE  (P68175)
Alkaloid biosynthesis; ajmaline biosynthesis
     PNAE_RAUSE  (Q9SE93)
****
And builds the relations: AC --> Reaction and Reaction --> AC 
It also builds an extract file controlled by the parameter --o ValidACs    which contains a list of the ACs that were output 
  **** This means that if all files need to be generated, the first step must be run first ****

At this point,  it generates the unipathway_pathways file and it can complete here.
However, it has the option to generate also a file with relations: Reaction --> Uniref50 and 90
If so,  it proceeds to read the Uniref50 and Uniref90 files in the same fashion as ReadSwisport.py does (Unzip the 50, 90 files,
glue them) and treats, like in the case of ReadSwissprot.py,  the AC table, as a transaction file and runs a Transaction vs.
Master process (AC Table vs. U5090 file of 80 million recs)  and this way updates the U50 and U90 for the particular AC and generates 
the extract:  Reaction, AC{s}, U50{s},U90{s}

To run the program:  
 python Build_mapping_Pathways_Uniprot.py --i /n/huttenhower_lab_nobackup/downloads/uniprot_pathways/2014_10/pathway.txt \
 --uniref50gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz \
 --uniref90gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef90.dat.gz \
 --oPathwaysACs  unipathway_pathways \
 --oValidACs  ../list_of_ACs \
 --oPathwaysUniref5090 PathwaysUniref5090   
  
  
The input Uniprot files can be downloaded from the following sites:
  http://www.uniprot.org/downloads
  http://www.uniprot.org/docs/pathway
  ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/
  ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/
 
2.2 Create the second file unipathway_uniprots.uniref by running ReadSwissprot.py
	
unipathway_uniprots.uniref is generated in the program ReadSwissprot.py by reading the uniprot_sprot.dat to gather the AC and all
ECs related to it,  and by gathering the Uniref50 and Uniref90 cluster names associated to it 
from the map_uniprot_UniRef50.dat.gz and map_uniprot_UniRef90.dat.gz files.
The logic of the program consists of building a table of ACs where each AC entry
contains all ECs related to it - this is gathered from uniprot_sprot.dat.
Once this table is built, we proceed to process Uniref50 and 90 files as following:
   a. We create a temporary directory 
   b Unzip the U50 and U90 files
   c. Because they contain exactly the same number of records (One record per AC, containing is Uniref cluster)
      we literally glue them together and generate a U50,90 singe file
   d. Because the U5090 file and the AC table that we gathered before from Swissprot are 
      sorted in the same sequence (AC)  we treat the AC table as a Transaction File (~250,00 records)
      which is processed against a Master (U5090 80 million records).
      When there is a match in the keys, we build the record: AC, EC{s}, U50,U90
      where the AC, EC(s) are taken from the table and the U50,90 are taken from the matched U5090 rec.

To run the job:	  
$ python ReadSwissprot.py \
  --i   /n/huttenhower_lab/data/uniprot/2014-09/uniprot_sprot.dat \
  --o unipathway_uniprots.uniref \
  --uniref50gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz \
  --uniref90gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef90.dat.gz
  
The input Uniprot files can be downloaded from the following sites:
  http://www.uniprot.org/downloads
  ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
  ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/
  ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/
  
  
  
How to run other scripts:  

1.  Instructions to run the ReadUniref.py   
#********************************************************************************************
#    Read Uniref Program                                                                    *
#    This program reads the mappings uniprot --> Uniref50                                   *
#    and uniprot --> Uniref90                                                               *
#    that currently reside in: /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz*
#    and /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef90.dat.gz                      *
#  -----------------------------------------------------------------------------------------*
#  Invoking the program:                                                                    *
#  ---------------------                                                                    *
#   python ReadUniref.py  map_uniprot_UniRef50.dat.gz map_uniprot_UniRef90.dat.gz  mcc outx *
#   Where:                                                                                  *
#   The first two files are the input mappings                                              *
#   The third file is the input mcc                                                         *
#   The fourth file is the output converted mcc file                                        *
#                                                                                           *
#                                                                                           *
#    Map to Uniref50 looks as follows:                                                      *
#   A0A008IWE9      UniRef50_I0C253                                                         *
#   A0A008IWF0      UniRef50_B0FFT6                                                         *
#   A0A008IWF1      UniRef50_Q5HRF6                                                         *
#   A0A008IWF2      UniRef50_O30875                                                         *
#                                                                                           * 
#    Map to Uniref90    looks as follows:                                                   *
#   A0A008IWE9      UniRef90_I0C253                                                         *
#   A0A008IWF0      UniRef90_X5E5F6                                                         *
#   A0A008IWF1      UniRef90_Q2G0I9                                                         *
#   A0A008IWF2      UniRef90_A5IQZ6                                                         * 
#                                                                                           *
#  The objective of the program is to translate the mcc file so that instead of             *
#  entries showing the uniprot id for a pathway,  we show the uniref50 and uniref90         *
#                                                                                           *
#  For example:   The following record                                                      *
#  BLASTICIDIN-S-DEAMINASE-RXN	EC-3.5.4.23	P33967	P78986                                  *      
#  will be converted to:                                                                    *
#  BLASTICIDIN-S-DEAMINASE-RXN EC-3.5.4.23 UniRef50_P33967 UniRef90_P33967                  *
#  Explanation:                                                                             *
#  P33967  maps in the translation files to UniRef50_P33967 P33967  UniRef90_P33967         *
#  and                                                                                      *
#  P78986 does not have a translation - thus no Uniref50 and Uniref90 entries were created  *
#  for P78986                                                                               *
#  LOGIC:                                                                                   *
#  The program uploads the Uniref translation files into a dictionary  and then reads       *
#  sequentially the mcc file and tries to convert each one of the Uniprot mappings to       *
#  the corresponding Uniref entry                                                           *
#  Logic of the Load of the table:                                                          *
#  The program creates a temporary directory and gunzips the input translation files        *
#  (sys.argv[1] and sys.argv[2])                                                            *
#  It then pastes the two of them together (In the temporary directory) and reads the       *
#  pasted file, uploading each of the Uniprot IDs as a key to the dictionary  with the      *
#  corresponding Uniref50 and 90 entries as a list in that dictionary entry                 *
#                                                                                           *
#  After the dictionary is built,  the temporary directory is deleted and we read the       *
#  mcc file looking for the uniprot id key in the dictionary and retrieving the             *
#  corresponding Uniref50 and 90 translations and posting that record into the output       *
#                                                                                           *
#   Written by George Weingart - george.weingart@gmail.com   8/28/2014                      *  
#********************************************************************************************


