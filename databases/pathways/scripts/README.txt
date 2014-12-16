HUMAnN2 Test Scripts

Steps to create a HUMAnN2 pathways database set of files:

1. Using the UniProt pathways database
TBD

2. Using the MetaCyc pathways database

2.1. Download the meta.tar.gz of flat-files from MetaCyc.
A description of the files along with download instructions 
can be found at http://bioinformatics.ai.sri.com/ptools/flatfile-format.html

2.2. Decompress the download.
$ tar zxvf meta.tar.gz

NOTE: Instructions that follow refer to the metacyc directory as $METACYC.
If you are running on hutlab3, $METACYC=/n/huttenhower_lab_nobackup/downloads/metacyc/

2.3. Create the humann2/data/metacyc_reactions.uniref data file using Uniprot EC mapping.

2.3.1 Download and decompress the UniProtKB SwissProt text file.
$ wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz
$ gunzip uniprot_sprot.dat.gz

NOTE: If you are running on hutlab3, this file can be found at /n/huttenhower_lab/data/uniprot/2014-09/ .

2.3.2 Download and create the UniProt to UniRef50/90 mappings from the full set.
$ ./UniProt_mapping.py -i UniRef50 -o map_uniprot_UniRef50.dat.gz
$ ./UniProt_mapping.py -i UniRef90 -o map_uniprot_UniRef90.dat.gz

NOTE: If you are running on hutlab3, these files can be found at /n/huttenhower_lab/data/idmapping/ .

2.3.3 Create the metacyc_reations.uniref file (pathways database file1).
$ python Reaction_to_Uniref5090.py --i_reactions $METACYC/18.1/data/reactions.dat  --i_sprot uniprot_sprot.dat  --uniref50gz map_uniprot_UniRef50.dat.gz --uniref90gz map_uniprot_UniRef90.dat.gz  --o metacyc_reactions.uniref

2.4. Create the humann2/data/metacyc_pathways data file (pathways database file2).
$ ./metacyc2mcpc.py < $METACYC/18.1/data/pathways.dat > metacyc_pathways

 *********************************************************************************
2.5  ********      GENERATION OF THE PATHWAYS FILES        ********
 *********************************************************************************
There are two important pathway files whose source is Swissprot:

	a. unipathway_uniprots.uniref
	b. unipathway_pathways
	
FIRST FILE: unipathway_uniprots.uniref  (ReadSwissprot.py)
	
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
  
  
  
SECOND FILE: unipathway_pathways (Build_mapping_Pathways_Uniprot)

The objective of this program is to map  Swissprot Pathways to Uniprot ACs  and Uniref50 and 90    

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
Then it proceeds to read the Uniref50 and Uniref90 files in the same fashion as ReadSwisport.py does (Unzip the 50, 90 files,
glue them) and treats, like in the case of ReadSwissprot.py,  the AC table, as a transaction file and runs a Transaction vs.
Master process (AC Table vs. U5090 file of 80 million recs)  and this way updates the U50 and U90 for the particular AC and generates 
the extract:  Reaction, AC{s}, U50{s},U90{s}


To run the program:  
 python Build_mapping_Pathways_Uniprot.py --i /n/huttenhower_lab_nobackup/downloads/uniprot_pathways/2014_10/pathway.txt \
 --uniref50gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz \
 --uniref90gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef90.dat.gz \
 --oPathwaysACs  unipathway_pathways \
 --oPathwaysUniref5090 PathwaysUniref5090   
  
  
The input Uniprot files can be downloaded from the following sites:
  http://www.uniprot.org/downloads
  http://www.uniprot.org/docs/pathway
  ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/
  ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/
  

 *********************************************************************************
2.5  ********      End doc GENERATION OF THE PATHWAYS FILES        ********
 *********************************************************************************

2.6  Instructions to run the ReadUniref.py   
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

2.7  Instructions to run the  Build_mapping_Pathways_Uniprot.py

#********************************************************************************************
#    Map Pathways to Uniprot IDs  and Uniref50 and 90                                       *
#                                                                                           *
#    The objective of this program is to map  Swissprot Pathways to Uniprot ACs             *
#    and Uniref50 and 90                                                                    *
#                                                                                           *
#  -----------------------------------------------------------------------------------------*
#  Invoking the program:                                                                    *
#  ---------------------                                                                    *
#  python Build_mapping_Pathways_Uniprot.py --i /n/huttenhower_lab_nobackup/downloads/uniprot_pathways/2014_10/pathway.txt \
# --uniref50gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz\
# --uniref90gz /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef90.dat.gz\
# --oPathwaysACs  unipathway_uniprots.uniref \
# --oPathwaysUniref5090 PathwaysUniref5090                                            
#                                                                                           *
#   Where:                                                                                  *
#    --i_reactions, is the pathways  file, which is currently located at                    *
#    /n/huttenhower_lab_nobackup/downloads/uniprot_pathways/2014_10/pathway.txt             *
#  and it was downloaded from  the site:  http://www.uniprot.org/help/pathway               *
#                                                                                           *
#   --uniref50gz and --uniref90gz are the Uniref50 and 90 mappings,                         *
#      Currently located at /n/huttenhower_lab/data/idmapping/map_uniprot_UniRef50.dat.gz and 90
#                                                                                           *
#   --oPathwaysACs  is the Output file containing the Pathways --> ACs relations            *
#                                                                                           *
#  --oPathwaysUniref5090 is the output file containing the Pathways --> Uniref50/90 relations
#     ****NOTE****  If this parameter is not supplied,  this file is not created            *
#  
#   Written by George Weingart  Oct. 20, 2014   george.weingart@gmail.com                   *
#********************************************************************************************

2.8  Instructions to run the  Reaction_to_Uniref5090.py


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







