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
