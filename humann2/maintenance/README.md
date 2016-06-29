## HUMAnN2 Maintenance Scripts ##

These scripts are used to build/maintain HUMAnN2 data files.
The typical user will not need to work with these.

### Steps to create a set of HUMAnN2 pathways database files ###

#### Create HUMANn2 MetaCyc pathways database files ####

1. Download and decompress the meta.tar.gz of flat-files from MetaCyc.
    * A description of the files along with download instructions can be found at http://bioinformatics.ai.sri.com/ptools/flatfile-format.html
    * To decompress: `` $ tar zxvf meta.tar.gz ``
    * NOTE: Instructions that follow refer to the metacyc directory as ``$METACYC``.
    * If you are running on hutlab machines, ``$METACYC=/n/huttenhower_lab_nobackup/downloads/metacyc/`` .

2. Create the humann2/data/pathways/metacyc_reactions_level4ec_only.uniref data file using Uniprot EC mapping.
    1. Download and decompress the UniProtKB SwissProt text file.
        * `` $ wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz ``
        * `` $ gunzip uniprot_sprot.dat.gz ``
        * If you are running on hutlab machines, this file can be found at /n/huttenhower_lab/data/uniprot/2014-09/ .

    2. Download and decompress the TrEMBL database and then concat with SwissProt for steps that follow.
        * `` $ wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz ``
        * `` $ gunzip uniprot_trembl.dat.gz ``
        * `` $ cat uniprot_sprot.dat uniprot_trembl.dat > uniprot_sprot_trembl.dat ``
        * This is a large file approx 130GB.

    3. Download and create the UniProt to UniRef50/90 mappings from the full set.
        * `` $ ./UniProt_mapping.py -i UniRef50 -o map_uniprot_UniRef50.dat.gz ``
        * `` $ ./UniProt_mapping.py -i UniRef90 -o map_uniprot_UniRef90.dat.gz ``
        * If you are running on hutlab machines, these files can be found at /n/huttenhower_lab/data/idmapping/ .

    4. Create the humann2/data/pathways/metacyc_reactions_level4ec_only.uniref file (default pathways database file1).
        * Create a reactions.dat file that only includes those with level 4 ECs and that also includes all reactions to AC mappings
            * `` $ ./map_reactions_to_uniprot.py --input-reactions $METACYC/19.1/data/reactions.dat --input-enzrxn $METACYC/19.1/data/enzrxns.dat --input-proteins $METACYC/19.1/data/proteins.dat --input-gene-links $METACYC/19.1/data/gene-links.dat --output reactions_level4ec_only.dat ``
        * Create the uniref to reactions file
            * `` $ python Reaction_to_Uniref5090.py --i_reactions reactions_level4ec_only.dat  --i_sprot uniprot_sprot_trembl.dat  --uniref50gz map_uniprot_UniRef50.dat.gz --uniref90gz map_uniprot_UniRef90.dat.gz  --o metacyc_reations_level4ec_only.uniref ``

3. Create the structured, filtered humann2/data/pathways/metacyc_pathways_structured_filtered file (default pathways database file2)
    * Create a structured pathways file
        * `` $ ./create_metacyc_structured_pathways_database.py --input $METACYC/19.1/data/pathways.dat --output metacyc_pathways_structured ``
    * Filter the structured pathways file
        * `` $ ./filter_pathways.py --input-pathways metacyc_pathways_structured --input-reactions metacyc_reactions_level4ec_only.uniref --output metacyc_structured_pathways_filtered ``

4. (Optional) Create the unstructured humann2/data/pathways/metacyc_pathways data file (optional pathways database file2).
    * `` $ ./metacyc2mcpc.py < $METACYC/19.1/data/pathways.dat > metacyc_pathways ``

#### Create HUMANn2 UniPathways pathways database files ####

1. Create the file humann2/data/pathways/unipathway_pathways (optional pathway database file2)
    * `` $ python Build_mapping_Pathways_Uniprot.py --i $UNIPROT/2015_05/pathway.txt --uniref50gz map_uniprot_UniRef50.dat.gz --uniref90gz map_uniprot_UniRef90.dat.gz --oPathwaysACs  unipathway_pathways --oValidACs  list_of_ACs --oPathwaysUniref5090 PathwaysUniref5090 ``   
        * The command above, and those that follow, ``$UNIPROT`` should be replaced will the full path to the uniprot_pathways folder.
        * For detailed documentation about what this script does, please refer to the script header.
    * The input Uniprot files can be downloaded from the following sites:
        * http://www.uniprot.org/downloads
        * http://www.uniprot.org/docs/pathway
        * ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/
        * ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/
    * If runnning on hutlab machines, the uniprot pathways download can be found at /n/huttenhower_lab_nobackup/downloads/uniprot_pathways/2014_10/ .

2. Create the file humann2/data/pathways/unipathway_uniprots.uniref (optional pathway database file1)
    * `` $ python ReadSwissprot.py --i uniprot_sprot.dat --o unipathway_uniprots.uniref --uniref50gz map_uniprot_UniRef50.dat.gz --uniref90gz map_uniprot_UniRef90.dat.gz ``
        * For detailed documentation about what this script does, please refer to the script header.

### Steps to create a diamond formatted UniRef database ###

#### Creating a full UniRef Database ####

1. Create the file available for download named uniref50_annotated_diamond.tar.gz (a default database for translated search)
    * `` $ ./create_uniref_database.py --input uniref50.fasta --output uniref50_annotated --format-database diamond --verbose ``
    * The uniref50.fasta can be downloaded from Uniprot from the following url: http://www.uniprot.org/downloads
    * This script will format the database for diamond and also annotate the fasta ids with the gene length.
    * To create a UniRef90 database, provide the file uniref90.fasta as input.

#### Create a EC filtered UniRef Database ####

1. Create a list of UniRef50 ids which have ECs of level 4
    * `` $ ./filter_uniref_ec_list.py --input map_EC_to_triplet_AC_U50_U90_Swissprot_and_Trembl.txt --uniref 50 --output uniref50_ec4_list.txt --min-ec-level 4 ``
    * The file map_EC_to_triplet_AC_U50_U90_Swissprot_and_Trembl.txt is generated by the script map_ECs_to_UniRef50_90.py .
    * To create a list of UniRef90 ids, use the option "--uniref 90". 

2. Create the file available for download named uniref50_EC4_filtered_annotated_diamond.tar.gz (a default database for translated search)
    * `` $ ./create_uniref_database.py --input uniref50.fasta --output uniref50_filtered_annotated --include-list uniref50_ec4_list.txt --format-database diamond --verbose ``
    * To create a UniRef90 database, provide the file uniref90.fasta as input along with the file uniref90_ec4_list.txt.
