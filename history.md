
# HUMAnN2 History #

## v0.3.0 08-20-2015 ##

### New Features ###

* Updated databases to MetaCyc v19.1
* Added UniProt TrEMBL data to reactions mappings
* Added gap filling
* Added automated install of dependencies (bowtie2 and diamond)
* Added alert for user when running with demo databases and non-demo input file
* Added new option "--names {ko,metacyc-pwy,uniref50}" to humann2_rename_table
* Added new option "--groups {uniref50_ko,uniref50_go,uniref50_ec,uniref50_rxn}" to humann2_regroup_table
* Added new option "--output-max-decimals" to set the number of decimals written to the output files

### Interface Changes ###

* The default mapping file format for humann2_regroup_table has been reversed (it is now features to groups)
* The humann2_rename_table option for a custom mapping file was changed from "--names" to "--custom"
* The humann2_regroup_table option for a custom mapping file was changed from "--groups" to "--custom" 

### Bug Fixes ###

* Updated fastq to fasta function to allow for "@" (special fastq character) at beginning of quality score lines

### Cross Platform Compatibility Updates ###

* Updated humann2_join_table to ignore dot files (ie .DS_Store on Apple OS X)

## v0.2.2 07-27-2015 ##

### Cross Platform Compatibility Updates ###

* Added to install process the setting execute permissions on glpsol binaries
    * On some platforms when glpsol is installed by setuptools it does not keep its execute permissions

## v0.2.1 07-23-2015 ##

### New Features ###

* Added "--metaphlan-options"
    * Allows user to set any options for MetaPhlAn2 including database locations
    * Example: --metaphlan-options="-t rel_ab"
* Error message reporting has been updated to include messages from external tools
    * Example: "metaphlan2.py: error: unrecognized arguments: --stat_e 1.0"
* Updated MetaCyc structured pathways database
    * Removed reactions with only non-specific ECs
    * Set unmappable reactions (those without links to gene families) as optional reactions
    * Filtered out small pathways

### Interface Changes ###

* Changed option names to be more general for users with custom flows
    * Changed "--chocophlan" to "--nucleotide-database"
    * Changed "--uniref" to "--protein-database"

### Bug Fixes ###

* Updated function which determines file format to allow for gene tables with scientific notation

## v0.2.0 07-01-2015 ##

### New Features ###

* Added MetaCyc structured pathways database
* Added structured pathways computations for abundance and coverage

## v0.1.10 06-30-2015 ##

### New Features ###

* Added pathways and gene family names to output files
* Changed from e-value to number of matches for alignment scoring 

### Bug Fixes ###

* Allow for non-comment, empty lines in gene tables

## v0.1.9 05-01-2015 ##

### New Features ###

* Updated download script to pull new ChocoPhlAn version
* Added humann2_reduce_table script

### Performance Updates ###

* Reduced memory and increased speed for humann2_join_tables script

## v0.1.8 04-30-2015 ##

### New Features ###

* Updated pathways to UniProt 2015_05
* Allow for a directory, bowtie2 index basename, or bowtie2 index file as input

## v0.1.7 04-24-2015 ##

### Bug Fixes ###

* Updated function which finds executable in $PATH to return first instance instead of last

## v0.1.6 04-20-2015 ##

### Performance Updates ###

* Reduced memory use for humann2_split_table

### Other ###

* Added python version check to install

## v0.1.5 04-09-2015 ##

### Performance Updates ###

* Reduced memory use of humann2 (and added --memory-use option)

### New Features ###

* Allow for input files from Picrust metagenome_contributions.py

