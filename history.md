
# HUMAnN History #

## v3.6.1 (03-02-2023) ##

* Added a new utility script to expand a UniRef90 gene family to all those included in its associated UniRef50 clusters.
* Improve error message that checks MetaPhlAn version.

## v3.6 (09-26-2022) ##

* Updated to require diamond v2.0+.

## v3.5 (08-22-2022) ##

* Updated to run with MetaPhlAn v4.0 (includes mapping of SGBs to species plus uses additional species column).

## v3.1.1 (07-27-2022) ##

* Small fix to config to check for protein database.

## v3.1.0 (07-26-2022) ##

* Update to sync with the latest MetaPhlAn database (v3.1).
  - Detection of MetaPhlAn v3+ versions.
  - ChocoPhlAn updated to include latest pangenomes.

## v3.0.1 (11-18-2021) ##

* Update default to point to MetaCyc v24 pathways database.
* Reduce the length of the legend and fix two typos in humann_barplot script.
* Fix humann case and version in humann_databases script stdout.
* Update humann_barplot options in readme.

## v3.0.0 (04-22-2021) ##

* ChocoPhlAn updated to include an additional 600 pangenomes.
* Updated to support the latest diamond version v0.+ (including default diamond install version plus databases).
* Utility mapping databases have been expanded to include additional EC mappings.
* Custom nucleotide database selection has been updated to only include pangenomes based on strict matching.
* Updates to the barplot script including support for matplotlib v3. 

## v3.0.0.alpha.5 (02-02-2021) ##

* Update pathways from v19 to Metacyc v24.
* Allow for shell in benchmark subprocess to allow for more commands.
* Use string instead of set to reduce memory usage by blastx_coverage.

## v3.0.0.alpha.4 (07-09-2020) ##

* Add option for user to provide database location (merged from v2.8 branch).

## v3.0.0.alpha.3 (06-24-2020) ##

* Change bowtie2 version check from error to warning to allow for conda bowtie2 package without version information.

## v3.0.0.alpha.2 (05-19-2020) ##

* Fix read in MetaPhlAn v3 abundances to allow for scientific notation.

## v3.0.0.alpha.1 (04-27-2020) ##

* Change in metaphlan taxonomic profile version tag to stay in sync.

## v3.0.0.alpha (04-23-2020) ##

* Add flexibility to allow for additional MetaPhlAn output column.
* Split "--identity-threshold" into two options, one for each search type with defaults of 0 (nucleotide) and 80/50 (translated).
* Update diamond default options to "--id 80.0 --top 1" for UniRef90 runs and "--id 50.0 --top 1" for UniRef50 runs.
* Updated databases to use the latest ChocoPhlAn and corresponding UniRef50/90 (release 01/2019).
* Added diamond options.
* Added bowtie2 options.
* Added gene based filtering to nucleotide search using the same method as translated search. To revert back to prior nucleotide filtering mode, run with the options "--nucleotide-subject-coverage-threshold 0 --nucleotide-query-coverage-threshold 0".
* Add flexibility to MetaPhlAn version check to allow for warning messages for numpy and biom-format.

## v2.8.2 04-03-2020 ##

* Add an option for the user to provide the location of the database.

## v2.8.1 07-12-2019 ##

* Updated MetaPhlAn2 commands to use new legacy database name (to be compatible with MetaPhlAn2 v2.9.13+).

## v2.8.0 07-01-2019 ##

* Updated MetaPhlAn2 commands to check for software and database versions required.

## v0.11.2 10-12-2018 ##

* Updated MetaPhlAn2 command to be compatible with the latest version.
* Fix executable search to allow for folders in $PATH to have the same names as executables (ie, bowtie2).
* Fix error in header format for humann2_unpack_pathways tool.

## v0.11.1 04-13-2017 ##

* Fix the write table function to allow the tools with optional output arguments (renorm, rename, regroup, infer_taxonomy) to write to stdout if an output path is not provided.

## v0.11.0 03-29-2017 ##

* Add a check and remove of spaces if present in fastq/fasta input files to address new illumina casava v1.8 format.

## v0.10.0 02-22-2017 ##

* The query threshold filtering computation was updated so strand directionality does not affect results (ie reverse was two percent less than forward equivalent).

## v0.9.10 02-21-2017 ##

* The query threshold filtering computation was updated to allow for query starts that are larger than query end indexes.

## v0.9.9 12-15-2016 ##

### Other Changes ###

* HUMAnN submitted as of this release.
* Renamed executable humann2_merge_abundance_tables to humann2_unpack_pathways.
* Removed executable blastx_coverage because only module is currently used.
* humann2_associate is now more robust to user-provided/misspecified/missing metadata.
* Descriptions of HUMAnN utility scripts and "tutorials" in the manual have been standardized as "Guides to HUMAnN utility scripts" and "Other HUMAnN guides." These include several new guides for previously available scripts, e.g. humann2_barplot.

## v0.9.8 11-29-2016 ##

### Other Changes ###

* Help messages were added to the humann2_benchmark tool.

## v0.9.7 11-22-2016 ##

### New Features ###

* A new script, humann2_benchmark, was added which captures MaxRSS and elapsed time for any command and all subprocesses spawned.

### Other Changes ###

* The gap fill option default was set to on.

## v0.9.6 11-10-2016 ##

### New Features ###

* Viral pangenomes were added to the full chocophlan database download.

### Other Changes ###

* The biom output format now uses the biom API to allow for empty matrices.
* All sequences that do not pass filtering are considered unaligned (even if they have alignments to the reference).

## v0.9.5 10-24-2016 ##

### New Features ###

* The script humann2_infer_taxonomy has been updated to enable assignment of approximate taxonomic annotations to a greater proportion of unclassified UniRef90 and UniRef50 stratifications. To use the updated script, please also update your HUMAnN utility mapping files (humann2_databases --download utility_mapping full $DIR).
* The script split stratified table has been updated to be compatible with gzip, bzip2, and biom formats.

## v0.9.4 10-04-2016 ##

### New Features ###

* Added option "--remove-column-description-output" which will remove the description from the output file columns.
* Added build simple method for diamond install for users without cmake installed.

## v0.9.3 09-21-2016 ##

### New Features ###

* Added biom input and output compatibility to the following scripts: humann2_infer_taxonomy, humann2_regroup_table, humann2_rename_table, humann2_renorm_table, humann2_rna_dna_norm, and humann2_strain_profiler.
* Added gzip and bz2 input file compatibility to the following scripts: humann2_join_tables, humann2_merge_abundance, and humann2_reduce_table.

### Other Changes ###

* Added functional end to end and also tools tests for biom compatibility features.
* Sorted the reactions and pathways provided to MinPath. This removes stochasticity seen when running with Python3.
* Changed diamond version check option to "--version" as this is backwards compatible with the diamond v0.7 series.

## v0.9.2 09-13-2016 ##

### New Features ###

* Software and databases were upgraded to be compatible with the latest Diamond (v0.8.22).
* Added new humann2_regroup_table mappings between EggNOG and UniRef50/90.
* Added new humann2_rename_table mappings for Pfam, EggNOG, GO, and Informative GO.
* Database files for humann2_infer_taxonomy are now discoverable after downloading the utility mapping dataset.

### Bug Fix ###

* Informative GO now uses bare GO identifiers to avoid a conflict with HUMAnN's stratification syntax. Names can be attached to Informative GO identifiers using humann2_rename_table.
* Remove requirement of future package for python2.

## v0.9.1 08-26-2016 ##

### New Features ###

* HUMAnN is now wheel compatible. 

## v0.9.0 08-23-2016 ##

### New Features ###

* HUMAnN is now python3 compatible. 

### Bug Fix ###

* The cigar string calculation has been updated to allow for multiple M fields along with I/D fields.

## v0.8.2 08-09-2016 ##

### Bug Fix ###

* Synced up the names of three rename files with those selected in the humann2_rename_table script.

## v0.8.1 08-05-2016 ##

### New Features ###

* A new visualization tool named "humann2_barplot" has been added.
* A new association tool named "humann2_associate" has been added.

## v0.8.0 07-26-2016 ##

### New Features ###

* A new tool named "humann2_infer_taxonomy" has been added.
* A new option, "--gap-fill <on/off>", was added. The default is set to off.
* A new option, "--taxonomy_level <Genus>", was added to humann2_split_tables. The default is "Genus". This option allows the user to select the level of taxonomy printed from input files that were created with PICRUSt predict_metagenomes.py. 
* A new minimal demo set ( input files, demo chocophlan, and demo uniref90 database) was added. This includes more unclassified hits.
* A new database download set ( utility_mapping ) has been added to the download databases script. This download includes the large rename, regroup, and infer taxonomy data files. After downloading this data set with the humann2_databases command, additional options are available for the rename and regroup scripts.
* Rename and regroup files for UniRef90 have been added. Also regroup files for pfam are now available.
* The uniref to metacyc mapping file has been updated to use centriods so that it is consistant with the latest regroup mapping files. 

## v0.7.1 04-27-2016 ##

### New Features ###

* A new option, "--search-mode <uniref50/uniref90>", was added. The default is set based on the translated database selected. This option sets the percent identity, chocophlan annotation index, and the translated search parameters.
* A new tool, humann2_gene_families_genus_level, and tutorial were added to create genus level gene families and pathways files.
* A new tool was added, humann2_split_stratified_table, to split output files into stratified and unstratified tables.

## v0.7.0 03-22-2016 ##

### New Features ###

* A new translated search filtering feature has been added. This by default filters alignments with less than 90% query coverage. This option is named "--translated-query-coverage-threshold <90.0>". The previous coverage option "--coverage-threshold <40.0>" has been renamed to "--translated-subject-coverage-threshold <50.0>" for clarity with the default increased from 40% to 50%. 
* New options have been added to the humann2_rename_table script. An option to switch between community and levelwise normalization was added with the default set to community normalization. Special features like UNMAPPED categories can now be excluded when normalizing. Also the flag "--norm" has been renamed "--units" to clarify this option only changes the units reported not the style of normalization.
* Four new translated search diamond formatted databases are now available. Two of the databases contain the full set of protein sequences (one for UniRef50 and one for UniRef90). The other databases are filtered to only include those proteins with ECs of level 4 or that are included in a MetaCyc pathway. The full databases allow users to identify unclassified proteins in their data set while the filtered databases require less memory and run time.

## v0.6.2 02-24-2016 ##

### New Features ###

* HUMAnN is now pip installable. The documentation has been updated to reflect this new feature.
* New options were added, to add names to kegg pathways and modules, to the humann2_rename_table script.
* A new tool was added to build a custom diamond database with taxonomic limitation. Adding this tool expands the Kegg legacy flow tutorial to start with quality controlled metagenome files and to include MetaPhlAn2 output.

## v0.6.1 02-10-2016 ##

### New Features ###

* The UniRef50_unknown and UniRef90_unknown values are now included in the gene families file.
    * For information on these values, please see the documentation section on the gene families file.
* More functional tests have been added
    * 12 functional end-to-end humann2 tests were added which run in about 19 minutes
    * 21 functional tool tests were added
* The tool humann2_join_tables now has the option to search sub-directories of the input folder for files to join.

### Bug Fixes ###

* If "--remove-stratified-output" is selected, the stratified output for the UNINTEGRATED pathways values are no longer printed
* The name of the UniRef50 to RXN mapping file, for script humann2_regroup_table, was updated
    * This resolves the error seen when running humann2_regroup_table with option "--groups uniref50_rxn"
* The step to build the fasta file for the custom ChocoPhlAn database has been split into multiple steps
    * This resolves the max arguments error seen when running on a MacOS with "--bypass-prescreen"

## v0.6.0 01-15-2016 ##

### New Features ###

* To reduce the influence of spurious hits in the translated search step, raw translated search results are now subjected to an initial protein coverage filter (default >50%). A read's weight is only divided over proteins that meet this coverage threshold, which (in addition to removing spurious hits) has the added benefit of substantially reducing the size of the genefamilies.tsv output file. To change the default coverage threshold, use the option "--coverage-threshold <50.0>".
* Unmapped reads are now included in the gene families abundance output file as a new "UNMAPPED" feature. This value represents the total unmapped reads after both searches, nucleotide and translated. For more information on this computation, see the documentation section about the gene families abundance file. 
* Unmapped read abundance is similarly carried through to the pathway abundance file, along with a new stratified feature, "UNINTEGRATED," which reflects the total abundance of genes that did not contribute to a metabolic pathway. These features are included in the pathway coverage file to maintain pathway abundance/coverage concordance (see below). For information on how these values are calculated, see the documentation section for each file. 

### Other Changes ###

* The included mapping from level-4 EC groups to UniRef50 clusters has been expanded based on annotations from TrEMBL. ~10x more UniRef50s can now be mapped to a level-4 EC group.
* The regroup_table script has been modified to include an "UNGROUPED" group, which captures the abundance of features that failed to map to another, non-trivial group. The "UNMAPPED" feature (see above) will always carry through to a regrouped table.
* The format of the pathway abundance and coverage output files has been updated so both files include the same order of pathways and species. The pathways are ordered by decreasing abundance. Pathways with zero abundance are not included.

## v0.5.0 10-22-2015 ##

### New Features ###

* Added default read length alignment normalization
    * This replaces the optional read length normalization using "--average-read-length L"
* Added new script humann2_merge_abundance_tables
    * This script merges the tables to stratify gene families with respect to pathways

### Other Changes ###

* Updated gap fill method to fill single lowest score
* Join script updated to allow for input files containing multiple samples
* Required Diamond version updated to v0.7.9+

## v0.4.0 09-16-2015 ##

### New Features ###

* Improvements to MetaCyc pathway quantification
    * Refined gene families to MetaCyc reactions mapping (update ignores ECs with less than 4 levels of classification)

### Cross Platform Compatibility Updates ###

* Modification to glpk/glpsol update version (resolves glpsol write permission errors and glpsol links to build directory)

## v0.3.1 08-24-2015 ##

### Cross Platform Compatibility Updates ###

* Added new option "--build-diamond" to the install (the default for Mac OS is to build diamond)

### Bug Fixes ###

* Updated the fastq to fasta function to allow for quality score lines that are similar to sequence lines (ie all alphabetical characters)

## v0.3.0 08-20-2015 ##

### New Features ###

* Improvements to structured pathway quantification identify more pathways
    * Updated databases to MetaCyc v19.1
    * Added UniProt TrEMBL data to reactions mappings
    * Added gap filling
* Added automated install of dependencies (bowtie2 and diamond)
* Added alert for user when running with demo databases and non-demo input file
* Added new option "--names {ko,metacyc-pwy,uniref50}" to humann2_rename_table
* Added new option "--groups {uniref50_ko,uniref50_go,uniref50_ec,uniref50_rxn}" to humann2_regroup_table
* Added new option "--output-max-decimals" to set the number of decimals written to the output files

### Interface Changes ###

* The default mapping file format for humann2_regroup_table has been reversed (it is now groups to features)
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

