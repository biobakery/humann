# ***ATTENTION***

Before opening a new issue here, please check the appropriate help channel on the bioBakery Support Forum (https://forum.biobakery.org) and consider opening or commenting on a thread there.

----

# HUMAnN 3.0 User Manual 

HUMAnN 3.0 is the next generation of HUMAnN (HMP Unified Metabolic Analysis Network).

----


**If you use HUMAnN 3.0 in your work, please cite the HUMAnN 3.0 paper:**

Francesco Beghini<sup>1</sup> ,Lauren J McIver<sup>2</sup> ,Aitor Blanco-Mìguez<sup>1</sup> ,Leonard Dubois<sup>1</sup> ,Francesco Asnicar<sup>1</sup> ,Sagun Maharjan<sup>2,3</sup> ,Ana Mailyan<sup>2,3</sup> ,Andrew Maltez Thomas<sup>1</sup> ,Paolo Manghi<sup>1</sup> ,Mireia Valles-Colomer<sup>1</sup> ,George Weingart<sup>2,3</sup> ,Yancong Zhang<sup>2,3</sup> ,Moreno Zolfo<sup>1</sup> ,Curtis Huttenhower<sup>2,3</sup> ,Eric A Franzosa<sup>2,3</sup> ,Nicola Segata<sup>1,4</sup>

[Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3](https://doi.org/10.7554/eLife.65088)

[eLife 2021;10:e65088](https://doi.org/10.7554/eLife.65088)
```
1 Department CIBIO, University of Trento, Italy
2 Harvard T. H. Chan School of Public Health, Boston, MA, USA
3 The Broad Institute of MIT and Harvard, Cambridge, MA, USA
4 IEO, European Institute of Oncology IRCCS, Milan, Italy
```
**And feel free to link to HUMAnN 3.0 in your Methods:**

[http://huttenhower.sph.harvard.edu/humann](http://huttenhower.sph.harvard.edu/humann)

----


**For additional information, read the [HUMAnN Tutorial](https://github.com/biobakery/biobakery/wiki/humann)**

HUMAnN is a pipeline for efficiently and accurately profiling the presence/absence and abundance of microbial pathways in a community from metagenomic or metatranscriptomic sequencing data (typically millions of short DNA/RNA reads). This process, referred to as functional profiling, aims to describe the metabolic potential of a microbial community and its members. More generally, functional profiling answers the question "What are the microbes in my community-of-interest doing (or capable of doing)?"

----

## Contents ##

* [Features](#features)
* [Workflows](#workflows)
    * [Main workflow](#main-workflow)
    * [Workflow by input file type](#workflow-by-input-file-type)
    * [Workflow by bypass mode](#workflow-by-bypass-mode)
    * [Workflow of the resume option](#workflow-of-the-resume-option)
* [Requirements](#requirements)
    * [Software](#software)
    * [Other](#other)
* [Initial Installation](#initial-installation)
    1. [Download HUMAnN](#1-download-humann)
    2. [Install HUMAnN](#2-install-humann)
    3. [Test the install](#3-test-the-install)
    4. [Try out a demo run](#4-try-out-a-demo-run)
    5. [Download the databases](#5-download-the-databases)
        * [Download the ChocoPhlAn database](#download-the-chocophlan-database)
        * [Download a translated search database](#download-a-translated-search-database)
* [Installation Update](#installation-update)
* [How to run](#how-to-run)
    * [Basic usage](#basic-usage)
    * [Demo runs](#demo-runs)
    * [Standard workflow](#standard-workflow)
* [Output files](#output-files)
    1. [Gene families file](#1-gene-families-file)
    2. [Pathway abundance file](#2-pathway-abundance-file)
    3. [Pathway coverage file](#3-pathway-coverage-file)
    4. [Intermediate temp output files](#4-intermediate-temp-output-files)
        1. [Bowtie2 alignment results](#1-bowtie2-alignment-results)
        2. [Bowtie2 reduced alignment results](#2-bowtie2-reduced-alignment-results)
        3. [Bowtie2 index files](#3-bowtie2-index-files)
        4. [Unaligned reads after Bowtie2](#4-unaligned-reads-after-bowtie2)
        5. [Custom ChocoPhlAn database](#5-custom-chocophlan-database)
        6. [MetaPhlAn Bowtie2 output](#6-metaphlan-bowtie2-output)
        7. [MetaPhlAn bugs list](#7-metaphlan-bugs-list)
        8. [Translated alignment results](#8-translated-alignment-results)
        9. [Translated alignment unaligned reads](#9-translated-alignment-unaligned-reads)
        10. [Log](#10-log)    
* [Databases](#databases)
* [Configuration](#configuration)
* [Guides to HUMAnN 3.0 utility scripts](#guides-to-humann-utility-scripts)
	* [humann_barplot](#humann_barplot)
    * [humann_config](#humann_config)
	* [humann_databases](#humann_databases)
	* [humann_gene_families_genus_level](#humann_gene_families_genus_level)
	* [humann_infer_taxonomy](#humann_infer_taxonomy)
    * [humann_join_tables](#humann_join_tables)
    * [humann_reduce_table](#humann_reduce_table)
    * [humann_regroup_table](#humann_regroup_table)
    * [humann_rename_table](#humann_rename_table)
    * [humann_renorm_table](#humann_renorm_table)
	* [humann_split_stratified_table](#humann_split_stratified_table)
	* [humann_split_table](#humann_split_table)
    * [humann_test](#humann_test)
    * [humann_unpack_pathways](#humann_unpack_pathways)
* [Other HUMAnN 3.0 guides](#other-humann-guides)
    * [Selecting a level of gene family resolution](#selecting-a-level-of-gene-family-resolution)
    * [Selecting a scope for translated search](#selecting-a-scope-for-translated-search)
    * [Paired-end reads](#humann-and-paired-end-sequencing-data)
    * [PICRUSt output](#picrust-output)
    * [Joint taxonomic profile](#joint-taxonomic-profile)
    * [Custom taxonomic profile](#custom-taxonomic-profile)
    * [Custom nucleotide reference database](#custom-nucleotide-reference-database)
    * [Custom protein reference database](#custom-protein-reference-database)
    * [Custom reference database annotations](#custom-reference-database-annotations)
    * [Custom pathways database](#custom-pathways-database)
    * [Core diversity analysis with QIIME](#core-diversity-analysis-with-qiime)
	* [Metatranscriptome analysis with HUMAnN 3.0](#metatranscriptome-analysis-with-humann)
    * [Genus level gene families and pathways](#genus-level-gene-families-and-pathways)
* [FAQs](#faqs)
* [Complete option list](#complete-option-list)

----

## Features ##

1. Community functional profiles stratified by known and unclassified organisms

    * [MetaPhlAn](http://huttenhower.sph.harvard.edu/metaphlan) and ChocoPhlAn pangenome database are used to facilitate fast, accurate, and organism-specific functional profiling
    * Organisms included are Archaea, Bacteria, Eukaryotes, and Viruses

2. Considerably expanded databases of genomes, genes, and pathways

    * [UniRef](http://www.uniprot.org/help/uniref) database provides gene family definitions
    * [MetaCyc](http://metacyc.org/) provides pathway definitions by gene family
    * [MinPath](http://omics.informatics.indiana.edu/MinPath/) is run to identify the set of minimum pathways

3. A simple user interface (single command driven flow)

    * The user only needs to provide a quality-controlled metagenome or metatranscriptome

4. Accelerated mapping of reads to reference databases (including run-time generated databases tailored to the input)

    * [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is run for accelerated nucleotide-level searches
    * [Diamond](http://ab.inf.uni-tuebingen.de/software/diamond/) is run for accelerated translated searches

----
	
## Workflows ##

### Main workflow ###

![](http://huttenhower.sph.harvard.edu/sites/default/files/humann_diamond_500x500.jpg)

### Workflow by input file type ###

There are four different types of files that can be provided as input to HUMAnN 3.0 . By default HUMAnN 3.0 will determine the type of the file. As shown in the figure below, the type of input file will determine where HUMAnN 3.0 will start the workflow. Files of type 2, 3, and 4 will begin the workflow after the alignment steps.

![](http://huttenhower.sph.harvard.edu/sites/default/files/humann_flow_by_file_type_reduced.png)

File Types:

*   File Type 1 (a quality-controlled metagenome or metatranscriptome)
    *   fastq (fastq.gz)
    *   fasta (fasta.gz)
*   File Type 2 (alignment results type 1)
    *   sam
    *   bam
*   File Type 3 (alignment results type 2)
    *   blast-like tsv
*   File Type 4 (gene table)
    *   tsv
    *   biom

----	
	
### Workflow by bypass mode ###

There are multiple bypass options that will allow you to adjust the standard workflow.

![](http://huttenhower.sph.harvard.edu/sites/default/files/humann_flow_bypass_modes_updated.png)

Bypass options:

*   --bypass-translated-search 
    *   runs all of the alignment steps except the translated search
*   --bypass-nucleotide-search 
    *   bypasses all of the alignment steps before the translated search
*   --bypass-prescreen 
    *   bypasses the taxomonic profiling step and uses the full ChocoPhlAn database
*   --taxonomic-profile bugs_list.tsv 
    *   bypasses the taxomonic profiling step and creates a custom ChocoPhlAn database of the species included in the list provided
*   --bypass-nucleotide-index
    *   starts the workflow with the nucleotide alignment step using the indexed database from "--nucleotide-database $DIR"

----	
	
### Workflow of the resume option ###

HUMAnN 3.0 includes a "--resume" option which will allow you to bypass alignment steps which have already been completed. For example, if you originally ran with a bypass option you can run just the step you bypassed with "--resume". This will only run the alignment step you bypassed and then recompute the gene families and pathways.

![](http://huttenhower.sph.harvard.edu/sites/default/files/humann_flow_resume_option_no_text.png)

When using the "--resume" option, the following steps will be bypassed if they have already been completed:

1.  Taxomonic profiling step
2.  Nucleotide alignment step
3.  Custom ChocoPhlAn database creation (merge and index)
4.  Translated alignment step

----

## Requirements ##

### Software ###

1. [MetaPhlAn 3.0](https://bitbucket.org/biobakery/metaphlan)
2. [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) (version >= 2.2) (automatically installed)
3. [Diamond](http://ab.inf.uni-tuebingen.de/software/diamond/) (version >= 0.9.24) (automatically installed)
4. [Python](http://www.python.org/) (version >= 3.7)
5. [MinPath](http://omics.informatics.indiana.edu/MinPath/) (automatically installed)
6. [Xipe](https://edwards.sdsu.edu/cgi-bin/xipe.cgi) (optional / included)
7. [Rapsearch2](http://omics.informatics.indiana.edu/mg/RAPSearch2/) (version >= 2.21) (only required if using rapsearch2 for translated search)
8. [Usearch](http://www.drive5.com/usearch/) (version >= 7.0) (only required if using usearch for translated search)
9. [SAMtools](http://samtools.sourceforge.net/) (only required if bam input files are provided)
10. [Biom-format](http://biom-format.org/) (only required if input or output files are in biom format)

Please install the required software in a location in your `$PATH` or provide the location with an optional argument to HUMAnN 3.0. 
For example, the location of the Bowtie2 install ($BOWTIE2_DIR) can be provided with "--bowtie2 $BOWTIE2_DIR".

Some of the required dependencies are installed automatically when installing HUMAnN 3.0 with pip. If these dependencies do not appear to be installed after installing HUMAnN 3.0 with pip, it might be that your environment is setup to use wheels instead of installing from source. HUMAnN 3.0 must be installed from source for it to also be able to install dependencies. To force pip to install HUMAnN 3.0 from source add one of the following options to your install command, "--no-use-wheel" or "--no-binary :all:".

If you always run with input files of type 2, 3, and 4 (for information on input file types, see section [Workflow by input file type](#markdown-header-workflow-by-input-file-type)),
MetaPhlAn2, Bowtie2, and Diamond are not required. Also if you always run with one or more bypass options (for information on bypass options, see section [Workflow by bypass mode](#markdown-header-workflow-by-bypass-mode)), 
the software required for the steps you bypass does not need to be installed.

### Other ###

1. Memory (>= 16 Gb)
2. Disk space (>= 15 Gb [to accommodate comprehensive sequence databases])
3. Operating system (Linux or Mac)

If always running with files of type 2, 3, and 4 (for information on file types, see section [Workflow by input file type](#markdown-header-workflow-by-input-file-type)),
less disk space is required. 

----

## Initial Installation ##

### 1. Download HUMAnN 3.0 ###
You can download the latest HUMAnN 3.0 release or the development version. The source contains example files. If installing with pip, it is optional to first download the HUMAnN 3.0 source.

Option 1: Latest Release (Recommended)

* [humann.tar.gz](https://pypi.python.org/pypi/humann) and unpack the latest release of HUMAnN 3.0

Option 2: Development Version

* Create a clone of the repository: 
    
	``$ git clone https://github.com/biobakery/humann ``

	Note: Creating a clone of the repository requires [Github](https://github.com/) to be installed. 

----	
	
### 2. Install HUMAnN 3.0 ###

#### Installing with pip ####

1. Install HUMAnN 3.0
    * `` $ pip install humann ``
    * This command will automatically install MinPath (and a new version of glpk) along with Bowtie2 and Diamond (if they are not already installed).
    * To bypass the install of Bowtie2 and Diamond, add the option "--install-option='--bypass-dependencies-install'" to the install command.
    * To build Diamond from source during the install, add the option "--install-option='--build-diamond'" to the install command.
    * To overwite existing installs of Bowtie2 and Diamond, add the option "--install-option='--replace-dependencies-install'" to the install command.
    * If you do not have write permissions to '/usr/lib/', then add the option "--user" to the HUMAnN install command. This will install the python package into subdirectories of '\~/.local' on Linux. Please note when using the "--user" install option on some platforms, you might need to add '~\/.local/bin/' to your $PATH as it might not be included by default. You will know if it needs to be added if you see the following message `humann: command not found` when trying to run HUMAnN after installing with the "--user" option.

----	
	
#### Installing from source ####

1. Move to the HUMAnN 3.0 directory

    * ``$ cd $HUMAnN_PATH `` 

2. Install HUMAnN 3.0

    * ``$ python setup.py install ``
    * This command will automatically install MinPath (and new version of glpk) along with Bowtie2 and Diamond (if they are not already installed).
    * To bypass the install of Bowtie2 and Diamond, add the option "--bypass-dependencies-install" to the install command. 
    * If you do not have write permissions to '/usr/lib/', then add the option ``--user`` to the install command. This will install the python package into subdirectories of '\~/.local'. Please note when using the "--user" install option on some platforms, you might need to add '\~/.local/bin/' to your $PATH as it might not be included by default. You will know if it needs to be added if you see the following message `humann: command not found` when trying to run HUMAnN after installing with the "--user" option.

----	
	
### 3. Test the install ###

1. Test out the install with unit and functional tests
     * `` $ humann_test``
     * To also run tool tests, add the option "--run-functional-tests-tools".
     * To also run end-to-end tests, add the option "-run-functional-tests-end-to-end". Please note these tests take about 20 minutes to run. Also they require all dependencies of HUMAnN 3.0 be installed in your PATH.

---- 
	 
### 4. Try out a demo run ###

With HUMAnN 3.0 installed you can try out a demo run using reduced versions of the databases. If installing from pip, please also download the source as it contains the demo examples.

``$ humann --input examples/demo.fastq --output $OUTPUT_DIR ``

Output from this demo run will be written to the folder $OUTPUT_DIR.

Please continue with the install directions to download the full databases before running with your sequencing data.

----

### 5. Download the databases ###

Downloading the databases is a required step if your input is a filtered shotgun sequencing metagenome file (fastq, fastq.gz, fasta, or fasta.gz format). If your input files will always be mapping results files (sam, bam or blastm8 format) or gene tables (tsv or biom format), you do not need to download the ChocoPhlAn and translated search databases. 

----

#### Download the ChocoPhlAn database ####

Download the ChocoPhlAn database providing $INSTALL_LOCATION as the location to install the database (approximate size = 16.4 GB).

`` $ humann_databases --download chocophlan full $INSTALL_LOCATION ``

NOTE: The humann config file will be updated to point to this location for the default chocophlan database. If you move this database, please use the "humann_config" command to update the default location of this database. Alternatively you can always provide the location of the chocophlan database you would like to use with the "--nucleotide-database <chocophlan>" option to humann.

----

#### Download a translated search database ####

While there is only one ChocoPhlAn database, users have a choice of translated search databases which offer trade-offs between resolution, coverage, and performance. For help selecting a translated search database, see the following guides:

* [Selecting a level of gene family resolution](#markdown-header-selecting-a-level-of-gene-family-resolution)
* [Selecting a scope for translated search](#markdown-header-selecting-a-scope-for-translated-search)

Download a translated search database providing $INSTALL_LOCATION as the location to install the database:

* To download the full UniRef90 database (20.7GB, **recommended**):
    * ``$ humann_databases --download uniref uniref90_diamond $INSTALL_LOCATION``

* To download the EC-filtered UniRef90 database (0.9GB):
    * ``$ humann_databases --download uniref uniref90_ec_filtered_diamond $INSTALL_LOCATION``

* To download the full UniRef50 database (6.9GB):
    * ``$ humann_databases --download uniref uniref50_diamond $INSTALL_LOCATION``

* To download the EC-filtered UniRef50 database (0.3GB):
    * ``$ humann_databases --download uniref uniref50_ec_filtered_diamond $INSTALL_LOCATION``
   
**NOTE:** The humann config file will be updated to point to this location for the default uniref database. If you move this database, please use the `humann_config` command to update the default location of this database. Alternatively you can always provide the location of the uniref database you would like to use with the `--protein-database <uniref>` option to humann.

**Note:** Please make sure that all the databases in the directory are of same version. Different versions of database must be separated in different directory. 

**NOTE:** By default HUMAnN 3.0 runs DIAMOND for translated alignment. If you would like to use RAPSearch2 for translated alignment, first download the RAPSearch2 formatted database by running this command with the rapsearch2 formatted database selected. It is suggested that you install both databases in the same folder so this folder can be the default uniref database location. This will allow you to switch between alignment software without having to specify a different location for the database.

----

## Installation Update ##

If you have already installed HUMAnN 3.0, using the [Initial Installation](#markdown-header-initial-installation) steps, and would like to upgrade your installed version to the latest version, please follow these steps.

1. [Download HUMAnN 3.0](#markdown-header-1-download-humann)
2. [Install HUMAnN 3.0](#markdown-header-2-install-humann)

Since you have already downloaded the databases in the initial installation, you do not need to download the databases again unless there are new versions available. However, you will want to update your latest HUMAnN 3.0 install to point to the databases you have downloaded as by default the new install configuration will point to the demo databases.

To update your HUMAnN 3.0 configuration file to include the locations of your downloaded databases, please use the following steps.

1. Update the location of the ChocoPhlAn database (replacing $DIR with the full path to the directory containing the ChocoPhlAn database)

    ``$ humann_config --update database_folders nucleotide $DIR ``

2. Update the location of the UniRef database (replacing $DIR with the full path to the directory containing the UniRef database)

    ``$ humann_config --update database_folders protein $DIR ``

Please note, after a new installation, all of the settings in the configuration file, like the database folders, will be reset to the defaults. If you have any additional settings that differ from the defaults, please update them at this time.

For more information on the HUMAnN 3.0 configuration file, please see the [Configuration](#markdown-header-configuration) section.

----

## How to run ##

### Basic usage ###

```
$ humann --input $SAMPLE --output $OUTPUT_DIR
```

`$SAMPLE` = a single file that is one of the following types:

1.  filtered shotgun sequencing metagenome file (fastq, fastq.gz, fasta, or fasta.gz format)
2.  alignment file (sam, bam or blastm8 format)
3.  gene table file (tsv or biom format)

`$OUTPUT_DIR` = the output directory

**Three output files will be created:**

1. `$OUTPUT_DIR/$SAMPLENAME_genefamilies.tsv`
2. `$OUTPUT_DIR/$SAMPLENAME_pathcoverage.tsv`
3. `$OUTPUT_DIR/$SAMPLENAME_pathabundance.tsv`

where `$SAMPLENAME` is the basename of `$SAMPLE`

*The gene families file will not be created if the input file type is a gene table.

**Intermediate temp files will also be created:**

1. `$DIR/$SAMPLENAME_bowtie2_aligned.sam`
	* the full alignment output from bowtie2 
2. `$DIR/$SAMPLENAME_bowtie2_aligned.tsv`
	* only the reduced aligned data from the bowtie2 output
3. `$DIR/$SAMPLENAME_bowtie2_index`
	* bowtie2 index files created from the custom chochophlan database
4. `$DIR/$SAMPLENAME_bowtie2_unaligned.fa`
	* a fasta file of unaligned reads after the bowtie2 step
5. `$DIR/$SAMPLENAME_custom_chocophlan_database.ffn` 
	* a custom chocophlan database of fasta sequences
6. `$DIR/$SAMPLENAME_metaphlan_bowtie2.txt`
	* the bowtie2 output from metaphlan
7. `$DIR/$SAMPLENAME_metaphlan_bugs_list.tsv` 
	* the bugs list output from metaphlan
8. `$DIR/$SAMPLENAME_$TRANSLATEDALIGN_aligned.tsv` 
	* the alignment results from the translated alignment step
9. `$DIR/$SAMPLENAME_$TRANSLATEDALIGN_unaligned.fa`
	* a fasta file of unaligned reads after the translated alignment step
10. `$DIR/$SAMPLENAME.log`
	* a log of the run
	
* `$DIR=$OUTPUT_DIR/$SAMPLENAME_humann_temp/`
* `$SAMPLENAME` is the basename of the fastq/fasta input file
* `$TRANSLATEDALIGN` is the translated alignment software selected (default: DIAMOND)

**NOTE:** `$SAMPLENAME` can be set by the user with the option `--output-basename <$NEWNAME>`. 

----

### Demo runs ###

The examples folder contains four demo example input files. These files are of fasta, fastq, sam, and blastm8 format. Blastm8 format is created by the following software: rapsearch2, usearch, and blast.

To run the fasta demo:

`` $ humann --input examples/demo.fasta --output $OUTPUT_DIR``

To run the fastq demo:

`` $ humann --input examples/demo.fastq --output $OUTPUT_DIR``

To run the sam demo:

`` $ humann --input examples/demo.sam --output $OUTPUT_DIR``

To run the blastm8 demo:

`` $ humann --input examples/demo.m8 --output $OUTPUT_DIR``

`$OUTPUT_DIR` is the output directory

Since sam and blastm8 are mapping results, using these files as input to HUMAnN 3.0 will bypass both the nucleotide and translated mapping portions of the flow.

----

### Standard workflow ###

The standard workflow involves running HUMAnN 3.0 on each filtered shotgun sequencing metagenome file, normalizing, and then merging the output files. 

Prior to running the workflow, filter your shotgun sequencing metagenome file. We recommend using [KneadData](http://huttenhower.sph.harvard.edu/kneaddata) for metagenome quality control including removal of host reads.

To run the standard workflow, follow these steps:

1. Run HUMAnN 3.0 on each of your filtered fastq files in `$INPUT_DIR` placing the results in `$OUTPUT_DIR`
    * for `$SAMPLE.fastq` in `$INPUT_DIR`
        * `` $ humann --input $SAMPLE.fastq --output $OUTPUT_DIR``
            * Replace `$SAMPLE.fastq` with the name of the fastq input file
            * Replace `$OUTPUT_DIR` with the full path to the folder to write output
        * If you have paired-end reads, please see the guide [Paired-end reads](#markdown-header-humann-and-paired-end-sequencing-data)
    * The results will be three main output files for each input file named `$SAMPLE_genefamilies.tsv`, `$SAMPLE_pathabundance.tsv`, and `$SAMPLE_pathcoverage.tsv`. 

2. Normalize the abundance output files
    * We recommend normalizing the abundance data based on the statistical tests you will perform.
    * Prior to nomalization, select the scheme to use (copies per million or relative abundance). For example, if using [MaAsLin](http://huttenhower.sph.harvard.edu/maaslin), select relative abundance.
    * Use the HUMAnN 3.0 tool renorm table, to compute the normalized abundances (relative abundance is selected in the example command below)
        * for `$SAMPLE_genefamilies.tsv` in `$OUTPUT_DIR`
            * `` $ humann_renorm_table --input $SAMPLE_genefamilies.tsv --output $SAMPLE_genefamilies_relab.tsv --units relab ``
    * Please note, gene family abundance is reported in RPK (reads per kilobase). This is computed as the sum of the scores for all alignments for a gene family. An alignment score is based on the number of matches to the reference gene for a specific sequence. It is divided by the length of the reference gene in kilobases to normalize for gene length. Each alignment score is also normalized to account for alignments for a single sequence to multiple reference genes. Alignments are not considered if they do not pass the e-value, identity, and coverage thresholds.
    * If you would like to normalize using the number of reads aligned per input file, this count along with the total number of reads and the percent unaligned reads after each alignment step is included in the log file. For more information on what is included in the log file, see the Intermediate temp output file section [Log](#markdown-header-10-log).
    * Alternatively, gene families can be regrouped to different functional categories prior to normalization. See the guide to [humann_regroup_table](#markdown-header-humann_regroup_table) for detailed information. 
    
3. Join the output files (gene families, coverage, and abundance) from the HUMAnN 3.0 runs from all samples into three files
    * `` $ humann_join_tables --input $OUTPUT_DIR --output humann_genefamilies.tsv --file_name genefamilies_relab ``
    * `` $ humann_join_tables --input $OUTPUT_DIR --output humann_pathcoverage.tsv --file_name pathcoverage ``
    * `` $ humann_join_tables --input $OUTPUT_DIR --output humann_pathabundance.tsv --file_name pathabundance_relab ``
    * For each command, replace `$OUTPUT_DIR` with the full path to the folder containing the HUMAnN 3.0 output files.
    * The resulting files from these commands are named `humann_genefamilies.tsv`, `humann_pathabundance.tsv`, and `humann_pathcoverage.tsv`.

----	
	
## Output files ##

When HUMAnN 3.0 is run, three main output files will be created (where `` $SAMPLENAME = the basename of $SAMPLE ``):

### 1. Gene families file ###

``` 
# Gene Family	$SAMPLENAME_Abundance-RPKs
UNMAPPED        187.0
UniRef50_unknown        150.0
UniRef50_unknown|g__Bacteroides.s__Bacteroides_fragilis 150.0
UniRef50_A6L0N6: Conserved protein found in conjugate transposon	67.0
UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_fragilis	57.0
UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_finegoldii	5.0
UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_stercoris	4.0
UniRef50_A6L0N6: Conserved protein found in conjugate transposon|unclassified	1.0
UniRef50_O83668: Fructose-bisphosphate aldolase	60.0
UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_vulgatus	31.0
UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_thetaiotaomicron	22.0
UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_stercoris	7.0
```
*   File name: `` $OUTPUT_DIR/$SAMPLENAME_genefamilies.tsv ``
*   This file details the abundance of each gene family in the community. Gene families are groups of evolutionarily-related protein-coding sequences that often perform similar functions.
*   Gene family abundance at the community level is stratified to show the contributions from known and unknown species. Individual species' abundance contributions sum to the community total abundance.
*   HUMAnN 3.0 uses the MetaPhlAn2 software along with the ChocoPhlAn database and translated search database for this computation.
*   Gene family abundance is reported in RPK (reads per kilobase) units to normalize for gene length; RPK units reflect relative gene (or transcript) copy number in the community. RPK values can be further sum-normalized to adjust for differences in sequencing depth across samples.
*   Please note the gene families file will not be created if the input file type is a gene table.
*   The "UNMAPPED" value is the total number of reads which remain unmapped after both alignment steps (nucleotide and translated search). Since other gene features in the table are quantified in RPK units, "UNMAPPED" can be interpreted as a single unknown gene of length 1 kilobase recruiting all reads that failed to map to known sequences.
* The UniRef50_unknown values represent the total abundance of reads which map to ChocoPhlAn nucleotide sequences which do not have a UniRef50 annotation.

----

### 2. Pathway abundance file ###

```
# Pathway	$SAMPLENAME_Abundance
UNMAPPED	140.0
UNINTEGRATED	87.0
UNINTEGRATED|g__Bacteroides.s__Bacteroides_caccae	23.0
UNINTEGRATED|g__Bacteroides.s__Bacteroides_finegoldii	20.0
UNINTEGRATED|unclassified	12.0
PWY0-1301: melibiose degradation	57.5
PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_caccae	32.5
PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_finegoldii	4.5
PWY0-1301: melibiose degradation|unclassified	3.0
PWY-5484: glycolysis II (from fructose-6P)	54.7
PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_caccae	16.7
PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_finegoldii	8.0
PWY-5484: glycolysis II (from fructose-6P)|unclassified	6.0
```
         
*   File name: `` $OUTPUT_DIR/$SAMPLENAME_pathabundance.tsv ``
*   This file details the abundance of each pathway in the community as a function of the abundances of the pathway's component reactions, with each reaction's abundance computed as the sum over abundances of genes catalyzing the reaction. 
*   Pathway abundance is computed once at the community level and again for each species (plus the "unclassified" stratum) using community- and species-level gene abundances along with the structure of the pathway.
*   The pathways are ordered by decreasing abundance with pathways for each species also sorted by decreasing abundance. Pathways with zero abundance are not included in the file.
*   Pathway abundance is proportional to the number of complete "copies" of the pathway in the community. Thus, for a simple linear pathway RXN1â†’RXN2â†’RXN3â†’RXN4, if RXN1 is 10 times as abundant as RXNs 2-4, the pathway abundance will be driven by the abundances of RXNs 2-4.
*   Unlike gene abundance, **a pathway's community-level abundance is not necessarily the sum of its stratified abundance values**. For example, continuing with the simple linear pathway example introduced above, if the abundances of RXNs 1-4 are [5, 5, 10, 10] in Species_A and [10, 10, 5, 5] in Species_B, HUMAnN 3.0 would report that Species_A and Species_B each contribute 5 complete copies of the pathway. However, at the community level, the reaction totals are [15, 15, 15, 15], and thus HUMAnN 3.0 would report 15 complete copies.
*   In greater detail, the abundance for each pathway is a recursive computation of abundances of sub-pathways with paths resolved to abundances based on the relationships and abundances of the reactions contained in each. Each path, the smallest portion of a pathway or sub-pathway which can't be broken down into sub-pathways, has an abundance that is the max or harmonic mean of the reaction abundances depending on the relationships of these reactions. Optional reactions are only added to the overall abundance if their abundance is greater than the harmonic mean of the required reactions.
*   Gap filling allows for a single required reaction to have a zero abundance. For all pathways, the required reaction with the lowest abundance is replaced with the abundance of the required reaction with the second lowest abundance.
*   By default, HUMAnN 3.0 uses MetaCyc pathway definitions and MinPath to identify a parsimonious set of pathways which explain observed reactions in the community.
*   The user has the option to provide a custom pathways database to HUMAnN 3.0 and to use all pathways instead of the minimal pathways computed by MinPath.
*   To account for non-linearity in the conversion of gene copy number to pathway copy number, we define a "compression constant" (*k*) equal to the total pathway abundance divided by the total abundance of genes that contributed to pathways. The "UNMAPPED" value reported in the pathway abundance table is equal to the total number of unmapped reads scaled by *k* (making it more comparable with pathway abundance values). Similarly, we define an "UNINTEGRATED" abundance for 1) the community, 2) each identified species, and 3) the "unclassified" stratum equal to the total abundance of genes in that level that did not contribute to pathways (scaled by *k*).
*   "UNINTEGRATED" does not appear for stratifications with no detected pathways.

----

### 3. Pathway coverage file ###

``` 
# Pathway	$SAMPLENAME_Coverage
UNMAPPED	1.0
UNINTEGRATED	1.0
UNINTEGRATED|g__Bacteroides.s__Bacteroides_caccae	1.0
UNINTEGRATED|g__Bacteroides.s__Bacteroides_finegoldii	1.0
UNINTEGRATED|unclassified	1.0
PWY0-1301: melibiose degradation	1.0
PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_caccae	1.0
PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_finegoldii	1.0
PWY0-1301: melibiose degradation|unclassified	1.0
PWY-5484: glycolysis II (from fructose-6P)	1.0
PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_caccae	0.7
PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_finegoldii	0.7
PWY-5484: glycolysis II (from fructose-6P)|unclassified	0.3
```

*   File name: `` $OUTPUT_DIR/$SAMPLENAME_pathcoverage.tsv ``
*   Pathway coverage provides an alternative description of the presence (1) and absence (0) of pathways in a community, independent of their quantitative abundance.
*   More specifically, HUMAnN 3.0 assigns a confidence score to each reaction detected in the community. Reactions with abundance greater than the median reaction abundance are considered to be more confidently detected than those below the median abundance.
*   HUMAnN 3.0 then computes pathway coverage using the same algorithms described above in the context of pathway abundance, but substituting reaction confidence for reaction abundance.
*   A pathway with coverage = 1 is considered to be confidently detected (independent of its abundance), as this implies that all of its member reactions were also confidently detected. A pathway with coverage = 0 is considered to less confidently detected (independent of its abundance), as this implies that some of its member reactions were not confidently detected.
*   Like pathway abundance, pathway coverage is computed for the community as a whole, as well as for each detected species and the unclassified stratum.
*   Much as community-level pathway abundance is not the strict sum of species-level contributions, **it is possible for a pathway to be confidently covered at the community level but never confidently detected from any single species**.
*   Pathway coverage is reported for any non-zero pathway abundance computed at the community-level or for an individual stratum (species or "unclassified").
*   The pathway coverage file follows the same order for pathways and species as the abundance file. Entries for "UNMAPPED" and "UNINTEGRATED" are included and set to 1.0 to further maintain this ordering, although the "coverage" of these features is not meaningful.

----

### 4. Intermediate temp output files ###

Ten intermediate temp output files will be created where:

```
$DIR = $OUTPUT_DIR/$SAMPLENAME_humann_temp/
$SAMPLENAME = basename of the fastq/fasta input file named $SAMPLE
$TRANSLATEDALIGN = translated alignment software selected (diamond, rapsearch2 or usearch)
```

**NOTE:** `$SAMPLENAME` can be set by the user with the option `--output-basename <$NEWNAME>`

#### 1. Bowtie2 alignment results ####

```
@HD	VN:1.0	SO:unsorted
@SQ	SN:g__Ruminococcus.s__Ruminococcus_bromii|UniRef90_D4L6K4|UniRef50_R6U703	LN:540
r99491	0	g__Bacteroides.s__Bacteroides_stercoris|UniRef90_R6B629|UniRef50_R5RCC8	1015	42	151M	*	0	0	CTGACCGATTTATATGAGGA	zzzzzzzzzzzzzzzzzzzz
r99526	0	g__Parabacteroides.s__Parabacteroides_merdae|UniRef90_unknown|UniRef50_D9RX34	155	42	151M	*	0	0	AATTTTCTTCAAAAAATATA	zzzzzzzzzzzzzzzzzzzz
r99581	16	g__Bacteroides.s__Bacteroides_stercoris|UniRef90_unknown|UniRef50_R6SXR7	2503	42	151M	*	0	0	TCTTTTATGCAGGGGATATG	zzzzzzzzzzzzzzzzzzzz
```

*   File name: `` $DIR/$SAMPLENAME_bowtie2_aligned.sam `` 
*   This file has the full alignment output from bowtie2.
*   This file is in SAM format. See the [SAM Format Specification](https://samtools.github.io/hts-specs/SAMv1.pdf) for more information.
*   The columns in this file (not including headers which start with ``@``) are as follows:
    * Column 1: Query sequence name
    * Column 2: bitwise flag
    * Column 3: Reference sequence name
    * Column 4: 1-based leftmost mapping position
    * Column 5: Mapping quality
    * Column 6: Cigar string
    * Column 7: Reference name of the mate/next read
    * Column 8: Position of the mate/next read
    * Column 9: Observed template length
    * Column 10: Segment sequence
    * Column 11: Quality scores

----	
	
#### 2. Bowtie2 reduced alignment results ####

``` 
r93	gi|423245752|ref|NZ_JH724135.1|:381614-382081|357276|g__Bacteroides.s__Bacteroides_dorei|UniRef90_A6L5K0|UniRef50_A6L5K0|468    99.0    100.0                                                   0
r113	gi|423245752|ref|NZ_JH724135.1|:381614-382081|357276|g__Bacteroides.s__Bacteroides_dorei|UniRef90_A6L5K0|UniRef50_A6L5K0|468    99.0    100.0                                                   0	
r704	gi|423245752|ref|NZ_JH724135.1|:381614-382081|357276|g__Bacteroides.s__Bacteroides_dorei|UniRef90_A6L5K0|UniRef50_A6L5K0|468    99.0    100.0                                                   0	
r663	gi|423245752|ref|NZ_JH724135.1|:381614-382081|357276|g__Bacteroides.s__Bacteroides_dorei|UniRef90_A6L5K0|UniRef50_A6L5K0|468    99.0    100.0                                                   0	
r940	gi|423245752|ref|NZ_JH724135.1|:381614-382081|357276|g__Bacteroides.s__Bacteroides_dorei|UniRef90_A6L5K0|UniRef50_A6L5K0|468    99.0    100.0                                                   0	 
```

*   File name: `` $DIR/$SAMPLENAME_bowtie2_aligned.tsv ``
*   This file contains the minimal amount of alignment results from the Bowtie2 sam file.
*   This file is in blast m8 format with some columns always empty to reduce the file size.
*   The columns in this file are as follows:
    * Column 1: Query sequence name
    * Column 2: Reference sequence name
    * Column 3: Percent identity
    * Column 4: Alignment length
    * Column 5-7: Empty or always 0 column
    * Column 8: Alignment length - 1
    * Column 9: Subject start index
    * Column 10: Subject end index
    * Column 11-12: Empty or always 0 column

----	
	
#### 3. Bowtie2 index files ####


*   Example not included as files are binary.
*   File name: `` $DIR/$SAMPLENAME_bowtie2_index* ``
*   These are a set of files containing the Bowtie2 index created from the custom ChocoPhlAn database.

----

#### 4. Unaligned reads after Bowtie2 ####

```
>r4370
GGCGGACGATCTTGTCGCCCAGCCTGTAGCCTTTCTGGTACACCGTGATGACGGTGCCGCTCTCCTGCCCGTCCGTGGCGGGGATCTGCTGG
>r4398
TGCCCGGACAGGATCTTCTCTTTCGTACCGGGCATCATCTGCTCCATGATCTCCACGCCTCGCATGAACTTTTCAGAACGGGCAACGTAGGA
```

*   File name: `` $DIR/$SAMPLENAME_bowtie2_unaligned.fa ``
*   This is a fasta file of unaligned reads after the Bowtie2 step.
*   These are the reads that will be provided as input in the translated alignment step.

----

#### 5. Custom ChocoPhlAn database ####

```
>gi|479150083|ref|NC_021013.1|:976220-976759|40518|g__Ruminococcus.s__Ruminococcus_bromii|UniRef90_D4L6K4|UniRef50_R6U703
ATGTTCTATGTATTTCTTGCAGAAGGCTTTGAAGAAACAGAGGCGCTTGCCCCCGTTGATGTAATGCGCAGGGCAAAGCT
TGATGTTAAAACAGTCGGTGTAACAGGCGAATGTGTTACAAGCTCACACGGTGTGCCTGTAAAAGCCGATATCACAATTG
ACAATATTGACCTTGACGATGTTCAGGGTGTTGTACTCCCCGGTGGTATGCCCGGAACTCTCAATCTTGAGGCAAACAAA
AAGGTTCTTGAGGCTGTTAAGTATAGCTGTGAAAACGGCAAAATCGTTGCCGCAATCTGTGCCGCTCCGTCAATTCTCGG
```

*   File name: `` $DIR/$SAMPLENAME_custom_chocophlan_database.ffn ``
*   This file is a custom ChocoPhlAn database of fasta sequences.

----

#### 6. MetaPhlAn Bowtie2 output ####

```
r113	gi|224485636|ref|NZ_EQ973490.1|:c728571-728107
r559	gi|479185170|ref|NC_021030.1|:c1678719-1677127
r663	gi|512436175|ref|NZ_KE159463.1|:c142391-139122
r704	gi|423310881|ref|NZ_JH724270.1|:c220428-218656
r1086	gi|238922432|ref|NC_012781.1|:c1988048-1987140 
```

*   File name: `` $DIR/$SAMPLENAME_metaphlan_bowtie2.txt ``
*   This file is the Bowtie2 output from MetaPhlAn.

----

#### 7. MetaPhlAn bugs list ####

```
k__Bacteria	100.0
k__Bacteria|p__Bacteroidetes	73.86802
k__Bacteria|p__Firmicutes	26.13198
k__Bacteria|p__Bacteroidetes|c__Bacteroidia	73.86802
k__Bacteria|p__Firmicutes|c__Clostridia	15.60912
k__Bacteria|p__Firmicutes|c__Negativicutes	10.52286
k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales	73.86802
k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales	15.60912
k__Bacteria|p__Firmicutes|c__Negativicutes|o__Selenomonadales	10.52286
k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae	51.32768
```

*   File name: `` $DIR/$SAMPLENAME_metaphlan_bugs_list.tsv ``
*   This file is the bugs list output from MetaPhlAn.

----

#### 8. Translated alignment results ####

```
r2805	UniRef50_E2ZJD8|627	37.50	48	30	0	147	4	152	199	5e-06	40.0
r2805	UniRef50_E2ZHU1|963	37.50	48	30	0	147	4	152	199	7e-06	39.7
r2805	UniRef50_K1UMQ9|612	35.42	48	31	0	147	4	35	82	2e-05	38.5
r3036	UniRef50_K1TCN4|540	100.00	22	0	0	150	85	148	169	8e-10	52.8
r3036	UniRef50_UPI00046A4B12|696	35.42	48	30	1	149	6	88	134	1e-05	38.9
```

*   File name: `` $DIR/$SAMPLENAME_$TRANSLATEDALIGN_aligned.tsv ``
*   This file is the alignment results from the translated alignment step.
*   This file is formatted as a tab-delimited blast-like results file (blast m8 format).
*   The columns in this file are as follows:
    * Column 1: Query sequence name
    * Column 2: Reference sequence name
    * Column 3: Percent identity
    * Column 4: Alignment length
    * Column 5: Mismatches
    * Column 6: Gap opening
    * Column 7: Query start
    * Column 8: Query end
    * Column 9: Reference start
    * Column 10: Reference end
    * Column 11: E-value
    * Column 12: Bit score

----	
	
#### 9. Translated alignment unaligned reads ####

```
>r4370
GGCGGACGATCTTGTCGCCCAGCCTGTAGCCTTTCTGGTACACCGTGATGACGGTGCCGCTCTCCTGCCCGTCCGTGGCGGGGATCTGCTGG
>r4398
TGCCCGGACAGGATCTTCTCTTTCGTACCGGGCATCATCTGCTCCATGATCTCCACGCCTCGCATGAACTTTTCAGAACGGGCAACGTAGGA
```
    
*   File name: `` $DIR/$SAMPLENAME_$TRANSLATEDALIGN_unaligned.fa ``
*   This is a fasta file of the unaligned reads after the translated alignment step

----

#### 10. Log ####

```
03/16/2015 01:09:52 PM - humann.utilities - INFO: File ( demo.fastq ) is of format:  fastq
03/16/2015 01:09:52 PM - humann.config - INFO: Run config settings:
DATABASE SETTINGS
nucleotide database folder = data/chocophlan_DEMO
protein database folder = data/uniref_DEM
```

*   File name: `` $DIR/$SAMPLENAME.log ``
*   This file is a log of the run.
*   Timestamps for each step in the flow are benchmarked in the log. Look for these as lines containing "TIMESTAMP".
*   The percent unaligned reads after each alignment step is included in the log. Look for these as lines containing the phrase "Unaligned reads after". 
*   The total number of reads for the input file is contained in the alignment statistics output from the nucleotide alignment step recorded in the log.

----

## Databases ##

HUMAnN 3.0 uses two databases for alignment, ChocoPhlAn and UniRef. There are different formats of these databases. The demo formats of both are included in the HUMAnN 3.0 install.

To see which databases are currently available, run (with example output included below):
```sh
$ humann_databases --available
HUMANnN Databases ( database : build = location )
chocophlan : DEMO = http://huttenhower.sph.harvard.edu/humann_data/chocophlan/DEMO_chocophlan.tar.gz
chocophlan : full = http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.tar.gz
uniref : DEMO_diamond = http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref50_GO_filtered/uniref50_DEMO_diamond.tar.gz
uniref : diamond = http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref50_GO_filtered/uniref50_GO_filtered_diamond.tar.gz
uniref : rapsearch2 = http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref50_GO_filtered/uniref50_GO_filtered_rapsearch2.tar.gz
```
To download a database, run:
```sh
$ humann_databases --download $DATABASE $BUILD $INSTALL_LOCATION
```

This will automatically update the HUMAnN 3.0 configuration. If you do not have permissions to edit the config file run with the option `--update-config no` and then include the locations of the databases when running HUMAnN.

----

## Configuration ##

HUMAnN 3.0 uses a configuration file to store user configuration settings. This configuration file is automatically updated when a database is installed.

To view the current settings of the configuration file, run (with example output included below):
```sh
$ humann_config --print
HUMAnN Configuration ( Section : Name = Value )
output_format : remove_stratified_output = False
output_format : output_max_decimals = 10
alignment_settings : prescreen_threshold = 0.01
alignment_settings : evalue_threshold = 1.0
alignment_settings : identity_threshold = 50.0
database_folders : nucleotide = data/chocophlan_DEMO
database_folders : protein = data/uniref_DEMO
run_modes : bypass_nucleotide_search = False
run_modes : verbose = False
run_modes : bypass_nucleotide_index = False
run_modes : bypass_translated_search = False
run_modes : resume = False
run_modes : threads = 1
run_modes : bypass_prescreen = False
```
To change a value in the configuration file, run:
```sh
$ humann_config --update $SECTION $NAME $VALUE 
```

The settings included in the configuration file plus others selected at run-time are printed to the beginning of the log file for each run.

Example log file configuration settings section:
```sh
03/16/2015 01:09:52 PM - humann.config - INFO: 
Run config settings: 

DATABASE SETTINGS
nucleotide database folder = data/chocophlan_DEMO
protein database folder = data/uniref_DEMO
pathways database file 1 = data/pathways/metacyc_reactions.uniref
pathways database file 2 = data/pathways/metacyc_pathways

RUN MODES
resume = False
verbose = True
bypass prescreen = False
bypass nucleotide index = False
bypass nucleotide search = False
bypass translated search = False
translated search = diamond
pick frames = off
threads = 4

ALIGNMENT SETTINGS
evalue threshold = 1.0
prescreen threshold = 0.01
identity threshold = 50.0

PATHWAYS SETTINGS
minpath = on
xipe = off

INPUT AND OUTPUT FORMATS
input file format = fastq
output file format = tsv
output max decimals = 10
remove stratified output = False
log level = DEBUG
```

----

## Guides to HUMAnN 3.0 utility scripts ##

### humann_barplot ###

* **Basic usage:** `$ humann_barplot --input $TABLE.tsv --focal-feature $FEATURE --outfile $FIGURE`
* `$TABLE.tsv` = a stratified HUMAnN 3.0 output file
* `$FEATURE` = Feature from the table to plot (defaults to first feature)
* `$FIGURE` = Where to save the figure
* Run with `-h` to see additional command line options

`humann_barplot` produces plots of stratified HUMAnN 3.0 features. **Note**: unlike many other HUMAnN 3.0 utilities, `humann_barplot` requires the Python scientific stack (notably `matplotlib`) to operate. `humann_barplot` includes many options for sorting and scaling data, as detailed in the `--help` menu. `humann_barplot` is also explored in greater detail in the [HUMAnN 3.0](https://github.com/biobakery/biobakery/wiki/humann).

Here is an example of a HUMAnN 3.0 barplot for a pathway (denitrification) that was preferentially enriched in Human Microbiome Project oral samples relative to other body sites. This figure uses many options from `humann_barplot`, including regrouping by genus, pseudolog scaling, and sorting samples by similarity and metadata:

![page_DENITRIFICATION-PWY.png](https://bitbucket.org/repo/bAAy6o/images/731303924-page_DENITRIFICATION-PWY.png)

----

### humann_config ###

Used for viewing and setting HUMAnN 3.0 defaults. See the [Configuration](#markdown-header-configuration) section above for details.

----

### humann_databases ###

Used for downloading ChocoPhlAn, translated search databases, and HUMAnN 3.0 utility mapping files. See the [Databases](#markdown-header-databases) section above for details.

----

### humann_gene_families_genus_level ###

See the [Genus level gene families and pathways](#markdown-header-genus-level-gene-families-and-pathways) guide below.

----

### humann_infer_taxonomy ###

* **Basic usage:** `$ humann_infer_taxonomy --input $INPUT_GENES.tsv --level $LEVEL --resolution $RES --output $OUTPUT.tsv`
* `$INPUT_GENES.tsv` = a file containing gene families abundance data with unclassified strata (tsv format)
* `$LEVEL` = desired taxonomic level for inference (default: Family)
* `$RES` = UniRef resolution of the input file
* `$OUTPUT.tsv` = the file to write the new merged abundance table (tsv format)
* Run with `-h` to see additional command line options

Relative to HUMAnN 3.0's pangenome mapping, we tend to be less confident about the taxonomic assignments made during translated search. For this reason, all translated search results are given as "taxonony=unclassified" by default (allowing the user to focus on functional assignments within these results). However, thanks to the structure of the UniRef50 and UniRef90 databases, it is possible to infer some degree of taxonomic information for unclassified hits when performing translated search against these databases.

Each UniRef family is associated with one or more protein sequences of a given degree of sequence homology, out of which a single representative is selected. Consider the following UniRef50 family represented by sequence O67749:

```
>UniRef50_O67749 GTPase Der n=10 Tax=Aquificaceae RepID=DER_AQUAE
```

This protein family contains 10 sequences (the `n=` field). All sequences derive from species belonging to the family-level clade Aquificaceae (the `tax=` field): the species' lowest common ancestor (LCA). If we observe translated hits to the UniRef50\_O67749 family, we can infer that the associated reads derive from species in Aquificaceae. This is not a guarantee, as the protein family may be distributed outside of the Aquificaceae, however no such sequences are known to UniRef. This method would not allow us to infer a _genus_ for reads mapping to UniRef50\_O67749 by translated search, as the 10 member sequences are heterogeneous at the genus level. However, clades above the LCA (e.g. kingdom=Bacteria) _can_ be inferred.

Based on this method, the HUMAnN 3.0 utility `humann_infer_taxonomy` can be used to assign approximate taxonomic information to translated search results as summarized in HUMAnN 3.0's `genefamilies.tsv` output files. By default, assignments are attempted at the taxonomic family level, but other levels (kingdom through genus) can be selected. If an assignment cannot be made at the target level, the taxonomy is left as "unclassified." By default, assignments are made to the gene family totals. Alternative modes allow the user to carry through known taxonomic information from the pangenome search ("stratified" mode) or to focus on unclassified hits only ("unclassified" mode). Consider the following examples based on this input file:

```
# genefamilies.tsv
UniRef50_O67749                                  50
UniRef50_O67749|g__Aquifex.s__Aquifex_pyrophilus 25
UniRef50_O67749|unclassified                     25
```

The command `humann_infer_taxonomy -i genefamilies.tsv -r uniref50` yields:

```
UniRef50_O67749                                  50
UniRef50_O67749|f__Aquificaceae                  50
```

The command `humann_infer_taxonomy -i genefamilies.tsv -r uniref50 -l Genus` yields:

```
UniRef50_O67749                                  50
UniRef50_O67749|unclassified                     50
```

The command `humann_infer_taxonomy -i genefamilies.tsv -r uniref50 -l Genus -m stratified` yields:

```
UniRef50_O67749                                  50
UniRef50_O67749|g__Aquifex                       25
UniRef50_O67749|unclassified                     25
```

In the last example, the LCA for the UniRef50\_O67749 family is insufficient to infer a genus, however the genus-level assignment from the pangenome search (i.e. to g\_\_Aquifex.s\_\_Aquifex_pyrophilus) is carried through.

The modified gene families output files can then be reprocessed through HUMAnN 3.0 to compute pathway abundance/coverage using the inferred taxonomic stratifications.

**Note:** ``humann_infer_taxonomy`` requires a data file to link UniRef families to LCAs and infer taxonomic relationships. Users can download these files using the HUMAnN 3.0 databases script:

```
$ humann_databases --download utility_mapping full $DIR
```

Where ``$DIR`` is the location where you'd like to store the mapping files. In addition to downloading the mapping files, this command also registers the mapping files with ``humann_infer_taxonomy``.

----

### humann_join_tables ###

* **Basic usage:** ``$ humann_join_tables --input $INPUT_DIR --output $TABLE``
* `$INPUT_DIR` = a directory containing gene/pathway tables (tsv or biom format)
* `$TABLE` = the file to write the new single gene table (biom format if input is biom format)
* Optional: ``--file_name $STR`` will only join gene tables with $STR in file name
* Run with `-h` to see additional command line options

This will join (merge) multiple single-sample output files into a single table with multiple samples. It can also be used to join multiple multi-sample output files into a single table.

----

### humann_reduce_table ###

Used for collapsing joined MetaPhlAn taxonomic profiles to a single joint profile. See the [Joint taxonomic profile](#markdown-header-joint-taxonomic-profile) guide below.

----

### humann_regroup_table ###

* **Basic usage:** `$ humann_regroup_table --input $TABLE --groups $GROUPS --output $TABLE2`
* `$TABLE` = gene/pathway table (tsv format)
* `$GROUPS` = options for regrouping table features (ex. uniref50 to metacyc reaction)
* `$TABLE2` = regrouped gene/pathway table
* Run with `-h` to see additional command line options

HUMAnN 3.0 `genefamilies.tsv` output can contain a very large number of features depending on the complexity of your underlying sample. One way to explore this information in a simplified manner is via HUMAnN 3.0's own pathway coverage and abundance, which summarize the values of their member genes. However, this approach does not apply to gene families that are not associated with metabolic pathways.

To further simplify the exploration of gene family abundance data, users can regroup gene families into other functional categories using `humann_regroup_table`. This script takes as arguments a gene family abundance table and a mapping (groups) file that indicates which gene families belong to which groups. Out of the box, HUMAnN 3.0 can regroup gene families to MetaCyc reactions (a step which is also used internally as part of MetaCyc pathway quantification). Users can download additional mapping files using the HUMAnN 3.0 databases script:

```
$ humann_databases --download utility_mapping full $DIR
```

Where ``$DIR`` is the location where you'd like to store the mapping files. In addition to downloading the mapping files, this command also registers the mapping files with `humann_regroup_table` (running the script with the `-h` flag will reveal an expanded list of regrouping options).

Mappings are available for both UniRef90 and UniRef50 gene families to the following systems:

* MetaCyc Reactions
* KEGG Orthogroups (KOs)
* Pfam domains
* Level-4 enzyme commission (EC) categories
* EggNOG (including COGs)
* Gene Ontology (GO)
* Informative GO

In most cases, mappings are directly inferred from the annotation of the corresponding UniRef centroid sequence in UniProt.

One exception to this are the "informative GO" (`infogo1000`) maps: These are informative subsets of GO computed from UniProt's annotations and the structure of the GO hierarchy specifically for HUMAnN 3.0 (each informative GO term has >1,000 UniRef centroids annotated to it, but none of its progeny terms have >1,000 centroids so annotated).

Users are free to create and use additional mapping files and pass them to `humann_regroup_table` via the `--custom` mapping flag. The format of a mapping file is:

```
group1  uniref1  uniref2  uniref3  ...
group2  uniref1  uniref5  ...
```

Where spaces between items above denote `TABS`. Mapping files can be gzip- or bzip2-compressed and will still be readable by `humann_regroup_table`. By default, feature abundances (such as gene families) are summed to produce group abundances.

If the "UNMAPPED" gene abundance feature is included in a user's input to ``humann_regroup_table``, it will automatically be carried forward to the final output. In addition, genes that do not group with a non-trivial feature are combined as an "UNGROUPED" group. By default, UNGROUPED reflects the total abundance of genes that did not belong to another group (similar in spirit to the "UNINTEGRATED" value reported in the pathway abundance file).

Some groups are not associated by default with human-readable names. To attach names to a regrouped table, use the `humann_rename_table` script. (The "GO" name map can be used for both raw GO and informative GO.)

----

### humann_rename_table ###

* **Basic usage:** ``$ humann_rename_table --input $TABLE --names $NAMES --output $TABLE2 ``
* `$TABLE` = gene/pathway table (tsv format)
* `$NAMES` = feature type to be renamed (run with "-h" to see included options)
* `$TABLE2` = gene/pathway table with new names attached
* Run with `-h` to see additional command line options

By default, larger HUMAnN 3.0 outputs do not attach names (glosses) to individual feature IDs to keep file sizes down. This script is used to attach those names. The most common use for this script is attaching gene names to UniRef90/UniRef50 features in the default `genefamilies.tsv` output files.

Features without names (common for many UniRef90/50 families) are assigned the default name: `NO_NAME`.

----

### humann_renorm_table ###

* ***Basic usage:*** ``$ humann_renorm_table --input $TABLE --units $CHOICE --output $TABLE2``
* `$TABLE` = gene/pathway table (tsv format)
* `$CHOICE` = "relab" (relative abundance) or "cpm" (copies per million)
* `$TABLE2` = normalized gene/pathway table
* Run with `-h` to see additional command line options

HUMAnN 3.0 quantifies genes and pathways in units of RPKs (reads per kilobase). These account for gene length but not sample sequencing depth. While there are some applications, e.g. strain profiling, where RPK units are superior to depth-normalized units, most of the time a user will renormalize their samples prior to downstream analysis.

This script provides the choice to normalize to relative abundance or copies per million (CPM) units. Both of these represent "total sum scaling (TSS)"-style normalization: in the former case, each sample is constrained to sum to 1, whereas in the latter case (CPMs) samples are constrained to sum to 1 million. Units out of 1 million are often more convenient for tables with many, many features (such as `genefamilies.tsv` tables).

**Note:** CPM as used here does *not* refer to unnormalized COUNTS per million, but rather copies per million. CPMs as used here are a generic analog of the TPM (transcript per million) unit in RNA-seq. You may wish to use the abbreviation CoPM for added clarity.

By default, this script normalizes all stratification levels to the sum of all community feature totals, but other options (such as level-wise normalization) are supported. "Special" features (such as UNMAPPED) can be included or excluded in the normalization process.

----

### humann_split_stratified_table ###

* **Basic usage:** `` $ humann_split_stratified_table --input $TABLE --output $OUTPUT_DIR ``
* `$TABLE` = gene/pathway table (tsv or biom format)
* `$OUTPUT_DIR` = the directory to write two (one stratified and one unstratified) new gene/pathway tables (in biom format if input is biom format)
* Run with `-h` to see additional command line options

This utility will split a table into two files (one stratified and one unstratified).

----

### humann_split_table ###

* **Basic usage:** `` $ humann_split_table --input $TABLE --output $OUTPUT_DIR ``
* `$TABLE` = gene/pathway table (tsv or biom format)
* `$OUTPUT_DIR` = the directory to write new gene/pathway tables (one per sample, in biom format if input is biom format)
* Run with `-h` to see additional command line options

This utility will split a merged feature table (multiple samples) into one file per sample. Some analyses can only accept one sample at a time.

----

### humann_test ###

Utility for executing HUMAnN 3.0 tests. See [Test the install](#markdown-header-3-test-the-install) above for further details.

----

### humann_unpack_pathways ###

* **Basic usage:** `$ humann_unpack_pathways --input-genes $INPUT_GENES.tsv --input-pathways $INPUT_PATHWAYS.tsv --output $OUTPUT.tsv`
* `$INPUT_GENES.tsv` = a file containing the gene families (or EC) abundance table (tsv format)
* `$INPUT_PATHWAYS.tsv` = a file containing the pathways abundance table (tsv format)
* `$OUTPUT.tsv` = the file to write the new abundance table (tsv format)
* Optional: `--remove-taxonomy` remove the taxonomy from the output file
* Run with `-h` to see additional command line options

This utility will unpack the pathways to show the genes for each. It adds another level of stratification to the pathway abundance table by including the gene families (or EC) abundances.

----

## Other HUMAnN 3.0 guides ##

### Selecting a level of gene family resolution ###

HUMAnN 3.0 uses UniRef protein clusters as a gene families system. UniRef clusters are constructed by clustering proteins from UniProt to remove redundancy (UniRef100), further clustering non-redundant proteins at 90% identity and selecting representative sequences (UniRef90), and further clustering UniRef90 representative sequences at 50% identity to produce broader clusters (UniRef50). The representative of a given UniRef cluster is generally the best-annotated member of the cluster (which may or may not be the true centroid of the cluster). Additional information about UniRef can be found [at the UniRef website](http://www.uniprot.org/help/uniref) and in the [original UniRef publication](http://www.ncbi.nlm.nih.gov/pubmed/17379688).

HUMAnN 3.0 can be configured to output gene family abundance at the resolution of UniRef90 clusters or UniRef50 clusters. The mode can be manually configured via the ``--search-mode`` flag. By default, the search mode is set based on your current translated search database. In addition to changing the resolution of the gene families output, the search mode optimizes HUMAnN 3.0 for searching against UniRef90 versus UniRef50 during translated search (e.g. toggling the minimum allowed amino acid mapping identity between 90% and 50%, respectively).

#### Should I pick UniRef90 or UniRef50? ####

For most applications, **UniRef90 is a good default choice** because its clusters are more conserved and therefore more likely to be isofunctional. UniRef50 clusters can be very broad, and thus pose the risk that the cluster representative might not reflect the function of the homologous sequence(s) found in your sample. One drawback of UniRef90 (relative to UniRef50) is that a slightly higher fraction of coding sequences in pangenomes are not assignable to a UniRef90 cluster, resulting in an increase in abundance for the ``UniRef90_unknown`` feature.

UniRef50 is preferable when dealing with very poorly characterized microbiomes. In such cases, requiring reads to map at 90% identity to UniRef90 may be too stringent, and so mapping at 50% identity to UniRef50 is expected to explain a larger portion of sample reads (albeit at reduced functional resolution).

Notably, performance is similar when mapping against UniRef90 versus UniRef50: while the former database is larger, the associated 90% percent identity alignment threshold results in fewer spurious seed events compared to aligning at 50% identity, which increases overall performance.

----

### Selecting a scope for translated search ###

HUMAnN 3.0 falls back to translated search for reads that failed to align to a known pangenome. The scope of this translated search can be tuned, resulting in a balance between the fraction of unclassified reads mapped and performance. There are three possible modes:

* **Bypass translated search** is selected via the ``--bypass-translated-search`` flag. In this mode, reads that failed to map to a pangenome are saved, but not otherwise searched at the protein level. There will be no ``unclassified`` strata in your output. 

* **Filtered translated search** is selected by using a EC-filtered protein database. These databases (which come in UniRef90 and UniRef50 flavors) have been filtered to only include proteins associated with a level-4 enzyme commission (EC) category. These are the only proteins that would be able to contribute to downstream pathway reconstruction. Thus, filtered translated search will provide a comprehensive metabolic reconstruction for the ``unclassified`` stratum, but the number of ``unclassified`` protein families in the output will be limited to those with EC annotation.

* **Comprehensive translated search** is selected by using a comprehensive protein database. These databases (which come in UniRef90 and UniRef50 flavors) have not been filtered, and are quite large (millions of sequences). Comprehensive search aims to explain as many sample reads as possible using reference-based approaches. Comprehensive search takes ~5x longer than filtered search [a difference of 5 versus 25 CPU-hours for a 10M read sample (using DIAMOND and UniRef50)]. Memory requirements are also larger for the larger database, though the increase is sublinear.

See [above](#markdown-header-download-a-translated-search-database) for instructions on how to acquire a translated search database.

#### Which translated search scope should I pick? ####

For many users, **filtered translated search** will serve as a good default option. This will still provide 1,000s of gene-level features and an uncompromised metabolic network for the ``unclassified`` stratum. For extremely well-characterized samples, bypassing translated search is a reasonable option. For example, in the healthy human gut, the majority of sample reads (~80%) that *can* be mapped to a reference are mapped prior to translated search. 

----

### HUMAnN 3.0 and paired-end sequencing data ###

End-pairing relationships are currently not taken into account during HUMAnN 3.0's alignment steps. This is due to the fact that HUMAnN 3.0 strictly aligns reads to isolated coding sequences: either at the nucleotide level or through translated search. As such, it will frequently be the case that one read (`READ1`) will map inside a given coding sequence while its mate-pair (`READ2`) will not.

Example:
```
GENEGENEGENE
     READ1-------READ2
```

Penalizing such cases would be overly strict: in the absence of a the gene's genomic context, this looks like a perfectly reasonable alignment (`READ2` may fall in a non-coding region and not align, or it may align to another [isolated] coding sequence). *As a result, the best way to use paired-end sequencing data with HUMAnN 3.0 is simply to concatenate all reads into a single FASTA or FASTQ file.*

If you have paired-end Illumina sequencing data, processed by CASAVA v1.8+ the sequence identifiers for the pairs have the same id after truncating by space. To keep track of the individual reads throughout the HUMAnN 3.0 workflow, the software will initially remove the spaces from the identifiers. This prevents read pairs from having the same id after being processed by bowtie2 and diamond. If this is the case for your read set, you will see an informative message printed to the screen and the log file indicating the removal of spaces from the sequence identifiers.

----

### PICRUSt output ###

You can run HUMAnN 3.0 with [PICRUSt](http://picrust.github.io/picrust/) output from predict_metagenomes.py or metagenome_contributions.py as input. Output from metagenome_contributions.py can include taxonomy information which will be used by HUMAnN 3.0. The steps that follow are the same for output files from either PICRUSt script.

For information on the PICRUSt software, please see the [project website](http://picrust.github.io/picrust/). For questions about PICRUSt, please contact the [PICRUSt user group](https://groups.google.com/forum/#!forum/picrust). The instructions that follow are for PICRUSt version 1.0.0-dev.

If you are running HUMAnN 3.0 with [PICRUSt](http://picrust.github.io/picrust/) output as input, please follow these steps:

1. Download the legacy kegg databases included in [HUMAnN 1.0](https://github.com/biobakery/humann_legacy/archive/0.99b.tar.gz)
    * The databases will be referred to in steps that follow with the path "humann1/data/*".

2. Split the picrust output file (picrust.biom or picrust.tsv) into a single file per sample (written to $OUTPUT_DIR)
    * `` $ humann_split_table --input picrust.biom --output $OUTPUT_DIR ``
    * The option `` --taxonomy_index -1 `` can be added if taxonomy information is included in the biom input file with column -1 associated with K0s.
    * If using biom input files, biom version 2.1+ must be installed.

3. Run HUMAnN 3.0 on each of the new files in $OUTPUT_DIR placing the results in $OUTPUT_DIR2
    * for $SAMPLE.biom in $OUTPUT_DIR
        * `` $ humann --input $SAMPLE.biom --output $OUTPUT_DIR2 --pathways-database humann1/data/keggc ``
    * The option ``--remove-stratified-output`` can be added if you do not want the data stratified by bug.
    * The option ``--output-format biom`` can be added if you want the output to be in biom format.
    * To run with the kegg modules instead of kegg pathways provide the modules file ``--pathways-database humann1/data/modulec``.
    * The input file can be in biom or tsv format.

4. Join the pathways data (coverage and abundance) files from the HUMAnN 3.0 runs from all samples into two files
    * `` $ humann_join_tables --input $OUTPUT_DIR2 --output humann_pathcoverage.tsv --file_name pathcoverage ``
    * `` $ humann_join_tables --input $OUTPUT_DIR2 --output humann_pathabundance.tsv --file_name pathabundance ``
    * The resulting files from these commands are named humann_pathcoverage.tsv and humann_pathabundance.tsv .
    * If the files being joined in this step are biom format, the ouput file will also be in biom format.

Please note the flag ``--verbose`` can be added to all commands.

----

### Joint taxonomic profile ###

A joint taxonomic profile can be created from all of the samples in your set. To create this file and use it for your HUMAnN 3.0 runs, please use the steps that follow.

1. Create taxonomic profiles for each of the samples in your set with [MetaPhlAn](https://bitbucket.org/biobakery/metaphlan)

2. Join all of the taxonomic profiles, located in directory $DIR, into a table of taxonomic profiles for all samples (joined_taxonomic_profile.tsv)
    * `` $ humann_join_tables --input $DIR --output joined_taxonomic_profile.tsv ``

3. Reduce this file into a taxonomic profile that represents the maximum abundances from all of the samples in your set
    * `` $ humann_reduce_table --input joined_taxonomic_profile.tsv --output max_taxonomic_profile.tsv --function max --sort-by level ``

4. Run HUMAnN 3.0 on all of the samples in your set, providing the max taxonomic profile
    * for $SAMPLE.fastq in samples
        * `` $ humann --input $SAMPLE.fastq --output $OUTPUT_DIR --taxonomic-profile max_taxonomic_profile.tsv ``

An alterative to step 4, which will save computing time, is to first run a single sample with the taxonomic profile. The HUMAnN 3.0 temp output folder for this sample will contain the bowtie2 indexed custom ChocoPhlAn database that can be provided when running your remaining samples. This will save compute time as this database will only be created once. Please see the steps below for the alternative to step 4.
    
1. Run HUMAnN 3.0 on one of your samples ($SAMPLE_1.fastq) providing the max taxonomic profile to create the custom indexed ChocoPhlAn database
    * `` $ humann --input $SAMPLE_1.fastq --output $OUTPUT_DIR --taxonomic-profile max_taxonomic_profile.tsv ``
    * The folder $OUTPUT_DIR/$SAMPLE_1_humann_temp/ will contain the custom indexed ChocoPhlAn database files
    
2. Run HUMAnN 3.0 on the rest of your samples providing the custom indexed ChocoPhlAn database ($OUTPUT_DIR/$SAMPLE_1_humann_temp/)
    * for $SAMPLE.fastq in samples
        * `` $ humann --input $SAMPLE.fastq --output $OUTPUT_DIR --nucleotide-database $OUTPUT_DIR/$SAMPLE_1_humann_temp/ --bypass-nucleotide-index ``

----

### Custom taxonomic profile ###

A custom taxonomic profile can be created to specify the taxa included in your samples. This file is used by HUMAnN 3.0 to create the custom ChocoPhlAn database for your samples.

The custom taxonomic profile must be in a tab-demilited format and contain two columns (taxon, with both genus and species, and percent abundance). An example follows:
```
g__Bacteroides|s__Bacteroides_thetaiotaomicron	12.16326
g__Bacteroides|s__Bacteroides_cellulosilyticus	12.02768
g__Bacteroides|s__Bacteroides_caccae	11.43394
g__Dialister|s__Dialister_invisus	10.52286
g__Bacteroides|s__Bacteroides_stercoris	10.42227
```

HUMAnN 3.0 uses the taxonomic profile to select pangenomes for the custom ChocoPhlAn database from the full ChocoPhlAn database. For example, the first line in the example above will add the centriods from Bacteroides thetaiotaomicron to the custom ChocoPhlAn database. Please note the taxa in the custom taxonomic profile must match the naming convention used by the full ChocoPhlAn database (ignoring case). 

To run HUMAnN 3.0 with the custom taxonomic profile ($FILE), use the option "--taxonomic-profile $FILE". This will bypass running MetaPhlAn, which creates a taxonomic profile, and instead will use the custom taxonomic profile provided. From the custom taxonomic profile, only those taxa that have a percent abundance greater than the default prescreen threshold will be considered. To change the default setting for the prescreen threshold to $THRESHOLD, use the option "--prescreen-threshold $THRESHOLD".

----

### Custom nucleotide reference database ###

A custom nucleotide reference database can be provided to HUMAnN 3.0 

This custom database must be formatted as a bowtie2 index. 

Please see the [Custom reference database annotations](#markdown-header-custom-reference-database-annotations) section for information on database annotations. Also please note, only alignments to genes included in the pathways databases will be considered in the pathways computations. The pathways databases included with HUMAnN 3.0 are for alignments to UniRef gene families. If you would like to create custom pathways databases for a different set of gene families, please see the [Custom pathways database](#markdown-header-custom-pathways-database) section for more information.

To run HUMAnN 3.0 with your custom nucleotide reference database (located in $DIR), use the option "--bypass-nucleotide-index" and provide the custom database as the ChocoPhlAn option with "--nucleotide-database $DIR". If you would like to bypass the translated alignment portion of HUMAnN 3.0, add the option "--bypass-translated-search". 

----

### Custom protein reference database ###

A custom protein reference database can be provided to HUMAnN 3.0 

This custom database must be formatted to be used by the translated alignment software selected (the default is Diamond). 

Please see the [Custom reference database annotations](#markdown-header-custom-reference-database-annotations) section for information on database annotations. Also please note, only alignments to genes included in the pathways databases will be considered in the pathways computations. The pathways databases included with HUMAnN 3.0 are for alignments to UniRef gene families. If you would like to create custom pathways databases for a different set of gene families, please see the [Custom pathways database](#markdown-header-custom-pathways-database) section for more information.

To run HUMAnN 3.0 with your custom protein reference database (located in $DIR), provide the custom database as the UniRef option with "--protein-database $DIR". Please note, HUMAnN 3.0 will run on all of the databases in this folder ($DIR) which have been formatted to be used by the translated alignment software selected. Also if you would like to bypass the nucleotide alignment portion of HUMAnN 3.0, add the option "--bypass-nucleotide-search".  

**Note:** HUMAnN can work with a protein database that's been split into chunks for computational efficiency (mapping the reads to them serially). Therefore HUMAnN models a protein database as a folder with one or more database chunks in it (derived from the same input set of protein sequences). If you mix unrelated chunks in the same folder HUMAnN will (correctly) complain at you. Hence if you want to have two separate DBs, even if each is only represented by a single chunk, they should be in separate folders.

----

### Custom reference database annotations ###

The annotations for sequences in a custom (nucleotide or protein) reference database can be as follows (delimited by "|"):

1. gene_family
2. gene_family|gene_length
3. gene_family|gene_length|taxonomy
4. identifier

For options 1 and 2, HUMAnN 3.0 defaults will be used. The default gene length is 1,000 bases and the default taxonomy is "unclassified".

Option 4 should be used along with a custom reference database annotation file. The custom reference database annotation file maps the identifiers to annotations. This file must be in a tab-delimited format and contain at least two columns (identifier and gene family). At most four columns of information can be included to describe each reference sequence. These columns should be organized as identifier, gene family, gene length, and taxonomy. An example follows:
```
256402719	UniRef50_C9LQU5	147	g__Dialister.s__Dialister_invisus
479150083	UniRef50_R6U703	540	g__Ruminococcus.s__Ruminococcus_bromii
423220654	UniRef50_I8UUJ6	1218	g__Bacteroides.s__Bacteroides_caccae
```

The first line of the example will use the gene family UniRef50_C9LQU5, gene length 147, and taxon ``g__Dialister.s__Dialister_invisus`` for any sequences in your reference databases with the identifier 256402719.

To run HUMAnN 3.0 with the custom reference database annotations ($FILE), use the option "--id-mapping $FILE". 

----

### Custom pathways database ###

The pathways databases included with HUMAnN 3.0 (from MetaCyc and UniProt) have been created to be used with alignments to UniRef gene families. A custom pathways database can be provided to HUMAnN 3.0, specifically made to work with your custom reference database(s).

One or two pathways database files (in a comma-delimited list) can be provided to HUMAnN 3.0 with the option "--pathways-database $FILE". If two files are provided, the first file provides a tab-delimited mapping while the second file provides the pathways mapping. For example, the first file could provide a mapping of gene families to reactions while the second file maps reactions to pathways. 

The first file, which is optional, should be organized to include at least two columns per line. The first column is the item to be mapped to (ie reactions from the example) while the remaining columns in the row are the gene families which can be mapped to the item in the first column. An example follows:
```
RXN-123	UniRef50_A0B6Z6	UniRef50_A3CRP6
RXN-456	UniRef50_A2RVM0	UniRef50_A4IGM4	UniRef50_A6NKP2	UniRef50_B8H806
```

The pathways file is a tab-delimited file with the first column the name of the pathway. It can be in a structured or unstructured format. In a structured format, the second column includes the items (space-delimited) contained in the pathway. The structure follows the same definition as that for Kegg modules. Each structure is a list (space-delimited) of items with a comma to indicate alternatives, a plus to indicate a complex, and a minus sign at the beginning of an item to indicate this is not a key item in the pathway. These items are gene families if only one pathways file is provided. If two files are provided, as in the example, these items would be reactions. In an unstructured format, the second column and any remaining columns in the row are the items that map to the pathway. For example, if a pathway contained four genes the line for this pathway in an unstructured format file would contain five columns. The first column is the pathway name and the remaining columns would each contain a single gene.

An example of a structured pathways file follows:
```
PWY-1	A B ( C , D )
PWY-2	( ( A + B ) , ( C + D ) ) E
```

An example of an unstructured pathways file follows:
```
PWY-3	A	B	C	D
PWY-4	A	B	C	D	E
```

----

### Core diversity analysis with QIIME ###

HUMAnN 3.0 output files can be provided to [QIIME](http://qiime.org/) as input to run core diversity analysis. For information on the QIIME software, please see the [project website](http://qiime.org/). For questions about QIIME, please contact the [QIIME user group](https://groups.google.com/forum/#!forum/qiime-forum). The instructions that follow are for QIIME version 1.9.1.


To run this analysis, run the following steps:

1. Run HUMAnN 3.0 on each of the samples (replacing $OUTPUT_DIR with the full path to the folder to write the output)
    * For each $SAMPLE.fastq in the set of all samples
        * `` $ humann --input $SAMPLE.fastq --output $OUTPUT_DIR --output-format biom --remove-stratified-output --output-max-decimals 0 ``
        * Each sample will have three main output files ($SAMPLE_genefamilies.tsv, $SAMPLE_pathabundance.tsv, and $SAMPLE_pathcoverage.tsv) in biom format.

2. Merge the three output files for each sample into three output files for all samples using [QIIME's merge_otu_tables.py](http://qiime.org/scripts/merge_otu_tables.html)
    * For this example, assume there are 3 samples named $SAMPLE1.fastq, $SAMPLE2.fastq, and $SAMPLE3.fastq
        * `` $ merge_otu_tables.py -i $SAMPLE1_genefamilies.biom,$SAMPLE2_genefamilies.biom,$SAMPLE3_genefamilies.biom -o genefamilies_all.biom ``
        * `` $ merge_otu_tables.py -i $SAMPLE1_pathabundance.biom,$SAMPLE2_pathabundance.biom,$SAMPLE3_pathabundance.biom -o pathabundance_all.biom ``
        * `` $ merge_otu_tables.py -i $SAMPLE1_pathcoverage.biom,$SAMPLE2_pathcoverage.biom,$SAMPLE3_pathcoverage.biom -o pathcoverage_all.biom ``

3. For each of the three merged biom files, run [QIIME's biom summarize-table](http://biom-format.org/documentation/summarizing_biom_tables.html) to obtain the sampling depth (e-value) required as input for the next step.
    * `` $ biom summarize-table -i genefamilies_all.biom -o genefamilies_summary.txt ``
    * `` $ biom summarize-table -i pathabundance_all.biom -o pathabundance_summary.txt ``
    * `` $ biom summarize-table -i pathcoverage_all.biom -o pathcoverage_summary.txt ``

4. Next run each merged biom file through [QIIME's core_diversity_analyses.py](http://qiime.org/scripts/core_diversity_analyses.html), providing the e-value from the prior step (replacing $EVALUE for each input file) and the [QIIME mapping file](http://qiime.org/documentation/file_formats.html#mapping-file-overview) you created (replacing $MAPPING_FILE with the full path to the mapping file).
    * `` $ core_diversity_analyses.py -i genefamilies_all.biom -o core_diversity_genefamilies -m $MAPPING_FILE --nonphylogenetic_diversity -e $EVALUE --suppress_taxa_summary ``
    * `` $ core_diversity_analyses.py -i pathabundance_all.biom -o core_diversity_pathabundance -m $MAPPING_FILE --nonphylogenetic_diversity -e $EVALUE --suppress_taxa_summary ``
    * `` $ core_diversity_analyses.py -i pathcoverage_all.biom -o core_diversity_pathcoverage -m $MAPPING_FILE --nonphylogenetic_diversity -e $EVALUE --suppress_taxa_summary ``

The output of [QIIME's core_diversity_analyses.py](http://qiime.org/scripts/core_diversity_analyses.html) will include PCoA plots of beta diversity analyses which can be viewed in the beta diversity folder, in the bray_curtis_emperor_pcoa_plot folder. Click on index.html, which will open a web browser with the plot. Using the tools on the right side, it is possible to change the colors by which the samples are labeled. 

Under the color tab, there are two variables that can be changed: 

1. the actual colors (Classic QIIME Colors is the default)
2. the mapping file category by which the samples are labeled. 

The more categories in your mapping file, the more options you have to see if your samples are separating based on that feature.

----

### Metatranscriptome analysis with HUMAnN 3.0 ###

The recommended HUMAnN 3.0 metatranscriptome (RNA) analysis protocol differs depending on whether or not you have metagenomic (DNA) reads available from the same biological sample.

**Analyzing a metatranscriptome with a paired metagenome.** In this case, we recommend analyzing the metagenome first. Then, rather than constructing a taxonomic profile directly from the metatranscriptome, use the taxonomic profile of the corresponding metagenome as an additional input to HUMAnN 3.0 via the ``--taxonomic-profile`` flag. This will guarantee that RNA reads are mapped to any species' pangenomes detected in the metagenome. RNA reads are otherwise provided as input to HUMAnN 3.0 just as DNA reads are. HUMAnN 3.0 RNA-level outputs (e.g. transcript family abundance) can then be normalized by corresponding DNA-level outputs to quantify microbial expression independent of gene copy number. **CAVEAT:** For low-abundance species, random sampling may lead to detection of transcripts for undetected genes. In these cases, we recommend smoothing DNA-level features to avoid divide-by-zero errors during normalization.

**Analyzing an isolated metatranscriptome (without a paired metagenome).** In this case, analyze RNA read data just as you would DNA data (provided as a fasta/fastq file). **CAVEAT 1:** Note that species are quantified in HUMAnN 3.0 based on recruitment of reads to species-specific marker genes. While each genome copy is assumed to donate ~1 copy of each marker to metagenome (DNA) data, the same assumption cannot be made for RNA data (markers may be more or less transcribed within a species compared to the species average). As long as a non-trivial fraction of a species' markers are expressed, HUMAnN 3.0 will still detect that species in the transcript pool. However, species relative abundance estimates from the taxonomic profile must be interpretted carefully: these values reflect species' relative contributions to the pool of species-specific transcripts, and not the overall transcript pool. **CAVEAT 2:** Transcript abundance inferred from a lone metatranscriptome is confounded with underlying gene copy number. For example, transcript X may be more abundant in sample A relative to sample B because (i) the same number of underlying X genes are more highly expressed in sample A relative to sample B or (ii) there are more copies of gene X in sample A relative to sample B (all of which are equally expressed). This is a general challenge in analyzing isolated metatranscriptomes (not specific to HUMAnN 3.0).

----

### Genus level gene families and pathways ###

By default, the gene families and pathways output files from HUMAnN 3.0 are species level. To obtain genus level gene families and pathways, follow these steps.

1. Create a genus level gene families file
    * `` $ humann_gene_families_genus_level --input $SAMPLE_genefamilies.tsv --output $SAMPLE_genefamilies_genus_level.tsv ``
    * In this command, replace ``$SAMPLE_genefamilies.tsv`` with the species level gene families file created by default by HUMAnN 3.0 and ``$SAMPLE_genefamilies_genus_level.tsv`` with the name of the gene families genus level file that will be created.

2. Run HUMAnN 3.0, with the genus level gene families file as input, to get genus level pathways output files
    * `` $ humann --input $SAMPLE_genefamilies_genus_level.tsv --output humann_genus_level_output ``
    * This run will be much faster and require less memory than the original run as HUMAnN 3.0 is provided gene family abundances so it only needs to compute the pathways.

----

## FAQs ##

HUMAnN 3.0 frequently asked questions:

1.  Is there a way to print more information to stdout during the run?
    *   Yes, add the ``--verbose`` flag
2.  How do I make use of multiple cores on the same machine?
    *   Add the ``--threads $CORES`` option
3.  How do I remove the intermediate temp output files?
    *   Add the ``--remove-temp-output`` flag
4.  Can I provide an alternative location for the ChocoPhlAn database?
    *   Yes, use the ``--nucleotide-database $DIR`` option
5.  Can I provide an alternative location for the UniRef database?
    *   Yes, use the ``--protein-database $DIR`` option
6.  I already have MetaPhlAn output. Can I start HUMAnN 3.0 with the MetaPhlAn output?
    *   Yes, use the ``--taxonomic-profile bugs_list.tsv`` option
7.  Is there a way to change $SAMPLENAME in the output file names?
    *   Yes, use the ``--output-basename $NAME`` option
8.  How do I remove the stratification by bug from the output files?
    *   Add the ``--remove-stratified-output`` flag
9.  How do I run with the unipathways databases?
    *   Add the ``--pathways unipathway`` option
10.  Is there a way to output files in biom format?
    *   Yes, use the ``--output-format biom`` option
11.  Can I change the identity threshold for alignments?
    *   Yes, use the ``--identity-threshold <50.0>`` option
12.  Can I change the options provided to MetaPhlAn?
    *   Yes, use the ``--metaphlan-options="-t rel_ab"`` option
    *   Please note the special formatting required for this option

----

## Complete option list ##

```
usage: humann [-h] -i <input.fastq> -o <output> [--threads <1>] [--version]
              [-r] [--bypass-nucleotide-index] [--bypass-nucleotide-search]
              [--bypass-prescreen] [--bypass-translated-search]
              [--taxonomic-profile <taxonomic_profile.tsv>]
              [--memory-use {minimum,maximum}]
              [--input-format {fastq,fastq.gz,fasta,fasta.gz,sam,bam,blastm8,genetable,biom}]
              [--search-mode {uniref50,uniref90}] [-v]
              [--metaphlan <metaphlan>]
              [--metaphlan-options <metaphlan_options>]
              [--prescreen-threshold <0.01>] [--bowtie2 <bowtie2>]
              [--bowtie-options <bowtie_options>]
              [--nucleotide-database <nucleotide_database>]
              [--nucleotide-identity-threshold <0.0>]
              [--nucleotide-query-coverage-threshold <90.0>]
              [--nucleotide-subject-coverage-threshold <50.0>]
              [--diamond <diamond>] [--diamond-options <diamond_options>]
              [--evalue <1.0>] [--pick-frames {on,off}]
              [--protein-database <protein_database>]
              [--rapsearch <rapsearch>]
              [--translated-alignment {usearch,rapsearch,diamond}]
              [--translated-identity-threshold <Automatically: 50.0 or 80.0, Custom: 0.0-100.0>]
              [--translated-query-coverage-threshold <90.0>]
              [--translated-subject-coverage-threshold <50.0>]
              [--usearch <usearch>] [--gap-fill {on,off}] [--minpath {on,off}]
              [--pathways {metacyc,unipathway}]
              [--pathways-database <pathways_database.tsv>] [--xipe {on,off}]
              [--annotation-gene-index <3>] [--id-mapping <id_mapping.tsv>]
              [--remove-temp-output]
              [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
              [--o-log <sample.log>] [--output-basename <sample_name>]
              [--output-format {tsv,biom}] [--output-max-decimals <10>]
              [--remove-column-description-output]
              [--remove-stratified-output]

HUMAnN : HMP Unified Metabolic Analysis Network 3

optional arguments:
  -h, --help            show this help message and exit

[0] Common settings:
  -i <input.fastq>, --input <input.fastq>
                        input file of type {fastq,fastq.gz,fasta,fasta.gz,sam,bam,blastm8,genetable,biom} 
                        [REQUIRED]
  -o <output>, --output <output>
                        directory to write output files
                        [REQUIRED]
  --threads <1>         number of threads/processes
                        [DEFAULT: 1]
  --version             show program's version number and exit

[1] Workflow refinement:
  -r, --resume          bypass commands if the output files exist
  --bypass-nucleotide-index
                        bypass the nucleotide index step and run on the indexed ChocoPhlAn database
  --bypass-nucleotide-search
                        bypass the nucleotide search steps
  --bypass-prescreen    bypass the prescreen step and run on the full ChocoPhlAn database
  --bypass-translated-search
                        bypass the translated search step
  --taxonomic-profile <taxonomic_profile.tsv>
                        a taxonomic profile (the output file created by metaphlan)
                        [DEFAULT: file will be created]
  --memory-use {minimum,maximum}
                        the amount of memory to use
                        [DEFAULT: minimum]
  --input-format {fastq,fastq.gz,fasta,fasta.gz,sam,bam,blastm8,genetable,biom}
                        the format of the input file
                        [DEFAULT: format identified by software]
  --search-mode {uniref50,uniref90}
                        search for uniref50 or uniref90 gene families
                        [DEFAULT: based on translated database selected]
  -v, --verbose         additional output is printed

[2] Configure tier 1: prescreen:
  --metaphlan <metaphlan>
                        directory containing the MetaPhlAn software
                        [DEFAULT: $PATH]
  --metaphlan-options <metaphlan_options>
                        options to be provided to the MetaPhlAn software
                        [DEFAULT: "-t rel_ab"]
  --prescreen-threshold <0.01>
                        minimum percentage of reads matching a species
                        [DEFAULT: 0.01]

[3] Configure tier 2: nucleotide search:
  --bowtie2 <bowtie2>   directory containing the bowtie2 executable
                        [DEFAULT: $PATH]
  --bowtie-options <bowtie_options>
                        options to be provided to the bowtie software
                        [DEFAULT: "--very-sensitive"]
  --nucleotide-database <nucleotide_database>
                        directory containing the nucleotide database
                        [DEFAULT: /home/ljmciver/miniconda3/envs/humann3/lib/python3.6/site-packages/humann/data/chocophlan_DEMO]
  --nucleotide-identity-threshold <0.0>
                        identity threshold for nuclotide alignments
                        [DEFAULT: 0.0]
  --nucleotide-query-coverage-threshold <90.0>
                        query coverage threshold for nucleotide alignments
                        [DEFAULT: 90.0]
  --nucleotide-subject-coverage-threshold <50.0>
                        subject coverage threshold for nucleotide alignments
                        [DEFAULT: 50.0]

[3] Configure tier 2: translated search:
  --diamond <diamond>   directory containing the diamond executable
                        [DEFAULT: $PATH]
  --diamond-options <diamond_options>
                        options to be provided to the diamond software
                        [DEFAULT: "--top 1 --outfmt 6"]
  --evalue <1.0>        the evalue threshold to use with the translated search
                        [DEFAULT: 1.0]
  --protein-database <protein_database>
                        directory containing the protein database
                        [DEFAULT: /home/ljmciver/miniconda3/envs/humann3/lib/python3.6/site-packages/humann/data/uniref_DEMO]
  --rapsearch <rapsearch>
                        directory containing the rapsearch executable
                        [DEFAULT: $PATH]
  --translated-alignment {usearch,rapsearch,diamond}
                        software to use for translated alignment
                        [DEFAULT: diamond]
  --translated-identity-threshold <Automatically: 50.0 or 80.0, Custom: 0.0-100.0>
                        identity threshold for translated alignments
                        [DEFAULT: Tuned automatically (based on uniref mode) unless a custom value is specified]
  --translated-query-coverage-threshold <90.0>
                        query coverage threshold for translated alignments
                        [DEFAULT: 90.0]
  --translated-subject-coverage-threshold <50.0>
                        subject coverage threshold for translated alignments
                        [DEFAULT: 50.0]
  --usearch <usearch>   directory containing the usearch executable
                        [DEFAULT: $PATH]

[5] Gene and pathway quantification:
  --gap-fill {on,off}   turn on/off the gap fill computation
                        [DEFAULT: on]
  --minpath {on,off}    turn on/off the minpath computation
                        [DEFAULT: on]
  --pathways {metacyc,unipathway}
                        the database to use for pathway computations
                        [DEFAULT: metacyc]
  --pathways-database <pathways_database.tsv>
                        mapping file (or files, at most two in a comma-delimited list) to use for pathway computations
                        [DEFAULT: metacyc database ]
  --xipe {on,off}       turn on/off the xipe computation
                        [DEFAULT: off]
  --annotation-gene-index <3>
                        the index of the gene in the sequence annotation
                        [DEFAULT: 3]
  --id-mapping <id_mapping.tsv>
                        id mapping file for alignments
                        [DEFAULT: alignment reference used]

[6] More output configuration:
  --remove-temp-output  remove temp output files
                        [DEFAULT: temp files are not removed]
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        level of messages to display in log
                        [DEFAULT: DEBUG]
  --o-log <sample.log>  log file
                        [DEFAULT: temp/sample.log]
  --output-basename <sample_name>
                        the basename for the output files
                        [DEFAULT: input file basename]
  --output-format {tsv,biom}
                        the format of the output files
                        [DEFAULT: tsv]
  --output-max-decimals <10>
                        the number of decimals to output
                        [DEFAULT: 10]
  --remove-column-description-output
                        remove the description in the output column
                        [DEFAULT: output column includes description]
  --remove-stratified-output
                        remove stratification from output
                        [DEFAULT: output is stratified]
```
----
## Contributions ##
Thanks go to these wonderful people:
	
## Support ##
 Please sign up for the [HUMAnN category in bioBakery Forum](https://forum.biobakery.org/c/Microbial-community-profiling/HUMAnN) if any questions or concerns.   
 
 Additionally, Google user group: <humann-users@googlegroups.com> (Read only) is available.  
