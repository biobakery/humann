# HUMAnN2 User Manual #

HUMAnN2 is the next generation of HUMAnN (HMP Unified Metabolic Analysis Network).

**If you use the HUMAnN2 software, please cite our manuscript: TBD**

HUMAnN is a pipeline for efficiently and accurately profiling the presence/absence and abundance of microbial pathways in a community from metagenomic or metatranscriptomic sequencing data (typically millions of short DNA/RNA reads). This process, referred to as functional profiling, aims to describe the metabolic potential of a microbial community and its members. More generally, functional profiling answers the question "What are the microbes in my community-of-interest doing (or capable of doing)?"

## Contents ##

* [Features](#markdown-header-features)
* [Workflows](#markdown-header-workflows)
    * [Main workflow](#markdown-header-main-workflow)
    * [Workflow by input file type](#markdown-header-workflow-by-input-file-type)
    * [Workflow by bypass mode](#markdown-header-workflow-by-bypass-mode)
    * [Workflow of the resume option](#markdown-header-workflow-of-the-resume-option)
* [Requirements](#markdown-header-requirements)
    * [Software](#markdown-header-software)
    * [Other](#markdown-header-other)
* [Installation](#markdown-header-installation)
    1. [Downloading HUMAnN2](#markdown-header-1-downloading-humann2)
    2. [Installing HUMAnN2](#markdown-header-2-installing-humann2)
    3. [Test the install](#markdown-header-3-test-the-install)
    4. [Try out a demo run](#markdown-header-4-try-out-a-demo-run)
    5. [Download the databases](#markdown-header-5-download-the-databases)
        * [Download the ChocoPhlAn database](#markdown-header-download-the-chocophlan-database)
        * [Download the UniRef50 database](#markdown-header-download-the-uniref50-database)
* [How to run](#markdown-header-how-to-run)
    * [Basic usage](#markdown-header-basic-usage)
    * [Demo runs](#markdown-header-demo-runs)
* [Output files](#markdown-header-output-files)
    1. [Gene families file](#markdown-header-1-gene-families-file)
    2. [Pathway coverage file](#markdown-header-2-pathway-coverage-file)
    3. [Pathway abundance file](#markdown-header-3-pathway-abundance-file)
    4. [Intermediate temp output files](#markdown-header-4-intermediate-temp-output-files)
        1. [Bowtie2 alignment results](#markdown-header-1-bowtie2-alignment-results)
        2. [Bowtie2 reduced alignment results](#markdown-header-2-bowtie2-reduced-alignment-results)
        3. [Bowtie2 index files](#markdown-header-3-bowtie2-index-files)
        4. [Unaligned reads after Bowtie2](#markdown-header-4-unaligned-reads-after-bowtie2)
        5. [Custom ChocoPhlAn database](#markdown-header-5-custom-chocophlan-database)
        6. [MetaPhlAn2 Bowtie2 output](#markdown-header-6-metaphlan2-bowtie2-output)
        7. [MetaPhlAn2 bugs list](#markdown-header-7-metaphlan2-bugs-list)
        8. [Translated alignment results](#markdown-header-8-translated-alignment-results)
        9. [Translated alignment unaligned reads](#markdown-header-9-translated-alignment-unaligned-reads)
        10. [Log](#markdown-header-10-log)    
* [Databases](#markdown-header-databases)
* [Configuration](#markdown-header-configuration)
* [Tools](#markdown-header-tools)
    * [Tools for tables](#markdown-header-tools-for-tables)
        1. [Split a table](#markdown-header-1-split-a-table)
        2. [Join tables](#markdown-header-2-join-tables)
        3. [Rename table features](#markdown-header-3-rename-table-features)
        4. [Renormalize table](#markdown-header-4-renormalize-table)
        5. [Regroup table features](#markdown-header-5-regroup-table-features)
        6. [Combine metagenomic and metatranscriptomic sequencing data](#markdown-header-6-combine-metagenomic-and-metatranscriptomic-sequencing-data)        
        7. [Strain-level functional profiling](#markdown-header-7-strain-level-functional-profiling)
        8. [Reduce table](#markdown-header-8-reduce-table)
* [Tutorials](#markdown-header-tutorials)
    * [PICRUSt output](#markdown-header-picrust-output)
    * [Legacy databases](#markdown-header-legacy-databases)
    * [Joint taxonomic profile](#markdown-header-joint-taxonomic-profile)
    * [Custom taxonomic profile](#markdown-header-custom-taxonomic-profile)
    * [Custom nucleotide reference database](#markdown-header-custom-nucleotide-reference-database)
    * [Custom protein reference database](#markdown-header-custom-protein-reference-database)
    * [Custom reference database annotations](#markdown-header-custom-reference-database-annotations)
    * [Custom pathways database](#markdown-header-custom-pathways-database)
    * [Analyzing metatranscriptomes](#markdown-header-analyzing-metatranscriptomes)
    * [Strain-level functional profiling](#markdown-header-strain-level-functional-profiling)
* [FAQs](#markdown-header-faqs)
* [Complete option list](#markdown-header-complete-option-list)


## Features ##


1. Community functional profiles stratified by known and unclassified organisms

    * [MetaPhlAn2](http://huttenhower.sph.harvard.edu/metaphlan2) and [ChocoPhlAn pangenome database](http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/chocophlan.tar.gz) are used to facilitate fast, accurate, and organism-specific functional profiling
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


## Workflows ##

### Main workflow ###

![](http://huttenhower.sph.harvard.edu/sites/default/files/humann2_diamond_500x500.jpg)


### Workflow by input file type ###

There are four different types of files that can be provided as input to HUMAnN2 . By default HUMAnN2 will determine the type of the file. As shown in the figure below, the type of input file will determine where HUMAnN2 will start the workflow. Files of type #2, #3, and #4 will begin the workflow after the alignment steps.

![](http://huttenhower.sph.harvard.edu/sites/default/files/humann2_flow_by_file_type_reduced.png)

File Types:

*   File Type #1 (a quality-controlled metagenome or metatranscriptome)
    *   fastq (fastq.gz)
    *   fasta (fasta.gz)
*   File Type #2 (alignment results type #1)
    *   sam
    *   bam
*   File Type #3 (alignment results type #2)
    *   blast-like tsv
*   File Type #4 (gene table)
    *   tsv
    *   biom


### Workflow by bypass mode ###

There are multiple bypass options that will allow you to adjust the standard workflow.

![](http://huttenhower.sph.harvard.edu/sites/default/files/humann2_flow_bypass_modes_updated.png)

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
    *   starts the workflow with the nucleotide alignment step using the indexed database from "--chocophlan $DIR"


### Workflow of the resume option ###

HUMAnN2 includes a "--resume" option which will allow you to bypass alignment steps which have already been completed. For example, if you originally ran with a bypass option you can run just the step you bypassed with "--resume". This will only run the alignment step you bypassed and then recompute the gene families and pathways.

![](http://huttenhower.sph.harvard.edu/sites/default/files/humann2_flow_resume_option_no_text.png)

When using the "--resume" option, the following steps will be bypassed if they have already been completed:

1.  Taxomonic profiling step
2.  Nucleotide alignment step
3.  Custom ChocoPhlAn database creation (merge and index)
4.  Translated alignment step


## Requirements ##

### Software ###

1. [MetaPhlAn2](https://bitbucket.org/biobakery/metaphlan2/)
2. [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) (version >= 2.2)
3. [Diamond](http://ab.inf.uni-tuebingen.de/software/diamond/) (version >= 0.7.3)
4. [Python](http://www.python.org/) (version >= 2.7)
5. [MinPath](http://omics.informatics.indiana.edu/MinPath/) (automatically downloaded/installed)
6. [Xipe](https://edwards.sdsu.edu/cgi-bin/xipe.cgi) (optional / included)
7. [Rapsearch2](http://omics.informatics.indiana.edu/mg/RAPSearch2/) (version >= 2.21) (only required if using rapsearch2 for translated search)
8. [Usearch](http://www.drive5.com/usearch/) (version >= 7.0) (only required if using usearch for translated search)
9. [SAMtools](http://samtools.sourceforge.net/) (only required if bam input files are provided)
10. [Biom-format](http://biom-format.org/) (only required if input or output files are in biom format)

Please install the required software in a location in your $PATH or provide the location with an optional argument to HUMAnN2. 
For example, the location of the Bowtie2 install ($BOWTIE2_DIR) can be provided with "--bowtie2 $BOWTIE2_DIR".

If you always run with input files of type #2, #3, and #4 (for information on input file types, see section [Workflow by input file type](#markdown-header-workflow-by-input-file-type)),
MetaPhlAn2, Bowtie2, and Diamond are not required. Also if you always run with one or more bypass options (for information on bypass options, see section [Workflow by bypass mode](#markdown-header-workflow-by-bypass-mode)), 
the software required for the steps you bypass does not need to be installed.

### Other ###

1. Memory (>= 16 Gb)
2. Disk space (>= 10 Gb [to accommodate comprehensive sequence databases])
3. Operating system (Linux or Mac)

If always running with files of type #2, #3, and #4 (for information on file types, see section [Workflow by input file type](#markdown-header-workflow-by-input-file-type)),
less disk space is required. 

## Installation ##

### 1. Downloading HUMAnN2 ###
You can download the latest HUMAnN2 release or the development version.

Option 1: Latest Release (Recommended)

* [Download](https://bitbucket.org/biobakery/humann2/downloads/humann2_v0.1.9.tar.gz) and unpack the latest release of HUMAnN2.

Option 2: Development Version

* Create a clone of the repository: 
    
	``hg clone https://bitbucket.org/biobakery/humann2 ``

	Note: Creating a clone of the repository requires [Mercurial](http://mercurial.selenic.com/) to be installed. 


### 2. Installing HUMAnN2 ###

1. Move to the HUMAnN2 directory

    * ``$ cd $HUMAnN2_PATH `` 

2. Install MinPath

    * ``$ python setup.py minpath ``
    * If you would like to update the glpk required by MinPath, add the option ``--update-glpk`` to the MinPath install command.
    * The glpk update is required if you are running on Mac OS and it will require gcc and make.
    
3. Install HUMAnN2

    * ``$ python setup.py install ``
    * If you do not have write permissions to '/usr/lib/', then add the option ``--user`` to the install command. This will install the python package into subdirectories of '~/.local'. Please note when using the "--user" install option on some platforms, you might need to add '~/.local/bin/' to your $PATH as it might not be included by default. You will know if it needs to be added if you see the following message `humann2: command not found` when trying to run HUMAnN2 after installing with the "--user" option.


### 3. Test the install ###

Test out the install of HUMAnN2 by running the unit tests.

``$ python setup.py test ``

### 4. Try out a demo run ###

With HUMAnN2 installed you can try out a demo run using reduced versions of the databases.

``$ humann2 --input examples/demo.fastq --output $OUTPUT_DIR ``

Output from this demo run will be written to the folder $OUTPUT_DIR.

Please continue with the install directions to download the full databases before running with your sequencing data.


### 5. Download the databases ###

Downloading the databases is a required step if your input is a filtered shotgun sequencing metagenome file (fastq, fastq.gz, fasta, or fasta.gz format). If your input files will always be mapping results files (sam, bam or blastm8 format) or gene tables (tsv or biom format), you do not need to download the ChocoPhlAn and UniRef50 databases. 

#### Download the ChocoPhlAn database ####

Download the ChocoPhlAn database providing $INSTALL_LOCATION as the location to install the database (approximate size = 5.6 GB).

`` $ humann2_databases --download chocophlan full $INSTALL_LOCATION ``

NOTE: The humann2 config file will be updated to point to this location for the default chocophlan database. If you move this database, please use the "humann2_config" command to update the default location of this database. Alternatively you can always provide the location of the chocophlan database you would like to use with the "--chocophlan <chocophlan>" option to humann2.


#### Download the UniRef50 database ####

Download the UniRef50 database providing $INSTALL_LOCATION as the location to install the database (approximate size = 2.8 GB).

`` $ humann2_databases --download uniref diamond $INSTALL_LOCATION ``

NOTE: The humann2 config file will be updated to point to this location for the default uniref database. If you move this database, please use the "humann2_config" command to update the default location of this database. Alternatively you can always provide the location of the uniref database you would like to use with the "--uniref <uniref>" option to humann2.

NOTE: By default HUMAnN2 runs diamond for translated alignment. If you would like to use rapsearch2 for translated alignment, first download the rapsearch2 formatted database by running this command with the rapsearch2 formatted database selected. It is suggested that you install both databases in the same folder so this folder can be the default uniref database location. This will allow you to switch between alignment software without having to specify a different location for the database.


## How to run ##

### Basic usage ###

```
$ humann2 --input $SAMPLE --output $OUTPUT_DIR
```


$SAMPLE = a single file that is one of the following types:

1.  filtered shotgun sequencing metagenome file (fastq, fastq.gz, fasta, or fasta.gz format)
2.  alignment file (sam, bam or blastm8 format)
3.  gene table file (tsv or biom format)

$OUTPUT_DIR = the output directory

**Three output files will be created:**

1. $OUTPUT_DIR/$SAMPLENAME_genefamilies.tsv*
2. $OUTPUT_DIR/$SAMPLENAME_pathcoverage.tsv
3. $OUTPUT_DIR/$SAMPLENAME_pathabundance.tsv

where $SAMPLENAME is the basename of $SAMPLE

*The gene families file will not be created if the input file type is a gene table.

**Intermediate temp files will also be created:**

1. $DIR/$SAMPLENAME_bowtie2_aligned.sam
	* the full alignment output from bowtie2 
2. $DIR/$SAMPLENAME_bowtie2_aligned.tsv
	* only the reduced aligned data from the bowtie2 output
3. $DIR/$SAMPLENAME_bowtie2_index*
	* bowtie2 index files created from the custom chochophlan database
4. $DIR/$SAMPLENAME_bowtie2_unaligned.fa 
	* a fasta file of unaligned reads after the bowtie2 step
5. $DIR/$SAMPLENAME_custom_chocophlan_database.ffn 
	* a custom chocophlan database of fasta sequences
6. $DIR/$SAMPLENAME_metaphlan_bowtie2.txt 
	* the bowtie2 output from metaphlan
7. $DIR/$SAMPLENAME_metaphlan_bugs_list.tsv 
	* the bugs list output from metaphlan
8. $DIR/$SAMPLENAME_$TRANSLATEDALIGN_aligned.tsv 
	* the alignment results from the translated alignment step
9. $DIR/$SAMPLENAME_$TRANSLATEDALIGN_unaligned.fa 
	* a fasta file of unaligned reads after the translated alignment step
10. $DIR/$SAMPLENAME.log 
	* a log of the run
	
* $DIR=$OUTPUT_DIR/$SAMPLENAME_humann2_temp/
* $SAMPLENAME is the basename of the fastq/fasta input file
* $TRANSLATEDALIGN is the translated alignment software selected (rapsearch2 or usearch)

NOTE: $SAMPLENAME can be set by the user with the option "--output-basename <$NEWNAME>". 

### Demo runs ###

The examples folder contains four demo example input files. These files are of fasta, fastq, sam, and blastm8 format. Blastm8 format is created by the following software: rapsearch2, usearch, and blast.


To run the fasta demo:

`` $ humann2 --input examples/demo.fasta --output $OUTPUT_DIR``

To run the fastq demo:

`` $ humann2 --input examples/demo.fastq --output $OUTPUT_DIR``

To run the sam demo:

`` $ humann2 --input examples/demo.sam --output $OUTPUT_DIR``

To run the blastm8 demo:

`` $ humann2 --input examples/demo.m8 --output $OUTPUT_DIR``

$OUTPUT_DIR is the output directory

Since sam and blastm8 are mapping results, using these files as input to HUMAnN2 will bypass both the nucleotide and translated mapping portions of the flow.

## Output files ##

When HUMAnN2 is run, three main output files will be created (where `` $SAMPLENAME = the basename of $SAMPLE ``):

### 1. Gene families file ###

``` 
# Gene Family   $SAMPLENAME
UniRef50_A6L0N6	67.0
UniRef50_A6L0N6|g__Bacteroides.s__Bacteroides_fragilis	8.0
UniRef50_A6L0N6|g__Bacteroides.s__Bacteroides_finegoldii	5.0
UniRef50_A6L0N6|g__Bacteroides.s__Bacteroides_stercoris	4.0
UniRef50_A6L0N6|unclassified	1.0
UniRef50_G9S1V7	60.0
UniRef50_G9S1V7|g__Bacteroides.s__Bacteroides_vulgatus	31.0
UniRef50_G9S1V7|g__Bacteroides.s__Bacteroides_thetaiotaomicron	22.0
UniRef50_G9S1V7|g__Bacteroides.s__Bacteroides_stercoris	7.0
```

*   File name: `` $OUTPUT_DIR/$SAMPLENAME_genefamilies.tsv ``
*   This file quantifies the abundance of each gene family in the community. Gene families are groups of evolutionarily-related protein-coding sequences that typically perform similar functions. Abundance is reported in RPK (reads per kilobase) units to normalize for gene length; RPK units reflect relative gene (or transcript) copy number in the community.
*   In addition to community-wide gene family abundance totals (as reported by HUMAnN), this file is stratified to indicate abundance contributions of known and unclassified organisms represented in the sample.
*   Please note the gene families file will not be created if the input file type is a gene table.
        
### 2. Pathway coverage file ###

``` 
# Pathway       $SAMPLENAME
PWY0-1301	1.0
PWY0-1301|g__Bacteroides.s__Bacteroides_caccae	1.0
PWY0-1301|g__Bacteroides.s__Bacteroides_finegoldii	1.0
PWY0-1301|unclassified	1.0
PWY-7134	1.0
PWY-7134|g__Bacteroides.s__Bacteroides_vulgatus	0.7
PWY-7134|g__Bacteroides.s__Bacteroides_thetaiotaomicron	0.7
PWY-7134|unclassified	0.3
PWY-7134|g__Parabacteroides.s__Parabacteroides_merdae	0.3
```

*   File name: `` $OUTPUT_DIR/$SAMPLENAME_pathcoverage.tsv ``
*   This file details the presence/absence of each pathway in the community. HUMAnN refers to pathway presence/absence as "coverage" and defines a pathway as a set of two or more gene families.
*   In addition to community-wide pathway coverage (as reported by HUMAnN), this file is stratified to indicate the coverage of the pathway by genomes of known and unclassified organisms represented in the sample.

### 3. Pathway abundance file ###

```
# Pathway       $SAMPLENAME
PWY-1921	57.5
PWY-1921|unclassified	32.5
PWY-1921|g__Bacteroides.s__Bacteroides_ovatus	4.5
PWY-1921|g__Alistipes.s__Alistipes_putredinis	3.0
PWY-1921|g__Bacteroides.s__Bacteroides_caccae	2.5
PWY0-1301	54.7
PWY0-1301|unclassified	16.7
PWY0-1301|g__Parabacteroides.s__Parabacteroides_merdae	8.0
PWY0-1301|g__Bacteroides.s__Bacteroides_caccae	6.0
```
         
*   File name: `` $OUTPUT_DIR/$SAMPLENAME_pathabundance.tsv ``
*   This file quantifies the abundance of each pathway in the community as a function of the abundance of its member gene families.
*   In addition to community-wide pathway abundance (as reported by HUMAnN), this file is stratified to indicate abundance contributions of known and unclassified organisms represented in the sample.

### 4. Intermediate temp output files ###

Ten intermediate temp output files will be created where:

```
$DIR = $OUTPUT_DIR/$SAMPLENAME_humann2_temp/
$SAMPLENAME = basename of the fastq/fasta input file named $SAMPLE
$TRANSLATEDALIGN = translated alignment software selected (diamond, rapsearch2 or usearch)
```

NOTE: $SAMPLENAME can be set by the user with the option --output-basename <$NEWNAME>

#### 1. Bowtie2 alignment results ####

```
@HD	VN:1.0	SO:unsorted
@SQ	SN:g__Ruminococcus.s__Ruminococcus_bromii|UniRef90_D4L6K4|UniRef50_R6U703	LN:540
r99491	0	g__Bacteroides.s__Bacteroides_stercoris|UniRef90_R6B629|UniRef50_R5RCC8	1015	42	151M	*	0	0	$SEQ	$QUAL
r99526	0	g__Parabacteroides.s__Parabacteroides_merdae|UniRef90_unknown|UniRef50_D9RX34	155	42	151M	*	0	0	$SEQ	$QUAL
r99581	16	g__Bacteroides.s__Bacteroides_stercoris|UniRef90_unknown|UniRef50_R6SXR7	2503	42	151M	*	0	0	$SEQ	$QUAL
```

*   In example above `` $SEQ = sequence and $QUAL = quality scores `` to fit in page.
*   File name: `` $DIR/$SAMPLENAME_bowtie2_aligned.sam `` 
*   This file has the full alignment output from bowtie2.

#### 2. Bowtie2 reduced alignment results ####

``` 
r93	g__Bacteroides.s__Bacteroides_cellulosilyticus|UniRef90_E2NEW2|UniRef50_E2NEW2	6.3095734448e-05
r113	g__Bacteroides.s__Bacteroides_cellulosilyticus|UniRef90_R6KNZ3|UniRef50_R6KNZ3	6.3095734448e-05	
r704	g__Bacteroides.s__Bacteroides_uniformis|UniRef90_unknown|UniRef50_E6STE9		0.794328234724	
r663	g__Bacteroides.s__Bacteroides_thetaiotaomicron|UniRef90_R7KKH7|UniRef50_R7KKH7		6.3095734448e-05	
r940	g__Ruminococcus.s__Ruminococcus_bromii|UniRef90_unknown|UniRef50_unknown		6.3095734448e-0	 
```

*   File name: `` $DIR/$SAMPLENAME_bowtie2_aligned.tsv ``
*   This file contains the minimal amount of alignment results from Bowtie2.

#### 3. Bowtie2 index files ####


*   Example not included as files are binary.
*   File name: `` $DIR/$SAMPLENAME_bowtie2_index* ``
*   These are a set of files containing the Bowtie2 index created from the custom ChocoPhlAn database.

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

#### 6. MetaPhlAn2 Bowtie2 output ####

```
r113	gi|224485636|ref|NZ_EQ973490.1|:c728571-728107
r559	gi|479185170|ref|NC_021030.1|:c1678719-1677127
r663	gi|512436175|ref|NZ_KE159463.1|:c142391-139122
r704	gi|423310881|ref|NZ_JH724270.1|:c220428-218656
r1086	gi|238922432|ref|NC_012781.1|:c1988048-1987140 
```

*   File name: `` $DIR/$SAMPLENAME_metaphlan_bowtie2.txt ``
*   This file is the Bowtie2 output from MetaPhlAn2.

#### 7. MetaPhlAn2 bugs list ####

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
*   This file is the bugs list output from MetaPhlAn2.

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
*   This file is formatted as tab-delimited blast-like results.

#### 9. Translated alignment unaligned reads ####

```
>r4370
GGCGGACGATCTTGTCGCCCAGCCTGTAGCCTTTCTGGTACACCGTGATGACGGTGCCGCTCTCCTGCCCGTCCGTGGCGGGGATCTGCTGG
>r4398
TGCCCGGACAGGATCTTCTCTTTCGTACCGGGCATCATCTGCTCCATGATCTCCACGCCTCGCATGAACTTTTCAGAACGGGCAACGTAGGA
```
    
*   File name: `` $DIR/$SAMPLENAME_$TRANSLATEDALIGN_unaligned.fa ``
*   This is a fasta file of the unaligned reads after the translated alignment step

#### 10. Log ####

```
03/16/2015 01:09:52 PM - humann2.utilities - INFO: File ( demo.fastq ) is of format:  fastq
03/16/2015 01:09:52 PM - humann2.config - INFO: Run config settings:
DATABASE SETTINGS
chocophlan database folder = data/chocophlan_DEMO
uniref database folder = data/uniref_DEM
```

*   File name: `` $DIR/$SAMPLENAME.log ``
*   This file is a log of the run.


## Databases ##

HUMAnN2 uses two databases for alignment, ChocoPhlAn and UniRef. There are different formats of these databases. The demo formats of both are included in the HUMAnN2 install.

To see which databases are currently available, run (with example output included below):
```sh
$ humann2_databases --available
HUMANnN2 Databases ( database : build = location )
chocophlan : DEMO = http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/DEMO_chocophlan.tar.gz
chocophlan : full = http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/full_chocophlan.tar.gz
uniref : DEMO_diamond = http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref50_GO_filtered/uniref50_DEMO_diamond.tar.gz
uniref : diamond = http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref50_GO_filtered/uniref50_GO_filtered_diamond.tar.gz
uniref : rapsearch2 = http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref50_GO_filtered/uniref50_GO_filtered_rapsearch2.tar.gz
```
To download a database, run:
```sh
$ humann2_databases --download $DATABASE $BUILD $INSTALL_LOCATION
```

This will automatically update the HUMAnN2 configuration. 


## Configuration ##

HUMAnN2 uses a configuration file to store user configuration settings. This configuration file is automatically updated when a database is installed.

To view the current settings of the configuration file, run (with example output included below):
```sh
$ humann2_config --print
HUMAnN2 Configuration ( Section : Name = Value )
output_format : remove_stratified_output = False
output_format : output_max_decimals = 10
alignment_settings : average_read_length = 1
alignment_settings : prescreen_threshold = 0.01
alignment_settings : evalue_threshold = 1.0
alignment_settings : identity_threshold = 0.4
database_folders : chocophlan = data/chocophlan_DEMO
database_folders : uniref = data/uniref_DEMO
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
$ humann2_config --update $SECTION $NAME $VALUE 
```

The settings included in the configuration file plus others selected at run-time are printed to the beginning of the log file for each run.

Example log file configuration settings section:
```sh
03/16/2015 01:09:52 PM - humann2.config - INFO: 
Run config settings: 

DATABASE SETTINGS
chocophlan database folder = data/chocophlan_DEMO
uniref database folder = data/uniref_DEMO
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
average read length = 1
prescreen threshold = 0.01
identity threshold = 0.4

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


## Tools ##


### Tools for tables ###

HUMAnN2 includes tools to be used with gene, pathway, and taxonomic profile table files.

#### 1.  Split a table ####

`` $ humann2_split_table --input $TABLE --output $OUTPUT_DIR ``

*   $TABLE = gene/pathway table (tsv or biom format)
*   $OUTPUT_DIR = the directory to write new gene/pathway tables (one per sample, in biom format if input is biom format)

#### 2.  Join tables ####

`` $ humann2_join_tables --input $INPUT_DIR --output $TABLE ``

*   $INPUT_DIR = a directory containing gene/pathway tables (tsv or biom format)
*   $TABLE = the file to write the new single gene table (biom format if input is biom format)
*   Optional: ``--file_name $STR`` will only join gene tables with $STR in file name

#### 3.  Rename table features ####

`` $ humann2_rename_table --input $TABLE --names $NAMES --output $TABLE2 ``

*   $TABLE = gene/pathway table (tsv format)
*   $NAMES = mapping of feature IDs to english names (tsv format)
*   $TABLE2 = gene/pathway table with new names attached
*   Note: A mapping of UniRef50 IDs to english names is provided with HUMAnN2

#### 4.  Renormalize table ####

`` $ humann2_renorm_table --input $TABLE --norm $CHOICE --output $TABLE2 ``

*   $TABLE = gene/pathway table (tsv format)
*   $CHOICE = "relab" (relative abundance) or "cpm" (copies per million)
*   $TABLE2 = normalized gene/pathway table

#### 5.  Regroup table features ####

`` $ humann2_regroup_table --input $TABLE --groups $GROUPS --output $TABLE2 ``

*   $TABLE = gene/pathway table (tsv format)
*   $GROUPS = mapping of features to superfeatures (.tsv or .tsv.gz format)
*   $TABLE2 = regrouped gene/pathway table

#### 6.  Combine metagenomic and metatranscriptomic sequencing data ####

`` $ humann2_rna_dna_norm --input_dna $TABLE_DNA --input_rna --groups $TABLE_DNA --output_basename $OUTPUT_BASENAME ``

*   $TABLE_DNA = gene families output for sample DNA (tsv format)
*   $TABLE_RNA = gene families output for sample RNA (tsv format)
*   $OUTPUT_BASENAME = Basename for the three output files (smoothed DNA, smoothed RNA, normalized RNA)
*   See the [tutorial](#markdown-header-analyzing-metatranscriptomes) below for assistance.

#### 7.  Strain-level functional profiling ####

`` $ humann2_strain_profiler --input $TABLE ``

*   $TABLE = merged gene families output for two or more samples (tsv format)
*   See the [tutorial](#markdown-header-strain-level-functional-profiling) below for assistance.

#### 8.  Reduce table ####

`` $ humann2_reduce_table --input $INPUT.tsv --output $OUTPUT.tsv ``

*   $INPUT.tsv = a file containing the table (tsv format)
*   $OUTPUT.tsv = the file to write the new reduced table (tsv format) 
*   Optional: ``--function {min|max|mean}`` the function to apply (default is max)
*   Optional: ``--sort-by {level|name|value}`` to indicate how the output should be sorted (default is original order)

## Tutorials ##

### PICRUSt output ###

You can run HUMAnN2 with [PICRUSt](http://picrust.github.io/picrust/) output from predict_metagenomes.py or metagenome_contributions.py as input. Output from 
metagenome_contributions.py can include taxonomy information which will be used by HUMAnN2. The steps that follow are the same for output files from either
PICRUSt script.

If you are running HUMAnN2 with [PICRUSt](http://picrust.github.io/picrust/) output as input, please follow these steps:

1. Download the legacy kegg databases included in [HUMAnN](https://bitbucket.org/biobakery/humann/downloads/humann-v0.99.tar.gz)

    * The databases will be referred to in steps that follow with the path "humann1/data/*".

2. Split the picrust output file (picrust.biom) into a single file per sample (written to $OUTPUT_DIR)

    * `` $ humann2_split_table --input picrust.biom --output $OUTPUT_DIR ``
    * The option `` --taxonomy_index -1 `` can be added if taxonomy information is included in the biom file with column -1 associated with K0s.

3. Run HUMAnN2 on each of the new files in $OUTPUT_DIR placing the results in $OUTPUT_DIR2

    * for $SAMPLE.biom in $OUTPUT_DIR
        * `` $ humann2 --input $SAMPLE.biom --output $OUTPUT_DIR2 --pathways-database humann1/data/keggc ``
    * The option ``--remove-stratified-output`` can be added if you do not want the data stratified by bug.
    * The option ``--output-format biom`` can be added if you want the output to be in biom format.
    
4. Join the pathways data (coverage and abundance) files from the HUMAnN2 runs from all samples into two files

    * `` $ humann2_join_tables --input $OUTPUT_DIR2 --output humann2_pathcoverage.tsv --file_name pathcoverage ``
    * `` $ humann2_join_tables --input $OUTPUT_DIR2 --output humann2_pathabundance.tsv --file_name pathabundance ``
    * The resulting files from these commands are named humann2_pathcoverage.tsv and humann2_pathabundance.tsv .
    * If the files being joined in this step are biom format, the ouput file will also be in biom format.

Please note the flag ``--verbose`` can be added to all commands.

### Legacy databases ###

The original version of HUMAnN used [Kegg](http://www.genome.jp/kegg/) databases. You can run with the legacy Kegg databases following these steps:

1. Download the legacy kegg databases included in [HUMAnN](https://bitbucket.org/biobakery/humann/downloads/humann-v0.99.tar.gz)

    * The databases will be referred to in steps that follow with the path "humann1/data/*".
    
2. Create an idmapping file formatted for HUMAnN2 using the legacy kegg databases and adding full names for the Kegg organisms

    * `` $ humann2_humann1_kegg --ikoc humann1/data/koc --igenels humann1/data/genels --o legacy_kegg_idmapping.tsv ``
    
3. Run HUMAnN2 on your kegg alignment file, $SAMPLE, with results written to $OUTPUT_DIR

    * `` $ humann2 --input $SAMPLE --output $OUTPUT_DIR --id-mapping legacy_kegg_idmapping.tsv --pathways-database humann1/data/keggc ``
    * To run with the kegg modules instead of kegg pathways provide the file ``humann1/data/modulec``.
    * For a demo run, provide the demo input file included with the original version of HUMAnN ``humann1/input/mock_even_lc.tsv``.

### Joint taxonomic profile ###

A joint taxonomic profile can be created from all of the samples in your set. To create this file and use it for your HUMAnN2 runs, 
please use the steps that follow.

1. Create taxonomic profiles for each of the samples in your set with [MetaPhlAn2](https://bitbucket.org/biobakery/metaphlan2/)
2. Join all of the taxonomic profiles, located in directory $DIR, into a table of taxonomic profiles for all samples (joined_taxonomic_profile.tsv)

    * `` $ humann2_join_tables --input $DIR --output joined_taxonomic_profile.tsv ``

3. Reduce this file into a taxonomic profile that represents the maximum abundances from all of the samples in your set

    * `` $ humann2_reduce_table --input joined_taxonomic_profile.tsv --output max_taxonomic_profile.tsv --function max --sort-by level ``

4. Run HUMAnN2 on all of the samples in your set, providing the max taxonomic profile

    * for $SAMPLE.fastq in samples
        * `` $ humann2 --input $SAMPLE.fastq --output $OUTPUT_DIR --taxonomic-profile max_taxonomic_profile.tsv ``

An alterative to step #4, which will save computing time, is to first run a single sample with the taxonomic profile. The HUMAnN2 temp output
folder for this sample will contain the bowtie2 indexed custom ChocoPhlAn database that can be provided when running your remaining samples. 
This will save compute time as this database will only be created once. Please see the steps below for the alternative to step #4.
    
1. Run HUMAnN2 on one of your samples ($SAMPLE_1.fastq) providing the max taxonomic profile to create the custom indexed ChocoPhlAn database

    * `` $ humann2 --input $SAMPLE_1.fastq --output $OUTPUT_DIR --taxonomic-profile max_taxonomic_profile.tsv ``
    * The folder $OUTPUT_DIR/$SAMPLE_1_humann2_temp/ will contain the custom indexed ChocoPhlAn database files
    
2. Run HUMAnN2 on the rest of your samples providing the custom indexed ChocoPhlAn database ($OUTPUT_DIR/$SAMPLE_1_humann2_temp/)

    * for $SAMPLE.fastq in samples
        * `` $ humann2 --input $SAMPLE.fastq --output $OUTPUT_DIR --chocophlan $OUTPUT_DIR/$SAMPLE_1_humann2_temp/ --bypass-nucleotide-index ``

### Custom taxonomic profile ###

A custom taxonomic profile can be created to specify the taxa included in your samples. This file is used by HUMAnN2 to create the custom ChocoPhlAn database for your samples.

The custom taxonomic profile must be in a tab-demilited format and contain two columns (taxon and percent abundance). An example follows:
```
g__Bacteroides|s__Bacteroides_thetaiotaomicron	12.16326
g__Bacteroides|s__Bacteroides_cellulosilyticus	12.02768
g__Bacteroides|s__Bacteroides_caccae	11.43394
g__Dialister|s__Dialister_invisus	10.52286
g__Bacteroides|s__Bacteroides_stercoris	10.42227
```

HUMAnN2 uses the taxonomic profile to select pangenomes for the custom ChocoPhlAn database from the full ChocoPhlAn database. For example, the first line in the example above will add the centriods from Bacteroides thetaiotaomicron to the custom ChocoPhlAn database. Please note the taxa in the custom taxonomic profile must match the naming convention used by the full ChocoPhlAn database (ignoring case). 

To run HUMAnN2 with the custom taxonomic profile ($FILE), use the option "--taxonomic-profile $FILE". This will bypass running MetaPhlAn2, which creates a taxonomic profile, and instead will use the custom taxonomic profile provided. From the custom taxonomic profile, only those taxa that have a percent abundance greater than the default prescreen threshold will be considered. To change the default setting for the prescreen threshold to $THRESHOLD, use the option "--prescreen-threshold $THRESHOLD".

### Custom nucleotide reference database ###

A custom nucleotide reference database can be provided to HUMAnN2. 

This custom database must be formatted as a bowtie2 index. 

Please see the [Custom reference database annotations](#markdown-header-custom-reference-database-annotations) section for information on database annotations. Also please note, only alignments to genes included in the pathways databases will be considered in the pathways computations. The
pathways databases included with HUMAnN2 are for UniRef genes. If you would like to create custom pathways databases for a different gene set, please see the [Custom pathways database](#markdown-header-custom-pathways-database) section for more information.

To run HUMAnN2 with your custom nucleotide reference database (located in $DIR), use the option "--bypass-nucleotide-index" and provide the custom database as the ChocoPhlAn option with "--chocophlan $DIR". If you would like to bypass the translated alignment portion of HUMAnN2, add the option "--bypass-translated-search". 

### Custom protein reference database ###

A custom protein reference database can be provided to HUMAnN2. 

This custom database must be formatted to be used by the translated alignment software selected (the default is Diamond). 

Please see the [Custom reference database annotations](#markdown-header-custom-reference-database-annotations) section for information on database annotations. Also please note, only alignments to genes included in the pathways databases will be considered in the pathways computations. The
pathways databases included with HUMAnN2 are for UniRef genes. If you would like to create custom pathways databases for a different gene set, please see the [Custom pathways database](#markdown-header-custom-pathways-database) section for more information.

To run HUMAnN2 with your custom protein reference database (located in $DIR), provide the custom database as the UniRef option with "--uniref $DIR". Please note, HUMAnN2 will run on all of the databases in this folder ($DIR) which have been formatted to be used by the translated alignment software selected. Also if you would like to bypass the nucleotide alignment portion of HUMAnN2, add the option "--bypass-nucleotide-search".  

### Custom reference database annotations ###

The annotations for sequences in a custom (nucleotide or protein) reference database can be as follows (delimited by "|"):

1. gene
2. gene|gene_length
3. gene|gene_length|taxonomy
4. identifier

For options #1 and #2, HUMAnN2 defaults will be used. The default gene length is 1,000 bases and the default taxonomy is "unclassified".

Option #4 should be used along with a custom reference database annotation file. The custom reference database annotation file maps the identifiers to annotations. This file must be in a tab-delimited format and contain at least two columns (identifer and gene). At most four columns of information can be 
included to describe each reference sequence. These columns should be organized as identifier, gene, gene length, and taxonomy. An example follows:
```
256402719	UniRef50_C9LQU5	147	g__Dialister.s__Dialister_invisus
479150083	UniRef50_R6U703	540	g__Ruminococcus.s__Ruminococcus_bromii
423220654	UniRef50_I8UUJ6	1218	g__Bacteroides.s__Bacteroides_caccae
```

The first line of the example will use the gene UniRef50_C9LQU5, gene length 147, and taxon ``g__Dialister.s__Dialister_invisus`` for any sequences in your
reference databases with the identifier 256402719.

To run HUMAnN2 with the custom reference database annotations ($FILE), use the option "--id-mapping $FILE". 

### Custom pathways database ###

The pathways databases included with HUMAnN2 (from MetaCyc and UniProt) have been created to be used with alignments to UniRef genes. A custom pathways
database can be provided to HUMAnN2, specifically made for your set of genes.

One or two pathways database files (in a comma-delimited list) can be provided to HUMAnN2 with the option "--pathways-database $FILE". If two files are
provided, the first file provides a tab-delimited mapping while the second file provides the pathways mapping. For example, the first file could provide a mapping of genes to reactions while the second file maps reactions to pathways. 

The first file, which is optional, should be organized to include at least three columns per line. The first column is the item to be mapped to (ie reactions from the example) while the second column (which can be blank) includes additional information about the item to be mapped to (ie EC number from the example). The remaining columns in the row are the genes which can be mapped to the item included in the first column.

The pathways file is a tab-delimited file with the first column the name of the pathway. The second column includes the items contained in the pathway. These items are genes if only one pathways file is provided. If two files are provided, as in the example, these items would be reactions.

### Analyzing Metatranscriptomes ###

The recommended HUMAnN2 metatranscriptome (RNA) analysis protocol differs depending on whether or not you have metagenomic (DNA) reads available from the same biological sample.

**Analyzing a metatranscriptome with a paired metagenome.** In this case, we recommend analyzing the metagenome first. Then, rather than constructing a taxonomic profile directly from the metatranscriptome, use the taxonomic profile of the corresponding metagenome as an additional input to HUMAnN2 via the ``--taxonomic-profile`` flag. This will guarantee that RNA reads are mapped to any species' pangenomes detected in the metagenome. RNA reads are otherwise provided as input to HUMAnN2 just as DNA reads are. HUMAnN2 RNA-level outputs (e.g. transcript family abundance) can then be normalized by corresponding DNA-level outputs to quantify microbial expression independent of gene copy number. **CAVEAT:** For low-abundance species, random sampling may lead to detection of transcripts for undetected genes. In these cases, we recommend smoothing DNA-level features to avoid divide-by-zero errors during normalization. The HUMAnN2 tool script ``humann2_rna_dna_norm`` will Witten-Bell smooth paired RNA and DNA output files and then return the smoothed DNA, smoothed RNA, and relative expression. The relative expression of gene *i* in sample *j* is the smoothed value of RNA(*i*,*j*) divided by the smoothed value of DNA(*i*,*j*) and adjusted for differences in sequencing depth. The ``humann2_rna_dna_norm`` script assumes: (i) that units in the table behave like counts (this is true for ``*_genefamilies.tsv`` files, which are measured in RPK) and (ii) that sample columns are in the same order.

**Analyzing an isolated metatranscriptome (without a paired metagenome).** In this case, analyze RNA read data just as you would DNA data (provided as a fasta/fastq file). **CAVEAT 1:** Note that species are quantified in HUMAnN2 based on recruitment of reads to species-specific marker genes. While each genome copy is assumed to donate ~1 copy of each marker to metagenome (DNA) data, the same assumption cannot be made for RNA data (markers may be more or less transcribed within a species compared to the species average). As long as a non-trivial fraction of a species' markers are expressed, HUMAnN2 will still detect that species in the transcript pool. However, species relative abundance estimates from the taxonomic profile must be interpretted carefully: these values reflect species' relative contributions to the pool of species-specific transcripts, and not the overall transcript pool. **CAVEAT 2:** Transcript abundance inferred from a lone metatranscriptome is confounded with underlying gene copy number. For example, transcript X may be more abundant in sample A relative to sample B because (i) the same number of underlying X genes are more highly expressed in sample A relative to sample B or (ii) there are more copies of gene X in sample A relative to sample B (all of which are equally expressed). This is a general challenge in analyzing isolated metatranscriptomes (not specific to HUMAnN2).

### Strain-level Functional Profiling ###

The HUMAnN2 script ``humann2_strain_profiler`` can help explore strain-level variation in your data. This approach assumes you have run HUMAnN2 on a series of samples and then merged the resulting ``*_genefamilies.tsv`` tables with ``humann2_merge_tables``. Cases will arise in which the same species was detected in two or more samples, but gene families within that species were not consistently present across samples. For example, four samples may contain the species *Dialister invisus*, but only two samples contain the gene family ``UniRef50_Q5WII6`` within *Dialister invisus*. This is a form of strain-level variation in the *Dialister invisus* species: one which we can connect directly to function based on annotations of the ``UniRef50_Q5WII6`` gene family.

``humann2_strain_profiler`` first looks for (species, sample) pairs where (i) a large number of gene families within the species were identified (default: 500) and (ii) the mean abundance of detected genes was high (default: mean > 10 RPK). For species that meet these criteria, we can infer that absent gene families are likely to be truly absent, as opposed to undersampled. Simulations suggest that the cutoff of 10 RPK results in a false negative rate below 0.001 (i.e. for every 1000 genes identified as absent, at most one would be present but missed due to undersampling). For a given species, if at least two samples pass these criteria, the species and passing samples are sliced from the merged table and saved as a strain profile.

Strain profiles can be additionally restricted to a subset of species (e.g. those from a particular genus) or to gene families with a high level of variability in the population (e.g. present in fewer than 80% of samples but more than 20% of samples). Additional thresholds (e.g. the minimum non-zero mean) can be configured with command line parameters.

## FAQs ##

HUMAnN2 frequently asked questions:

1.  Is there a way to print more information to stdout during the run?
    *   Yes, add the ``--verbose`` flag
2.  How do I make use of multiple cores on the same machine?
    *   Add the ``--threads $CORES`` option
3.  How do I remove the intermediate temp output files?
    *   Add the ``--remove-temp-output`` flag
4.  Can I provide an alternative location for the ChocoPhlAn database?
    *   Yes, use the ``--chocophlan $DIR`` option
5.  Can I provide an alternative location for the UniRef database?
    *   Yes, use the ``--uniref $DIR`` option
6.  I already have MetaPhlAn2 output. Can I start HUMAnN2 with the MetaPhlAn2 output?
    *   Yes, use the ``--taxonomic-profile bugs_list.tsv`` option
7.  Is there a way to change $SAMPLENAME in the output file names?
    *   Yes, use the ``--output-basename $NAME`` option
8.  How do I remove the stratification by bug from the output files?
    *   Add the ``--remove-stratified-output`` flag
9.  How do I run with the unipathways databases?
    *   Add the ``--pathways unipathway`` option
10.  Is there a way to output files in biom format?
    *   Yes, use the ``--output-format biom`` option
11.  Can I change the e-value threshold for alignments?
    *   Yes, use the ``--evalue <1.0>`` option

## Complete option list ##

```
usage: humann2 [-h] [--version] [-v] [-r] [--bypass-prescreen]
               [--bypass-nucleotide-index] [--bypass-translated-search]
               [--bypass-nucleotide-search] -i <input.fastq> -o <output>
               [-c <chocophlan>] [--chocophlan-gene-index <-1>] [-u <uniref>]
               [--average-read-length <1>] [--evalue <1.0>]
               [--metaphlan <metaphlan>] [--o-log <sample.log>]
               [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
               [--remove-temp-output] [--bowtie2 <bowtie2>] [--threads <1>]
               [--prescreen-threshold <0.01>] [--identity-threshold <0.4>]
               [--usearch <usearch>] [--rapsearch <rapsearch>]
               [--diamond <diamond>]
               [--taxonomic-profile <taxonomic_profile.tsv>]
               [--id-mapping <id_mapping.tsv>]
               [--translated-alignment {usearch,rapsearch,diamond}]
               [--xipe {on,off}] [--minpath {on,off}] [--pick-frames {on,off}]
               [--output-format {tsv,biom}] [--output-basename <sample_name>]
               [--remove-stratified-output]
               [--input-format {fastq,fastq.gz,fasta,fasta.gz,sam,bam,blastm8,genetable,biom}]
               [--pathways-database <pathways_database.tsv>]
               [--pathways {metacyc,unipathway}]
               [--memory-use {minimum,maximum}]

HUMAnN2 : HMP Unified Metabolic Analysis Network 2

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -v, --verbose         additional output is printed
  -r, --resume          bypass commands if the output files exist
  --bypass-prescreen    bypass the prescreen step and run on the full ChocoPhlAn database
  --bypass-nucleotide-index
                        bypass the nucleotide index step and run on the indexed ChocoPhlAn database
  --bypass-translated-search
                        bypass the translated search step
  --bypass-nucleotide-search
                        bypass the nucleotide search steps
  -i <input.fastq>, --input <input.fastq>
                        input file of type {fastq,fastq.gz,fasta,fasta.gz,sam,bam,blastm8,genetable,biom} 
                        [REQUIRED]
  -o <output>, --output <output>
                        directory to write output files
                        [REQUIRED]
  -c <chocophlan>, --chocophlan <chocophlan>
                        directory containing the ChocoPhlAn database
                        [DEFAULT: data/chocophlan_DEMO]
  --chocophlan-gene-index <-1>
                        the index of the gene in the sequence annotation
                        [DEFAULT: -1]
  -u <uniref>, --uniref <uniref>
                        directory containing the UniRef database
                        [DEFAULT: data/uniref_DEMO]
  --average-read-length <1>
                        the average length of the reads
                        [DEFAULT: 1]
  --evalue <1.0>        the evalue threshold
                        [DEFAULT: 1.0]
  --metaphlan <metaphlan>
                        directory containing the MetaPhlAn software
                        [DEFAULT: $PATH]
  --o-log <sample.log>  log file
                        [DEFAULT: temp/sample.log]
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        level of messages to display in log
                        [DEFAULT: DEBUG]
  --remove-temp-output  remove temp output files
                        [DEFAULT: temp files are not removed]
  --bowtie2 <bowtie2>   directory of the bowtie2 executable
                        [DEFAULT: $PATH]
  --threads <1>         number of threads/processes
                        [DEFAULT: 1]
  --prescreen-threshold <0.01>
                        minimum percentage of reads matching a species
                        [DEFAULT: 0.01]
  --identity-threshold <0.4>
                        identity threshold to use with the translated search
                        [DEFAULT: 0.4]
  --usearch <usearch>   directory containing the usearch executable
                        [DEFAULT: $PATH]
  --rapsearch <rapsearch>
                        directory containing the rapsearch executable
                        [DEFAULT: $PATH]
  --diamond <diamond>   directory containing the diamond executable
                        [DEFAULT: $PATH]
  --taxonomic-profile <taxonomic_profile.tsv>
                        a taxonomic profile (the output file created by metaphlan)
                        [DEFAULT: file will be created]
  --id-mapping <id_mapping.tsv>
                        id mapping file for alignments
                        [DEFAULT: alignment reference used]
  --translated-alignment {usearch,rapsearch,diamond}
                        software to use for translated alignment
                        [DEFAULT: diamond]
  --xipe {on,off}       turn on/off the xipe computation
                        [DEFAULT: off]
  --minpath {on,off}    turn on/off the minpath computation
                        [DEFAULT: on]
  --pick-frames {on,off}
                        turn on/off the pick_frames computation
                        [DEFAULT: off]
  --output-format {tsv,biom}
                        the format of the output files
                        [DEFAULT: tsv]
  --output-basename <sample_name>
                        the basename for the output files
                        [DEFAULT: input file basename]
  --remove-stratified-output
                        remove stratification from output
                        [DEFAULT: output is stratified]
  --input-format {fastq,fastq.gz,fasta,fasta.gz,sam,bam,blastm8,genetable,biom}
                        the format of the input file
                        [DEFAULT: format identified by software]
  --pathways-database <pathways_database.tsv>
                        mapping file (or files, at most two in a comma-delimited list) to use for pathway computations
                        [DEFAULT: metacyc database ]
  --pathways {metacyc,unipathway}
                        the database to use for pathway computations
                        [DEFAULT: metacyc]
  --memory-use {minimum,maximum}
                        the amount of memory to use
                        [DEFAULT: minimum]
```

