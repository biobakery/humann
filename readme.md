# HUMAnN2: HMP Unified Metabolic Analysis Network 2 #

HUMAnN is a pipeline for efficiently and accurately determining the presence/absence and abundance of microbial pathways in a community from metagenomic or metatranscriptomic sequencing data (typically millions of short DNA/RNA reads). 

HUMAnN2 is the next generation of HUMAnN. HUMAnN2 incorporates several new features including an expanded database of microbial genomes (>4x the size of the database included in HUMAnN), a simple user interface (single command driven flow), and bug-specific output files. 

The HUMAnN2 pipeline is a single command driven flow requiring the user to only provide a filtered fastq file to produce gene and pathway summary files ready for analysis. The pipeline converts sequence reads into coverage and abundance tables summarizing the gene families and pathways in one or more microbial communities. 

## Contents ##
* [Requirements](#markdown-header-requirements)
    * [Software](#markdown-header-software)
    * [Other](#markdown-header-other)
* [Installation](#markdown-header-installation)
    * [Downloading HUMAnN2](#markdown-header-downloading-humann2)
    * [Installing HUMAnN2](#markdown-header-installing-humann2)
    * [Test the install](#markdown-header-test-the-install)
    * [Try out a demo run](#markdown-header-try-out-a-demo-run)
    * [Download the databases](#markdown-header-download-the-databases)
        * [Download the ChocoPhlAn database](#markdown-header-download-the-chocophlan-database)
        * [Download the UniRef50 database](#markdown-header-download-the-uniref50-database)
* [How to run](#markdown-header-how-to-run)
    * [Basic usage](#markdown-header-basic-usage)
    * [Demo runs](#markdown-header-demo-runs)
* [Output files](#markdown-header-output-files)
    * [Gene families](#markdown-header-gene-families)
    * [Pathway coverage](#markdown-header-pathway-coverage)
    * [Pathway abundance](#markdown-header-pathway-abundance)
* [Complete option list](#markdown-header-complete-option-list)

## Requirements ##

### Software ###

1. [MetaPhlAn2](https://bitbucket.org/biobakery/metaphlan2/)
2. [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) (version >= 2.2)
3. [diamond](http://ab.inf.uni-tuebingen.de/software/diamond/) (version >= 0.7.3)
4. [Python](http://www.python.org/) (version >= 2.7)
5. [MinPath](http://omics.informatics.indiana.edu/MinPath/) (automatically downloaded/installed)
6. [Xipe](https://edwards.sdsu.edu/cgi-bin/xipe.cgi) (optional / included)
7. [rapsearch2](http://omics.informatics.indiana.edu/mg/RAPSearch2/) (version >= 2.21) (optional)
8. [usearch](http://www.drive5.com/usearch/) (version >= 7.0) (optional)

If MetaPhlAn2, bowtie2, and diamond (or rapsearch2, usearch) are not installed in a location in your $PATH,
then add them to your $PATH or use the HUMAnN2 parameters to indicate the locations of their
directories (--metaphlan $METAPHLAN/, --bowtie2 $BOWTIE2/, --diamond $DIAMOND/ (or --rapsearch $RAPSEARCH, --usearch $USEARCH/)). By default diamond is run for the translated alignment but this can be changed to rapsearch2 or usearch by setting "--translated-alignment {rapsearch|usearch}".

### Other ###
1. Memory (>= 10 Gb)
1. Disk space (>= 10 Gb)
1. Operating system (Linux or Mac)

## Installation ##

### Downloading HUMAnN2 ###
You can download the latest HUMAnN2 release or the development version.

Option 1: Latest Release (Recommended)

* [Download](https://bitbucket.org/biobakery/humann2/get/0.1.2.tar.gz) and unpack the latest release of HUMAnN2.

Option 2: Development Version

* Create a clone of the repository: 
	
	``hg clone https://bitbucket.org/biobakery/humann2 ``

	Note: Creating a clone of the repository requires [Mercurial](http://mercurial.selenic.com/) to be installed. Once the repository has been cloned upgrading to the latest development version of HUMAnN2 is simple. Just type ``hg pull -u`` from within the repository.

For the steps that follow, $HUMAnN2_PATH is the location of the HUMAnN2 directory (ie $HUMAnN2_PATH=/home/user/humann2/ with this readme file found in this folder).

### Installing HUMAnN2 ###

1. Move to the HUMAnN2 directory : ``$ cd $HUMAnN2_PATH ``
1. Install MinPath : ``$ python setup.py minpath ``
1. Install HUMAnN2 : ``$ python setup.py install ``

Note: If you do not have write permissions to '/usr/lib/', then add the option "--user" to the install command. This will install the python package into subdirectories of '~/.local'.

### Test the install ###

Test out the install of HUMAnN2 by running the unit tests.

``$ python setup.py test ``

### Try out a demo run ###

With HUMAnN2 installed you can try out a demo run using reduced versions of the databases.

``$ humann2 --input examples/demo.fastq --output $OUTPUT_DIR ``

Output from this demo run will be written to the folder $OUTPUT_DIR.

Please continue with the install directions to download the full databases before running with your sequencing data.


### Download the databases ###

Downloading the databases is a required step if your input is a filtered shotgun sequencing metagenome file (fastq, fastq.gz, fasta, or fasta.gz format). If your input files will always be mapping results files (sam, bam or blastm8 format) or gene tables (tsv or biom format), you do not need to download the ChocoPhlAn and UniRef50 databases. 

#### Download the ChocoPhlAn database ####

Download the ChocoPhlAn database providing $INSTALL_LOCATION as the location to install the database.

```
$ humann2_databases --download chocophlan full $INSTALL_LOCATION
```

NOTE: The humann2 config file will be updated to point to this location for the default chocophlan database. If you move this database, please use the "humann2_config" command to update the default location of this database. Alternatively you can always provide the location of the chocophlan database you would like to use with the "--chocophlan <chocoplan>" option to humann2.


#### Download the UniRef50 database ####

Download the UniRef50 database providing $INSTALL_LOCATION as the location to install the database.

```
$ humann2_databases --download uniref diamond $INSTALL_LOCATION
```

NOTE: The humann2 config file will be updated to point to this location for the default uniref database. If you move this database, please use the "humann2_config" command to update the default location of this database. Alternatively you can always provide the location of the uniref database you would like to use with the "--uniref <uniref>" option to humann2.

NOTE: By default HUMAnN2 runs diamond for translated alignment. If you would like to use rapsearch2 for translated alignment, first download the rapsearch2 formatted database by running this command with the rapsearch2 formatted database selected. It is suggested that you install both databases in the same folder so this folder can be the default uniref database location. This will allow you to switch between alignment software without having to specify a different location for the database.


## How to Run ##

### Basic usage ###

Type the command:

`` humann2 --input $SAMPLE --output $OUTPUT_DIR``

$SAMPLE = a single file that is one of the following types:

1. filtered shotgun sequencing metagenome file (fastq, fastq.gz, fasta, or fasta.gz format)
2. mapping results file (sam, bam or blastm8 format)
3. gene table file (tsv or biom format)

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


To run the fasta demo type the command:

`` humann2 --input examples/demo.fasta --output $OUTPUT_DIR``

To run the fastq demo type the command:

`` humann2 --input examples/demo.fastq --output $OUTPUT_DIR``

To run the sam demo type the command:

`` humann2 --input examples/demo.sam --output $OUTPUT_DIR``

To run the blastm8 demo type the command:

`` humann2 --input examples/demo.m8 --output $OUTPUT_DIR``

$OUTPUT_DIR is the output directory

Since sam and blastm8 are mapping results, using these files as input to HUMAnN2 will bypass both the nucleotide and translated mapping portions of the flow.

## Output files ##

HUMAnN2 produces three output files which by default are tab-delimited text. There is an option to print out files in biom format. 

### Gene Families ###

```
# Gene Family	$SAMPLENAME
UniRef50_A6L0N6	67.0
UniRef50_A6L0N6|s__Bacteroides_fragilis	8.0
UniRef50_A6L0N6|s__Bacteroides_finegoldii	5.0
UniRef50_A6L0N6|s__Bacteroides_stercoris	4.0
UniRef50_A6L0N6|unclassified	1.0
UniRef50_G9S1V7	60.0
UniRef50_G9S1V7|s__Bacteroides_vulgatus	31.0
UniRef50_G9S1V7|s__Bacteroides_thetaiotaomicron	22.0
UniRef50_G9S1V7|s__Bacteroides_stercoris	7.0
```

* This file includes the abundance of each orthologous gene family in the community organized by bug. Orthologous families are groups of genes that perform roughly the same biological roles. 
* HUMAnN2 uses the MetaPhlAn2 software along with the ChocoPhlAn database and UniRef for this computation.

### Pathway Coverage ###

```
# Pathway	$SAMPLENAME
PWY0-1301	1.0
PWY0-1301|s__Bacteroides_caccae	1.0
PWY0-1301|s__Bacteroides_finegoldii	1.0
PWY0-1301|unclassified	1.0
PWY-7134	1.0
PWY-7134|s__Bacteroides_vulgatus	0.666666666667
PWY-7134|s__Bacteroides_thetaiotaomicron	0.666666666667
PWY-7134|unclassified	0.333333333333
PWY-7134|s__Parabacteroides_merdae	0.333333333333
```

* This file includes the presence/absence of each pathway in the community grouped by bug. HUMAnN refers to pathway presence/absence as "coverage" and defines a pathway as a set of two or more genes. 
* HUMAnN2 uses MetaCyc pathways along with MinPath for this computation. 
* The user has the option to provide a custom pathways database to HUMAnN2 and to use all pathways instead of the minimal pathways computed by MinPath.

### Pathway Abundance ###

```
# Pathway	$SAMPLENAME
PWY-1921	57.0136768635
PWY-1921|unclassified	32.2636768635
PWY-1921|s__Bacteroides_ovatus	4.5
PWY-1921|s__Alistipes_putredinis	3.0
PWY-1921|s__Bacteroides_caccae	2.25
PWY0-1301	54.9996450867
PWY0-1301|unclassified	16.9996450867
PWY0-1301|s__Parabacteroides_merdae	8.0
PWY0-1301|s__Bacteroides_caccae	6.0
```

* This file includes the abundance of each pathway in the community grouped by bug. This is the total number of “copies” of the pathways present. 
* HUMAnN2 uses MetaCyc pathways along with MinPath for this computation. 
* The user has the option to provide a custom pathways database to HUMAnN2 and to use all pathways instead of the minimal pathways computed by MinPath.

## Complete option list ##
```
usage: humann2 [-h] [--version] [-v] [-r] [--bypass-prescreen]
               [--bypass-nucleotide-index] [--bypass-translated-search]
               [--bypass-nucleotide-search] -i <input.fastq> -o <output>
               [-c <chocophlan>] [--chocophlan-gene-index <-1>] [-u <uniref>]
               [--average-read-length <1>] [--evalue <1.0>]
               [--metaphlan <metaplhan>] [--o-log <sample.log>]
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
  --metaphlan <metaplhan>
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
```