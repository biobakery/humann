[TOC]

# HUMAnN2: HMP Unified Metabolic Analysis Network 2 #

## Description ##
HUMAnN is a pipeline for efficiently and accurately determining the presence/absence and abundance of microbial pathways in a community from metagenomic or metatranscriptomic sequencing data (typically millions of short DNA/RNA reads). 

HUMAnN2 is the next generation of HUMAnN. HUMAnN2 incorporates several new features including an expanded database of microbial genomes (>4x the size of the database included in HUMAnN), a simple user interface (single command driven flow), and bug-specific output files. 

The HUMAnN2 pipeline is a single command driven flow requiring the user to only provide a filtered fastq file to produce gene and pathway summary files ready for analysis. The pipeline converts sequence reads into coverage and abundance tables summarizing the gene families and pathways in one or more microbial communities. 

## Requirements ##

### Software ###

1. [MetaPhlAn2](https://bitbucket.org/biobakery/metaphlan2/)
1. [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) (version >= 2.1)
1. [rapsearch2](http://omics.informatics.indiana.edu/mg/RAPSearch2/)
1. [Python](http://www.python.org/) (version >= 2.7)
1. [MinPath](http://omics.informatics.indiana.edu/MinPath/) (automatically downloaded/installed)
1. [Xipe](https://edwards.sdsu.edu/cgi-bin/xipe.cgi) (optional / included)
1. [usearch](http://www.drive5.com/usearch/) (version = v7.0.1001) (optional)

If MetaPhlAn2, bowtie2, and rapsearch2 (or usearch) are not installed in a location in your $PATH,
then add them to your $PATH or use the HUMAnN2 parameters to indicate the locations of their
directories (--metaphlan $METAPHLAN/, --bowtie2 $BOWTIE2/, --rapsearch $RAPSEARCH/ (or --usearch $USEARCH/)). By default rapsearch is run for the translated alignment but this can be changed to usearch by setting "--translated_alignment usearch".

### Other ###
1. Memory (>= 10 Gb)
1. Disk space (>= 40 Gb)
1. Operating system (Linux or Mac)

## Installation ##

### Downloading HUMAnN2 ###
HUMAnN2 can be downloaded in two ways:

* [Download](https://bitbucket.org/biobakery/humann2/downloads) a compressed set of files.
* Create a clone of the repository on your computer with the command: 
	
	``hg clone https://bitbucket.org/biobakery/humann2 ``

Note: Creating a clone of the repository requires [Mercurial](http://mercurial.selenic.com/) to be installed. Once the repository has been cloned upgrading to the latest release of HUMAnN2 is simple. Just type ``hg pull -u`` from within the repository which will download the latest release.

For the steps that follow, $HUMANn2_PATH is the location that HUMAnN2 was download (ie $HUMAnN2_PATH=/home/user/humann2/ with the file "humann2.py" found in this folder).


### Downloading the databases ###

#### Downloading the [ChocoPhlAn database](http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/chocophlan.tar.gz) ####

```
$ cd $HUMANn2_PATH/databases
$ wget http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/chocophlan.tar.gz
$ tar zxvf chocophlan.tar.gz 
$ rm chocophlan.tar.gz
```

$HUMANn2_PATH = the full path to the HUMAnN2 download

NOTE: These steps download the ChocoPhlAn database to the humann2/databases folder. This folder is the default location that HUMAnN2 will look for this database. If you place it in another location (ie $DIR), provide this to HUMAnN2 with the "--chocophlan $DIR" option.



#### Downloading the [UniRef50 database](http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref50_rapsearch/uniref50_rapsearch.tar.gz) ####

```
$ cd $HUMANn2_PATH/databases
$ wget http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref50_rapsearch/uniref50_rapsearch.tar.gz
$ tar zxvf uniref50_rapsearch.tar.gz
$ rm uniref50_rapsearch.tar.gz
```

$HUMANn2_PATH = the full path to the HUMAnN2 download


NOTE: These steps download the UniRef50 database formatted for rapsearch2 to the humann2/databases folder. This folder is the default location that HUMAnN2 will look for this database. If you place it in another location (ie $DIR), provide this to HUMAnN2 with the "--uniref $DIR" option.

NOTE: By default HUMAnN2 runs rapsearch2 for translated alignment. If usearch is selected for translated alignment, provide a database that has been formatted for usearch.

### Updating the environment ###
To update the environment, add the path to the HUMAnN2 download directory ($HUMAnN2_PATH) to your $PATH.

1. Add this line to the .bashrc file located in your home directory ($HOME) : `` export PATH=$PATH:$HUMAnN2_PATH ``
1.  Run this command to update your current environment: ``$ source $HOME/.bashrc ``
1. HUMAnN2 can now be run without providing the location of the install directory: `` $ humann2.py --help ``


The update environment step is optional if the path to the HUMAnN2 executable is always provided. 
For example, calling HUMAnN2 as follows does not require updates to the environment: `` $ $HUMAnN2_PATH/humann2.py --help ``


## How to Run ##

### Basic usage ###

Type the command:

`` humann2.py --input $SAMPLE --output $OUTPUT_DIR``

$SAMPLE = a single file that is one of the following types:

1. filtered shotgun sequencing metagenome file (fastq, fastq.gz, fasta, or fasta.gz format)
1. mapping results file (sam, bam or blastm8 format)
1. gene table file (tsv or biom format)

$OUTPUT_DIR = the output directory

**Three output files will be created:**

1. $OUTPUT_DIR/$SAMPLENAME_genefamilies.tsv*
1. $OUTPUT_DIR/$SAMPLENAME_pathcoverage.tsv
1. $OUTPUT_DIR/$SAMPLENAME_pathabundance.tsv

where $SAMPLENAME is the basename of $SAMPLE

*The gene families file will not be created if the input file type is a gene table.

**Intermediate temp files will also be created:**

1. $DIR/$SAMPLENAME_bowtie2_aligned.sam
	* the full alignment output from bowtie2 
1. $DIR/$SAMPLENAME_bowtie2_aligned.tsv
	* only the reduced aligned data from the bowtie2 output
1. $DIR/$SAMPLENAME_bowtie2_index*
	* bowtie2 index files created from the custom chochophlan database
1. $DIR/$SAMPLENAME_bowtie2_unaligned.fa 
	* a fasta file of unaligned reads after the bowtie2 step
1. $DIR/$SAMPLENAME_custom_chocophlan_database.ffn 
	* a custom chocophlan database of fasta sequences
1. $DIR/$SAMPLENAME_metaphlan_bowtie2.txt 
	* the bowtie2 output from metaphlan
1. $DIR/$SAMPLENAME_metaphlan_bugs_list.tsv 
	* the bugs list output from metaphlan
1. $DIR/$SAMPLENAME_$TRANSLATEDALIGN_aligned.tsv 
	* the alignment results from the translated alignment step
1. $DIR/$SAMPLENAME_$TRANSLATEDALIGN_unaligned.fa 
	* a fasta file of unaligned reads after the translated alignment step
1. $DIR/$SAMPLENAME.log 
	* a log of the run
	
* $DIR=$OUTPUT_DIR/$SAMPLENAME_humann2_temp/
* $SAMPLENAME is the basename of the fastq/fasta input file
* $TRANSLATEDALIGN is the translated alignment software selected (rapsearch2 or usearch)

NOTE: $SAMPLENAME can be set by the user with the option "--output_basename <$NEWNAME>". 

### Demo runs ###

The input folder contains four demo input files. These files are of fasta, fastq, sam, and blastm8 format. Blastm8 format is created by the following software: rapsearch2, usearch, and blast.


To run the fasta demo type the command:

`` humann2.py --input input/demo.fasta --output $OUTPUT_DIR``

To run the fastq demo type the command:

`` humann2.py --input input/demo.fastq --output $OUTPUT_DIR``

To run the sam demo type the command:

`` humann2.py --input input/demo.sam --output $OUTPUT_DIR``

To run the blastm8 demo type the command:

`` humann2.py --input input/demo.m8 --output $OUTPUT_DIR``

$OUTPUT_DIR is the output directory

Since sam and blastm8 are mapping results, using these files as input to HUMAnN2 will bypass both the nucleotide and translated mapping portions of the flow.

### Output files ###

HUMAnN2 produces three output files which by default are tab-delimited text. There is an option to print out files in biom format. 

#### Gene Families ####

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

#### Pathway Coverage ####

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

#### Pathway Abundance ####

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

### Additional ways to run ####

1. To run with additional output printed to stdout: add the ``--verbose`` flag
1. To run using multiple cores: add the ``--threads $CORES`` option
1. To remove the intermediate temp output files: add the ``--remove_temp_output`` flag
1. To bypass the MetaPhlAn prescreen step: add the ``--bypass_prescreen`` flag
1. To bypass the prescreen and the nucleotide alignment index step and start with the bowtie2 alignment step: add the ``--bypass_nucleotide_index`` flag
	* If using this flag provide the bowtie2 index as the input to the chocophlan parameter instead of the chocoplan directory. For example, run with "--bypass_nucleotide_index --chocophlan chocophlan_dir/chocophlan_bowtie2_index"  if the first bowtie2 index file is located at chocophlan_dir/chocophlan_bowtie2_index.1.bt2 .
1. To provide the location of the ChocoPhlAn database: add the ``--chocophlan $DIR `` option
1. To provide the location of the UniRef database: add the ``--uniref $DIR`` option

### Complete option list ###
```
usage: humann2.py [-h] [-v] [-r] [--bypass_prescreen]
                  [--bypass_nucleotide_index] [--bypass_translated_search]
                  [--bypass_nucleotide_search] -i <input.fastq> -o <output>
                  [-c <chocophlan>] [--chocophlan_gene_index <-1>]
                  [-u <uniref>] [--average_read_length <1>]
                  [--metaphlan <metaplhan>] [--o_log <sample.log>]
                  [--log_level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
                  [--remove_temp_output] [--bowtie2 <bowtie2>] [--threads <1>]
                  [--prescreen_threshold <0.01>] [--identity_threshold <0.4>]
                  [--usearch <usearch>] [--rapsearch <rapsearch>]
                  [--metaphlan_output <bugs_list.tsv>]
                  [--id_mapping <id_mapping.tsv>]
                  [--translated_alignment {usearch,rapsearch}]
                  [--xipe {on,off}] [--minpath {on,off}]
                  [--pick_frames {on,off}] [--output_format {tsv,biom}]
                  [--output_basename <sample_name>]
                  [--remove_stratified_output]
                  [--input_format {fastq,fastq.gz,fasta,fasta.gz,sam,bam,blastm8,genetable,biom}]
                  [--pathways_databases <pathways_database_part1.tsv> <pathways_database_part2.tsv>]

HUMAnN2 : HMP Unified Metabolic Analysis Network 2

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         additional output is printed
  -r, --resume          bypass commands if the output files exist
  --bypass_prescreen    bypass the prescreen step and run on the full ChocoPhlAn database
  --bypass_nucleotide_index
                        bypass the nucleotide index step and run on the indexed ChocoPhlAn database
  --bypass_translated_search
                        bypass the translated search step
  --bypass_nucleotide_search
                        bypass the nucleotide search steps
  -i <input.fastq>, --input <input.fastq>
                        input file of type {fastq,fastq.gz,fasta,fasta.gz,sam,bam,blastm8,genetable,biom} 
                        [REQUIRED]
  -o <output>, --output <output>
                        directory to write output files
                        [REQUIRED]
  -c <chocophlan>, --chocophlan <chocophlan>
                        directory containing the ChocoPhlAn database
                        [DEFAULT: databases/chocophlan/]
  --chocophlan_gene_index <-1>
                        the index of the gene in the sequence annotation
                        [DEFAULT: -1]
  -u <uniref>, --uniref <uniref>
                        directory containing the UniRef database
                        [DEFAULT: databases/uniref/]
  --average_read_length <1>
                        the average length of the reads
                        [DEFAULT: 1]
  --metaphlan <metaplhan>
                        directory containing the MetaPhlAn software
                        [DEFAULT: $PATH]
  --o_log <sample.log>  log file
                        [DEFAULT: temp/sample.log]
  --log_level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        level of messages to display in log
                        [DEFAULT: DEBUG]
  --remove_temp_output  remove temp output files
                        [DEFAULT: temp files are not removed]
  --bowtie2 <bowtie2>   directory of the bowtie2 executable
                        [DEFAULT: $PATH]
  --threads <1>         number of threads/processes
                        [DEFAULT: 1]
  --prescreen_threshold <0.01>
                        minimum percentage of reads matching a species
                        [DEFAULT: 0.01]
  --identity_threshold <0.4>
                        identity threshold to use with the translated search
                        [DEFAULT: 0.4]
  --usearch <usearch>   directory containing the usearch executable
                        [DEFAULT: $PATH]
  --rapsearch <rapsearch>
                        directory containing the rapsearch executable
                        [DEFAULT: $PATH]
  --metaphlan_output <bugs_list.tsv>
                        output file created by metaphlan
                        [DEFAULT: file will be created]
  --id_mapping <id_mapping.tsv>
                        id mapping file for alignments
                        [DEFAULT: alignment reference used]
  --translated_alignment {usearch,rapsearch}
                        software to use for translated alignment
                        [DEFAULT: rapsearch]
  --xipe {on,off}       turn on/off the xipe computation
                        [DEFAULT: off]
  --minpath {on,off}    turn on/off the minpath computation
                        [DEFAULT: on]
  --pick_frames {on,off}
                        turn on/off the pick_frames computation
                        [DEFAULT: off]
  --output_format {tsv,biom}
                        the format of the output files
                        [DEFAULT: tsv]
  --output_basename <sample_name>
                        the basename for the output files
                        [DEFAULT: input file basename]
  --remove_stratified_output
                        remove stratification from output
                        [DEFAULT: output is stratified]
  --input_format {fastq,fastq.gz,fasta,fasta.gz,sam,bam,blastm8,genetable,biom}
                        the format of the input file
                        [DEFAULT: format identified by software]
  --pathways_databases <pathways_database_part1.tsv> <pathways_database_part2.tsv>
                        the two mapping files to use for pathway computations
                        [DEFAULT: databases/pathways/metacyc_reactions.uniref , databases/pathways/metacyc_pathways]
```