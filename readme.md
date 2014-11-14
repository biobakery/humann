[TOC]

# HUMAnN2: HMP Unified Metabolic Analysis Network 2 #

## Description ##
HUMAnN is a pipeline for efficiently and accurately determining the presence/absence and abundance of microbial pathways in a community from metagenomic sequencing data (typically millions of short DNA/RNA reads). 

HUMAnN2 incorporates multiple upgrades including an expanded database of microbial genomes (>4x the size of the database included in HUMAnN), a simple user interface (single command driven flow), and bug-specific output files. 

HUMAnN2 produces three output files:

1. Gene Families
	* This file includes the abundance of each orthologous gene family in the community organized by bug. Orthologous families are groups of genes that perform roughly the same biological roles. 
	* HUMAnN2 uses the MetaPhlAn2 software along with the ChocoPhlAn database and UniRef for this computation.

2. Pathway Coverage
	* This file includes the presence/absence of each pathway in the community grouped by bug. HUMAnN refers to pathway presence/absence as "coverage" and defines a pathway as a set of two or more genes. 
	* HUMAnN2 uses MetaCyc pathways along with MinPath for this computation. 
	* The user has the option to provide a custom pathways database to HUMAnN2 and to use all pathways instead of the minimal pathways computed by MinPath.

3. Pathway Abundance
	* This file includes the abundance of each pathway in the community grouped by bug. This is the total number of “copies” of the pathways present. 
	 * HUMAnN2 uses MetaCyc pathways along with MinPath for this computation. 
	 * The user has the option to provide a custom pathways database to HUMAnN2 and to use all pathways instead of the minimal pathways computed by MinPath.

The HUMAnN2 pipeline is a single command driven flow requiring the user to only provide a filtered fastq file to produce gene and pathway summary files ready for analysis. The pipeline converts sequence reads into coverage and abundance tables summarizing the gene families and pathways in one or more microbial communities. 

## Requirements ##

### Software ###

1. [MetaPhlAn](https://bitbucket.org/biobakery/metaphlan2/)
1. [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) (version >= 2.1)
1. [rapsearch2](http://omics.informatics.indiana.edu/mg/RAPSearch2/)
1. [Python](http://www.python.org/) (version >= 2.7)
1. [MinPath](http://omics.informatics.indiana.edu/MinPath/) (automatically downloaded)
1. [Xipe](https://edwards.sdsu.edu/cgi-bin/xipe.cgi) (optional / included)
1. [usearch](http://www.drive5.com/usearch/) (version = v7.0.1001) (optional)

If MetaPhlAn, bowtie2, and rapsearch (or usearch) are not installed in a location in your $PATH,
then add them to your $PATH or use the HUMAnN2 parameters to indicate the locations of their
directories (--metaphlan $METAPHLAN/, --bowtie2 $BOWTIE2/, --rapsearch $RAPSEARCH/ (or --usearch $USEARCH/)). By default rapsearch is run for the translated alignment but this can be changed to usearch by setting "--translated_alignment usearch".

### Other ###
1. Memory (>= 10 Gb)
1. Disk space (>= 40 Gb)
1. Operating system (Linux or Mac)

## Installation ##

### For the impatient ###

```
$ hg clone https://bitbucket.org/biobakery/humann2
$ cd humann2
$ python setup.py install # or 'develop'
$ python setup.py download
```

### Downloading HUMAnN2 ###
HUMAnN2 can be downloaded in two ways:

* [Download](https://bitbucket.org/biobakery/humann2/downloads) a compressed set of files.
* Create a clone of the repository on your computer with the command: 
	
	``hg clone https://bitbucket.org/biobakery/humann2 ``

Note: Creating a clone of the repository requires [Mercurial](http://mercurial.selenic.com/) to be installed. Once the repository has been cloned upgrading to the latest release of HUMAnN2 is simple. Just type ``hg -u pull`` from within the repository which will download the latest release.

For the steps that follow, $HUMANn2_PATH is the location that HUMAnN2 was download (ie $HUMAnN2_PATH=/home/user/humann2/ with the directory "humann2lib" and file "readme.md" found in this folder).


### Installing HUMAnN2 scripts and libraries ###
Install HUMAnN2 using the following command:

```
$ cd $HUMANn2_PATH/
$ python setup.py install
```

### Installing HUMAnN2 scripts and libraries in editable mode ###
To install the HUMAnN2 scripts such that edits made to the files in `$HUMAnN2_PATH` take effect immediately, use the `develop` command as follows:

```
$ cd $HUMANn2_PATH/
$ python setup.py develop
```


### Downloading the databases ###

#### Downloading the [ChocoPhlAn](http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/) and [UniRef50](http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref50_rapsearch/) databases####

```
$ cd $HUMANn2_PATH/
$ python setup.py download
```

$HUMANn2_PATH = the full path to the HUMAnN2 download

NOTE: By default, this step download the databases to the $HUMAnN2_PATH/databases directory. To download to a different directory use the `--to` option. For example `python setup.py download --to /opt/humann2` will download the databses to the `/opt/humann2` directory. Remember to provide this directory to HUMAnN2 with the "--chocophlan $DIR" option.

NOTE: By default HUMAnN2 runs rapsearch2 for translated alignment. Thus, the UniRef50 database downloaded in the step above will be formatted for rapsearch2. If usearch is selected for translated alignment, you will need to provide a database that has been formatted for usearch yourself.

## How to Run ##

### Basic usage ###

Type the command:

`` humann2.py --input $SAMPLE --output $OUTPUT_DIR``

where $SAMPLE is your filtered fasta or fastq file (ie metagenome.fastq) and $OUTPUT_DIR is the output directory

Three output files will be created:

1. $OUTPUT_DIR/$SAMPLENAME_genefamilies.tsv
1. $OUTPUT_DIR/$SAMPLENAME_pathcoverage.tsv
1. $OUTPUT_DIR/$SAMPLENAME_pathabundance.tsv

where $SAMPLENAME is the basename of $SAMPLE

NOTE: To keep all of the intermediate temp files use the "--temp" flag.

### Additional ways to run ####

1. To run with additional output printed to stdout: add the ``--verbose`` flag
1. To run using multiple cores: add the ``--threads $CORES`` option
1. To keep the intermediate output files: add the ``--temp`` flag
	* The intermediate files are:
		1. $DIR/$SAMPLENAME_bowtie2_aligned.sam
			* the full alignment output from bowtie2 
		1. $DIR/$SAMPLENAME_bowtie2_aligned.tsv
			* only the aligned data from the bowtie2 output
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
	* $DIR=$OUTPUT_DIR/$SAMPLENAME_HUMAnN2_temp/
	* $SAMPLENAME is the basename of the fastq/fasta input file
	* $TRANSLATEDALIGN is the translated alignment software selected (rapsearch2 or usearch)

1. To bypass the MetaPhlAn prescreen step: add the ``--bypass_prescreen`` flag
1. To bypass the prescreen and the nucleotide alignment index step and start with the bowtie2 alignment step: add the ``--bypass_nucleotide_index`` flag
	* If using this flag provide the bowtie2 index as the input to the chocophlan parameter instead of the chocoplan directory. For example, run with "--bypass_nucleotide_index --chocophlan chocophlan_dir/chocophlan_bowtie2_index"  if the first bowtie2 index file is located at chocophlan_dir/chocophlan_bowtie2_index.1.bt2 .
1. To provide the location of the ChocoPhlAn database: add the ``--chocophlan $DIR `` option
1. To provide the location of the UniRef database: add the ``--uniref $DIR`` option

### Complete option list ###
```
usage: humann2.py [-h] [-v] [-r] [--bypass_prescreen]
                  [--bypass_nucleotide_index] -i <input.fastq> -o <output>
                  [-c <chocophlan>] [-u <uniref>] [--metaphlan <metaplhan>]
                  [--o_log <sample.log>]
                  [--log_level {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--temp]
                  [--bowtie2 <bowtie2>] [--threads <1>]
                  [--prescreen_threshold <0.01>] [--identity_threshold <0.4>]
                  [--usearch <usearch>] [--rapsearch <rapsearch>]
                  [--metaphlan_output <bugs_list.tsv>]
                  [--translated_alignment {usearch,rapsearch}]
                  [--xipe {on,off}] [--minpath {on,off}]
                  [--output_format {tsv,biom}]
                  [--pathways_databases <pathways_database_part1.tsv> <pathways_database_part2.tsv>]

HUMAnN2 : HMP Unified Metabolic Analysis Network 2

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         additional output is printed
  -r, --resume          bypass commands if the output files exist
  --bypass_prescreen    bypass the prescreen step and run on the full ChocoPhlAn database
  --bypass_nucleotide_index
                        bypass the nucleotide index step and run on the indexed ChocoPhlAn database
  -i <input.fastq>, --input <input.fastq>
                        fastq/fasta input file
                        [REQUIRED]
  -o <output>, --output <output>
                        directory to write output files
                        [REQUIRED]
  -c <chocophlan>, --chocophlan <chocophlan>
                        directory containing the ChocoPhlAn database
                        [DEFAULT: databases/chocophlan/ ]
  -u <uniref>, --uniref <uniref>
                        directory containing the UniRef database
                        [DEFAULT: databases/uniref/ ]
  --metaphlan <metaplhan>
                        directory containing the MetaPhlAn software
                        [DEFAULT: $PATH]
  --o_log <sample.log>  log file
                        [DEFAULT: temp/sample.log]
  --log_level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        level of messages to display in log
                        [DEFAULT: DEBUG ]
  --temp                keep temp output files
                        [DEFAULT: temp files are removed]
  --bowtie2 <bowtie2>   directory of the bowtie2 executable
                        [DEFAULT: $PATH]
  --threads <1>         number of threads/processes
                        [DEFAULT: 1 ]
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
  --translated_alignment {usearch,rapsearch}
                        software to use for translated alignment
                        [DEFAULT: rapsearch]
  --xipe {on,off}       turn on/off the xipe computation
                        [DEFAULT: off ]
  --minpath {on,off}    turn on/off the minpath computation
                        [DEFAULT: on ]
  --output_format {tsv,biom}
                        the format of the output files
                        [DEFAULT: tsv ]
  --pathways_databases <pathways_database_part1.tsv> <pathways_database_part2.tsv>
                        the two mapping files to use for pathway computations
                        [DEFAULT: databases/pathways/metacyc_reactions.uniref , databases/pathways/metacyc_pathways ]
```