## HUMAnN2 User Manual 

HUMAnN2 is the next generation of HUMAnN (HMP Unified Metabolic Analysis Network).

**If you use the HUMAnN2 software, please cite our manuscript: TBD**

HUMAnN is a pipeline for efficiently and accurately profiling the presence/absence and abundance of microbial pathways in a community from metagenomic or metatranscriptomic sequencing data (typically millions of short DNA/RNA reads). This process, referred to as functional profiling, aims to describe the metabolic potential of a microbial community and its members. More generally, functional profiling answers the question "What are the microbes in my community-of-interest doing (or capable of doing)?"

**Table of Contents**

[TOC]

### Requirements

1.  [MetaPhlAn2](http://huttenhower.sph.harvard.edu/metaphlan2)
2.  [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (version >= 2.1)
3.  [Diamond](http://ab.inf.uni-tuebingen.de/software/diamond/)
4.  [Python](http://www.python.org/) (version >= 2.7)
5.  Memory (>= 10 GB)
6.  Disk space (>= 10 GB [to accommodate comprehensive sequence databases])
7.  Operating system (Linux or Mac)

If always running with files of type #2, #3, and #4 (for information on file types, see section *Workflow by input file type*), the requirements are reduced. MetaPhlAn2, Bowtie2, Diamond, and the amount of disk space listed are not required. Also if you always run with one or more bypass options (for information on bypass options, see section *Workflow by bypass mode*), the requirements might also be reduced.

Please note there are additional requirements if you are using input files of type sam or biom. The [SAMtools](http://samtools.sourceforge.net/) software is required for bam files and the [biom-format](http://biom-format.org/) software is required for biom files.

### Installation

1. Download and unpack the latest release of the [HUMAnN2 software](https://bitbucket.org/biobakery/humann2/get/0.1.tar.gz)
2. Install [MinPath](http://omics.informatics.indiana.edu/MinPath/) and the HUMAnN2 software (see notes 1 and 2)
 
    ```$ python setup.py minpath
	   $ python setup.py install
	```
    
3. Test the HUMAnN2 install (Optional)
 
     `` $ python setup.py test``

4. Try out a HUMAnN2 demo run (Optional)

    `` $ humann2 --input humann2/examples/demo.fastq --output $OUTPUT_DIR ``

5. Download the ChocoPhlAn database with $INSTALL_LOCATION = the location you have selected to install the database (approx. size = 5.6 GB) (see note 3)

    ``$ humann2_databases --download chocophlan full $INSTALL_LOCATION``
    
6. Download the UniRef database with $INSTALL_LOCATION = the location you have selected to install the database (approx. size = 2.8 GB) (see note 3)

    ``$ humann2_databases --download uniref diamond $INSTALL_LOCATION``

NOTE 1: If you would like to update the glpk required by MinPath, add the option "--update-glpk" to the MinPath install command.

NOTE 2: If you do not have write permissions to '/usr/lib/', then add the option "--user" to the HUMAnN2 install command. This will install the python package into subdirectories of '~/.local'.

NOTE 3: Downloading and installing the ChocoPhlAn and UniRef databases is not required if always running with files of type #2, #3, and #4 (for information on file types, see section *Workflow by input file type*). It is also not required to run the demo which runs on demo versions of the ChocoPhlAn and UniRef databases included as part of the HUMAnN2 install.

### How to run

To run HUMAnN2:
```
$ humann2 --input $SAMPLE --output $OUTPUT_DIR
```


$SAMPLE = a single file that is one of the following types:

1.  filtered shotgun sequencing metagenome file (fastq, fastq.gz, fasta, or fasta.gz format)
2.  alignment file (sam, bam or blastm8 format)
3.  gene table file (tsv or biom format)

$OUTPUT_DIR = the output directory

### Output files

When HUMAnN2 is run, three main output files will be created (where `` $SAMPLENAME = the basename of $SAMPLE ``):

#### Gene families file

``` 
# Gene Family	$SAMPLE
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

*   File name: `` $OUTPUT_DIR/$SAMPLENAME_genefamilies.tsv ``
*   This file quantifies the abundance of each gene family in the community. Gene families are groups of evolutionarily-related protein-coding sequences that typically perform similar functions. Abundance is reported in RPK (reads per kilobase) units to normalize for gene length; RPK units reflect relative gene (or transcript) copy number in the community.
*   In addition to community-wide gene family abundance totals (as reported by HUMAnN), this file is stratified to indicate abundance contributions of known and unclassified organisms represented in the sample.
* Please note the gene families file will not be created if the input file type is a gene table.
        
#### Pathway coverage file

``` 
# Pathway	$SAMPLE
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

*   File name: `` $OUTPUT_DIR/$SAMPLENAME_pathcoverage.tsv ``
*   This file details the presence/absence of each pathway in the community. HUMAnN refers to pathway presence/absence as "coverage" and defines a pathway as a set of two or more gene families.
*   In addition to community-wide pathway coverage (as reported by HUMAnN), this file is stratified to indicate the coverage of the pathway by genomes of known and unclassified organisms represented in the sample.

#### Pathway abundance file

```
# Pathway	$SAMPLE
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
         
*   File name: `` $OUTPUT_DIR/$SAMPLENAME_pathabundance.tsv ``
*   This file quantifies the abundance of each pathway in the community as a function of the abundance of its member gene families.
*   In addition to community-wide pathway abundance (as reported by HUMAnN), this file is stratified to indicate abundance contributions of known and unclassified organisms represented in the sample.

#### Intermediate temp output files

Ten intermediate temp output files will be created where:

$DIR = $OUTPUT_DIR/$SAMPLENAME_humann2_temp/
$SAMPLENAME = basename of the fastq/fasta input file
$TRANSLATEDALIGN = translated alignment software selected (diamond, rapsearch2 or usearch)

NOTE: $SAMPLENAME can be set by the user with the option --output-basename <$NEWNAME>

##### Bowtie2 alignment results

> @HD	VN:1.0	SO:unsorted
> @SQ	SN:g__Ruminococcus.s__Ruminococcus_bromii|UniRef90_D4L6K4|UniRef50_R6U703	LN:540
> r99491	0	g__Bacteroides.s__Bacteroides_stercoris|UniRef90_R6B629|UniRef50_R5RCC8	1015	42	151M	*	0	0	$SEQ	$QUAL
> r99526	0	g__Parabacteroides.s__Parabacteroides_merdae|UniRef90_unknown|UniRef50_D9RX34	155	42	151M	*	0	0	$SEQ	$QUAL
> r99581	16	g__Bacteroides.s__Bacteroides_stercoris|UniRef90_unknown|UniRef50_R6SXR7	2503	42	151M	*	0	0	$SEQ	$QUAL

*   File name: `` $DIR/$SAMPLENAME_bowtie2_aligned.sam `` 
*   This file has the full alignment output from bowtie2.
*   In example above `` $SEQ = sequence and $QUAL = quality scores `` to fit in page.

##### Bowtie2 reduced alignment results

``` 
r93	g__Bacteroides.s__Bacteroides_cellulosilyticus|UniRef90_E2NEW2|UniRef50_E2NEW2	6.3095734448e-05
r113	g__Bacteroides.s__Bacteroides_cellulosilyticus|UniRef90_R6KNZ3|UniRef50_R6KNZ3	6.3095734448e-05	
r704	g__Bacteroides.s__Bacteroides_uniformis|UniRef90_unknown|UniRef50_E6STE9		0.794328234724	
r663	g__Bacteroides.s__Bacteroides_thetaiotaomicron|UniRef90_R7KKH7|UniRef50_R7KKH7		6.3095734448e-05	
r940	g__Ruminococcus.s__Ruminococcus_bromii|UniRef90_unknown|UniRef50_unknown		6.3095734448e-0	 
```

*   File name: `` $DIR/$SAMPLENAME_bowtie2_aligned.tsv ``
*   This file contains the minimal amount of alignment results from Bowtie2.

##### Bowtie2 index files


*   Example not included as files are binary.
*   File name: `` $DIR/$SAMPLENAME_bowtie2_index* ``
*   These are a set of files containing the Bowtie2 index created from the custom ChocoPhlAn database.

##### Unaligned reads after Bowtie2

```
\>r4370
GGCGGACGATCTTGTCGCCCAGCCTGTAGCCTTTCTGGTACACCGTGATGACGGTGCCGCTCTCCTGCCCGTCCGTGGCGGGGATCTGCTGG
\>r4398
TGCCCGGACAGGATCTTCTCTTTCGTACCGGGCATCATCTGCTCCATGATCTCCACGCCTCGCATGAACTTTTCAGAACGGGCAACGTAGGA
```

*   File name: `` $DIR/$SAMPLENAME_bowtie2_unaligned.fa ``
*   This is a fasta file of unaligned reads after the Bowtie2 step.
*   These are the reads that will be provided as input in the translated alignment step.

##### Custom ChocoPhlAn database

```
\>gi|479150083|ref|NC_021013.1|:976220-976759|40518|g__Ruminococcus.s__Ruminococcus_bromii|UniRef90_D4L6K4|UniRef50_R6U703
ATGTTCTATGTATTTCTTGCAGAAGGCTTTGAAGAAACAGAGGCGCTTGCCCCCGTTGATGTAATGCGCAGGGCAAAGCT
TGATGTTAAAACAGTCGGTGTAACAGGCGAATGTGTTACAAGCTCACACGGTGTGCCTGTAAAAGCCGATATCACAATTG
ACAATATTGACCTTGACGATGTTCAGGGTGTTGTACTCCCCGGTGGTATGCCCGGAACTCTCAATCTTGAGGCAAACAAA
AAGGTTCTTGAGGCTGTTAAGTATAGCTGTGAAAACGGCAAAATCGTTGCCGCAATCTGTGCCGCTCCGTCAATTCTCGG
```

*   File name: `` $DIR/$SAMPLENAME_custom_chocophlan_database.ffn ``
*   This file is a custom ChocoPhlAn database of fasta sequences.

##### MetaPhlAn2 Bowtie2 output

```
r113	gi|224485636|ref|NZ_EQ973490.1|:c728571-728107
r559	gi|479185170|ref|NC_021030.1|:c1678719-1677127
r663	gi|512436175|ref|NZ_KE159463.1|:c142391-139122
r704	gi|423310881|ref|NZ_JH724270.1|:c220428-218656
r1086	gi|238922432|ref|NC_012781.1|:c1988048-1987140 
```

*   File name: `` $DIR/$SAMPLENAME_metaphlan_bowtie2.txt ``
*   This file is the Bowtie2 output from MetaPhlAn2.

##### MetaPhlAn2 bugs list

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

##### Translated alignment results

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

##### Translated alignment unaligned reads

```
\>r4370
GGCGGACGATCTTGTCGCCCAGCCTGTAGCCTTTCTGGTACACCGTGATGACGGTGCCGCTCTCCTGCCCGTCCGTGGCGGGGATCTGCTGG
\>r4398
TGCCCGGACAGGATCTTCTCTTTCGTACCGGGCATCATCTGCTCCATGATCTCCACGCCTCGCATGAACTTTTCAGAACGGGCAACGTAGGA
```
    
*   File name: `` $DIR/$SAMPLENAME_$TRANSLATEDALIGN_unaligned.fa ``
*   This is a fasta file of the unaligned reads after the translated alignment step


##### Log

```
03/16/2015 01:09:52 PM - humann2.utilities - INFO: File ( demo.fastq ) is of format:  fastq
03/16/2015 01:09:52 PM - humann2.config - INFO: Run config settings:
DATABASE SETTINGS
chocophlan database folder = data/chocophlan_DEMO
uniref database folder = data/uniref_DEM
```

*   File name: `` $DIR/$SAMPLENAME.log ``
*   This file is a log of the run.


### Workflows


#### Workflow by input file type

There are four different types of files that can be provided as input to HUMAnN2\. By default HUMAnN2 will determine the type of the file. As shown in the figure below, the type of input file will determine where HUMAnN2 will start the workflow. Files of type #2, #3, and #4 will begin the workflow after the alignment steps.

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


#### Workflow by bypass mode

There are multiple bypass options that will allow you to adjust the standard workflow.

![](http://huttenhower.sph.harvard.edu/sites/default/files/humann2_flow_bypass_modes.png)

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
    *   starts the workflow with the nucleotide alignment step using the indexed database provided with "--chocophlan $DIR/bowtie2_index"


#### Workflow of the resume option

HUMAnN2 includes a "--resume" option which will allow you to bypass alignment steps which have already been completed. For example, if you originally ran with a bypass option you can run just the step you bypassed with "--resume". This will only run the alignment step you bypassed and then recompute the gene families and pathways.

![](http://huttenhower.sph.harvard.edu/sites/default/files/humann2_flow_resume_option_no_text.png)

When using the "--resume" option, the following steps will be bypassed if they have already been completed:

1.  Taxomonic profiling step
2.  Nucleotide alignment step
3.  Custom ChocoPhlAn database creation (merge and index)
4.  Translated alignment step


### Databases

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


### Configuration

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


### Tools


#### Tools for tables

HUMAnN2 includes tools to be used with gene or pathway table files.

1.  Split a table
    *   $TABLE = gene/pathway table (tsv or biom format)
    *   $OUTPUT_DIR = the directory to write new gene/pathway tables (one per sample, in biom format if input is biom format)
2.  Join tables
    *   $INPUT_DIR = a directory containing gene/pathway tables (tsv or biom format)
    *   $TABLE = the file to write the new single gene table (biom format if input is biom format)
    *   Optional: ``--file_name $STR`` will only join gene tables with $STR in file name
3.  Rename table feature entries
    *   $TABLE = gene/pathway table (tsv format)
    *   $NAMES = mapping of feature IDs to english names (tsv format)
    *   $TABLE2 = gene/pathway table with new names attached
    *   Note: A mapping of UniRef50 IDs to english names is provided with HUMAnN2
4.  Normalize sample columns
    *   $TABLE = gene/pathway table (tsv format)
    *   $CHOICE = "relab" (relative abundance) or "cpm" (copies per million)
    *   $TABLE2 = normalized gene/pathway table
    *   Note: Can be combined with renaming


### FAQs

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

