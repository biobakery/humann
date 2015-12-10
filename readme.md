# HUMAnN2 #

----

 * Download the HUMAnN2 software ( [humann2.tar.gz](https://bitbucket.org/biobakery/humann2/downloads/humann2_v0.5.0.tar.gz) ) then follow the [steps to install and run](#markdown-header-getting-started-with-humann2).

 * For additional information, please see the [HUMAnN2 User Manual](http://huttenhower.sph.harvard.edu/humann2/manual).

 * Please direct questions to the [HUMAnN2 google group](https://groups.google.com/forum/#!forum/humann-users) (subscribe to receive HUMAnN2 news).

 * If you use the HUMAnN2 software, please cite our manuscript: TBD

----

HUMAnN2 is the next generation of [HUMAnN (HMP Unified Metabolic Analysis Network)](http://huttenhower.sph.harvard.edu/humann).

HUMAnN is a pipeline for efficiently and accurately profiling the presence/absence and abundance of microbial pathways in a community from metagenomic or metatranscriptomic sequencing data (typically millions of short DNA/RNA reads). This process, referred to as functional profiling, aims to describe the metabolic potential of a microbial community and its members. More generally, functional profiling answers the question "What are the microbes in my community-of-interest doing (or capable of doing)?"


## HUMAnN2 Features ##

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


## HUMAnN2 Workflow ##

![](http://huttenhower.sph.harvard.edu/sites/default/files/humann2_diamond_500x500.jpg)


## Getting Started with HUMAnN2 ##

### Requirements ###

1.  [MetaPhlAn2](http://huttenhower.sph.harvard.edu/metaphlan2)
2.  [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (version >= 2.1) (see NOTE)
3.  [Diamond](http://ab.inf.uni-tuebingen.de/software/diamond/) (version >= 0.7.3) (see NOTE)
4.  [MinPath](http://omics.informatics.indiana.edu/MinPath/) (see NOTE)
5.  [Python](http://www.python.org/) (version >= 2.7)
6.  Memory (>= 16 GB)
7.  Disk space (>= 10 GB [to accommodate comprehensive sequence databases])
8.  Operating system (Linux or Mac)

NOTE: Bowtie2, Diamond, and MinPath are automatically installed when installing HUMAnN2. 


### Installation ###

Before installing HUMAnN2, please download and install [MetaPhlAn2](http://huttenhower.sph.harvard.edu/metaphlan2). You can then add the MetaPhlAn2 folder to your $PATH or you can provide its location when running HUMAnN2 with the option "--metaphlan $DIR" (replacing $DIR with the full path to the MetaPhlAn2 folder). Bowtie2, Diamond, and Minpath will be automatically installed when you install HUMAnN2. 

1. Download and unpack the HUMAnN2 software
    * Download the software: [humann2.tar.gz](https://bitbucket.org/biobakery/humann2/downloads/humann2_v0.5.0.tar.gz)
    * `` $ tar zxvf humann2.tar.gz ``
    * `` $ cd humann2 ``
3. Install the HUMAnN2 software
    * `` $ python setup.py install ``
    * This command will automatically install MinPath (and a new version of glpk) along with Bowtie2 and Diamond (if they are not already installed).
    * To bypass the install of Bowtie2 and Diamond, add the option "--bypass-dependencies-install" to the install command.
    * If you do not have write permissions to '/usr/lib/', then add the option "--user" to the HUMAnN2 install command. This will install the python package into subdirectories of '~/.local'. Please note when using the "--user" install option on some platforms, you might need to add '~/.local/bin/' to your $PATH as it might not be included by default. You will know if it needs to be added if you see the following message `humann2: command not found` when trying to run HUMAnN2 after installing with the "--user" option.
4. Test the HUMAnN2 install (Optional)
     * `` $ python setup.py test``
5. Try out a HUMAnN2 demo run (Optional)
    * `` $ humann2 --input humann2/examples/demo.fastq --output $OUTPUT_DIR ``
    * When running this command, $OUTPUT_DIR should be replaced with the full path to the directory you have selected to write the output from the HUMAnN2 demo run.
6. Download the ChocoPhlAn database to $DIR (approx. size = 5.6 GB)
    * ``$ humann2_databases --download chocophlan full $DIR``
    * When running this command, $DIR should be replaced with the full path to the directory you have selected to store the database.
    * This command will update the HUMAnN2 configuration file, storing the location you have selected for the ChocoPhlAn database. If you move this database and would like to change the configuration file, please see the [Configuration Section of the HUMAnN2 User Manual](http://huttenhower.sph.harvard.edu/humann2/manual#markdown-header-configuration). Alternatively, if you move this database, you can provide the location by adding the option "--nucleotide-database $DIR" when running HUMAnN2.
7. Download the UniRef database to $DIR (approx. size = 2.8 GB)
    * ``$ humann2_databases --download uniref diamond $DIR``
    * When running this command, $DIR should be replaced with the full path to the directory you have selected to store the database.
    * This command will update the HUMAnN2 configuration file, storing the location you have selected for the UniRef database. If you move this database and would like to change the configuration file, please see the [Configuration Section of the HUMAnN2 User Manual](http://huttenhower.sph.harvard.edu/humann2/manual#markdown-header-configuration). Alternatively, if you move this database, you can provide the location by adding the option "--protein-database $DIR" when running HUMAnN2.


### How to Run ###

#### Basic usage ####

`` $ humann2 --input $SAMPLE --output $OUTPUT_DIR``

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

#### Demo runs ####

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

### Standard workflow ###

The standard workflow involves running HUMAnN2 on each filtered shotgun sequencing metagenome file, normalizing, and then merging output files.

For detailed information on the standard workflow, please see the [HUMAnN2 User Manual Standard Workflow Section](http://huttenhower.sph.harvard.edu/humann2/manual/#markdown-header-standard-workflow).

### Output files ###

HUMAnN2 produces three output files which by default are tab-delimited text. There is an option to print out files in biom format. 

#### Gene Families ####

```
# Gene Family	$SAMPLENAME_Abundance
UniRef50_A6L0N6: Conserved protein found in conjugate transposon	67.0
UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_fragilis	8.0
UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_finegoldii	5.0
UniRef50_A6L0N6: Conserved protein found in conjugate transposon|g__Bacteroides.s__Bacteroides_stercoris	4.0
UniRef50_A6L0N6: Conserved protein found in conjugate transposon|unclassified	1.0
UniRef50_O83668: Fructose-bisphosphate aldolase	60.0
UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_vulgatus	31.0
UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_thetaiotaomicron	22.0
UniRef50_O83668: Fructose-bisphosphate aldolase|g__Bacteroides.s__Bacteroides_stercoris	7.0
```

* This file includes the abundance of each orthologous gene family in the community organized by bug. Orthologous families are groups of genes that perform roughly the same biological roles. 
* HUMAnN2 uses the MetaPhlAn2 software along with the ChocoPhlAn database and UniRef for this computation.
* Gene family abundance is reported in RPK (reads per kilobase) units to normalize for gene length; RPK units reflect relative gene (or transcript) copy number in the community. RPK values can be further sum-normalized to adjust for differences in sequencing depth across samples. For more information on these units and normalization, please see the [HUMAnN2 User Manual Standard Workflow Section](http://huttenhower.sph.harvard.edu/humann2/manual/#markdown-header-standard-workflow)

#### Pathway Coverage ####

```
# Pathway	$SAMPLENAME_Coverage
PWY0-1301: melibiose degradation	1.0
PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_caccae	1.0
PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_finegoldii	1.0
PWY0-1301: melibiose degradation|unclassified	1.0
PWY-5484: glycolysis II (from fructose-6P)	1.0
PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_vulgatus	0.7
PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_thetaiotaomicron	0.7
PWY-5484: glycolysis II (from fructose-6P)|unclassified	0.3
PWY-5484: glycolysis II (from fructose-6P)|g__Parabacteroides.s__Parabacteroides_merdae	0.3
```

* This file includes the presence/absence of each pathway in the community grouped by bug. HUMAnN refers to pathway presence/absence as "coverage" and defines a pathway as a set of two or more genes. 
* HUMAnN2 uses MetaCyc pathways along with MinPath for this computation. 
* The user has the option to provide a custom pathways database to HUMAnN2 and to use all pathways instead of the minimal pathways computed by MinPath.

#### Pathway Abundance ####

```
# Pathway	$SAMPLENAME_Abundance
PWY0-1301: melibiose degradation	57.5
PWY0-1301: melibiose degradation|unclassified	32.5
PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_ovatus	4.5
PWY0-1301: melibiose degradation|g__Alistipes.s__Alistipes_putredinis	3.0
PWY0-1301: melibiose degradation|g__Bacteroides.s__Bacteroides_caccae	2.5
PWY-5484: glycolysis II (from fructose-6P)	54.7
PWY-5484: glycolysis II (from fructose-6P)|unclassified	16.7
PWY-5484: glycolysis II (from fructose-6P)|g__Parabacteroides.s__Parabacteroides_merdae	8.0
PWY-5484: glycolysis II (from fructose-6P)|g__Bacteroides.s__Bacteroides_caccae	6.0
```

* This file includes the abundance of each pathway in the community grouped by bug. This is the total number of “copies” of the pathways present. 
* HUMAnN2 uses MetaCyc pathways along with MinPath for this computation. 
* The user has the option to provide a custom pathways database to HUMAnN2 and to use all pathways instead of the minimal pathways computed by MinPath.
