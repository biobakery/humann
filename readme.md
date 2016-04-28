# HUMAnN2 #

----

 * Install HUMAnN2 with `` $ pip install humann2 `` then follow the [steps to start running HUMAnN2](#markdown-header-getting-started-with-humann2).

 * Optionally, download the HUMAnN2 software source and demos ( [humann2.tar.gz](https://pypi.python.org/pypi/humann2) ).

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

When installing HUMAnN2, please also download and install [MetaPhlAn2](http://huttenhower.sph.harvard.edu/metaphlan2). You can then add the MetaPhlAn2 folder to your $PATH or you can provide its location when running HUMAnN2 with the option "--metaphlan $DIR" (replacing $DIR with the full path to the MetaPhlAn2 folder). Bowtie2, Diamond, and Minpath will be automatically installed when you install HUMAnN2. 

1. Install HUMAnN2
    * `` $ pip install humann2 ``
    * This command will automatically install MinPath (and a new version of glpk) along with Bowtie2 and Diamond (if they are not already installed).
    * To bypass the install of Bowtie2 and Diamond, add the option "--install-option='--bypass-dependencies-install'" to the install command.
    * To build Diamond from source during the install, add the option "--install-option='--build-diamond'" to the install command.
    * To overwite existing installs of Bowtie2 and Diamond, add the option "--install-option='--replace-dependencies-install'" to the install command.
    * If you do not have write permissions to '/usr/lib/', then add the option "--user" to the HUMAnN2 install command. This will install the python package into subdirectories of '~/.local' on Linux. Please note when using the "--user" install option on some platforms, you might need to add '~/.local/bin/' to your $PATH as it might not be included by default. You will know if it needs to be added if you see the following message `humann2: command not found` when trying to run HUMAnN2 after installing with the "--user" option.
2. Test the HUMAnN2 install (Optional)
     * `` $ humann2_test``
     * To also run tool tests, add the option "--run-functional-tests-tools".
     * To also run end-to-end tests, add the option "-run-functional-tests-end-to-end". Please note these tests take about 20 minutes to run. Also they require all dependencies of HUMAnN2 be installed in your PATH.
3. Try out a HUMAnN2 demo run (Optional)
    * Download the HUMAnN2 source with demos: [humann2.tar.gz](https://pypi.python.org/pypi/humann2)
    * `` $ tar zxvf humann2.tar.gz ``
    * `` $ humann2 --input humann2/examples/demo.fastq --output $OUTPUT_DIR ``
    * When running this command, $OUTPUT_DIR should be replaced with the full path to the directory you have selected to write the output from the HUMAnN2 demo run.
    * Other types of demo files are included in this folder and can be run with the exact same command.
    * Demo ChocoPhlAn and UniRef databases are also included in the download. The demo ChocoPhlAn database is located a humann2/data/chocophlan_DEMO and the demo UniRef database is located a humann2/data/uniref_DEMO. Until the full databases are downloaded HUMAnN2 will run with the demo database by default.
4. Download the ChocoPhlAn database (approx. size = 5.6 GB)
    * ``$ humann2_databases --download chocophlan full $DIR``
    * When running this command, $DIR should be replaced with the full path to the directory you have selected to store the database.
    * This command will update the HUMAnN2 configuration file, storing the location you have selected for the ChocoPhlAn database. If you move this database and would like to change the configuration file, please see the [Configuration Section of the HUMAnN2 User Manual](http://huttenhower.sph.harvard.edu/humann2/manual#markdown-header-configuration). Alternatively, if you move this database, you can provide the location by adding the option "--nucleotide-database $DIR" when running HUMAnN2.
5. Download a UniRef database (only download one database: UniRef50 full, UniRef50 EC filtered, UniRef90 full, or UniRef90 EC filtered)
    * Download one of the following databases (replacing $DIR with the location to store the database):
        * To download the UniRef90 EC filtered database (RECOMMENDED, approx. size = 846 MB): 
            * ``$ humann2_databases --download uniref uniref90_ec_filtered_diamond $DIR``
        * To download the full UniRef90 database (approx. size = 11 GB): 
            * ``$ humann2_databases --download uniref uniref90_diamond $DIR``
        * To download the UniRef50 EC filtered database (approx. size = 239 MB): 
            * ``$ humann2_databases --download uniref uniref50_ec_filtered_diamond $DIR``
        * To download the full UniRef50 database (approx. size = 4.6 GB): 
            * ``$ humann2_databases --download uniref uniref50_diamond $DIR``
    * Select a full database if you are interested in identifying uncharacterized proteins in your data set. Alternatively, select an EC filtered database if you have limited disk space and/or memory. For example, a run with 13 million reads (approximately 7 GB fastq file) passed as input to the translated search step, using a single core, ran in about 4 hours with a maximum of 6 GB of memory using the UniRef50 EC filtered database. The same input file using the UniRef50 full database ran in 25 hours, with a single core, with a maximum of 11 GB of memory.  
    * The download command will update the HUMAnN2 configuration file, storing the location you have selected for the UniRef database. If you move this database and would like to change the configuration file, please see the [Configuration Section of the HUMAnN2 User Manual](http://huttenhower.sph.harvard.edu/humann2/manual#markdown-header-configuration). Alternatively, if you move this database, you can provide the location by adding the option "--protein-database $DIR" when running HUMAnN2.


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

### Standard workflow ###

The standard workflow involves running HUMAnN2 on each filtered shotgun sequencing metagenome file, normalizing, and then merging output files.

For detailed information on the standard workflow, please see the [HUMAnN2 User Manual Standard Workflow Section](http://huttenhower.sph.harvard.edu/humann2/manual/#markdown-header-standard-workflow).

