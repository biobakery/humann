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
4. Download the ChocoPhlAn database to $DIR (approx. size = 5.6 GB)
    * ``$ humann2_databases --download chocophlan full $DIR``
    * When running this command, $DIR should be replaced with the full path to the directory you have selected to store the database.
    * This command will update the HUMAnN2 configuration file, storing the location you have selected for the ChocoPhlAn database. If you move this database and would like to change the configuration file, please see the [Configuration Section of the HUMAnN2 User Manual](http://huttenhower.sph.harvard.edu/humann2/manual#markdown-header-configuration). Alternatively, if you move this database, you can provide the location by adding the option "--nucleotide-database $DIR" when running HUMAnN2.
5. Download the UniRef database to $DIR (approx. size = 4.6 GB for the full database, 239 MB for the EC filtered database)
    * Download one database (full or EC filtered):
        * To download the full database (RECOMMENDED): ``$ humann2_databases --download uniref uniref50_diamond $DIR``
        * To download the EC filtered database: ``$ humann2_databases --download uniref uniref50_ec_filtered_diamond $DIR``
        * Select the full database if you are interested in identifying uncharacterized proteins in your data set. Alternatively, select the EC filtered database if you have limited disk space and/or memory. For example, a run with 13 million reads (approximately 7 GB fastq file) passed as input to the translated search step, using a single core, ran in about 4 hours with a maximum of 6 GB of memory using the EC filtered database. The same input file using the full database ran in 25 hours, with a single core, with a maximum of 11 GB of memory.  
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

The examples folder in the HUMAnN2 source download contains four demo example input files. These files are of fasta, fastq, sam, and blastm8 format. Blastm8 format is created by the following software: rapsearch2, usearch, and blast.


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
UNMAPPED        187.0
UniRef50_unknown	150.0
UniRef50_unknown|g__Bacteroides.s__Bacteroides_fragilis	150.0
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

* This file details the abundance of each gene family in the community. Gene families are groups of evolutionarily-related protein-coding sequences that often perform similar functions.
* Gene family abundance at the community level is stratified to show the contributions from known and unknown species. Individual species' abundance contributions sum to the community total abundance. 
* HUMAnN2 uses the MetaPhlAn2 software along with the ChocoPhlAn database and UniRef for this computation.
* Gene family abundance is reported in RPK (reads per kilobase) units to normalize for gene length; RPK units reflect relative gene (or transcript) copy number in the community. RPK values can be further sum-normalized to adjust for differences in sequencing depth across samples. For more information on these units and normalization, please see the [HUMAnN2 User Manual Standard Workflow Section](http://huttenhower.sph.harvard.edu/humann2/manual/#markdown-header-standard-workflow)
* The "UNMAPPED" value is the total number of reads which remain unmapped after both alignment steps (nucleotide and translated search). Since other gene features in the table are quantified in RPK units, "UNMAPPED" can be interpreted as a single unknown gene of length 1 kilobase recruiting all reads that failed to map to known sequences.
* The UniRef50_unknown values represent the total abundance of reads which map to ChocoPhlAn nucleotide sequences which do not have a UniRef50 annotation.

#### Pathway Coverage ####

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

* This file includes the presence/absence of each pathway in the community grouped by species. HUMAnN refers to pathway presence/absence as "coverage" and defines a pathway as a set of two or more genes.
*   Pathway coverage at the community level is stratified to show the contributions from known and unknown species. **A pathway's community-level coverage is not necessarily the sum of its stratified coverage values.** For example, in the two-gene pathway {A, B}, if species 1 contributes abundances {A=5, B=5} and species 2 contributes abundance {A=10, B=10}, the pathway has coverage=1.0 in species 1, species 2, and at the community level.
* HUMAnN2 uses MetaCyc pathways along with MinPath for this computation. 
* The user has the option to provide a custom pathways database to HUMAnN2 and to use all pathways instead of the minimal pathways computed by MinPath.
*   This file follows the same order for pathways and species as the abundance file. The values for UNMAPPED and UNINTEGRATED are set to 1.0 included so that this file will match the format of the abundance file exactly.

#### Pathway Abundance ####

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

* This file includes the abundance of each pathway in the community as a function of the abundances of the pathway's member genes. Pathway abundance is proportional to the total number of “copies” of the pathways present. 
*   Pathway abundance at the community level is stratified to show the contributions from known and unknown species. **A pathway's community-level abundance is not necessarily the sum of its stratified abundance values.** For example, in the two-gene pathway {A, B}, if species 1 contributes abundances {A=5, B=10} and species 2 contributes abundances {A=10, B=5}, species 1 and 2 each contribute 5 complete copies of the pathway, but at the community level there are 15 complete copies.
* HUMAnN2 uses MetaCyc pathways along with MinPath for this computation. 
* The user has the option to provide a custom pathways database to HUMAnN2 and to use all pathways instead of the minimal pathways computed by MinPath.
* To account for non-linearity in the conversion of gene copy number to pathway copy number, we define a “compression constant” (*k*) equal to the total pathway abundance divided by the total abundance of genes that contributed to pathways. The "UNMAPPED" value reported in the pathway abundance table is equal to the total number of unmapped reads scaled by *k* (making it more comparable with pathway abundance values). Similarly, we define an "UNINTEGRATED" abundance for 1) the community, 2) each identified species, and 3) unclassified species equal to the total abundance of genes in that level that did not contribute to pathways scaled by *k*.
* The pathways are ordered by decreasing abundance with pathways for each species also sorted by decreasing abundance. Pathways with zero abundance are not included in the file.
