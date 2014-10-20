
[TOC]

#** HUMAnN2: HMP Unified Metabolic Analysis Network 2 **#

##** Description **##
HUMAnN is a pipeline for efficiently and accurately determining the presence/absence and abundance of microbial pathways in a community from metagenomic sequencing data (typically millions of short DNA/RNA reads). HUMAnN2 incorporates multiple upgrades including an expanded database of microbial genomes (>4x the size of the database included in HUMAnN), a simple user interface (single command driven flow), and bug-specific output files. 

HUMAnN2 produces three output files:

1. "genefamilies.tsv": The abundance of each orthologous gene family in the community organized by bug. Orthologous families are groups of genes that perform roughly the same biological roles. HUMAnN2 uses the MetaPhlAn2 software along with the ChocoPhlAn database and UniRef for this computation.

2. "pathcoverage.tsv": The presence/absence of each pathway in the community grouped by bug. HUMAnN refers to pathway presence/absence as "coverage," and defines a pathway as a set of two or more genes. By deafult, HUMAnN2 uses MetaCyc along with MinPath for this computation. The user has the option to provide a custom pathways database to HUMAnN2 and to use all pathways instead of the minimal pathways computed by MinPath.

3. "pathabundance.tsv": The abundance of each pathway in the community grouped by bug. This is the total number of “copies” of the pathways present. By deafult, HUMAnN2 uses MetaCyc along with MinPath for this computation. The user has the option to provide a custom pathways database to HUMAnN2 and to use all pathways instead of the minimal pathways computed by MinPath.

The HUMAnN2 pipeline is a single command driven flow requiring the user to only provide a filtered fastq file to produce gene and pathway summary files ready for analysis. The pipeline converts sequence reads into coverage and abundance tables summarizing the gene families and pathways in one or more microbial communities. This lets you analyze a collection of metagenomes as a matrix of gene/pathway abundances, just like you might analyze a collection of microarrays.

##** Prerequisites **##

|Software|Source|
|-|-|
| MetaPhlAn | https://bitbucket.org/biobakery/metaphlan2/ |
| bowtie2 | http://bowtie-bio.sourceforge.net/bowtie2/ |
| rapsearch2 | http://omics.informatics.indiana.edu/mg/RAPSearch2/ |
| Python (>= version 2.7) | http://www.python.org/ |
| usearch v7.0.1001 (optional) | http://www.drive5.com/usearch/ |

|Database|Source|
|-|-|
| ChocoPhlAn | Compressed distribution (link TBD) |
| UniRef50 | ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz |

1. Download the ChocoPhlAn and UniRef50 databases from the links provided above.
2. Format the UniRef50 database for the translated alignment software you will use (rapsearch or usearch).
3. If running with rapsearch, the command to format is:
"prerapsearch -d uniref50.fasta -n processed-db-file"
4. If running with usearch, the command to format is:
"usearch -makeudb_usearch uniref50.fasta --output uniref50.udb" 
HUMAnN2 expects the UniRef50 database for usearch runs to be of the fasta format or files with 
the *.udb extension. Also to meet the maximum memory 32-bit usearch requirements it is suggested to break up
the UniRef50 database into files of at most 10,000 sequences before running the formatting step. 
HUMAnN2 can run with multiple database files for the UniRef50 database.
5. Install MetaPhlAn, bowtie2, and rapsearch (or usearch) from the links provided above. 
6. If MetaPhlAn, bowtie2, and rapsearch (or usearch) are not installed in a location in your $PATH,
then add them to your $PATH or use the HUMAnN2 parameters to indicate the locations of their
directories (--metaphlan $METAPHLAN/, --bowtie2 $BOWTIE2/, --rapsearch $RAPSEARCH/ (or --usearch $USEARCH/)). By default rapsearch is run for the translated alignment but this can be changed to usearch by setting "--translated_alignment usearch".

##** Installation **##

HUMAnN2 can be downloaded in two ways:

1. [Download](https://bitbucket.org/biobakery/humann2/downloads) a compressed set of files.
1. Create a clone of the repository on your computer with the command: 

* ``hg clone https://bitbucket.org/biobakery/humann2 ``

Note: Creating a clone of the repository requires [Mercurial](http://mercurial.selenic.com/) to be installed. Once the repository has been cloned upgrading to the latest release of HUMAnN2 is simple. Just type "hg -u pull" from within the repository which will download the latest release.

Once HUMAnN2 is downloaded, we now need to add the location of the code to the paths.
Type these commands at the prompt or include them in your .bashrc file where $HUMANn2_PATH is the location that HUMAnN2 was download (ie $HUMAnN2_PATH=/home/user/humann2/ with the file "humann2.py" found in this folder).

1. ``export PATH=$PATH:$HUMAnN2_PATH``
1. ``export PYTHONPATH=$PYTHONPATH:$HUMAnN2_PATH/src``

##** FAQS: **##

### If I run with "--temp $DIR", what are the contents of the files placed in $DIR? ###

The files placed in $DIR are temporary intermediate files. 

Their contents are as follows:

1. $DIR/$SAMPLE_bowtie2_aligned.sam : the full alignment output from bowtie2 
1. $DIR/$SAMPLE_bowtie2_aligned.tsv : only the aligned data from the bowtie2 output
1. $DIR/$SAMPLE_bowtie2_index* : bowtie2 index files created from the custom chochophlan database
1. $DIR/$SAMPLE_bowtie2_unaligned.fa : a fasta file of unaligned reads after the bowtie2 step
1. $DIR/$SAMPLE_custom_chocophlan_database.ffn : a custom chocophlan database of fasta sequences
1. $DIR/$SAMPLE_metaphlan_bowtie2.txt : the bowtie2 output from metaphlan
1. $DIR/$SAMPLE_metaphlan_bugs_list.tsv : the bugs list output from metaphlan
1. $DIR/$SAMPLE_$TRANSLATEDALIGN_aligned.tsv : the alignment results from the translated alignment step
1. $DIR/$SAMPLE_$TRANSLATEDALIGN_unaligned.fa : a fasta file of unaligned reads after the translated alignment step
1. $DIR/$SAMPLE.log : a log of the run for the $SAMPLE

* $DIR is the directory provided to store the temp files
* $SAMPLE is the basename of the fastq/fasta input file
* $TRANSLATEDALIGN is the translated alignment software selected (rapsearch2 or usearch)

### How do I bypass the MetaPhlAn prescreen step and run with the full ChocoPhlAn database? ###

Provide the "--bypass_prescreen" flag to run with the full ChocoPhlAn database.

### How do I bypass the prescreen and the nucleotide alignment index step to have HUMAnN2 start with the bowtie2 alignment step? ###

Provide the "--bypass_nucleotide_index" flag to start the HUMAnN2 run at the bowtie2 alignment step. If using this flag provide the bowtie2 index as the input to the chocophlan parameter instead of the chocoplan directory. For example, run with "--bypass_nucleotide_index --chocophlan chocophlan_dir/chocophlan_bowtie2_index"  if the first bowtie2 index file is located at chocophlan_dir/chocophlan_bowtie2_index.1.bt2 .
