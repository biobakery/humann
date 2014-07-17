#!/usr/bin/env python
"""
Index database, run alignment, find unused reads
"""

import os
import utilities

def alignment(custom_database, user_fastq, debug_dir, threads):
    """
    Index database and run alignment with bowtie2
    """
    # name the index
    custom_database_name = os.path.splitext(os.path.basename(custom_database))[0]
    index_name = os.path.join(debug_dir, custom_database_name + "_bowtie2_index")
  
    # if indexes are ever large (>4G) then use the --large-index flag
    exe="bowtie2-build"
    args=" -f " + custom_database + " " + index_name

    outfiles=[index_name + ".1.bt2"]

    # index the database
    print "\nRunning " + exe + " ........\n"
    #utilities.execute_command(exe,args,[custom_database],outfiles)

    # name the alignment file
    sample_name = os.path.splitext(os.path.basename(user_fastq))[0]
    alignment_file = os.path.join(debug_dir, custom_database_name + "_chocophlan_align.sam")


    # align user input to database
    exe="bowtie2"

    #determine input type as fastq or fasta
    input_type = utilities.fasta_or_fastq(user_fastq)

    #determine input type flag
    #default flag to fastq
    input_type_flag = " -q "
    if input_type == "fasta":
        input_type_flag=" -f "

    args=input_type_flag + " -x " + index_name + " -U " + user_fastq + " -S " + alignment_file
    
    #add threads
    if threads > 1:
        args+=" -p " + threads

    # run the bowtie2 alignment
    print "\nRunning " + exe + " ........\n"
    #utilities.execute_command(exe,args,[user_fastq],[alignment_file])

    return alignment_file

def unaligned_reads(alignment_file, debug_dir):
    """ 
    Return file of just the unaligned reads
    """

    # name the index
    alignment_file_name = os.path.splitext(os.path.basename(alignment_file))[0]
    unaligned_reads_file = os.path.join(debug_dir, alignment_file_name + "_unaligned_reads.sam")

    # pull out the unaligned reads from the sam file
    exe = "samtools"
    args = " view -Sf 4 " + alignment_file + " -o " + unaligned_reads_file

    print "\nRunning " + exe + " ........\n"
    utilities.execute_command(exe,args,[alignment_file],[unaligned_reads_file])

    return unaligned_reads_file
