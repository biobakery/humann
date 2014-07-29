"""
Run alignment, find unused reads
"""

import os
import utilities, config

def alignment(uniref, unaligned_reads_file_fastq, identity_threshold, 
    temp_dir, threads):
    """
    Run usearch for alignment
    """

    unaligned_name = os.path.splitext(os.path.basename(unaligned_reads_file_fastq))[0]
    alignment_file_base = os.path.join(temp_dir, 
        unaligned_name + config.translated_alignment_name)

    exe="usearch"
    opts=config.usearch_opts

    args=["-usearch_global",unaligned_reads_file_fastq,"-id",identity_threshold,
        "-userfields","query+target+id"]

    if threads > 1:
        args+=["-threads",threads]

    print "\nRunning muliple " + exe + " ........\n"

    args+=opts

    #run the search on each of the databases in the directory
    index=1
    temp_out_files=[]
    for database in os.listdir(uniref):
        input_database=os.path.join(uniref,database)
        full_args=args+["-db",input_database]

        # name temp output file
        ext= "." + str(index)+ config.usearch_name_temp
        temp_out_file=alignment_file_base + ext
        temp_out_files.append(temp_out_file)

        full_args+=["-userout",temp_out_file]

        utilities.execute_command(exe,full_args,[input_database],[temp_out_file],"")
        index+=1
      
    # merge the temp output files
    exe="cat"
    alignment_file=alignment_file_base + ".tsv"

    utilities.execute_command(exe,temp_out_files,temp_out_files,[alignment_file],
        alignment_file)  

    # remove the temp files which have been merged
    for temp_file in temp_out_files:
        utilities.remove_temp_file(temp_file)

    return alignment_file

def unaligned_reads(input_fastq, alignment_file, temp_dir):
    """
    Create a fasta/fastq file of the unaligned reads
    """

    #name the index
    alignment_file_name = os.path.splitext(os.path.basename(alignment_file))[0]
    unaligned_reads_file_fastq = os.path.join(temp_dir,
        alignment_file_name + config.unaligned_reads_name_no_ext)

    #determine the index to use for the fastq/fasta file
    #use the same as that that was used by the user for the input file
    original_extension = os.path.splitext(os.path.basename(input_fastq))[1]
    unaligned_reads_file_fastq+= original_extension

    utilities.unaligned_reads_from_tsv(input_fastq, alignment_file, 
        unaligned_reads_file_fastq)

    return unaligned_reads_file_fastq

