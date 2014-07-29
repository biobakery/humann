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
    alignment_file_base = os.path.join(temp_dir, unaligned_name + "_translated_aligned")

    exe="usearch"
    opts=config.usearch_opts

    args=["-usearch_global",unaligned_reads_file_fastq,"-id",identity_threshold,
        "-userfields","query+target+id"]

    if threads > 1:
        args+=["-threads",threads]

    print "\nRunning " + exe + " ........\n"

    args+=opts

    #run the search on each of the databases in the directory
    index=1
    temp_out_files=[]
    for database in os.listdir(uniref):
        input_database=os.path.join(uniref,database)
        full_args=args+["-db",input_database]

        # name temp output file
        ext= "." + str(index)+ ".usearch.tmp"
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
