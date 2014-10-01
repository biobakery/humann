"""
Run alignment, find unused reads
"""

import os
import tempfile
import re
import numbers

import utilities
import config
import store

def usearch_alignment(alignment_file, uniref, unaligned_reads_file_fastq):
    """
    Run usearch alignment with memory management
    Individual runs can be threaded
    """

    bypass=utilities.check_outfiles([alignment_file])

    exe="usearch"
    opts=config.usearch_opts

    args=["-id",config.identity_threshold]

    print("\nRunning " + exe + " ........\n")

    if not bypass:

        args+=opts

        #convert unaligned reads file to fasta
        temp_fasta_files=[]
        if utilities.fasta_or_fastq(unaligned_reads_file_fastq) == "fastq":
            utilities.unaligned_reads_file_fasta=fastq_to_fasta(unaligned_reads_file_fastq)
            temp_fasta_files.append(unaligned_reads_file_fasta)
        else:
            unaligned_reads_file_fasta=unaligned_reads_file_fastq

        #break up input file into smaller files for memory requirement
        temp_in_files=utilities.break_up_fasta_file(unaligned_reads_file_fasta,
            config.usearch_max_seqs)

       #run the search on each of the databases in the directory
        temp_out_files=[]
        command_args=[]
        for input_file in temp_in_files:
            for database in os.listdir(uniref):
                input_database=os.path.join(uniref,database)
                full_args=["-usearch_global",input_file]+args+["-db",input_database]

                # create temp output file
                file_out, temp_out_file=tempfile.mkstemp()
                os.close(file_out)
                temp_out_files.append(temp_out_file)

                full_args+=["-blast6out",temp_out_file]

                command_args.append([exe,full_args,[input_database],[],"",""])
                
        results=utilities.command_multiprocessing(config.threads,command_args)

        # merge the temp output files
        exe="cat"

        utilities.execute_command(exe,temp_out_files,temp_out_files,[alignment_file],
            alignment_file)

        # remove the temp files which have been merged
        for temp_file in temp_fasta_files + temp_in_files + temp_out_files:
            utilities.remove_file(temp_file)
    else:
        print("Bypass")


def rapsearch_alignment(alignment_file,uniref,
        unaligned_reads_file_fastq):
    """
    Run rapsearch alignment on database formatted for rapsearch
    """

    bypass=utilities.check_outfiles([alignment_file])

    exe="rapsearch"
    opts=config.rapsearch_opts

    args=["-q",unaligned_reads_file_fastq,"-b",0]

    if config.threads > 1:
        args+=["-z",config.threads]

    print("\nRunning " + exe + " ........\n")

    if not bypass:

        args+=opts

        temp_out_files=[]
        for database in os.listdir(uniref):
            # ignore the *.info database files
            if not database.endswith(config.rapsearch_database_extension):
                input_database=os.path.join(uniref,database)
                full_args=args+["-d",input_database]

                # create temp output file
                file_out, temp_out_file=tempfile.mkstemp()
                os.close(file_out)
                utilities.remove_file(temp_out_file)

                temp_out_files.append(temp_out_file+config.rapsearch_output_file_extension)

                full_args+=["-o",temp_out_file]

                utilities.execute_command(exe,full_args,[input_database],[])
        
        # merge the temp output files
        utilities.execute_command("cat",temp_out_files,temp_out_files,[alignment_file],
            alignment_file)

        # remove the temp files which have been merged
        for temp_file in temp_out_files:
            utilities.remove_file(temp_file)
    else:
        print("Bypass")



def alignment(uniref, unaligned_reads_file_fasta):
    """
    Run rapsearch2 or usearch for alignment
    """

    alignment_file = utilities.name_temp_file( 
        "_" + config.translated_alignment_selected 
        + config.translated_alignment_name)
    
    # Check that the file of reads to align is fasta
    temp_file=""
    if utilities.fasta_or_fastq(unaligned_reads_file_fasta) == "fastq":
        input_fasta=utilities.fastq_to_fasta(unaligned_reads_file_fasta)
        temp_file=input_fasta
    else:
        input_fasta=unaligned_reads_file_fasta

    if config.translated_alignment_selected == "usearch":
        usearch_alignment(alignment_file, uniref, input_fasta)
    else:
        rapsearch_alignment(alignment_file,uniref,
            input_fasta)
        
    # Remove the temp fasta file if exists
    if temp_file:
        utilities.remove_file(temp_file)

    return alignment_file

def unaligned_reads(unaligned_reads_store, alignment_file_tsv, alignments):
    """
    Create a fasta file of the unaligned reads
    Store the alignment results
    """

    #create a fasta file of unaligned reads
    unaligned_file_fasta= utilities.name_temp_file(
        "_" + config.translated_alignment_selected + 
        config.translated_unaligned_reads_name_no_ext + 
        config.fasta_extension)
    
    utilities.file_exists_readable(alignment_file_tsv)

    # read through the alignment file to identify ids
    # that correspond to aligned reads
    # all translated alignment files will be of the tabulated blast format
    file_handle=open(alignment_file_tsv,"r")
    line=file_handle.readline()

    aligned_ids=[]
    while line:
        if not re.search("^#",line):
            alignment_info=line.split(config.blast_delimiter)
            identity=float(alignment_info[config.blast_identity_index])
            
            if identity >= config.identity_threshold:
                # only store those alignments which meet the identity threshold
                
                referenceid=alignment_info[config.blast_reference_index]
                queryid=alignment_info[config.blast_query_index]
                aligned_length=int(alignment_info[config.blast_aligned_length_index])
                evalue=alignment_info[config.blast_evalue_index]
                
                # remove the id of the alignment from the unaligned reads store
                unaligned_reads_store.remove_id(queryid)                
                
                if config.translated_alignment_selected == "rapsearch":
                    try:
                        evalue=math.pow(10.0, float(evalue))
                    except:
                        evalue=1.0 
                else:
                    if not isinstance(evalue, numbers.Number):
                        evalue=1.0
            
                alignments.add(referenceid, queryid, evalue,"unclassified")
            
                aligned_ids+=[queryid]
        line=file_handle.readline()

    file_handle.close()

    # create unaligned file using list of remaining unaligned stored data
    file_handle_write=open(unaligned_file_fasta,"w")

    file_handle_write.write(unaligned_reads_store.get_fasta())

    file_handle_write.close()

    return unaligned_file_fasta

