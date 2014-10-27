"""
HUMAnN2: translated_search module
Run alignment, find unused reads

Copyright (c) 2014 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import os
import re
import numbers
import logging
import math
import traceback

import utilities
import config
import store

# name global logging instance
logger=logging.getLogger(__name__)

def usearch_alignment(alignment_file, uniref, unaligned_reads_file_fasta):
    """
    Run usearch alignment with memory management
    Individual runs can be threaded
    """

    bypass=utilities.check_outfiles([alignment_file])

    exe="usearch"
    opts=config.usearch_opts

    args=["-id",config.identity_threshold]

    message="Running " + exe + " ........"
    logger.info(message)
    print("\n"+message+"\n")

    if not bypass:

        args+=opts

        #break up input file into smaller files for memory requirement
        logger.debug("Break up fasta file into smaller files for memory management")
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
                temp_out_file=utilities.unnamed_temp_file()
                temp_out_files.append(temp_out_file)

                full_args+=["-blast6out",temp_out_file]

                command_args.append([exe,full_args,[input_database],[],"",""])
                
        results=utilities.command_multiprocessing(config.threads,command_args)

        # merge the temp output files
        exe="cat"

        utilities.execute_command(exe,temp_out_files,temp_out_files,[alignment_file],
            alignment_file)

    else:
        message="Bypass"
        logger.info(message)
        print(message)


def rapsearch_alignment(alignment_file,uniref, unaligned_reads_file_fasta):
    """
    Run rapsearch alignment on database formatted for rapsearch
    """

    bypass=utilities.check_outfiles([alignment_file])

    exe="rapsearch"
    opts=config.rapsearch_opts

    args=["-q",unaligned_reads_file_fasta,"-b",0]

    if config.threads > 1:
        args+=["-z",config.threads]

    message="Running " + exe + " ........"
    logger.info(message)
    print("\n"+message+"\n")

    if not bypass:

        args+=opts

        temp_out_files=[]
        for database in os.listdir(uniref):
            # ignore the *.info database files
            if not database.endswith(config.rapsearch_database_extension):
                input_database=os.path.join(uniref,database)
                full_args=args+["-d",input_database]

                # create temp output file
                temp_out_file=utilities.unnamed_temp_file()
                utilities.remove_file(temp_out_file)

                temp_out_files.append(temp_out_file+config.rapsearch_output_file_extension)

                full_args+=["-o",temp_out_file]

                utilities.execute_command(exe,full_args,[input_database],[])
        
        # merge the temp output files
        utilities.execute_command("cat",temp_out_files,temp_out_files,[alignment_file],
            alignment_file)

    else:
        message="Bypass"
        logger.info(message)
        print(message)



def alignment(uniref, unaligned_reads_file):
    """
    Run rapsearch2 or usearch for alignment
    """

    alignment_file = utilities.name_temp_file( 
        "_" + config.translated_alignment_selected 
        + config.translated_alignment_name)
    
    # Check that the file of reads to align is fasta
    temp_file=""
    if utilities.fasta_or_fastq(unaligned_reads_file) == "fastq":
        logger.debug("Convert unaligned reads fastq file to fasta")
        input_fasta=utilities.fastq_to_fasta(unaligned_reads_file)
        temp_file=input_fasta
    else:
        input_fasta=unaligned_reads_file

    if config.translated_alignment_selected == "usearch":
        usearch_alignment(alignment_file, uniref, input_fasta)
    else:
        rapsearch_alignment(alignment_file, uniref, input_fasta)
        
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
                    except ValueError:
                        logger.warning("Unable to convert rapsearch e-value: %s", evalue)
                        logger.warning("Traceback: \n" + traceback.format_exc())
                        evalue=1.0 
                else:
                    if not isinstance(evalue, numbers.Number):
                        logger.warning("Usearch e-value is not a number: %s", evalue)
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

