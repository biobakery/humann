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
                if database.endswith(config.usearch_database_extension):
                    input_database=os.path.join(uniref,database)
                    full_args=["-ublast",input_file]+args+["-db",input_database]
    
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
        
        # Find the rapsearch database files in the directory
        # These will be files of the same name as the *.info files
        files=os.listdir(uniref)
        rapsearch_databases=[]
        for file in files:
            if file.endswith(config.rapsearch_database_extension):
                # Check for the corresponding database file
                database_file=re.sub(config.rapsearch_database_extension+"$","",file)
                if database_file in files:
                    rapsearch_databases.append(database_file)
                
        for database in rapsearch_databases:
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


def diamond_alignment(alignment_file,uniref, unaligned_reads_file_fasta):
    """
    Run diamond alignment on database formatted for diamond
    """

    bypass=utilities.check_outfiles([alignment_file])

    exe="diamond"
    
    # Select the command based on a protein or nucleotide database search
    args=[]
    if config.pick_frames_toggle == "on":
        args=[config.diamond_cmmd_protein_search]
    else:
        args=[config.diamond_cmmd_nucleotide_search]
        
    opts=config.diamond_opts

    args+=["--query",unaligned_reads_file_fasta]

    if config.threads > 1:
        args+=["--threads",config.threads]

    message="Running " + exe + " ........"
    logger.info(message)
    print("\n"+message+"\n")

    if not bypass:

        args+=opts

        temp_out_files=[]
        for database in os.listdir(uniref):
            # ignore any files that are not the database files
            if database.endswith(config.diamond_database_extension):
                # Provide the database name without the extension
                input_database=os.path.join(uniref,database)
                input_database_extension_removed=re.sub(config.diamond_database_extension
                    +"$","",input_database)
                full_args=args+["--db",input_database_extension_removed]
    
                # create temp output file
                temp_out_file=utilities.unnamed_temp_file()
                utilities.remove_file(temp_out_file)
    
                temp_out_files.append(temp_out_file)
    
                full_args+=["--out",temp_out_file,"--tmpdir",os.path.dirname(temp_out_file)]
    
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
        # Convert file to fasta, also pick frames if selected
        if config.pick_frames_toggle == "on":
            logger.debug("Applying pick frames")
            input_fasta=utilities.fastq_to_fasta(unaligned_reads_file,
                apply_pick_frames=True)
        else:
            input_fasta=utilities.fastq_to_fasta(unaligned_reads_file)
        # set the file as a temp to be removed later
        temp_file=input_fasta
    elif config.bypass_nucleotide_search and config.pick_frames_toggle == "on":
        # Process the fasta file to pick frames
        logger.debug("Applying pick frames")
        input_fasta=utilities.pick_frames_from_fasta(unaligned_reads_file)
        # set the file as a temp to be removed later
        temp_file=input_fasta
    else:
        input_fasta=unaligned_reads_file

    if config.translated_alignment_selected == "usearch":
        usearch_alignment(alignment_file, uniref, input_fasta)
    elif config.translated_alignment_selected == "rapsearch":
        rapsearch_alignment(alignment_file, uniref, input_fasta)
    elif config.translated_alignment_selected == "diamond":
        diamond_alignment(alignment_file, uniref, input_fasta)
    else:
        sys.exit("CRITICAL ERROR: The translated alignment software selected is not"
            + " available: " + config.translated_alignment_selected )
        
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
    try:
        utilities.file_exists_readable(alignment_file_tsv,raise_IOError=True)
    except IOError:
        message="No alignment results found from translated search"
        logger.critical(message)
        print(message)
        return unaligned_file_fasta
        

    # read through the alignment file to identify ids
    # that correspond to aligned reads
    # all translated alignment files will be of the tabulated blast format
    file_handle=open(alignment_file_tsv,"r")
    line=file_handle.readline()

    aligned_ids=[]
    log_evalue=False
    while line:
        if re.search("^#",line):
            # Check for the rapsearch2 header to determine if these are log(e-value)
            if re.search(config.blast_delimiter,line):
                data=line.split(config.blast_delimiter)
                if len(data)>config.blast_evalue_index:
                    if re.search("log",data[config.blast_evalue_index]):
                        log_evalue=True
        else:
            alignment_info=line.split(config.blast_delimiter)
            
            # try to obtain the identity value to determine if threshold is met
            process_alignment=True
            try:
                identity=float(alignment_info[config.blast_identity_index])
                if identity < config.identity_threshold:
                    process_alignment=False
            except ValueError:
                if alignment_info[config.blast_identity_index]:
                    logger.debug("Unexpected value in identity field: %s",
                        alignment_info[config.blast_identity_index])
                    process_alignment=False
                else:
                    identity="NA"
            
            if process_alignment:
                    
                queryid=alignment_info[config.blast_query_index]
                evalue=alignment_info[config.blast_evalue_index]
                
                # remove the id of the alignment from the unaligned reads store
                unaligned_reads_store.remove_id(queryid)    
                
                # try converting evalue to float to check if it is a number
                try:
                    evalue=float(evalue)
                except ValueError:
                    logger.warning("E-value is not a number: %s", evalue)
                    evalue=1.0
                                
                # convert rapsearch evalue to blastm8 format if logged
                if log_evalue:
                    try:
                        evalue=math.pow(10.0, evalue)
                    except (ValueError, OverflowError):
                        logger.warning("Unable to convert rapsearch e-value: %s", evalue)
                        logger.warning("Traceback: \n" + traceback.format_exc())
                        evalue=1.0 
            
                # only store alignments with evalues less than threshold
                if evalue<1.0:
                    alignments.add_annotated(queryid, evalue, 
                        alignment_info[config.blast_reference_index])
                else:
                    logger.debug("Not including alignment based on large e-value: %s", evalue)
            
                aligned_ids+=[queryid]
        line=file_handle.readline()

    file_handle.close()

    # create unaligned file using list of remaining unaligned stored data
    file_handle_write=open(unaligned_file_fasta,"w")
    file_handle_write.write(unaligned_reads_store.get_fasta())
    file_handle_write.close()

    return unaligned_file_fasta

