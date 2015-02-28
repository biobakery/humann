"""
HUMAnN2: nucleotide_search module
Index database, run alignment, find unused reads

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
import math
import logging
import traceback

from humann2 import utilities
from humann2 import config
from humann2 import store
from humann2.search import pick_frames

# name global logging instance
logger=logging.getLogger(__name__)

def index(custom_database):
    """
    Index database and run alignment with bowtie2
    """
    # name the index
    index_name = utilities.name_temp_file( 
        config.bowtie2_index_name)
  
    exe="bowtie2-build"
    opts=config.bowtie2_build_opts

    args=["-f",custom_database,index_name]

    outfiles=[index_name + ext for ext in config.bowtie2_index_ext_list] 

    # if custom_database is large (>4G) then use the --large-index flag
    if os.path.getsize(custom_database) > config.bowtie2_large_index_threshold:
        args+=["--large-index"]
        outfiles=[index_name + config.bowtie2_large_index_ext]
        
    # index the database
    message="Running " + exe + " ........"
    logger.info(message)
    print("\n"+message+"\n")

    args+=opts
    
    # create temp file for stdout and stderr
    tmpfile=utilities.unnamed_temp_file()
    tmpfile2=utilities.unnamed_temp_file()
    
    utilities.execute_command(exe,args,[custom_database],outfiles,
        stdout_file=tmpfile, stderr_file=tmpfile2)

    return index_name

def alignment(user_fastq, index_name):
    """
    Run alignment with bowtie2
    """
    
    # name the alignment file
    alignment_file = utilities.name_temp_file(
        config.chocophlan_alignment_name)

    # align user input to database
    exe="bowtie2"
    opts=config.bowtie2_align_opts

    #determine input type as fastq or fasta
    input_type = utilities.fasta_or_fastq(user_fastq)

    logger.debug("Nucleotide input file is of type: %s", input_type)

    #determine input type flag
    #default flag to fastq
    input_type_flag = "-q"
    if input_type == "fasta":
        input_type_flag="-f"

    args=[input_type_flag,"-x",index_name,"-U",user_fastq,"-S",alignment_file]
    
    #add threads
    if config.threads > 1:
        args+=["-p",config.threads]

    # run the bowtie2 alignment
    message="Running " + exe + " ........"
    print("\n"+message+"\n")
    
    args+=opts

    utilities.execute_command(exe,args,[user_fastq],[alignment_file])

    return alignment_file



def unaligned_reads(sam_alignment_file, alignments, unaligned_reads_store, keep_sam=None):
    """ 
    Return file and data structure of the unaligned reads 
    Store the alignments and return
    """

    #for translated search create fasta unaligned reads file
    #even if original reads file is fastq
    unaligned_reads_file_fasta= utilities.name_temp_file(
        config.nucleotide_unaligned_reads_name_no_ext + config.fasta_extension)
    
    # if set to run frame picker, create named temp file
    write_picked_frames=False
    if config.pick_frames_toggle == "on":
        logger.debug("Creating picked frames file")
        unaligned_reads_file_picked_frames_fasta = utilities.name_temp_file( 
            config.nucleotide_unaligned_reads_picked_frames_name_no_ext + 
            config.fasta_extension)
        file_handle_write_unaligned_frames=open(unaligned_reads_file_picked_frames_fasta, "w")
        write_picked_frames=True

    #name the reduced aligned reads file with tsv extension
    reduced_aligned_reads_file=utilities.name_temp_file(
        config.nucleotide_aligned_reads_name_tsv)

  
    utilities.file_exists_readable(sam_alignment_file)
    file_handle_read=open(sam_alignment_file, "r")
    
    file_handle_write_unaligned=open(unaligned_reads_file_fasta, "w")
    file_handle_write_aligned=open(reduced_aligned_reads_file, "w")

    # read through the file line by line
    line = file_handle_read.readline()
    query_ids={}
    no_frames_found_count=0
    large_evalue_count=0
    while line:
        # ignore headers ^@ 
        if not re.search("^@",line):
            info=line.split(config.sam_delimiter)
            query_ids[info[config.blast_query_index]]=1
            # check flag to determine if unaligned
            if int(info[config.sam_flag_index]) & config.sam_unmapped_flag != 0:
                file_handle_write_unaligned.write(">"+
                    info[config.sam_read_name_index]+"\n")
                file_handle_write_unaligned.write(info[config.sam_read_index]+"\n")
                
                # find the frames for the sequence and write to file
                if write_picked_frames:
                    picked_frames=pick_frames.pick_frames(info[config.sam_read_index])
                    if not picked_frames:
                        no_frames_found_count+=1
                    for frame in picked_frames:
                        file_handle_write_unaligned_frames.write(">"+
                            info[config.sam_read_name_index]+"\n")
                        file_handle_write_unaligned_frames.write(frame+"\n")
                
                # store the unaligned reads data
                unaligned_reads_store.add(info[config.sam_read_name_index], 
                    info[config.sam_read_index])
            else:
                # convert the e-value from global to local
                try:
                    evalue=math.pow(10.0, float(info[config.sam_mapq_index])/-10.0)
                except ValueError:
                    logger.warning("Unable to convert bowtie2 e-value: %s", 
                        info[config.sam_mapq_index])
                    logger.warning("Traceback: \n" + traceback.format_exc())
                    evalue=1.0 
                
                query=info[config.sam_read_name_index]
                # write output to be blastm8-like
                new_info=[""] * config.blast_total_columns
                new_info[config.blast_query_index]=query
                new_info[config.blast_reference_index]=info[config.sam_reference_index]
                new_info[config.blast_evalue_index]=str(evalue)
                file_handle_write_aligned.write(config.blast_delimiter.join(new_info)+"\n")
                   
                # only store alignments with evalues less than threshold
                if evalue<config.evalue_threshold:
                    alignments.add_annotated(query,evalue,info[config.sam_reference_index])
                else:
                    large_evalue_count+=1
                    
        line=file_handle_read.readline()

    if write_picked_frames:
        logger.debug("Total sequences without frames found: " + str(no_frames_found_count))
    logger.debug("Total nucleotide alignments not included based on large e-value: " +
        str(large_evalue_count))
    
    file_handle_read.close()
    file_handle_write_unaligned.close()   
    file_handle_write_aligned.close()
    
    # set the total number of queries
    unaligned_reads_store.set_initial_read_count(len(query_ids))
    
    if write_picked_frames:
        file_handle_write_unaligned_frames.close()

    # remove the alignment file as it will be replaced by the two files created
    if not config.resume:
        if keep_sam:
            logger.debug("Keeping sam file")
        else:
            logger.debug("Remove sam file")
            utilities.remove_file(sam_alignment_file)

    # return the picked frames file if written
    return_list=[unaligned_reads_file_fasta, reduced_aligned_reads_file]
    if write_picked_frames:
        return_list=[unaligned_reads_file_picked_frames_fasta, reduced_aligned_reads_file]

    return return_list
