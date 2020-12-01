"""
HUMAnN: nucleotide_search module
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
import logging
import traceback
import sys

from .. import utilities
from .. import config
from .. import store
from ..search import pick_frames
from ..search import blastx_coverage

# name global logging instance
logger=logging.getLogger(__name__)

def find_index(directory):
    """
    Search through the directory for the name of the bowtie2 index files
    Or if a file name is provided check it is a bowtie2 index
    """
    
    index=""
    bowtie2_extensions=config.bowtie2_index_ext_list+[config.bowtie2_large_index_ext]
    
    if not os.path.isdir(directory):
        # check if this is the bowtie2 index file
        if os.path.isfile(directory):
            # check for the bowtie2 extension
            for ext in bowtie2_extensions:
                if re.search(ext+"$",directory):
                    index=directory.replace(ext,"")
                    break
        else:
            # check if this is the basename of the bowtie2 index files
            small_index=directory+config.bowtie2_index_ext_list[0]
            large_index=directory+config.bowtie2_large_index_ext
            if os.path.isfile(small_index) or os.path.isfile(large_index):
                index=directory
    else:
        # Search through the files to find one with the bowtie2 extension
        for file in os.listdir(directory):
            # Look for an extension for a standard and large bowtie2 indexed database
            for ext in [config.bowtie2_index_ext_list[-1],config.bowtie2_large_index_ext]:
                if re.search(ext+"$",file):
                    index=os.path.join(directory,file.replace(ext,""))
                    break
            if index:
                break
    
    if not index:
        sys.exit("CRITICAL ERROR: Unable to find bowtie2 index files in directory: " + directory)
    
    return index
            

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
    tmpfile=utilities.unnamed_temp_file("bowtie2_stdout_")
    tmpfile2=utilities.unnamed_temp_file("bowtie2_stderr_")
    
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

def calculate_percent_identity(cigar_string, md_field):
    """
    Calculate the percent identity using the cigar string and md field from the sam file
    Returns the percent identity and the alignment length
    """
    
    match_numbers=re.compile("\d+")
    match_non_numbers=re.compile("\D+")
    
    # find the sets of numbers and identifers from the cigar string
    cigar_numbers=match_numbers.findall(cigar_string)
    cigar_identifiers=match_non_numbers.findall(cigar_string)

    # find the index for all of the match/mismatch/insert/delete
    match_mismatch_indel_index = []
    reference_length_index = []
    for index, cigar_identifier in enumerate(cigar_identifiers):
        if cigar_identifier in config.sam_cigar_match_mismatch_indel_identifiers:
            match_mismatch_indel_index.append(index)
        if cigar_identifier in config.sam_cigar_add_to_reference_identifiers:
            reference_length_index.append(index)
    
    # get reference length
    try:
        reference_length = int(sum([float(cigar_numbers[index]) for index in reference_length_index]))
    except (IndexError, ValueError):
        reference_length = 0

    # identify the total number of match/mismatch/indel
    try:
        match_mismatch_indel_count=sum([float(cigar_numbers[index]) for index in match_mismatch_indel_index])
    except (IndexError, ValueError):
        match_mismatch_indel_count=0.0
    
    # remove the tag from the md field
    md_field=md_field.split(config.sam_md_field_identifier)[-1]
    
    # find the sets of numbers from the md field
    md_field_numbers=match_numbers.findall(md_field)
    
    # sum the md field numbers to get the total number of matches
    try:
        matches=sum(int(n) for n in md_field_numbers)
    except ValueError:
        matches=0.0
    
    percent_identity=0.0
    if match_mismatch_indel_count > 0.0:
        percent_identity = 100.0 * ( matches / ( match_mismatch_indel_count * 1.0 ) )
        
    return percent_identity, match_mismatch_indel_count, reference_length 
    
def find_md_field(info):
    """
    Using the array of data from an alignment line, find the md field
    """
    
    # Search the data, starting with the first optional column to find the md field
    md_field=""
    for data in info[config.sam_start_optional_index:]:
        if re.match(config.sam_md_field_identifier,data):
            md_field=data
            break
        
    return md_field

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
    file_handle_read=open(sam_alignment_file, "rt")
    
    file_handle_write_aligned=open(reduced_aligned_reads_file, "w")

    # read through the file line by line
    # generate blast-like output file of alignments
    line = file_handle_read.readline()
    while line:
        # ignore headers ^@ 
        unaligned_read=False
        if not re.search("^@",line):
            info=line.split(config.sam_delimiter)
            # check flag to determine if unaligned
            if int(info[config.sam_flag_index]) & config.sam_unmapped_flag != 0:
                unaligned_read=True
            else:
                # convert the cigar string and md field to percent identity
                cigar_string=info[config.sam_cigar_index]
                md_field=find_md_field(info)
                identity, alignment_length, reference_length=calculate_percent_identity(cigar_string, md_field)
                
                query=info[config.sam_read_name_index]
                # write output to be blastm8-like
                new_info=[""] * config.blast_total_columns
                new_info[config.blast_query_index]=query
                new_info[config.blast_reference_index]=info[config.sam_reference_index]
                new_info[config.blast_subject_start_index]=str(info[config.sam_pos_index])
                new_info[config.blast_subject_end_index]=str(int(info[config.sam_pos_index])+reference_length)
                new_info[config.blast_evalue_index]="0"
                new_info[config.blast_identity_index]=str(identity)
                new_info[config.blast_aligned_length_index]=str(alignment_length)
                new_info[config.blast_query_start_index]="0"
                new_info[config.blast_query_end_index]=str(alignment_length-1)
                file_handle_write_aligned.write(config.blast_delimiter.join(new_info)+"\n")
        line=file_handle_read.readline()
                   
    file_handle_read.close()
    file_handle_write_aligned.close()

    # process alignments to determine genes for filtering
    allowed_genes = blastx_coverage.blastx_coverage(reduced_aligned_reads_file,
        config.nucleotide_subject_coverage_threshold, alignments, log_messages=True, apply_filter=True,
        nucleotide=True, query_coverage_threshold=config.nucleotide_query_coverage_threshold,
        identity_threshold = config.nucleotide_identity_threshold)

    file_handle_read=open(sam_alignment_file, "rt")
    file_handle_write_unaligned=open(unaligned_reads_file_fasta, "w")

    # read through the file line by line
    # capture alignments and also write out unaligned reads for next step in processing
    line = file_handle_read.readline()
    query_ids=set()
    no_frames_found_count=0
    small_identity_count=0
    filtered_genes_count=0
    query_coverage_count=0
    while line:
        # ignore headers ^@ 
        unaligned_read=False
        if not re.search("^@",line):
            info=line.split(config.sam_delimiter)
            query_ids.add(info[config.blast_query_index])
            # check flag to determine if unaligned
            if int(info[config.sam_flag_index]) & config.sam_unmapped_flag != 0:
                unaligned_read=True
            else:
                # convert the cigar string and md field to percent identity
                cigar_string=info[config.sam_cigar_index]
                md_field=find_md_field(info)
                identity, alignment_length, reference_length=calculate_percent_identity(cigar_string, md_field)

                # only store alignments with identity greater than threshold
                # and with genes included in the filtered list
                query=info[config.sam_read_name_index]

                gene_name, gene_length, bug = alignments.process_reference_annotation(
                info[config.sam_reference_index])

                if not gene_name in allowed_genes:
                    filtered_genes_count+=1
                    unaligned_read=True

                query_start_index=0
                query_stop_index=alignment_length-1
                if utilities.filter_based_on_query_coverage(alignment_length, query_start_index, query_stop_index, config.nucleotide_query_coverage_threshold):
                    query_coverage_count+=1
                    unaligned_read=True

                if identity > config.identity_threshold:
                    matches=identity/100.0*alignment_length
                    if not unaligned_read:
                        alignments.add_annotated(query,matches,info[config.sam_reference_index],
                            alignment_length)
                else:
                    small_identity_count+=1
                    unaligned_read=True

            if unaligned_read:
                annotated_sam_read_name=utilities.add_length_annotation(info[config.sam_read_name_index],
                                                                    len(info[config.sam_read_index]))
                file_handle_write_unaligned.write(">"+annotated_sam_read_name+"\n")
                file_handle_write_unaligned.write(info[config.sam_read_index]+"\n")
                
                # find the frames for the sequence and write to file
                if write_picked_frames:
                    picked_frames=pick_frames.pick_frames(info[config.sam_read_index])
                    if not picked_frames:
                        no_frames_found_count+=1
                    for frame in picked_frames:
                        file_handle_write_unaligned_frames.write(">"+
                            annotated_sam_read_name+"\n")
                        file_handle_write_unaligned_frames.write(frame+"\n")
                
                # store the unaligned reads data
                unaligned_reads_store.add(info[config.sam_read_name_index], 
                    info[config.sam_read_index])
                    
        line=file_handle_read.readline()

    if write_picked_frames:
        logger.debug("Total sequences without frames found: " + str(no_frames_found_count))
    logger.debug("Total nucleotide alignments not included based on filtered genes: " +
        str(filtered_genes_count))
    logger.debug("Total nucleotide alignments not included based on small percent identities: " +
        str(small_identity_count))
    logger.debug("Total nucleotide alignments not included based on query coverage threshold: " +
        str(query_coverage_count))
    
    file_handle_read.close()
    file_handle_write_unaligned.close()   
    file_handle_write_aligned.close()
    
    # set the total number of queries
    unaligned_reads_store.set_initial_read_count(len(query_ids))
    
    # set the unaligned reads file to read sequences from
    unaligned_reads_store.set_file(unaligned_reads_file_fasta)
    
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
