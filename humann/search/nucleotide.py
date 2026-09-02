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

# The outcome of a single alignment for a read, ordered from the least to the most
# likely to be kept, so the outcome for a read with more than one alignment is the
# largest of the outcomes of its alignments
READ_NOT_MAPPED=0
READ_ALIGNMENT_FILTERED=1
READ_SUBJECT_NOT_COVERED=2
READ_ALIGNED=3

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

    args=["-f",custom_database,index_name,"--threads",config.threads]

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
        if line[0] != "@":
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

    # process alignments to determine the reference sequences for filtering
    # NOTE: these are the full reference annotations, not the gene family names, since
    # coverage is computed per reference sequence (see blastx_coverage)
    allowed_references = blastx_coverage.blastx_coverage(reduced_aligned_reads_file,
        config.nucleotide_subject_coverage_threshold, alignments, log_messages=True, apply_filter=True,
        nucleotide=True, query_coverage_threshold=config.nucleotide_query_coverage_threshold,
        identity_threshold = config.nucleotide_identity_threshold)

    file_handle_read=open(sam_alignment_file, "rt")

    # A read can have more than one alignment in the sam file (for example when the
    # aligner is run with "-k" or reports secondary alignments), so the aligned or
    # unaligned decision has to be made per read and not per alignment. Track the best
    # outcome seen for each read: a read is only unaligned when none of its alignments
    # is kept, and the status recorded for an unaligned read is the one closest to being
    # kept, which gives the breakdown of why reads are unaligned after this search.
    read_status={}
    primary_alignment_count=0
    no_frames_found_count=0
    small_identity_count=0
    filtered_genes_count=0
    query_coverage_count=0

    # read through the file line by line to capture the alignments
    line = file_handle_read.readline()
    while line:
        # ignore headers ^@
        if not re.search("^@",line):
            info=line.split(config.sam_delimiter)
            query=info[config.sam_read_name_index]
            flag=int(info[config.sam_flag_index])

            # every read has exactly one primary alignment, so counting them gives the
            # number of reads in the input even when reads have more than one alignment
            if flag & (config.sam_secondary_alignment_flag|config.sam_supplementary_alignment_flag) == 0:
                primary_alignment_count+=1

            # check flag to determine if unaligned
            if flag & config.sam_unmapped_flag != 0:
                status=READ_NOT_MAPPED
            else:
                # convert the cigar string and md field to percent identity
                cigar_string=info[config.sam_cigar_index]
                md_field=find_md_field(info)
                identity, alignment_length, reference_length=calculate_percent_identity(cigar_string, md_field)

                # only store alignments with identity greater than threshold
                # and with genes included in the filtered list
                status=READ_ALIGNED

                if not info[config.sam_reference_index] in allowed_references:
                    filtered_genes_count+=1
                    status=READ_SUBJECT_NOT_COVERED

                query_start_index=0
                query_stop_index=alignment_length-1
                if utilities.filter_based_on_query_coverage(alignment_length, query_start_index, query_stop_index, config.nucleotide_query_coverage_threshold):
                    query_coverage_count+=1
                    status=min(status,READ_ALIGNMENT_FILTERED)

                if identity > config.nucleotide_identity_threshold:
                    matches=identity/100.0*alignment_length
                    if status == READ_ALIGNED:
                        alignments.add_annotated(query,matches,info[config.sam_reference_index],
                            alignment_length)
                else:
                    small_identity_count+=1
                    status=min(status,READ_ALIGNMENT_FILTERED)

            read_status[query]=max(read_status.get(query,READ_NOT_MAPPED),status)

        line=file_handle_read.readline()

    file_handle_read.close()
    file_handle_write_aligned.close()

    # count the reads in each category before writing out the unaligned reads, since the
    # statuses are removed as they are written
    total_reads=len(read_status)

    # reads are tracked by name, so reads that share a name are counted as a single read
    # and only the sequence of the first is searched in the translated search
    if primary_alignment_count > total_reads:
        message="WARNING: The input includes "+str(primary_alignment_count-total_reads)+\
            " read(s) with a name that is not unique. Reads with the same name are "+\
            "counted as a single read, so the read counts reported will be smaller "+\
            "than the number of reads in the input."
        logger.warning(message)
        print("\n"+message+"\n")

    status_counts={status:0 for status in
        [READ_NOT_MAPPED, READ_ALIGNMENT_FILTERED, READ_SUBJECT_NOT_COVERED, READ_ALIGNED]}
    for status in read_status.values():
        status_counts[status]+=1
    total_unaligned=total_reads-status_counts[READ_ALIGNED]

    # read through the file a second time to write out each unaligned read once, for the
    # next step in processing
    file_handle_read=open(sam_alignment_file, "rt")
    file_handle_write_unaligned=open(unaligned_reads_file_fasta, "w")

    line = file_handle_read.readline()
    while line:
        # ignore headers ^@
        if not re.search("^@",line):
            info=line.split(config.sam_delimiter)
            query=info[config.sam_read_name_index]

            # remove the status so each read is only written once, which also releases
            # the statuses as the file is processed
            if read_status.pop(query,READ_ALIGNED) != READ_ALIGNED:
                annotated_sam_read_name=utilities.add_length_annotation(query,
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
                unaligned_reads_store.add(query, info[config.sam_read_index])

        line=file_handle_read.readline()

    if write_picked_frames:
        logger.debug("Total sequences without frames found: " + str(no_frames_found_count))
    logger.debug("Total nucleotide alignments not included based on filtered genes: " +
        str(filtered_genes_count))
    logger.debug("Total nucleotide alignments not included based on small percent identities: " +
        str(small_identity_count))
    logger.debug("Total nucleotide alignments not included based on query coverage threshold: " +
        str(query_coverage_count))

    # report the reads, rather than the alignments, that are unaligned after this search
    message="Total reads: "+str(total_reads)
    message+="\nTotal reads aligned after nucleotide search: "+str(status_counts[READ_ALIGNED])
    message+="\nTotal reads unaligned after nucleotide search: "+str(total_unaligned)
    message+="\n  never mapped to the database: "+str(status_counts[READ_NOT_MAPPED])
    message+="\n  mapped but all alignments filtered: "+str(status_counts[READ_ALIGNMENT_FILTERED])
    message+="\n  mapped but no reference sequence hit was well covered: "+str(status_counts[READ_SUBJECT_NOT_COVERED])
    logger.info(message)

    file_handle_read.close()
    file_handle_write_unaligned.close()

    # set the total number of queries
    unaligned_reads_store.set_initial_read_count(total_reads)

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
