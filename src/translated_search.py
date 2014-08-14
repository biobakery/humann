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
    alignment_file=alignment_file_base + ".tsv"

    if config.translated_alignment_selected == "usearch":
        utilities.usearch_alignment(alignment_file,threads,identity_threshold,
            uniref, unaligned_reads_file_fastq)
    else:
        utilities.rapsearch_alignment(alignment_file,threads,uniref,
            unaligned_reads_file_fastq)

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

