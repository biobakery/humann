#!/usr/bin/env python
"""
Identify initial list of bugs from user supplied fasta/fastq
"""

import os
import utilities

def run_alignment(metaphlan_dir, input, threads, debug_dir):
    """
    Runs metaphlan to identify initial list of bugs
    """
    
    exe="metaphlan2.py"
    opts="-t rel_ab"
    
    #determine input type as fastq or fasta
    input_type="multi" + utilities.fasta_or_fastq(input)
    
    # outfile name
    sample_name = os.path.splitext(os.path.basename(input))[0]
    bug_file = os.path.join(debug_dir, sample_name + "_bugs_list.tsv")
    bowtie2_out = os.path.join(debug_dir, sample_name + "_bowtie2_out.txt") 

 
    infiles=[input, os.path.join(metaphlan_dir, "db_v20/mpa_v20_m200.pkl")]
    
    # location of the index name to multiple files
    infiles_index=[os.path.join(metaphlan_dir, "db_v20/mpa_v20_m200")]
    
    outfiles=[bug_file, bowtie2_out]
    
    params=infiles[0] + " --bowtie2db " + infiles_index[0] + " --mpa_pkl " + \
        infiles[1] + " --input_type " + input_type + " -o " + outfiles[0] + \
        " --bowtie2out " + bowtie2_out
    
    if threads >1:
        params=params + " --nproc " + threads

    utilities.execute_software(exe, params + " " + opts, infiles, outfiles)
    
    return bug_file

def create_custom_database(chocophlan_dir, bug_file, debug_dir):
	"""
	Using ChocoPhlAn creates a custom database based on the bug_file
	"""

	custom_database=""

	return custom_database

