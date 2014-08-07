"""
Index database, run alignment, find unused reads
"""

import os
import utilities, config

def alignment(custom_database, user_fastq, temp_dir, threads):
    """
    Index database and run alignment with bowtie2
    """
    # name the index
    custom_database_name = os.path.splitext(os.path.basename(custom_database))[0]
    index_name = os.path.join(temp_dir, 
        custom_database_name + config.bowtie2_index_name)
  
    exe="bowtie2-build"
    opts=config.bowtie2_build_opts

    args=["-f",custom_database,index_name]

    outfiles=[index_name + ext for ext in config.bowtie2_index_ext_list] 

    # if custom_database is large (>4G) then use the --large-index flag
    if os.path.getsize(custom_database) > config.bowtie2_large_index_threshold:
        args+=["--large-index"]
        outfiles=[index_name + config.bowtie2_large_index_ext]
        
    # index the database
    print "\nRunning " + exe + " ........\n"

    args+=opts
    
    utilities.execute_command(exe,args,[custom_database],outfiles,"","")

    # name the alignment file
    sample_name = os.path.splitext(os.path.basename(user_fastq))[0]
    alignment_file = os.path.join(temp_dir, 
        custom_database_name + config.chocophlan_alignment_name)


    # align user input to database
    exe="bowtie2"
    opts=config.bowtie2_align_opts

    #determine input type as fastq or fasta
    input_type = utilities.fasta_or_fastq(user_fastq)

    #determine input type flag
    #default flag to fastq
    input_type_flag = "-q"
    if input_type == "fasta":
        input_type_flag="-f"

    args=[input_type_flag,"-x",index_name,"-U",user_fastq,"-S",alignment_file]
    
    #add threads
    if threads > 1:
        args+=["-p",threads]

    # run the bowtie2 alignment
    print "\nRunning " + exe + " ........\n"
    
    args+=opts

    utilities.execute_command(exe,args,[user_fastq],[alignment_file],"","")

    return alignment_file

def unaligned_reads(input_fastq, alignment_file, temp_dir):
    """ 
    Return file of just the unaligned reads
    """

    #name the index
    alignment_file_name = os.path.splitext(os.path.basename(alignment_file))[0]
    unaligned_reads_file_fastq = os.path.join(temp_dir, 
        alignment_file_name + config.unaligned_reads_name_no_ext)

    #determine the index to use for the fastq/fasta file
    #use the same as that that was used by the user for the input file
    original_extension = os.path.splitext(os.path.basename(input_fastq))[1]
    unaligned_reads_file_fastq+= original_extension

    #name the reduced aligned reads file with tsv extension
    reduced_aligned_reads_file=os.path.join(temp_dir,
        alignment_file_name + config.aligned_reads_name_tsv)

    #create two files of aligned and unaligned reads
    utilities.unaligned_reads_from_sam(alignment_file, utilities.fasta_or_fastq(input_fastq), 
        unaligned_reads_file_fastq, reduced_aligned_reads_file)

    return [ unaligned_reads_file_fastq, reduced_aligned_reads_file ]
