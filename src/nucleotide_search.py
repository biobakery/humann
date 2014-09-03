"""
Index database, run alignment, find unused reads
"""

import os, re, math
import utilities, config, store

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
    print "\nRunning " + exe + " ........\n"

    args+=opts
    
    utilities.execute_command(exe,args,[custom_database],outfiles,"","")

    return index_name

def alignment(user_fastq, threads, index_name):
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



def unaligned_reads(input_fastq, sam_alignment_file, alignments):
    """ 
    Return file of just the unaligned reads
    Store the alignments
    """

    #determine the index to use for the fastq/fasta file
    #use the same as that that was used by the user for the input file
    original_extension = os.path.splitext(os.path.basename(input_fastq))[1]

    #for translated search create fasta unaligned reads file
    #even if original reads file is fastq
    unaligned_reads_file_fasta= utilities.name_temp_file(
        config.nucleotide_unaligned_reads_name_no_ext + config.fasta_extension)

    #name the reduced aligned reads file with tsv extension
    reduced_aligned_reads_file=utilities.name_temp_file(
        config.nucleotide_aligned_reads_name_tsv)

  
    utilities.file_exists_readable(sam_alignment_file)

    file_handle_read=open(sam_alignment_file, "r")
    file_handle_write_unaligned=open(unaligned_reads_file_fasta, "w")
    file_handle_write_aligned=open(reduced_aligned_reads_file, "w")

    # read through the file line by line
    line = file_handle_read.readline()

    while line:
        # ignore headers ^@ 
        if not re.search("^@",line):
            info=line.split(config.sam_delimiter)
            # check flag to determine if unaligned
            if int(info[config.sam_flag_index]) & config.sam_unmapped_flag != 0:
                file_handle_write_unaligned.write(">"+
                    info[config.sam_read_name_index]+"\n")
                file_handle_write_unaligned.write(info[config.sam_read_index]+"\n")
            else:
                # convert the e-value from global to local
                try:
                    evalue=math.pow(10.0, float(info[config.sam_mapq_index])/-10.0)
                except:
                    evalue=1.0 
                reference_info=info[config.sam_reference_index].split(
                    config.chocophlan_delimiter)
                # identify bug and gene families
                bug=reference_info[config.chocophlan_bug_index]
                uniref=reference_info[config.chocophlan_uniref_index]
                query=info[config.sam_read_name_index]
                newline=("\t").join([query,"",info[config.sam_reference_index],
                    "",str(evalue)])
                file_handle_write_aligned.write(newline+"\n")
                   
                # store the alignment data
                alignments.add(uniref,query,evalue,1,1,bug)
                    
        line=file_handle_read.readline()

    file_handle_read.close()
    file_handle_write_unaligned.close()   
    file_handle_write_aligned.close()

    # remove the alignment file as it will be replaced by the two files created
    #if not config.debug:
    #    remove_file(sam_alignment_file)

    return [ unaligned_reads_file_fasta, reduced_aligned_reads_file ]
