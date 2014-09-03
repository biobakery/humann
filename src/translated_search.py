"""
Run alignment, find unused reads
"""

import os, tempfile
import utilities, config

def usearch_alignment(alignment_file,threads,identity_threshold,
        uniref, unaligned_reads_file_fastq):
    """
    Run usearch alignment with memory management
    Individual runs can be threaded
    """

    bypass=utilities.check_outfiles([alignment_file])

    exe="usearch"
    opts=config.usearch_opts

    args=["-id",identity_threshold]

    if threads > 1:
        args+=["-threads",threads]

    print "\nRunning muliple " + exe + " ........\n"

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
        for input_file in temp_in_files:
            for database in os.listdir(uniref):
                input_database=os.path.join(uniref,database)
                full_args=["-usearch_global",input_file]+args+["-db",input_database]

                # create temp output file
                file_out, temp_out_file=tempfile.mkstemp()
                os.close(file_out)
                temp_out_files.append(temp_out_file)

                full_args+=["-blast6out",temp_out_file]

                utilities.execute_command(exe,full_args,[input_database],[],"","")

        # merge the temp output files
        exe="cat"

        utilities.execute_command(exe,temp_out_files,temp_out_files,[alignment_file],
            alignment_file,"")

        # remove the temp files which have been merged
        for temp_file in temp_fasta_files + temp_in_files + temp_out_files:
            utilities.remove_file(temp_file)
    else:
        print "Bypass"


def rapsearch_alignment(alignment_file,threads, uniref,
        unaligned_reads_file_fastq):
    """
    Run rapsearch alignment on database formatted for rapsearch
    """

    bypass=utilities.check_outfiles([alignment_file])

    exe="rapsearch"
    opts=config.rapsearch_opts

    args=["-q",unaligned_reads_file_fastq,"-b",0]

    if threads > 1:
        args+=["-z",threads]

    print "\nRunning muliple " + exe + " ........\n"

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
                temp_out_files.append(temp_out_file)

                full_args+=["-o",temp_out_file]

                utilities.execute_command(exe,full_args,[input_database],[],"","")
        # merge the temp output files
        exe="cat"

        file_out, temp_merge_file=tempfile.mkstemp()
        os.close(file_out)

        utilities.execute_command(exe,temp_out_files,temp_out_files,[temp_merge_file],
            temp_merge_file,"")

        # reformat output file to fit blast6out format
        utilities.ConvertRapsearch2ToBlastM8Format(temp_merge_file,alignment_file)

        # remove the temp files which have been merged
        for temp_file in temp_out_files + [temp_merge_file]:
            utilities.remove_file(temp_file)
    else:
        print "Bypass"



def alignment(uniref, unaligned_reads_file_fasta, identity_threshold, 
    threads):
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
        usearch_alignment(alignment_file,threads,identity_threshold,
            uniref, input_fasta)
    else:
        rapsearch_alignment(alignment_file,threads,uniref,
            input_fasta)
        
    # Remove the temp fasta file if exists
    if temp_file:
        utilities.remove_file(temp_file)

    return alignment_file

def unaligned_reads(input_fastq, alignment_file):
    """
    Create a fasta/fastq file of the unaligned reads
    """

    #determine the index to use for the fastq/fasta file
    #use the same as that that was used by the user for the input file
    original_extension = os.path.splitext(os.path.basename(input_fastq))[1]
    unaligned_reads_file_fastq= utilities.name_temp_file(
        "_" + config.translated_alignment_selected + 
        config.translated_unaligned_reads_name_no_ext + original_extension)

    utilities.unaligned_reads_from_tsv(input_fastq, alignment_file, 
        unaligned_reads_file_fastq)

    return unaligned_reads_file_fastq

