"""
Utilities relating to executing third party software, file permissions,
and file formats
"""

import os, sys, subprocess, re

import config

def file_exists_readable(file):
    """
    Exit with error if file does not exist or is not readable
    """
    if not os.path.isfile(file):
        sys.exit("ERROR: Can not find file "+ file)
		
    if not os.access(file, os.R_OK):
        sys.exit("ERROR: Not able to read file " + file)

def add_exe_to_path(exe_dir):
    """ 
    Add path to executable to $PATH
    """
    os.environ["PATH"] += os.pathsep + exe_dir	

def find_exe_in_path(exe):
    """
    Check that an executable exists in $PATH
    """
    paths = os.environ["PATH"].split(os.pathsep)
    for path in paths:
        fullexe = os.path.join(path,exe)
        if os.path.exists(fullexe):
            if os.access(fullexe,os.X_OK):
                return True
    return False	

def return_exe_path(exe):
    """
    Return the location of the exe in $PATH
    """
    paths = os.environ["PATH"].split(os.pathsep)
    for path in paths:
        fullexe = os.path.join(path,exe)
        if os.path.exists(fullexe):
            if os.access(fullexe,os.X_OK):
                return path
    return "Error"	

def remove_if_exists(file):
    """
    If file exists, then remove
    """
    if os.path.isfile(file):
        os.unlink(file)

def remove_temp_file(file):
    """
    If file exists, then remove
    If debug mode, do not remove
    """
    if os.path.isfile(file):
        if not "debug" == config.run_mode:
            os.unlink(file)

def check_outfiles(outfiles):
    """
    If outfiles already_exist, then remove or bypass
    """

    bypass=False
    for file in outfiles:
        if os.path.isfile(file):
            if "debug" == config.run_mode:
                bypass=True
                break
            else:
                os.unlink(file)

    return bypass

def execute_command(exe, args, infiles, outfiles, stdout_file):
    """
    Execute third party software or shell command with files
    """
	
    # check that the executable can be found
    if not find_exe_in_path(exe):
        sys.exit("ERROR: Can not find executable " + exe)
	
    # check that the input files exist and are readable
    for file in infiles:
        file_exists_readable(file)
        
    # check if outfiles already exist
    bypass=check_outfiles(outfiles)

    # convert numbers to strings
    args=[str(i) for i in args]

    if not bypass:

        print "\n" + exe + " " + " ".join(args) + "\n"

        cmd=[exe]+args
	
        try:
            if stdout_file:
                p = subprocess.call(cmd, stdout=open(stdout_file,"w"))
            else:
                p = subprocess.call(cmd)
	
        except OSError as e:
            sys.exit("Error: Problem executing " + cmd + "\n" + e.strerror)
		
        # check that the output files exist and are readable
        for file in outfiles:
            file_exists_readable(file)
    else:
        print "Bypass: \n" + exe + " " + " ".join(args) + "\n"

def fasta_or_fastq(file):
    """
    Check to see if a file is of fasta or fastq format
	
    Fastq format short example:
    @SEQ_ID
    GATCTGG
    +
    !''****
	
    Fasta format short example:
    >SEQ_INFO
    GATCTGG
	
    Returns error if not of fasta or fastq format
    """
	
    format="error"
	
    # check file exists
    file_exists_readable(file)
	
    # read in first 4 lines of file to check format
    file_handle = open(file, "r")
	
    first_line = file_handle.readline()
    second_line = file_handle.readline()
	
    # check that second line is only nucleotides
    if re.search("^[aATtcCgGN]+$", second_line):
        # check first line to determine fasta or fastq format
        if re.search("^@",first_line):
            format="fastq"		
        if re.search("^>",first_line):
            format="fasta"
			
    file_handle.close()

    return format

def sam_aligned_reads(sam_alignment_file, output_type, unaligned_fastq_file, aligned_tsv_file):	
    """
    Create two files of aligned and unaligned reads in fasta/fastq and tsv format
    """
    sam_read_name_index=0
    sam_flag_index=1
    sam_reference_index=2
    sam_mapq_index=4
    sam_read_index=9
    sam_read_quality=10
	
    sam_unmapped_flag=0x4

    # check input and output files
    bypass=check_outfiles([unaligned_fastq_file, aligned_tsv_file])
  
    if not bypass:
        file_exists_readable(sam_alignment_file)

        file_handle_read=open(sam_alignment_file, "r")
        file_handle_write_unaligned=open(unaligned_fastq_file, "w")
        file_handle_write_aligned=open(aligned_tsv_file, "w")

        # read through the file line by line
        line = file_handle_read.readline()

        while line:
            # ignore headers ^@ 
            if not re.search("^@",line):
                info=line.split("\t")
                # check flag to determine if unaligned
                if int(info[sam_flag_index]) & sam_unmapped_flag != 0:
                    if output_type == "fastq":
                        file_handle_write_unaligned.write("@"+info[sam_read_name_index]+"\n")
                        file_handle_write_unaligned.write(info[sam_read_index]+"\n")
                        file_handle_write_unaligned.write("+\n")
                        file_handle_write_unaligned.write(info[sam_read_quality]+"\n")            
                    #default is fasta
                    else:
                        file_handle_write_unaligned.write(">"+info[sam_read_name_index]+"\n")
                        file_handle_write_unaligned.write(info[sam_read_index]+"\n")
                else:
                    newline=("\t").join([info[sam_read_name_index],info[sam_reference_index],
                        info[sam_mapq_index]])
                    file_handle_write_aligned.write(newline+"\n")
            line=file_handle_read.readline()

        file_handle_read.close()
        file_handle_write_unaligned.close()   
        file_handle_write_aligned.close()

        # remove the alignment file as it will be replaced by the two files created
        remove_temp_file(sam_alignment_file)
