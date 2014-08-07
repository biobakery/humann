"""
Utilities relating to executing third party software, file permissions,
and file formats
"""

import os, sys, subprocess, re, shutil, tempfile

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

def check_software_version(exe,version_flag,
    version_required):
    """
    Determine if the software is of the correct version
    """

    if return_exe_path(exe) == "Error":
        sys.exit("ERROR: Can not find software " + exe)
    else:
        try:
            p_out = subprocess.check_output([exe,version_flag])
        except OSError as e:
            sys.exit("Error: Can not call software version\n" + e.strerror)
        
    if not re.search(version_required, p_out):
        sys.exit("Error: Please update " + exe + " from version " +
            p_out.rstrip('\n') + " to version " + version_required )


def remove_file(file):
    """
    If file exists, then remove
    """
    if os.path.isfile(file):
        os.unlink(file)
        if config.verbose:
            print "Remove file: " + file

def check_outfiles(outfiles):
    """
    If outfiles already_exist, then remove or bypass
    """
    bypass=[]
    for file in outfiles:
        if os.path.isfile(file):
            if config.debug and os.path.getsize(file) > 0:
                bypass.append(True)
            else:
                bypass.append(False)
        else:
            bypass.append(False)

    if False in bypass or not bypass:
        # remove any existing files
        for file in outfiles:
            remove_file(file)
        return False
    else:
        return True

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

        if config.verbose:
            print "\n" + exe + " " + " ".join(args) + "\n"

        cmd=[exe]+args
	
        if stdout_file:
            try:
                p = subprocess.call(cmd, stdout=open(stdout_file,"w"))
            except OSError as e:
                sys.exit("Error: Problem executing " + " ".join(cmd) + "\n" + e.strerror)
        else:
            try:
                p_out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
                if config.verbose:
                    print p_out
            except OSError as e:
                sys.exit("Error: Problem executing " + " ".join(cmd) + "\n" + e.strerror)

        # check that the output files exist and are readable
        for file in outfiles:
            file_exists_readable(file)
    else:
        if config.verbose:
            print "Bypass: \n" + exe + " " + " ".join(args) + "\n"
        else:
            print "Bypass\n"

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
	
    # check that second line is only nucleotides or amino acids
    if re.search("^[A-Z|a-z]+$", second_line):
        # check first line to determine fasta or fastq format
        if re.search("^@",first_line):
            format="fastq"		
        if re.search("^>",first_line):
            format="fasta"
			
    file_handle.close()

    return format

def unaligned_reads_from_sam(sam_alignment_file, output_type, unaligned_fastq_file, aligned_tsv_file):	
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
        if not config.debug:
            remove_file(sam_alignment_file)

def estimate_unaligned_reads(input_fastq, unaligned_fastq):
    """
    Calculate an estimate of the percent of reads unaligned
    """

    # check files exist and are readable
    file_exists_readable(input_fastq)
    file_exists_readable(unaligned_fastq)

    aligned_size=os.path.getsize(input_fastq)
    unaligned_size=os.path.getsize(unaligned_fastq)

    percent=int(unaligned_size/float(aligned_size) * 100)

    return str(percent)

def unaligned_reads_from_tsv(input_fastq, alignment_file_tsv, unaligned_file_fastq):
    """
    Create a fasta/fastq file of the unaligned reads
    """

    # check input and output files
    bypass=check_outfiles([unaligned_file_fastq])

    if not bypass:
        file_exists_readable(input_fastq)
        file_exists_readable(alignment_file_tsv)

        # read through the alignment file to identify ids
        # that correspond to aligned reads

        file_handle=open(alignment_file_tsv,"r")

        line=file_handle.readline()

        aligned_ids=[]
        while line:
            # expected format has the query ids in column 1
            aligned_ids+=[line.split("\t")[0]]
            line=file_handle.readline()

        file_handle.close()

        # create unaligned file using list of aligned ids
        file_handle_read=open(input_fastq,"r")
        file_handle_write=open(unaligned_file_fastq,"w")

        line=file_handle_read.readline()

        fastq_type=fasta_or_fastq(input_fastq)

        id_indicator="@"
        if fastq_type == "fasta":
            id_indicator=">"

        print_flag=False
        while line:
            # check for id line which will start with "@" or ">"
            if re.search("^"+id_indicator,line):
                id=line.strip(id_indicator+"\n")

                if not id in aligned_ids:
                    print_flag=True
                    file_handle_write.write(line)
                else:
                    print_flag=False
            else:
                if print_flag:
                    file_handle_write.write(line)
            line=file_handle_read.readline()

        file_handle_write.close()
        file_handle_read.close()

def remove_directory(dir):
    """
    Remove directory if exists
    """
    if os.path.isdir(dir):
        try:
            if config.verbose:
                print "Remove directory: " + dir
            shutil.rmtree(dir)
        except OSError: 
            sys.exit("ERROR: Unable to delete directory " + dir)

def break_up_fasta_file(fasta_file, max_seqs):
    """
    Break up a fasta file into smaller fasta files with max_seqs
    """

    # check file exists
    file_exists_readable(fasta_file)
    
    file_handle_read=open(fasta_file,"r")

    line=file_handle_read.readline()

    file_out, new_file=tempfile.mkstemp()
    fasta_files=[new_file]

    current_seq=0
    while line:
        if not re.search("^>",line):
            os.write(file_out,line)
        else:
            if current_seq > max_seqs:
                # close current file and open new
                os.close(file_out)
                file_out, new_file=tempfile.mkstemp()
                fasta_files+=[new_file]
                os.write(file_out, line)
                current_seq=1
            else:
                current_seq+=1
                os.write(file_out, line)
        line=file_handle_read.readline()

    os.close(file_out)
    file_handle_read.close()

    return fasta_files

def fastq_to_fasta(file):
    """
    Convert fastq file to fasta
	
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
	
    # check file exists
    file_exists_readable(file)
	
    file_handle_read = open(file, "r")
	
    line = file_handle_read.readline()
	
    file_out, new_file=tempfile.mkstemp()

    while line:
        if re.search("^@",line):
            os.write(file_out,line.replace("@",">",1))
            os.write(file_out,file_handle_read.readline())
        line=file_handle_read.readline()
  	
    os.close(file_out)	
    file_handle_read.close()

    return new_file

def usearch_alignment(alignment_file,threads,identity_threshold,
        uniref, unaligned_reads_file_fastq):
    """
    Run usearch alignment with memory management
    Individual runs can be threaded
    """

    bypass=check_outfiles([alignment_file])

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
        if fasta_or_fastq(unaligned_reads_file_fastq) == "fastq":
            unaligned_reads_file_fasta=fastq_to_fasta(unaligned_reads_file_fastq)
            temp_fasta_files.append(unaligned_reads_file_fasta)
        else:
            unaligned_reads_file_fasta=unaligned_reads_file_fastq

        #break up input file into smaller files for memory requirement
        temp_in_files=break_up_fasta_file(unaligned_reads_file_fasta,
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

                execute_command(exe,full_args,[input_database],[],"")

        # merge the temp output files
        exe="cat"

        execute_command(exe,temp_out_files,temp_out_files,[alignment_file],
            alignment_file)

        # remove the temp files which have been merged
        for temp_file in temp_fasta_files + temp_in_files + temp_out_files:
            remove_file(temp_file)
    else:
        print "Bypass"

