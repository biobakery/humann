"""
Utilities relating to executing third party software, file permissions,
and file formats
"""

import os, sys, subprocess, re, shutil, tempfile
import fileinput
import math
from decimal import *

import config

def name_temp_file(file_name):
    """
    Return the full path to a new temp file 
    using the sample name and temp dir location
    """

    return os.path.join(config.temp_dir,
       config.file_basename + file_name) 

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

def execute_command(exe, args, infiles, outfiles, stdout_file, stdin_file):
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
            if stdin_file:
                try:
                    p = subprocess.call(cmd, stdin=open(stdin_file,"r"),stdout=open(stdout_file,"w"))
                except OSError as e:
                    sys.exit("Error: Problem executing " + " ".join(cmd) + "\n" + e.strerror)
            else:
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
                    newline=("\t").join([info[sam_read_name_index],"",info[sam_reference_index],
                        "",info[sam_mapq_index]])
                    file_handle_write_aligned.write(newline+"\n")
            line=file_handle_read.readline()

        file_handle_read.close()
        file_handle_write_unaligned.close()   
        file_handle_write_aligned.close()

        # remove the alignment file as it will be replaced by the two files created
        #if not config.debug:
        #    remove_file(sam_alignment_file)

def count_reads(file):
    """
    Count the total number of reads in a file
    """

    file_handle_read=open(file,"r")

    line=file_handle_read.readline()

    if fasta_or_fastq(file) == "fastq":
        read_token="@"
    else:
        read_token=">"

    sequence_count=0
    while line:
        if re.search("^"+read_token,line):
            sequence_count+=1
        line=file_handle_read.readline()

    file_handle_read.close()

    return sequence_count


def estimate_unaligned_reads(input_fastq, unaligned_fastq):
    """
    Calculate an estimate of the percent of reads unaligned
    """

    # check files exist and are readable
    file_exists_readable(input_fastq)
    file_exists_readable(unaligned_fastq)

    percent=int(count_reads(unaligned_fastq)/float(count_reads(input_fastq)) * 100)

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
            if current_seq == max_seqs:
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

	
#********************************************************************************************
#*   ConvertRapsearch2ToBlastM8Format                                                       *
#*  <<------------------------------------------------------------------------------->>     *    
#*  The following three subroutines:                                                        *
#*  1.  SciStr                                                                              *
#*  2. ProcessRecord(InputLine)                                                             *
#*  3. ConvertRapsearch2ToBlastM8Format(InputFileName,OutputFileName)                       *
#*  Work together.                                                                          *
#*  Their purpose is to take the output of rapsearch2 and convert it to blast m8 format     *
#*  Only two items are modified: rapsearch2 returns log of eValue so we convert it to       *
#*     eValue by calculating 10^^(log(eValue)) and formatting it to a string that contains  *
#*     its value in Scientific format                                                       *  
#*                                                                                          *
#*  The subroutines are invoked as follows:                                                 *
#*  RC = ConvertRapsearch2ToBlastM8Format(InputFileName,OutputFileName)                     *
#*  Where:                                                                                  *
#*  InputFileName = Rapsearch2 output file name  (Input to these modules)                   *   
#*  OutputFileName = Output of these modules:  Compatible to Blast m8 format                *
#*                                                                                          *
#*  If the conversion of the log(eValue) fails - we post spaces for that field              *                                       
#*                                                                                          *
#********************************************************************************************


#********************************************************************************************
#*    Convert eValue into a scientific string to comply with blast m8 format                *    
#*    the routine to convert rapsearch2 file format to blast m8 format                      *
#******************************************************************************************** 
def SciStr(dec):
	ScientificString = '%.1E' % Decimal(math.pow(10, dec))
	return  ScientificString 
 
#********************************************************************************************
#*    Process the Input Record - this is part of                                            *    
#*    the routine to convert rapsearch2 file format to blast m8 format                      *
#******************************************************************************************** 
def ProcessRecord(InputLine):
	SplitLine = InputLine.split('\t')
 
	try:
		SplitLine[10]  =  SciStr(float(SplitLine[10]))	#Convert the eValue into Scientific format
	except:
		SplitLine[10] = " "								#If it did not convert - set up to space
	ProcessedRecord =  '\t'.join(SplitLine)				#Rebuild the record
	return ProcessedRecord   
	
#********************************************************************************************
#*    Routine to convert rapsearch2 file format to blast m8 format                          *
#********************************************************************************************
def  ConvertRapsearch2ToBlastM8Format(	InputFileName,OutputFileName):
	OutputFile = open(OutputFileName,'w')
	
	for InputLine in fileinput.input(InputFileName):
		if  InputLine.startswith('#'):			# If it is one of the header lines
			if  InputLine.startswith('# Fields: Query'):    # If the record starts with #Fields 
				InputLine = InputLine.replace("log(e-value)", "e-value", 1)   #  Replace log(eValue) with eValue 
			OutputFile.write(InputLine) 		#     just write it
		else:
			ProcessedRecord = ProcessRecord(InputLine)	# Process the Input Record
			OutputFile.write(ProcessedRecord ) 	 		#    and write it

	OutputFile.close() 						# Finished - close the output file
	return 0
	

#********************************************************************************************
#*   End of ConvertRapsearch2ToBlastM8Format Subroutines                                    *
#*  <<------------------------------------------------------------------------------->>     *    
#********************************************************************************************		
		
		
		
