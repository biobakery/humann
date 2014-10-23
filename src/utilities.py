"""
HUMAnN2: utilities module
Utilities relating to third party software, file permissions, and file formats

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
import sys
import subprocess
import re
import shutil
import tempfile
import urllib
import tarfile
import multiprocessing
import logging
import traceback

import config

# name global logging instance
logger=logging.getLogger(__name__)

def unnamed_temp_file():
    """
    Return the full path to an unnamed temp file
    stored in the temp folder
    """
    
    file_out, new_file=tempfile.mkstemp(dir=config.unnamed_temp_dir)
    os.close(file_out)
    
    return(new_file)
    

def name_temp_file(file_name):
    """
    Return the full path to a new temp file 
    using the sample name and temp dir location
    """

    return os.path.join(config.temp_dir,
       config.file_basename + file_name) 

def file_exists_readable(file, raise_IOError=None):
    """
    Exit with error if file does not exist or is not readable
    Or raise an IOerror if selected
    """
    
    logger.debug("Check file exists and is readable: %s",file)
    
    if not os.path.isfile(file):
        message="Can not find file "+ file
        logger.critical(message)
        if raise_IOError:
            print("CRITICAL ERROR: " + message)   
            raise IOError 
        else:
            sys.exit("CRITICAL ERROR: " + message)
		
    if not os.access(file, os.R_OK):
        message="Not able to read file " + file
        logger.critical(message)
        if raise_IOError:
            print("CRITICAL ERROR: " + message)
            raise IOError
        else:
            sys.exit("CRITICAL ERROR: " + message)

def add_exe_to_path(exe_dir):
    """ 
    Add path to executable to $PATH
    """
    
    logger.debug("Add directory, %s, to path", exe_dir)
    
    os.environ["PATH"] += os.pathsep + exe_dir	

def find_exe_in_path(exe):
    """
    Check that an executable exists in $PATH
    """
    
    logger.debug("Find executable, %s, in path", exe)
    
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
    full_path=""
    for path in paths:
        fullexe = os.path.join(path,exe)
        if os.path.exists(fullexe):
            if os.access(fullexe,os.X_OK):
                full_path=path
    return full_path

def check_software_version(exe,version_flag,
    version_required):
    """
    Determine if the software is of the correct version
    """
    
    logger.debug("Check software, %s, for correct version, %s",exe,version_required)

    if not find_exe_in_path(exe):
        message="Can not find software " + exe
        logger.critical(message)
        sys.exit("CRITICAL ERROR: " + message)
    else:
        try:
            p_out = subprocess.check_output([exe,version_flag])
        except EnvironmentError:
            message="Can not call software version for " + exe
            logger.warning(message)
            print("WARNING: " + message + "\n")
        
    if not re.search(version_required, p_out):
        message=("Please update " + exe + " from version " +
            p_out.rstrip('\n') + " to version " + version_required)
        logger.critical(message) 
        sys.exit("CRITICAL ERROR: " + message)

def download_tar_and_extract(url, filename, folder):
    """
    Download the file at the url
    """
    
    message="Download URL:" + url
    logger.info(message)
    if config.verbose:
        print(message) 

    try:
        file, headers = urllib.urlretrieve(url,filename)
        message="Extracting:" + filename
        logger.info(message)
        if config.verbose:
            print(message)
        tarfile_handle=tarfile.open(filename,'r:gz')
        tarfile_handle.extractall(path=folder)
    except EnvironmentError:
        message="Unable to download and extract from URL: " + url
        logger.critical(message)
        logger.critical("Traceback: \n" + traceback.print_exc())
        sys.exit("CRITICAL ERROR: " + message)

def remove_file(file):
    """
    If file exists, then remove
    """
    
    try:
        if os.path.isfile(file):
            os.unlink(file)
            logger.debug("Remove file: %s", file)
    except OSError:
        message="Unable to remove file"
        logger.error(message)

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

def command_multiprocessing(threads, args, function=None):
    """
    Execute command in parallel, maximum of threads total
    """
    
    if not function:
        function=execute_command_args_convert
    
    logger.debug("Create %s processes for function %s", threads, function)
    
    pool=multiprocessing.Pool(threads)
    try:
        results=pool.map(function, args)
    except (EnvironmentError, ValueError, CalledProcessError):
        logger.critical("TRACEBACK: \n" + traceback.print_exc())
        message=("Error in one of the processes. See the log file for additional" + 
            " details including tracebacks.")
        logger.critical(message)
        sys.exit(message)
        
    return results
    
def execute_command_args_convert(args):
    """
    Convert the list of args to function arguments
    """
    
    return execute_command(*args)

def execute_command(exe, args, infiles, outfiles, stdout_file=None, stdin_file=None):
    """
    Execute third party software or shell command with files
    """
	
    # check that the executable can be found
    if not find_exe_in_path(exe):
        message="Can not find executable " + exe
        print("CRITICAL ERROR: " + message)
        logger.critical(message)
        raise EnvironmentError
	
    # check that the input files exist and are readable
    for file in infiles:
        logger.debug("Check input file exists and is readable: %s",file)
        file_exists_readable(file, raise_IOError=True)
        
    # check if outfiles already exist
    bypass=check_outfiles(outfiles)

    # convert numbers to strings
    args=[str(i) for i in args]

    if not bypass:
        
        cmd=[exe]+args

        message=" ".join(cmd)
        logger.info("Execute command: "+ message)
        if config.verbose:
            print("\n"+message+"\n")
	
        try:
            if stdout_file:
                if stdin_file:
                    p = subprocess.call(cmd, stdin=open(stdin_file,"r"),stdout=open(stdout_file,"w"))
                else:
                    p = subprocess.call(cmd, stdout=open(stdout_file,"w"))
            else:
                p_out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
                logger.debug(p_out)            
        except (EnvironmentError, CalledProcessError):
            message="Problem executing " + " ".join(cmd) + "\n"
            logger.critical(message)
            logger.critical("TRACEBACK: \n" + traceback.print_exc())
            print("CRITICAL ERROR: " + message)
            raise

        # check that the output files exist and are readable
        for file in outfiles:
            logger.debug("Check output file exists and is readable: %s",file)
            file_exists_readable(file, raise_IOError=True)
    
    else:
        if config.verbose:
            print("Bypass: \n" + exe + " " + " ".join(args) + "\n")
        else:
            print("Bypass\n")

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

def remove_directory(dir):
    """
    Remove directory if exists
    """
    if os.path.isdir(dir):
        try:
            shutil.rmtree(dir)
            logger.debug("Remove directory: " + dir)
        except EnvironmentError: 
            logger.error("Unable to remove directory: " + dir)
    else:
        logger.debug("Request to remove directory that does not exist: " + dir)

def break_up_fasta_file(fasta_file, max_seqs):
    """
    Break up a fasta file into smaller fasta files with max_seqs
    """

    # check file exists
    file_exists_readable(fasta_file)
    
    file_handle_read=open(fasta_file,"r")

    line=file_handle_read.readline()

    new_file=unnamed_temp_file()
    file_out=open(new_file,"w")
    
    fasta_files=[new_file]

    current_seq=0
    while line:
        if not re.search("^>",line):
            file_out.write(line)
        else:
            if current_seq == max_seqs:
                # close current file and open new
                file_out.close()
                     
                new_file=unnamed_temp_file()
                file_out=open(new_file,"w")
                
                fasta_files+=[new_file]
                file_out.write(line)
                current_seq=1
            else:
                current_seq+=1
                file_out.write(line)
        line=file_handle_read.readline()

    file_out.close()
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
	
    new_file=unnamed_temp_file()
    file_out=open(new_file,"w")

    while line:
        if re.search("^@",line):
            file_out.write(line.replace("@",">",1))
            line=file_handle_read.readline()
            while line:
                if re.search("^\\+",line):
                    break
                else:
                    file_out.write(line)
                line=file_handle_read.readline()
        line=file_handle_read.readline()
  	
    file_out.close()	
    file_handle_read.close()

    return new_file	
		
def install_minpath():
    """
    Download and install the minpath software
    """
    
    # Download the minpath software v1.2
    # Check to see if already downloaded
    
    fullpath_scripts=os.path.dirname(os.path.realpath(__file__))
    minpath_exe=os.path.join(fullpath_scripts,config.minpath_folder,
        config.minpath_script)

    # Check for write permission to the target folder
    if not os.access(fullpath_scripts, os.W_OK):
        sys.exit("ERROR: The directory set to install MinPath is not writeable: "+
            fullpath_scripts + " . Please modify the permissions.")

    if not os.path.isfile(minpath_exe):
        if config.verbose:
            print("Installing minpath ... ")
        download_tar_and_extract(config.minpath_url, 
            os.path.join(fullpath_scripts, config.minpath_file),fullpath_scripts)
    
def tsv_to_biom(tsv_file, biom_file):
    """
    Convert from a tsv to biom file
    """

    exe="biom"
    args=["convert","-i",tsv_file,"-o",biom_file,"--table-type","OTU table","--to-json"]
    
    # Remove output file if already exists
    if os.path.isfile(biom_file):
        remove_file(biom_file)
    
    execute_command(exe, args, [tsv_file], [biom_file])
    
    
    
		
