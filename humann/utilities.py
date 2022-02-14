"""
HUMAnN: utilities module
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

# try to import urllib.request.urlretrieve for python3
try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve

import tarfile
import logging
import traceback
import gzip
import shutil
import threading

# try to import the python2 module Queue
# if unable to import, try to import the python3 module queue
try:
    import Queue as queue
except ImportError:
    import queue

import datetime
import time
import math

from . import config
from .search import pick_frames

# name global logging instance
logger=logging.getLogger(__name__)

def determine_file_format(file):
    """
    Return the type of file based on the format
    
    Formats recognized are: fasta, fastq, sam, and blastm8
    
    Fastq format short example:
    @SEQ_ID
    GATCTGG
    +
    !''****
    
    Fasta format short example:
    >SEQ_INFO
    GATCTGG
    
    Sam format short example:
    @HD_optional_headers
    Qname    flag    rname    pos    mapq    cigar    rnext    pnext    tlen   \
    seq    qual optional_tags
    
    Blastm8 format short example:
    # optional comment line
    Qname    Sname    id    len    mismatched    gap    Qstart    Qend    Sstart \
    Send    e-value    bit_score
    
    Bam, biom, and gzipped files are recognized based on their extensions
    
    Error is return if the file is not of a known format
    """
    
    # check file exists
    file_exists_readable(file)
    
    format=""
    
    # read in the first 2 lines of the file to check format  
    gzipped=False
    try:      
        # check for gzipped files
        if file.endswith(".gz"):
            file_handle = gzip.open(file, "rt")
            gzipped=True
        else:
            file_handle = open(file, "rt")
        
        first_line = file_handle.readline().rstrip()
        while re.search("^#",first_line):
            first_line = file_handle.readline().rstrip()
        second_line = file_handle.readline().rstrip()
    except (EnvironmentError, UnicodeDecodeError):
        # if unable to open and read the file, set the format to unknown
        format="unknown"
        first_line = ""
        second_line = ""  
    finally:
        file_handle.close()

    # check for a bam file
    if file.endswith(".bam"):
        format="bam"
    elif file.endswith(".biom"):
        format="biom"
    # check that second line is only nucleotides or amino acids
    elif re.search("^[A-Z|a-z]+$", second_line):
        # check first line to determine fasta or fastq format
        if re.search("^@",first_line):
            format="fastq"        
        if re.search("^>",first_line):
            format="fasta"
    else:
        # check for sam format with header on first line
        if re.search("^@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$",
            first_line) or re.search("^@CO\t.*",first_line):
            format="sam"
        
        # check for formats that have tabs in the first line
        if re.search(("\t"),first_line) and not format:
            data=first_line.split("\t")
            if len(data)>config.sam_read_quality:
                # check for sam format
                if re.search("\*|[A-Za-z=.]+",data[config.sam_read_index]):
                    format="sam"
                # check for standard blastm8 format (blast, usearch, rapsearch)
                # this will have only numeric values in the column that for
                # a sam file would be the read sequence
                elif re.search("^[0-9]+$",data[config.sam_read_index]):
                    format="blastm8"
                else:
                    try:
                        # the location of the blastm8 evalue contains quality
                        # score information in the sam format
                        evalue=float(data[config.blast_evalue_index])
                        format="blastm8"
                    except ValueError:
                        if re.search("[!-~]+",data[config.blast_evalue_index]):
                            format="sam"
            # check for gene table for a single sample
            elif len(data)==config.gene_table_total_columns:
                # check that the data column is numerical
                if re.search("^[0-9E\-.]+$",data[config.gene_table_value_index]):
                    format="genetable"
    if not format:
        format="unknown"
    elif gzipped:
        format+=".gz"
                
    message="File ( " + file + " ) is of format:  " + format
    if config.verbose:
        print(message+"\n")
    
    logger.info(message)       
                
    return format

def space_in_identifier(file):
    """ Check if there are spaces in the fasta/fastq identifier by
    checking the first line of the file """
    
    space_found = False
    try:
        file_handle = open(file, "rt")
        line = file_handle.readline()
        if " " in line:
            space_found = True
        file_handle.close()
    except (EnvironmentError, UnicodeDecodeError):
        logger.info("Unable to check file for identifier format")
    
    return space_found

def remove_spaces_from_file(file):
    """ Remove any spaces in the file, creating a new file of the output """
    
    # create an unnamed temp file
    new_file=unnamed_temp_file()
    
    try:
        file_handle_read = open(file, "rt")
        file_handle_write = open(new_file, "wt")
        for line in file_handle_read:
            file_handle_write.write(line.replace(" ",""))
        file_handle_read.close()
        file_handle_write.close()
    except (EnvironmentError, UnicodeDecodeError):
        logger.info("Unable to write new file after removing spaces in identifier")
        new_file=""
        
    return new_file

def bam_to_sam(bam_file):
    """
    Convert from a bam to sam file
    """
    
    # create a unnamed temp file
    new_file=unnamed_temp_file()

    exe="samtools"
    args=["view","-h",bam_file,"-o",new_file]
    
    message="Converting bam file to sam format ..."
    print(message)
    logger.info(message)
    
    execute_command(exe, args, [bam_file], [new_file])
    
    return new_file

def gunzip_file(gzip_file):
    """
    Return a new copy of the file that is not gzipped
    The new file will be placed in the unnamed temp folder
    """
    
    message="Decompressing gzipped file ..."
    print(message+"\n")
    logger.info(message)    
    
    try:
        file_handle_gzip=gzip.open(gzip_file,"rt")
        
        # create a unnamed temp file
        new_file=unnamed_temp_file()
        
        # write the gunzipped file
        file_handle=open(new_file,"wt")
        shutil.copyfileobj(file_handle_gzip, file_handle)
        
    except EnvironmentError:
        print("Critical Error: Unable to unzip input file: " + gzip_file)
        new_file=""
    finally:
        file_handle.close()
        file_handle_gzip.close()
        
    return new_file

def double_sort(pathways_dictionary):
    """
    Return the keys to a dictionary sorted with top values first
    then for duplicate values sorted alphabetically by key
    """

    sorted_keys=[]
    prior_value=""
    store=[]
    for pathway in sorted(pathways_dictionary, key=pathways_dictionary.get, reverse=True):
        if prior_value == pathways_dictionary[pathway]:
            if not store:
                store.append(sorted_keys.pop())
            store.append(pathway)
        else:
            if store:
                sorted_keys+=sorted(store)
                store=[]
            prior_value=pathways_dictionary[pathway]
            sorted_keys.append(pathway)

    if store:
        sorted_keys+=sorted(store)
    return sorted_keys

def unnamed_temp_file(prefix=None):
    """
    Return the full path to an unnamed temp file
    stored in the unnamed temp folder
    """
    
    if not prefix:
        prefix="tmp"
        
    try:
        file_out, new_file=tempfile.mkstemp(dir=config.unnamed_temp_dir,prefix=prefix)
        os.close(file_out)
    except EnvironmentError:
        sys.exit("ERROR: Unable to create temp file in directory: " + config.unnamed_temp_dir)
    
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
    
    os.environ["PATH"] = exe_dir + os.pathsep + os.environ["PATH"]
    
def add_directory_to_pythonpath(dir):
    """
    Prepend the directory to the PYTHONPATH if not currently included
    """
    
    # check if the directory is already at the beginning of the list
    if not os.path.samefile(sys.path[0],dir):
        # add directory to beginning of path list
        logger.debug("Add directory, %s, to pythonpath", dir)
        sys.path.insert(1,dir)
        
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
        if os.path.exists(fullexe) and os.path.isfile(fullexe):
            if os.access(fullexe,os.X_OK):
                return path
    return ""

def return_module_path(module):
    """
    Return the full path to the python module
    """
    
    # if this is the full path to the module then return the parent directory
    if os.path.isabs(module):
        if os.path.exists(module):
            return os.path.abspath(os.path.join(module,os.pardir))
    else:
    # search for the path to the module
        for path in sys.path:
            full_module = os.path.join(path,module)
            if os.path.exists(full_module):
                return path
    return ""

def check_software_version(exe,version,warning=False):
    """
    Determine if the software is of the correct version
    """
    
    version_required=str(version["major"])+"."+str(version["minor"])
    if "second minor" in version:
        version_required+="."+str(version["second minor"])
        
    logger.debug("Check software, %s, for required version, %s",exe,version_required)

    if not find_exe_in_path(exe):
        message="Can not find software " + exe
        logger.critical(message)
        sys.exit("CRITICAL ERROR: " + message)

    try:
        process = subprocess.Popen([exe,version["flag"]],stderr=subprocess.STDOUT,
            stdout=subprocess.PIPE, universal_newlines=True)
        process_out=process.communicate()[0]
    except EnvironmentError:
        message="Error trying to call software version"
        logger.debug(message)
        
    try:
        # find the version string and remove a "v" and ":" if present
        version_line=list(filter(lambda x: x, process_out.split("\n")))[version["line"]]
        version_line_split=version_line.split(" ")
        version_string=version_line_split[version["column"]]
        current_version=version_string.replace("v","").replace(":","").split(".")
        current_major_version=int(current_version[0])
        if current_version[1].isdecimal():
            current_minor_version=int(current_version[1])
        elif len(current_version)> 2:
            current_minor_version=int(current_version[2])
        else:
            current_minor_version=int(current_version[1])
        if "second minor" in version:
            current_second_minor_version=int(current_version[2])
        current_version=str(current_major_version)+"."+str(current_minor_version)
        if "second minor" in version:
            current_version+="."+str(current_second_minor_version)
    except (NameError,KeyError,ValueError,IndexError):
        message="Can not call software version for " + exe
        if warning:
            logger.warning(message)
            print("WARNING: " + message + "\n") 
            current_major_version=0
            current_minor_version=0
            current_version="UNK"      
        else:
            logger.critical(message)
            sys.exit("CRITICAL ERROR: " + message + "\n")       
        
    prior_version=False
    if version["major"] > current_major_version:
        prior_version=True
    elif (version["major"] == current_major_version and version["minor"] > current_minor_version):
        prior_version=True
    elif "second minor" in version:
        if (version["major"] == current_major_version and version["minor"] == current_minor_version
            and version["second minor"] > current_second_minor_version):
            prior_version=True
        
    if prior_version and not current_version=="UNK": 
        message=("Please update " + exe + " from version " + current_version               
                 + " to version " + version_required)
        if warning:
            logger.warning(message)
            print("WARNING: " + message + "\n") 
        else:
            logger.critical(message) 
            sys.exit("CRITICAL ERROR: " + message + "\n")
        
    # log version of software
    message="Using " + exe + " version " + current_version
    logger.info(message)

        
class ReportHook():
    def __init__(self):
        self.start_time=time.time()
        
    def report(self, blocknum, block_size, total_size):
        """
        Print download progress message
        """
        
        if blocknum == 0:
            self.start_time=time.time()
            if total_size > 0:
                print("Downloading file of size: " + "{:.2f}".format(byte_to_gigabyte(total_size)) + " GB\n")
        else:
            total_downloaded=blocknum*block_size
            status = "{:3.2f} GB ".format(byte_to_gigabyte(total_downloaded))
                    
            if total_size > 0:
                percent_downloaded=total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stdout to overwrite stdout
                download_rate=total_downloaded/(time.time()-self.start_time)
                estimated_time=(total_size-total_downloaded)/download_rate
                estimated_minutes=int(estimated_time/60.0)
                estimated_seconds=estimated_time-estimated_minutes*60.0
                status +="{:3.2f}".format(percent_downloaded) + " %  " + \
                    "{:5.2f}".format(byte_to_megabyte(download_rate)) + " MB/sec " + \
                    "{:2.0f}".format(estimated_minutes) + " min " + \
                    "{:2.0f}".format(estimated_seconds) + " sec "
            status+="        \r"
            sys.stdout.write(status)
            

def download_tar_and_extract_with_progress_messages(url, filename, folder):
    """
    Download the file at the url
    """
    
    # check for local file
    local_file = False
    if os.path.isfile(url):
        local_file = True   

    if not local_file:
        print("Download URL: " + url)

    try:
        if not local_file:
            url_handle = urlretrieve(url, filename, reporthook=ReportHook().report)
        else:
            filename = url

        print("\nExtracting: " + filename)
        tarfile_handle=tarfile.open(filename)
        tarfile_handle.extractall(path=folder)
    except (EnvironmentError, tarfile.ReadError):
        if local_file:
            sys.exit("CRITICAL ERROR: Unable to extract from local file: " + url)
        else:
            sys.exit("CRITICAL ERROR: Unable to download and extract from URL: " + url)

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
            if config.resume and os.path.getsize(file) > 0:
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
    
class Worker(threading.Thread):
    """
    Create a worker class for threads with a work queue
    Threads also record exit codes to a locked dictionary
    """

    exit_codes_lock = threading.Lock()

    def __init__(self, work_queue, exit_codes):
        super(Worker, self).__init__()
        self.work_queue = work_queue
        self.exit_codes = exit_codes
        
    def run(self):
        """
        Get work from the queue and process
        """
        while True:
            try:
                id,command = self.work_queue.get()
                self.process(id,command)
            finally:
                self.work_queue.task_done()
                
    def process(self, id, command):
        """
        Execute the command for the task id
        Record the exit code
        """
        try:
            execute_command_args_convert(command)
            with self.exit_codes_lock:
                self.exit_codes[id]=0
        # Record the error in the exit codes              
        except (EnvironmentError, subprocess.CalledProcessError):
            with self.exit_codes_lock:
                self.exit_codes[id]=-1
        
def command_threading(threads,commands):
    """
    Process a set of commands using a set of worker threads
    """
    
    # Create a queue to hold work
    work_queue = queue.Queue()
    
    # Create a set of worker threads
    exit_codes={}
    for number in range(threads):
        worker = Worker(work_queue, exit_codes)
        worker.daemon = True
        worker.start()
    
    # Add the work to the queue and start the queue
    commands_by_id={}
    for id,command in enumerate(commands):
        work_queue.put((id,command))
        command=" ".join([command[0]]+[str(i) for i in command[1]])
        commands_by_id[id]=command
    work_queue.join()
    
    # Check for any errors in the threads
    error_commands=[]
    for id in exit_codes:
        if exit_codes[id] != 0:
            error_commands.append("Error message returned from command for thread task " 
                + str(id) + ": " + commands_by_id[id]+"\n")
            
    if error_commands:
        message="\nCRITICAL ERROR: Unable to process all thread commands.\n\n"
        message+="\n".join(error_commands)
        logger.critical(message)
        sys.exit(message)

def execute_command_args_convert(args):
    """
    Convert the list of args to function arguments
    """
    
    return execute_command(*args)


def execute_command(exe, args, infiles, outfiles, stdout_file=None, 
        stdin_file=None, raise_error=None, stderr_file=None):
    """
    Execute third party software or shell command with files
    """
	
    if exe == sys.executable:
        # check that the python module can be found
        module_path=return_module_path(args[0])
        if not module_path:
            message="Can not find python module " + args[0]
            logger.critical(message)
            if raise_error:
                raise EnvironmentError
            else:
                sys.exit("CRITICAL ERROR: " + message)
        # update the module to the full path if not already the full path
        elif not os.path.isabs(args[0]):
            args[0]=os.path.join(module_path,args[0])
            
        logger.debug("Using python module : " + args[0])
    else:
        # check that the executable can be found
        exe_path=return_exe_path(exe)
        if not exe_path:
            message="Can not find executable " + exe
            logger.critical(message)
            if raise_error:
                raise EnvironmentError
            else:
                sys.exit("CRITICAL ERROR: " + message)
        # update the executable to the full path
        else:
            exe=os.path.join(exe_path,exe)
	
        logger.debug("Using software: " + exe)
    
    # check that the input files exist and are readable
    for file in infiles:
        file_exists_readable(file, raise_IOError=raise_error)
        
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
            
        # Open the input and output files (stdin, stdout, stderr)
        stdin=None
        stdout=None
        stderr=None
        
        if stdin_file:
            try:
                stdin=open(stdin_file,"rt")
            except EnvironmentError:
                message="Unable to open file: " + stdin_file
                logger.critical(message)
                if raise_error:
                    raise EnvironmentError
                else:
                    sys.exit("CRITICAL ERROR: " + message)
        
        if stdout_file:
            # check for file open mode
            try:
                stdout_file_name, mode = stdout_file
            except ValueError:
                stdout_file_name = stdout_file
                mode = "w"
            
            try:
                stdout=open(stdout_file_name,mode)
            except EnvironmentError:
                message="Unable to open file: " + stdout_file_name
                logger.critical(message)
                if raise_error:
                    raise EnvironmentError
                else:
                    sys.exit("CRITICAL ERROR: " + message)
                    
        if stderr_file:
            try:
                stderr=open(stderr_file,"w")
            except EnvironmentError:
                message="Unable to open file: " + stderr_file
                logger.critical(message)
                if raise_error:
                    raise EnvironmentError
                else:
                    sys.exit("CRITICAL ERROR: " + message)
	
        try:
            if stdin_file or stdout_file or stderr_file:
                # run command, raise CalledProcessError if return code is non-zero
                p = subprocess.check_call(cmd, stdin=stdin, stdout=stdout, stderr=stderr)
            else:
                p_out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
                logger.debug(p_out)            
        except (EnvironmentError, subprocess.CalledProcessError) as e:
            message="Error executing: " + " ".join(cmd) + "\n"
            if hasattr(e, 'output') and e.output:
                message+="\nError message returned from " + os.path.basename(exe) + " :\n" + e.output.decode("utf-8")
            logger.critical(message)
            logger.critical("TRACEBACK: \n" + traceback.format_exc())
            log_system_status()
            if raise_error:
                raise
            else:
                sys.exit("CRITICAL ERROR: " + message)

        # check that the output files exist and are readable
        for file in outfiles:
            file_exists_readable(file, raise_IOError=raise_error)
    
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
	
    # read in first 2 lines of file to check format
    file_handle = open(file, "rt")
	
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

    file_handle_read=open(file,"rt")

    line=file_handle_read.readline()

    file_type=fasta_or_fastq(file)
    if file_type == "fastq":
        read_token="@"
    else:
        read_token=">"

    sequence_count=0
    lines_since_last_sequence_id=0
    while line:
        if re.search("^"+read_token,line):
            # check that this is not a quality score line
            if file_type == "fastq":
                if lines_since_last_sequence_id>2 or not sequence_count:
                    sequence_count+=1
                    lines_since_last_sequence_id=0
                else:
                    lines_since_last_sequence_id+=1
            else:
                sequence_count+=1
        elif file_type == "fastq":
            lines_since_last_sequence_id+=1
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

    percent=count_reads(unaligned_fastq)/float(count_reads(input_fastq)) * 100

    return format_float_to_string(percent)

def estimate_unaligned_reads_stored(input_fastq, unaligned_store):
    """
    Calculate an estimate of the percent of reads unaligned and stored
    """

    # check if the total number of reads from the input file is stored
    if not unaligned_store.get_initial_read_count():
        # check files exist and are readable
        file_exists_readable(input_fastq)
        unaligned_store.set_initial_read_count(count_reads(input_fastq))

    percent=unaligned_store.count_reads()/float(unaligned_store.get_initial_read_count()) * 100

    return format_float_to_string(percent)

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
    
    file_handle_read=open(fasta_file,"rt")

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

def fastq_to_fasta(file, apply_pick_frames=None, length_annotation=None):
    """
    Convert fastq file to fasta
    Also pick frames for sequences if set
	
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
	
    file_handle_read = open(file, "rt")
	
    line = file_handle_read.readline()
	
    new_file=unnamed_temp_file()
    file_out=open(new_file,"w")

    sequence=""
    sequence_id=""
    while line:
        if re.search("^@",line):
            # write previous sequence
            if sequence:
                if apply_pick_frames:
                    sequences=pick_frames.pick_frames(sequence)
                else:
                    sequences=[sequence]
                    
                if length_annotation:
                    sequence_id=add_length_annotation(sequence_id,len(sequence))
                    
                for sequence in sequences:
                    file_out.write(sequence_id+"\n")
                    file_out.write(sequence+"\n")
            
            sequence_id=line.replace("@",">",1).rstrip()
            sequence=""
        elif re.search("^[A|a|T|t|G|g|C|c|N|n]+$", line):
            sequence+=line.rstrip()
        line=file_handle_read.readline()
        
    # write out the last sequence
    if sequence:
        if apply_pick_frames:
            sequences=pick_frames.pick_frames(sequence)
        else:
            sequences=[sequence]
            
        if length_annotation:
            sequence_id=add_length_annotation(sequence_id,len(sequence))
        
        for sequence in sequences:
            file_out.write(sequence_id+"\n")
            file_out.write(sequence+"\n")    

    file_out.close()	
    file_handle_read.close()

    return new_file	

def pick_frames_from_fasta(file, length_annotation=None):
    """
    Convert fasta file to picked frames
    """
    
    # check file exists
    file_exists_readable(file)
    
    file_handle_read = open(file, "rt")
    
    line = file_handle_read.readline()
    
    new_file=unnamed_temp_file()
    file_out=open(new_file,"w")

    sequence=""
    while line:
        if not re.search("^>",line):
            sequence+=line.rstrip()
        else:
            # if a sequence has been read in then pick frames and write
            if sequence:
                sequences=pick_frames.pick_frames(sequence)
                
                if length_annotation:
                    sequence_id=add_length_annotation(sequence_id,len(sequence))
                    
                for sequence in sequences:
                    file_out.write(sequence_id+"\n")
                    file_out.write(sequence+"\n")
                sequence=""
            sequence_id=line.rstrip()
        line=file_handle_read.readline()
    # if a sequence has been read in then pick frames and write
    if sequence:
        sequences=pick_frames.pick_frames(sequence)

        if length_annotation:
            sequence_id=add_length_annotation(sequence_id,len(sequence))

        for sequence in sequences:
            file_out.write(sequence_id+"\n")
            file_out.write(sequence+"\n")
    file_out.close()    
    file_handle_read.close()

    return new_file

def length_annotate_fasta(file):
    """
    Add annotations of the lengths of the sequences to the fasta sequence ids
    """
    
    # check file exists
    file_exists_readable(file)
    
    file_handle_read = open(file, "rt")
    
    line = file_handle_read.readline()
    
    new_file=unnamed_temp_file()
    file_out=open(new_file,"w")

    sequence=""
    while line:
        if not re.search("^>",line):
            sequence+=line.rstrip()
        else:
            # if a sequence has been read in then annotate and write
            if sequence:
                sequence_id=add_length_annotation(sequence_id,len(sequence))
                file_out.write(sequence_id+"\n")
                file_out.write(sequence+"\n")
                sequence=""
            sequence_id=line.rstrip()
        line=file_handle_read.readline()
    # if a sequence has been read in then annotate and write
    if sequence:
        sequence_id=add_length_annotation(sequence_id,len(sequence))
        file_out.write(sequence_id+"\n")
        file_out.write(sequence+"\n")
        
    file_out.close()    
    file_handle_read.close()

    return new_file      
    
def tsv_to_biom(tsv_file, biom_file, table_type):
    """
    Convert from a tsv to biom file using the biom API
    """
        
    try:
        import biom
    except ImportError:
        sys.exit("Could not find the biom software."+
            " This software is required since the output file is a biom file.")
            
    try:
        import numpy
    except ImportError:
        sys.exit("Could not find the numpy software."+
            " This software is required since the output file is a biom file.")
            
    try:
        import h5py
    except ImportError:
        sys.exit("Could not find the h5py software."+
            " This software is required since the output file is a biom file.")
            
    # read the tsv file
    ids=[]
    data=[]
    with open(tsv_file) as file_handle:
        samples=file_handle.readline().rstrip().split("\t")[1:]
        for line in file_handle:
            row=line.rstrip().split("\t")
            ids.append(row[0])
            data.append(row[1:])
            
    # reformat the rows into a biom table
    table=biom.Table(numpy.array(data), ids, samples)
        
    # write a h5py biom table
    with h5py.File(biom_file, 'w') as file_handle:
        table.to_hdf5(file_handle, table_type)
    
def biom_to_tsv(biom_file):
    """
    Convert from a biom to tsv file
    """

    # create a unnamed temp file
    new_tsv_file=unnamed_temp_file()

    exe="biom"
    args=["convert","-i",biom_file,"-o",new_tsv_file,"--to-tsv"]
    
    message="Converting biom file to tsv ..."
    logger.info(message)
    
    execute_command(exe, args, [biom_file], [new_tsv_file])
    
    return new_tsv_file
    
def format_float_to_string(number):
    """
    Format a float to a string using the config max number of decimals
    """
    
    return "{:.{digits}f}".format(number, digits=config.output_max_decimals)

def byte_to_gigabyte(byte):
    """
    Convert byte value to gigabyte
    """
    
    return byte / (1024.0**3)

def byte_to_megabyte(byte):
    """
    Convert byte value to megabyte
    """
    
    return byte / (1024.0**2)

def byte_to_kilobyte(byte):
    """
    Convert byte value to kilobyte
    """
    
    return byte / 1024.0

def log_system_status():
    """
    Print the status of the system
    """
    
    module_available=True
    try:
        import psutil
    except ImportError:
        module_available=False
        
    if module_available:
        try:
            # record the memory used
            memory = psutil.virtual_memory()
            logger.info("Total memory = " + str(byte_to_gigabyte(memory.total)) + " GB")
            logger.info("Available memory = " + str(byte_to_gigabyte(memory.available)) + " GB")
            logger.info("Free memory = " + str(byte_to_gigabyte(memory.free)) + " GB")
            logger.info("Percent memory used = " + str(memory.percent) + " %")
    
            # record the cpu info
            logger.info("CPU percent = " + str(psutil.cpu_percent()) + " %")
            logger.info("Total cores count = " + str(psutil.cpu_count()))
            
            # record the disk usage
            disk = psutil.disk_usage('/')
            logger.info("Total disk = " + str(byte_to_gigabyte(disk.total)) + " GB")
            logger.info("Used disk = "+ str(byte_to_gigabyte(disk.used)) + " GB")
            logger.info("Percent disk used = " + str(disk.percent) + " %")

            # record information about this current process
            process=psutil.Process()
            process_memory=process.memory_info()
            process_create_time=datetime.datetime.fromtimestamp(
                process.create_time()).strftime("%Y-%m-%d %H:%M:%S")
            process_cpu_times=process.cpu_times()
            # two calls required to cpu percent for non-blocking as per documentation
            process_cpu_percent=process.cpu_percent()
            process_cpu_percent=process.cpu_percent()
            
            logger.info("Process create time = " + process_create_time)
            logger.info("Process user time = " + str(process_cpu_times.user) + " seconds")
            logger.info("Process system time = " + str(process_cpu_times.system) + " seconds")
            logger.info("Process CPU percent = " + str(process_cpu_percent) + " %")
            logger.info("Process memory RSS = " + str(byte_to_gigabyte(process_memory.rss)) + " GB")
            logger.info("Process memory VMS = " + str(byte_to_gigabyte(process_memory.vms)) + " GB")
            logger.info("Process memory percent = " + str(process.memory_percent()) + " %")
            
        except (AttributeError, OSError, TypeError, psutil.Error):
            pass
    
def add_length_annotation(id, length):
    """
    Add the length to the query id
    """
    
    # add the length and handle spaces as translated search will split on spaces
    return id.split(" ")[0]+config.query_length_annotation_delimiter+str(length)

def remove_length_annotation(id):
    """
    Remove the length from the query id
    """
    
    return config.query_length_annotation_delimiter.join(id.split(config.query_length_annotation_delimiter)[:-1])

def get_length_annotation(id):
    """
    Try to get the length annotation from the query id
    """
    
    # check for the annotation delimiter
    if config.query_length_annotation_delimiter in id:
        info=id.split(config.query_length_annotation_delimiter)
        try:
            # the last item is the length
            length=int(info.pop())
            # the first and remaining items are the id
            new_id=config.query_length_annotation_delimiter.join(info)
        except (ValueError, IndexError):
            length=1
            new_id=id
    else:
        # if not present, then return full id and default length
        new_id=id
        length=1

    return new_id, length    
   
def filter_based_on_query_coverage(query_length, query_start_index, query_stop_index,
    query_coverage_threshold):
    """
    Determine if read should be filtered based on query coverage threshold
    """

    if query_length > 1:
        query_coverage = ( ( abs(query_stop_index - query_start_index) + 1) / float(query_length) )* 100.0
    else:
        # if the query length is not provided, default coverage to greater than threshold
        query_coverage = query_coverage_threshold + 1.0
                
    filter = False
    if query_coverage < query_coverage_threshold:
        filter = True

    return filter
 
def get_filtered_translated_alignments(alignment_file_tsv, alignments, apply_filter=None,
                            log_filter=None, unaligned_reads_store=None,
                            query_coverage_threshold=config.translated_query_coverage_threshold, identity_threshold=None):
    """
    Read through the alignment file, yielding filtered alignments
    Filter based on identity threshold, evalue, and coverage threshold
    Remove from unaligned reads store if set
    """

    # if identity threshold is not set, use the config default
    if identity_threshold is None:
        identity_threshold = config.identity_threshold

    # read through the alignment file to identify ids
    # that correspond to aligned reads
    # all translated alignment files will be of the tabulated blast format
    file_handle=open(alignment_file_tsv,"rt")
    line=file_handle.readline()

    log_evalue=False
    large_evalue_count=0
    small_identity_count=0
    small_query_coverage_count=0
    percent_identity_convert_error=0
    alignment_length_convert_error=0
    evalue_convert_error=0
    rapsearch_evalue_convert_error=0
    while line:
        if re.search("^#",line):
            # Check for the rapsearch2 header to determine if these are log(e-value)
            if re.search(config.blast_delimiter,line):
                data=line.split(config.blast_delimiter)
                if len(data)>config.blast_evalue_index:
                    if re.search("log",data[config.blast_evalue_index]):
                        log_evalue=True
        else:
            alignment_info=line.split(config.blast_delimiter)
            
            # try to obtain the identity value to determine if threshold is met
            identity=alignment_info[config.blast_identity_index]
            try:
                identity=float(identity)
            except ValueError:
                percent_identity_convert_error+=1
                identity=0.0

            queryid=alignment_info[config.blast_query_index]
                
            # try converting the alignment length to a number
            alignment_length=alignment_info[config.blast_aligned_length_index]
            try:
                alignment_length=float(alignment_length)
            except ValueError:
                alignment_length_convert_error+=1
                alignment_length=0.0
                
            # try converting evalue to float to check if it is a number
            evalue=alignment_info[config.blast_evalue_index] 
            try:
                evalue=float(evalue)
            except ValueError:
                evalue_convert_error+=1
                evalue=1.0
                                
            # try to get the start and end positions for the query
            try:
                query_start_index = int(alignment_info[config.blast_query_start_index])
                query_stop_index = int(alignment_info[config.blast_query_end_index])
            except (ValueError, IndexError):
                query_start_index=0
                query_stop_index=0
                
            # check for query length annotation
            queryid, query_length = get_length_annotation(queryid)
                
            # try to get the start and end positions for the subject
            try:
                subject_start_index = int(alignment_info[config.blast_subject_start_index])
                subject_stop_index = int(alignment_info[config.blast_subject_end_index])
            except (ValueError, IndexError):
                subject_start_index=0
                subject_stop_index=0

            # convert rapsearch evalue to blastm8 format if logged
            if log_evalue:
                try:
                    evalue=math.pow(10.0, evalue)
                except (ValueError, OverflowError):
                    rapsearch_evalue_convert_error+=1
                    evalue=1.0 
            
            # compute the number of matches
            matches=identity/100.0*alignment_length
            
            # get the protein alignment information
            protein_name, gene_length, bug = alignments.process_reference_annotation(
                alignment_info[config.blast_reference_index])
            
            # check if percent identity is less then threshold
            filter=False
            if identity < identity_threshold:
                filter=True
                small_identity_count+=1
                
            # filter alignments with evalues greater than threshold
            if evalue > config.evalue_threshold:
                filter=True
                large_evalue_count+=1
            
            # filter alignments that do not meet query coverage threshold    
            if filter_based_on_query_coverage(query_length, query_start_index, query_stop_index, query_coverage_threshold):
                filter=True
                small_query_coverage_count+=1
                
            if apply_filter:
                if not filter:
                    yield ( protein_name, gene_length, queryid, matches, bug, 
                            alignment_length, subject_start_index, subject_stop_index )
                elif unaligned_reads_store:
                    # remove the read from the unaligned reads store
                    unaligned_reads_store.remove_id(queryid)
            else:
                yield ( protein_name, gene_length, queryid, matches, bug, 
                        alignment_length, subject_start_index, subject_stop_index )
            
        line=file_handle.readline()
        
    if log_filter:
        logger.debug("Total alignments where percent identity is not a number: " + str(percent_identity_convert_error))
        logger.debug("Total alignments where alignment length is not a number: " + str(alignment_length_convert_error))
        logger.debug("Total alignments where E-value is not a number: " + str(evalue_convert_error))
        if log_evalue:
            logger.debug("Total alignments unable to convert rapsearch e-value: " + str(rapsearch_evalue_convert_error))
        logger.debug("Total alignments not included based on large e-value: " + 
                     str(large_evalue_count))
        logger.debug("Total alignments not included based on small percent identity: " + 
                     str(small_identity_count))
        logger.debug("Total alignments not included based on small query coverage: " + 
                     str(small_query_coverage_count))
    
