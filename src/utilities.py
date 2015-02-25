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
import logging
import traceback
import gzip
import shutil
import threading
import Queue

import config
import pick_frames

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
            file_handle = gzip.open(file, "r")
            gzipped=True
        else:
            file_handle = open(file, "r")
        
        first_line = file_handle.readline()
        while re.search("^#",first_line):
            first_line = file_handle.readline()
        second_line = file_handle.readline()
    except EnvironmentError:
        # if unable to open and read the file, return unknown
        return "unknown"   
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
                if re.search("^[0-9.]+$",data[config.gene_table_value_index]):
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
        file_handle_gzip=gzip.open(gzip_file,"r")
        
        # create a unnamed temp file
        new_file=unnamed_temp_file()
        
        # write the gunzipped file
        file_handle=open(new_file,"w")
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

def unnamed_temp_file():
    """
    Return the full path to an unnamed temp file
    stored in the unnamed temp folder
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
    full_path=""
    for path in paths:
        fullexe = os.path.join(path,exe)
        if os.path.exists(fullexe):
            if os.access(fullexe,os.X_OK):
                full_path=path
    return full_path

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
        logger.critical("Traceback: \n" + traceback.format_exc())
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
    work_queue = Queue.Queue()
    
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
            error_commands.append("Error message returned from command for threading task " 
                + str(id) + ": " + commands_by_id[id])
            
    if error_commands:
        message="CRITICAL ERROR: Unable to process all threading commands.\n"
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
                stdin=open(stdin_file,"r")
            except EnvironmentError:
                message="Unable to open file: " + stdin_file
                logger.critical(message)
                if raise_error:
                    raise EnvironmentError
                else:
                    sys.exit("CRITICAL ERROR: " + message)
        
        if stdout_file:
            try:
                stdout=open(stdout_file,"w")
            except EnvironmentError:
                message="Unable to open file: " + stdout_file
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
                p = subprocess.call(cmd, stdin=stdin, stdout=stdout, stderr=stderr)
            else:
                p_out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
                logger.debug(p_out)            
        except (EnvironmentError, subprocess.CalledProcessError):
            message="Problem executing " + " ".join(cmd) + "\n"
            logger.critical(message)
            logger.critical("TRACEBACK: \n" + traceback.format_exc())
            log_system_status()
            if raise_error:
                raise subprocess.CalledProcessError
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

    percent=int(count_reads(unaligned_fastq)/float(count_reads(input_fastq)) * 100)

    return str(percent)

def estimate_unaligned_reads_stored(input_fastq, unaligned_store):
    """
    Calculate an estimate of the percent of reads unaligned and stored
    """

    # check if the total number of reads from the input file is stored
    if not unaligned_store.get_initial_read_count():
        # check files exist and are readable
        file_exists_readable(input_fastq)
        unaligned_store.set_initial_read_count(count_reads(input_fastq))

    percent=int(unaligned_store.count_reads()/float(unaligned_store.get_initial_read_count()) * 100)

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

def fastq_to_fasta(file, apply_pick_frames=None):
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
	
    file_handle_read = open(file, "r")
	
    line = file_handle_read.readline()
	
    new_file=unnamed_temp_file()
    file_out=open(new_file,"w")

    while line:
        if re.search("^@",line):
            sequence_id=line.replace("@",">",1)
            line=file_handle_read.readline()
            sequence=""
            while line:
                if re.search("^\\+",line):
                    sequences=[sequence]
                    if apply_pick_frames:
                        sequences=pick_frames.pick_frames(sequence)
                    for sequence in sequences:
                        file_out.write(sequence_id)
                        file_out.write(sequence+"\n")
                    break
                else:
                    sequence+=line.rstrip()
                line=file_handle_read.readline()
        line=file_handle_read.readline()
  	
    file_out.close()	
    file_handle_read.close()

    return new_file	

def pick_frames_from_fasta(file):
    """
    Convert fasta file to picked frames
    """
    
    # check file exists
    file_exists_readable(file)
    
    file_handle_read = open(file, "r")
    
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
                for sequence in sequences:
                    file_out.write(sequence_id)
                    file_out.write(sequence+"\n")
                sequence=""
            sequence_id=line
        line=file_handle_read.readline()
    # if a sequence has been read in then pick frames and write
    if sequence:
        sequences=pick_frames.pick_frames(sequence)
        for sequence in sequences:
            file_out.write(sequence_id)
            file_out.write(sequence+"\n")
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
        config.minpath_original_script)

    # Check for write permission to the target folder
    if not os.access(fullpath_scripts, os.W_OK):
        sys.exit("ERROR: The directory set to install MinPath is not writeable: "+
            fullpath_scripts + " . Please modify the permissions.")

    if not os.path.isfile(minpath_exe):
        if config.verbose:
            print("Installing minpath ... ")
        download_tar_and_extract(config.minpath_url, 
            os.path.join(fullpath_scripts, config.minpath_file),fullpath_scripts)
    
def tsv_to_biom(tsv_file, biom_file, table_type):
    """
    Convert from a tsv to biom file
    """

    exe="biom"
    args=["convert","-i",tsv_file,"-o",biom_file,"--table-type",table_type+" table","--to-hdf5"]
    
    # Remove output file if already exists
    if os.path.isfile(biom_file):
        remove_file(biom_file)
    
    execute_command(exe, args, [tsv_file], [biom_file])
    
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

def log_system_status():
    """
    Print the status of the system
    """
    
    try:
        import psutil
        
        # record the memory used
        memory = psutil.virtual_memory()
        logger.info("Total memory = " + str(memory.total))
        logger.info("Available memory = " + str(memory.available))
        logger.info("Used memory = " + str(memory.used))
        logger.info("Percent memory used = " + str(memory.percent) + " %")

        # record the cpu info
        logger.info("CPU percent = " + str(psutil.cpu_percent()))
        logger.info("Total cores count = " + str(psutil.cpu_count()))
        
        # record the disk usage
        disk = psutil.disk_usage('/')
        logger.info("Total disk = " + str(disk.total))
        logger.info("Used disk = "+ str(disk.used))
        logger.info("Percent disk used = " + str(disk.percent) + " %")

    except (ImportError, AttributeError, OSError):
        pass
    
    
    
    
