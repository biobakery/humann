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
import gzip
import shutil

import config

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

    # check that second line is only nucleotides or amino acids
    if re.search("^[A-Z|a-z]+$", second_line):
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
    if not format:
        format="unknown"
    elif gzipped:
        format+=".gz"
                
    message="File ( " + file + " ) is of format:  " + format
    if config.verbose:
        print(message+"\n")
    
    logger.info(message)       
                
    return format

def gunzip_file(gzip_file):
    """
    Return a new copy of the file that is not gzipped
    The new file will be placed in the unnamed temp folder
    """
    
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
    
    os.environ["PATH"] = exe_dir + os.pathsep + os.environ["PATH"]

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
    
class MultiprocessingWorker(multiprocessing.Process):
    """
    Class to get/put to worker queue
    """
    
    def __init__(self, work_queue, results_queue, function):
        super(MultiprocessingWorker, self).__init__()
        self.work_queue=work_queue
        self.results_queue=results_queue
        self.function=function
        
    def run(self):
        while True:
            try:
                args = self.work_queue.get()
                results = self.function(args)
                self.results_queue.put(results)
            finally:
                self.work_queue.task_done()
        

def command_multiprocessing(threads, args, function=None, lock=None):
    """
    Execute command in parallel, maximum of threads total
    """
    
    if not function:
        function=execute_command_args_convert
    
    try:
        if threads>1:
            logger.debug("Create %s processes for function %s", threads, function)
            if lock:
                work_queue = multiprocessing.JoinableQueue()
                results_queue = multiprocessing.Queue()
                # set up thread number of workers
                for i in range(threads):
                    worker=MultiprocessingWorker(work_queue, results_queue, function)
                    worker.daemon=True
                    worker.start()
                    
                # add the arguments for each run to the queue
                for arg in args:
                    work_queue.put(arg)
                work_queue.join()
                
                # collect the results
                results=[]
                for i in range(len(args)):
                    results.append(results_queue.get())
            else:
                pool=multiprocessing.Pool(threads)
                results=pool.map(function, args)
        else:
            results=[function(arg) for arg in args]
    except (EnvironmentError, ValueError, subprocess.CalledProcessError):
        logger.critical("TRACEBACK: \n" + traceback.format_exc())
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
        except (EnvironmentError, subprocess.CalledProcessError):
            message="Problem executing " + " ".join(cmd) + "\n"
            logger.critical(message)
            logger.critical("TRACEBACK: \n" + traceback.format_exc())
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
    
def process_chocophlan_length(location,uniref):
    """
    Return the length given the sequence location
    """
    
    try:
        if config.chocophlan_multiple_location_delimiter in location:
            locations=location.split(config.chocophlan_multiple_location_delimiter)
        else:
            locations=[location]
        length=0
        for location in locations:
            start, end = re.sub(config.chocophlan_location_extra_characters,
                '',location).split(config.chocophlan_location_delimiter)
            length=length+abs(int(end)-int(start))+1
    except (ValueError, IndexError):
        length=0
        logger.debug("Unable to compute length for gene: " + uniref)
    
    return length

def process_reference_annotation(reference):
    """
    Process the reference string for information on gene, gene length, and bug
    Allow for chocophlan annotations, gene|gene_length, gene_length|gene, and gene
    """
    
    reference_info=reference.split(config.chocophlan_delimiter)
    
    # identify bug and gene families
    location=""
    length=0
    uniref=reference_info[0]
    try:
        bug=reference_info[config.chocophlan_bug_index]
        uniref=reference_info[config.chocophlan_uniref_index]
        location=reference_info[config.chocophlan_location_index]
    except IndexError:
        # try to find gene length if present
        bug="unclassified"
        if len(reference_info)==2:
            if re.search("^[0-9]+$",reference_info[0]):
                length=int(reference_info[0])
                uniref=reference_info[1]
            elif re.search("^[0-9]+$",reference_info[1]):
                length=int(reference_info[1])
                uniref=reference_info[0]
                        
    # compute the length of the gene from the location provided
    if location:
        length=process_chocophlan_length(location, uniref)
        
    return [uniref,length,bug]
    
    
		
