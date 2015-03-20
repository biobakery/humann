"""
HUMAnN2: prescreen module
Identify initial list of bugs from user supplied fasta/fastq

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
import re
import sys
import logging

from .. import utilities
from .. import config

# name global logging instance
logger=logging.getLogger(__name__)

def alignment(input):
    """
    Runs metaphlan to identify initial list of bugs
    """
   
    exe="metaphlan2.py"
    opts=config.metaphlan_opts  

    # find the location of the metaphlan dir
    metaphlan_dir=utilities.return_exe_path(exe)
 
    #determine input type as fastq or fasta
    input_type="multi" + utilities.fasta_or_fastq(input)
    
    # outfile name
    bug_file = utilities.name_temp_file(config.bugs_list_name)
    bowtie2_out = utilities.name_temp_file(config.metaphlan_bowtie2_name) 

    infiles=[input, os.path.join(metaphlan_dir, config.metaphlan_pkl_file)]
    
    # location of the index name to multiple files
    infiles_index=[os.path.join(metaphlan_dir, config.metaphlan_mpa_index)]
    
    outfiles=[bug_file, bowtie2_out]
    
    args=[infiles[0],"--bowtie2db",infiles_index[0],
        "-o",outfiles[0],"--input_type",input_type,
        "--bowtie2out",bowtie2_out,
        "--mpa_pkl",infiles[1]]
    
    if config.threads >1:
        args+=["--nproc",config.threads]

    args+=opts

    message="Running " + exe + " ........"
    logger.info(message)
    print("\n"+message+"\n")
    utilities.execute_command(exe, args, infiles, outfiles)
    
    return bug_file

def create_custom_database(chocophlan_dir, bug_file):
    """
    Using ChocoPhlAn creates a custom database based on the bug_file
    """

    # outfile name
    custom_database = utilities.name_temp_file( 
        config.chocophlan_custom_database_name)
    
    species_found = []
    total_reads_covered = 0
    if bug_file != "Empty":
        # Identify the species that pass the threshold
        utilities.file_exists_readable(bug_file)
        file_handle = open(bug_file, "r")

        line = file_handle.readline()
        while line:

            # if we see taxon-level we are done processing
            if re.search("t__", line):
                break
        
            # search for the lines that have the species-level information
            if re.search("s__", line):
                # check threshold
                read_percent=float(line.split("\t")[1])
                if read_percent >= config.prescreen_threshold:
                    total_reads_covered += read_percent
                    organism_info=line.split("\t")[0]
                    # use the genus and species
                    species=organism_info.split("|")[-1]
                    genus=organism_info.split("|")[-2]
                    message=("Found " + genus + "|" + species + " : " +
                        "{:.2f}".format(read_percent) + "% of mapped reads")
                    logger.info(message)
                    print(message)
                    species_found.append(genus + "." + species)

            line = file_handle.readline()
    
    # compute total species found
    if not config.bypass_prescreen:
        message="Total species selected from prescreen: " + str(len(species_found))
        logger.info(message)
        print("\n"+message+"\n")
        message="Selected species explain " + "{:.2f}".format(total_reads_covered) + "% of predicted community composition"
        print(message+"\n")

    # identify the files to be used from the ChocoPhlAn database
    species_file_list = []
    if not config.bypass_prescreen:
        for species_file in os.listdir(chocophlan_dir):
            for species in species_found:
                if re.search(species.lower(), species_file.lower()): 
                    species_file_list.append(os.path.join(chocophlan_dir,species_file))
                    logger.debug("Adding file to database: " + species_file)   
    else:
        for species_file in os.listdir(chocophlan_dir):
            species_file_list.append(os.path.join(chocophlan_dir,species_file))
            logger.debug("Adding file to database: " + species_file)   


    # create new fasta file containing only those species found
    if not species_file_list:
        logger.debug("The custom ChocoPhlAn database is empty")   
        return "Empty"
    else:
        message="Creating custom ChocoPhlAn database ........"
        logger.info(message)
        print("\n"+message+"\n")   

        # determine if the files are compressed with gzip
        ext=os.path.splitext(species_file_list[0])[1]

        exe="cat"
        args=species_file_list
        if ext == ".gz":
            exe="gunzip"
            args=["-c"]+species_file_list
 
        utilities.execute_command(exe,args,species_file_list,[custom_database],custom_database)

        return custom_database

