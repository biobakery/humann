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

    # set the pkl file, if not set in options
    if not "--mpa_pkl" in opts:
        pkl_file=os.path.join(metaphlan_dir, config.metaphlan_pkl_file)
        # check that the file exists
        if os.path.isfile(pkl_file):
            opts+=["--mpa_pkl",pkl_file]
        else:
            sys.exit("The database file for MetaPhlAn does not exist at "
                + pkl_file + " . Please provide the location with --metaphlan-options .")
    
    # set the bowtie2 database file, if not set in options
    if not "--bowtie2db" in opts:
        mpa_index=os.path.join(metaphlan_dir, config.metaphlan_mpa_index)
        # check that the database exits, by adding the extension to the index name
        mpa_file=mpa_index+config.bowtie2_index_ext_list[0]
        if os.path.isfile(mpa_file):
            opts+=["--bowtie2db",mpa_index]
        else:
            sys.exit("The bowtie2 database file for MetaPhlAn does not exist at "
                + mpa_file + " . Please provide the location with --metaphlan-options ." )
    
    args=[input]+opts+["-o",bug_file,"--input_type",input_type, "--bowtie2out",bowtie2_out]
    
    if config.threads >1:
        args+=["--nproc",config.threads]

    message="Running " + exe + " ........"
    logger.info(message)
    print("\n"+message+"\n")
    utilities.execute_command(exe, args, [input], [bug_file, bowtie2_out])
    
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
        file_handle = open(bug_file, "rt")

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
                    try:
                        species=organism_info.split("|")[-1]
                        genus=organism_info.split("|")[-2]
                    except IndexError:
                        species=""
                        genus=""
                        logger.debug("Unable to process species: " + line)
                        
                    if species and genus:
                        message=("Found " + genus + "." + species + " : " +
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
        message="\n\n"
        if len(species_found) > 0:
            message+="None of the species selected from the prescreen were found in the ChocoPhlAn database.\n"
        elif config.bypass_prescreen:
            message+="The ChocoPhlAn database is empty.\n"
        else:
            message+="No species were selected from the prescreen.\n"
        message+="Because of this the custom ChocoPhlAn database is empty.\n"
        message+="This will result in zero species-specific gene families and pathways.\n\n"
        logger.debug(message)
        print(message)
        return "Empty"
    else:
        message="Creating custom ChocoPhlAn database ........"
        logger.info(message)
        print("\n"+message+"\n")   

        # determine if the files are compressed with gzip
        ext=os.path.splitext(species_file_list[0])[1]

        exe="cat"
        args=[]
        if ext == ".gz":
            exe="gunzip"
            args=["-c"]
        
        # check if set to bypass this step
        bypass=utilities.check_outfiles([custom_database])
        
        if not bypass:
            # run the command with chunks of input files to not exceed max arguments
            species_file_list_subsets=[species_file_list[i:i+config.max_arguments] for i in range(0, len(species_file_list), config.max_arguments)]
            
            # run the first subset to create the new custom database
            first_subset=species_file_list_subsets.pop(0)
            utilities.execute_command(exe,args+first_subset,first_subset,[],custom_database)
            
            # append the remaining subsets to the existing database
            for subset in species_file_list_subsets:
                utilities.execute_command(exe,args+subset,subset,[],[custom_database,"a"])

        return custom_database

