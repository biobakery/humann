"""
HUMAnN: prescreen module
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
   
    exe="metaphlan"
    opts=config.metaphlan_opts  

    # find the location of the metaphlan dir
    metaphlan_dir=utilities.return_exe_path(exe)
 
    #determine input type as fastq or fasta
    input_type=utilities.fasta_or_fastq(input)
    
    # outfile name
    bug_file = utilities.name_temp_file(config.bugs_list_name)
    bowtie2_out = utilities.name_temp_file(config.metaphlan_bowtie2_name) 

    args=[input]+opts+["-o",bug_file,"--input_type",input_type, "--bowtie2out",bowtie2_out]
    
    if config.threads >1:
        args+=["--nproc",config.threads]

    message="Running " + exe + " ........"
    logger.info(message)
    print("\n"+message+"\n")
    utilities.execute_command(exe, args, [input], [bug_file, bowtie2_out])
    
    return bug_file

def read_sgb_mapping_file():
    """
    Read in the file mapping SGB to species names
    """

    sgbs_to_species={}
    with open(config.sgb_to_species_file) as file_handle:
        header=file_handle.readline()

        for line in file_handle:
            data=line.rstrip().split("\t")
            sgb="t__"+data[0]
            species=["s__"+data[1]]

            # add additional species if included
            try:
                species+=[ "s__"+name for name in filter(lambda x: x, data[3].split(","))]
            except IndexError:
                pass

            if not sgb in sgbs_to_species:
                sgbs_to_species[sgb]=species
            else:
                sgbs_to_species[sgb]+=species
    return sgbs_to_species

def get_abundance(line):
    """
    Read in the abundance value from the taxonomy file
    """
    try:
        data=line.split("\t")
        if data[-1].replace(".","").replace("e-","").isdigit():
            read_percent=float(data[-1])
        else:
            read_percent=float(data[-2])
    except ValueError:
        message="The MetaPhlAn taxonomic profile provided was not generated with the expected database version. Please update your version of MetaPhlAn to at least v3.0."
        logger.error(message)
        sys.exit("\n\nERROR: "+message)

    return read_percent
                   
def get_species_name(line, sgb=False):
    """
    Parse the line for the genus.species name
    """
 
    offset=0
    if sgb:
        offset=-1

    try:
        species=line.split("|")[-1+offset]
        genus=line.split("|")[-2+offset]
        name=genus + "." + species
    except IndexError:
        name=""
        logger.debug("Unable to process species: " + line)

    return name                        

def create_custom_database(chocophlan_dir, bug_file):
    """
    Using ChocoPhlAn creates a custom database based on the bug_file
    """

    # outfile name
    custom_database = utilities.name_temp_file( 
        config.chocophlan_custom_database_name)
    
    species_found = []
    sgb_species_found = []
    abundances = {}
    sgb_abundances = {}
    additional_species_found = {}
    total_reads_covered = 0
    if bug_file != "Empty":
        # Identify the species that pass the threshold
        utilities.file_exists_readable(bug_file)
        file_handle = open(bug_file, "rt")

        # read in the sgb file
        sgb_to_species=read_sgb_mapping_file()

        line = file_handle.readline()
        version_found = False
        while line:


            if line.startswith("#") and (config.metaphlan_v3_db_version in line or config.metaphlan_v4_db_version in line):
                version_found = True

            # if we see strain-level, look for possible SGBs
            if re.search("t__", line) and re.search("s__", line):
                # check threshold
                read_percent=get_abundance(line)

                if read_percent >= config.prescreen_threshold:
                    organism_info=line.split("\t")[0]

                    try:
                        additional_species=list(filter(lambda x: "s__" in x,line.rstrip().split("\t")[-1].split(",")))
                    except AttributeError:
                        additional_species=[]


                    # get the species from the sgb mapping file
                    sgb=organism_info.split("|")[-1]
                    species=[]
                    try:
                        species=sgb_to_species[sgb]
                    except KeyError:
                        species=[]
                        logger.debug("Taxon not in mapping file: " + line)

                    # also include the genus and species listed in the abundance file
                    genus_species=get_species_name(organism_info, sgb=True)
                    additional_taxa=list(filter(None,[get_species_name(info) for info in additional_species]))

                    if species:
                        sgb_abundances[sgb]=read_percent
                        sgb_species_found+=species
                        additional_species_found[sgb]=additional_taxa
                    else:
                        sgb_abundances[genus_species]=read_percent
                        additional_species_found[genus_species]=additional_taxa
                    sgb_species_found+=[genus_species]
                    sgb_species_found+=additional_taxa

            # search for the lines that have the species-level information
            if re.search("s__", line):
                # check threshold
                read_percent=get_abundance(line)
                    
                if read_percent >= config.prescreen_threshold:
                    organism_info=line.split("\t")[0]
                    # use the genus and species
                    genus_species=get_species_name(organism_info)
                        
                    if genus_species:
                        abundances[genus_species]=read_percent
                        species_found.append(genus_species)

            line = file_handle.readline()
   
        if not version_found:
            message="The MetaPhlAn taxonomic profile provided was not generated with the database version "+\
                config.metaphlan_v3_db_version+" or "+config.metaphlan_v4_db_version+" . Please update your version of MetaPhlAn to at least v3.0 or if you are using MetaPhlAn v4 please use the database "+\
                config.metaphlan_v4_db_version+"."
            logger.error(message)
            sys.exit("\n\nERROR: "+message)
        
    # if sgbs are provided, use those instead of species
    if sgb_species_found:
        species_found=sgb_species_found
        for sgb in sgb_abundances:
            total_reads_covered += sgb_abundances[sgb]
            extra_species = sgb_to_species.get(sgb,[]) + additional_species_found[sgb]
            message=("Found " + sgb + " : " + "{:.2f}".format(sgb_abundances[sgb]) + "% of mapped reads ( "+",".join(extra_species)+" )")
                
            logger.info(message)
            print(message)
    else:
        for species_genus in abundances:
            total_reads_covered += abundances[species_genus]
            message=("Found " + species_genus + " : " + "{:.2f}".format(abundances[species_genus]) + "% of mapped reads")
            logger.info(message)
            print(message)

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
                # match the exact genus and species from the MetaPhlAn (or custom) list
                new_database_file=os.path.join(chocophlan_dir,species_file)
                if ( re.search(species.lower()+"\.", species_file.lower()) or re.search(species.lower()+"_group\.", species_file.lower()) ) and not new_database_file in species_file_list: 
                    species_file_list.append(new_database_file)
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

