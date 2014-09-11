"""
Identify initial list of bugs from user supplied fasta/fastq
"""

import os, re, sys
import utilities, config

def alignment(input, threads):
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
    
    if threads >1:
        args+=["--nproc",threads]

    args+=opts

    print "\nRunning " + exe + " ........\n"
    utilities.execute_command(exe, args, infiles, outfiles)
    
    return bug_file

def create_custom_database(chocophlan_dir, threshold, bug_file):
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
                if read_percent >= threshold:
                    total_reads_covered += read_percent
                    organism_info=line.split("\t")[0]
                    # use the genus and species
                    species=organism_info.split("|")[-1]
                    genus=organism_info.split("|")[-2]
                    print "Found " + genus + "|" + species + " : " + \
                        str(read_percent) + "% of mapped reads"
                    species_found.append(genus + "." + species)

            line = file_handle.readline()
    
    # compute total species found
    if not config.bypass_prescreen:
        print "\nTotal species indentified in prescreen: " + str(len(species_found)) + "\n"
        print "Species cover " + str(total_reads_covered) + "% of all mapped reads\n"

    # identify the files to be used from the ChocoPhlAn database
    species_file_list = []
    if not config.bypass_prescreen:
        for species_file in os.listdir(chocophlan_dir):
            for species in species_found:
                if re.search(species, species_file): 
                    species_file_list.append(os.path.join(chocophlan_dir,species_file))
                    if config.verbose:
                        print "Adding file to database: " + species_file   
    else:
        for species_file in os.listdir(chocophlan_dir):
            species_file_list.append(os.path.join(chocophlan_dir,species_file))
            if config.verbose:
                print "Adding file to database: " + species_file   


    # create new fasta file containing only those species found
    if not species_file_list:
        print "The custom ChocoPhlAn database is empty"   
        return "Empty"
    else:
        print "\nCreating custom ChocoPhlAn database ........\n"   

        # determine if the files are compressed with gzip
        ext=os.path.splitext(species_file_list[0])[1]

        exe="cat"
        if ext == ".gz":
            exe="zcat"
 
        utilities.execute_command(exe,species_file_list,species_file_list,[custom_database],custom_database)

        return custom_database

