"""
Identify initial list of bugs from user supplied fasta/fastq
"""

import os, re, sys
import utilities

def run_alignment(input, threads, temp_dir):
    """
    Runs metaphlan to identify initial list of bugs
    """
   
    exe="metaphlan2.py"
    opts=["-t","rel_ab"]   

    # find the location of the metaphlan dir
    metaphlan_dir=utilities.return_exe_path(exe)
 
    #determine input type as fastq or fasta
    input_type="multi" + utilities.fasta_or_fastq(input)
    
    # outfile name
    sample_name = os.path.splitext(os.path.basename(input))[0]
    bug_file = os.path.join(temp_dir, sample_name + "_bugs_list.tsv")
    bowtie2_out = os.path.join(temp_dir, sample_name + "_bowtie2_out.txt") 

 
    infiles=[input, os.path.join(metaphlan_dir, "db_v20/mpa_v20_m200.pkl")]
    
    # location of the index name to multiple files
    infiles_index=[os.path.join(metaphlan_dir, "db_v20/mpa_v20_m200")]
    
    outfiles=[bug_file, bowtie2_out]
    
    args=[infiles[0],"--bowtie2db",infiles_index[0],
        "-o",outfiles[0],"--input_type",input_type,
        "--bowtie2out",bowtie2_out,
        "--mpa_pkl",infiles[1]]
    
    if threads >1:
        args+=["--nproc",threads]

    args+=opts

    print "\nRunning " + exe + " ........\n"
    utilities.execute_command(exe,args, infiles, outfiles,"")
    
    return bug_file

def create_custom_database(chocophlan_dir, threshold, bug_file, temp_dir):
    """
    Using ChocoPhlAn creates a custom database based on the bug_file
    """

    # Identify the species that pass the threshold
    file_handle = open(bug_file, "r")

    # outfile name
    bug_sample_name = os.path.splitext(os.path.basename(bug_file))[0]
    custom_database = os.path.join(temp_dir, bug_sample_name + "_custom_database.ffn")
    
    species_found = []
    total_reads_covered = 0
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
                    str(read_percent) + "% of reads"
                species_found.append(genus + "." + species)

        line = file_handle.readline()
    
    # compute total species found
    print "\nTotal species indentified in prescreen: " + str(len(species_found)) + "\n"
    print "Species cover " + str(total_reads_covered) + "% of all reads in input\n"

    # identify the files to be used from the ChocoPhlAn database
    species_file_list = []
    for species_file in os.listdir(chocophlan_dir):
        for species in species_found:
            if re.search(species, species_file): 
                species_file_list.append(os.path.join(chocophlan_dir,species_file))
                print "Adding file to database: " + species_file   

    # create new fasta file containing only those species found
    if not species_file_list:
        sys.exit("ERROR: The custom ChocoPhlAn database is empty.\n")   
    
    print "\nCreating custom ChocoPhlAn database ........\n"   
 
    utilities.execute_command("cat",species_file_list,species_file_list,[custom_database],custom_database)

    return custom_database

