""" 
Generate pathway coverage and abundance
"""
import os
import utilities, config, shutil, tempfile

def reactions(hits_file):
    """
    Identify the reactions from the hits found
    """

    reactions_file=utilities.name_temp_file(
        config.reactions_name)

    gene_to_reactions=os.path.join(config.data_folder,
        config.metacyc_gene_to_reactions)
    utilities.execute_command(
        config.humann1_script_hits_to_reactions,
        [gene_to_reactions],[hits_file],[],
        reactions_file, hits_file)

    return reactions_file

def pathways(reactions_file):
    """
    Compute the pathways for the reactions found
    """

    # Download the minpath software v1.2
    # Check to see if already downloaded
    minpath_exe=os.path.join(config.data_folder,config.minpath_folder,
        config.humann1_script_minpath)

    if not os.path.isfile(minpath_exe):
        utilities.download_tar_and_extract(config.minpath_url, 
            os.path.join(config.data_folder, config.minpath_file))

        # Copy the updated executable
        shutil.copy(os.path.join(config.humann1_scripts,
            config.humann1_script_minpath), os.path.join(config.data_folder,
            config.minpath_folder))

    metacyc_datafile=os.path.join(config.data_folder,
        config.metacyc_reactions_to_pathways)

    pathways_file=utilities.name_temp_file(
        config.pathways_minpath)

    # Identify pathways
    utilities.execute_command(
        config.humann1_script_enzymes_to_pathways,
        [minpath_exe,metacyc_datafile],[reactions_file],[],
        pathways_file, reactions_file)

    return pathways_file

def pathways_abundance_and_coverage(pathways_file):
    """
    Compute the abundance and coverage of the pathways
    """
    
    metacyc_datafile=os.path.join(config.data_folder,
        config.metacyc_reactions_to_pathways)

    # Compute abundance
    utilities.execute_command(
        "cat",[],[pathways_file],[],
        config.pathabundance_file, pathways_file)

    # Compute coverage 
    utilities.execute_command(
        "cat",[],[pathways_file],[],
        config.pathcoverage_file, pathways_file)

    return config.pathabundance_file, config.pathcoverage_file
