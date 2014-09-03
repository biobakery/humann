""" 
Generate pathway coverage and abundance
"""
import os, shutil, tempfile
import utilities, config, store
import hits

def pickled_hits(alignments):
    """
    Create pickled hits files by bug
    Also create a file of all of the hits
    """
    hits_files={}
    
    file_all_out, all_file=tempfile.mkstemp()
    os.close(file_all_out)
    hits_files["all"]=all_file
    pHits_all = hits.CHits()
    
    for bug in alignments.bug_list():
        file_out, new_file=tempfile.mkstemp()
        os.close(file_out)
        hits_files[bug]=new_file
        
        # create a new pickled hits instance
        pHits = hits.CHits()
        reference_index=alignments.find_index("reference")
        query_index=alignments.find_index("query")
        evalue_index=alignments.find_index("evalue")
        identity_index=alignments.find_index("identity")
        coverage_index=alignments.find_index("coverage")
        
        for hit in alignments.hits_for_bug(bug):
            pHits.add(hit[reference_index],hit[query_index],
                hit[evalue_index],hit[identity_index],hit[coverage_index])
            pHits_all.add(hit[reference_index],hit[query_index],
                hit[evalue_index],hit[identity_index],hit[coverage_index])
            
        file_handle=open(new_file,"w")
        pHits.save(file_handle)
        file_handle.close()
    
    file_handle=open(all_file,"w")
    pHits_all.save(file_handle)
    file_handle.close()
        
    return hits_files

def reactions(threads, hits_files):
    """
    Identify the reactions from the hits found
    """
    
    gene_to_reactions=os.path.join(config.data_folder,
        config.metacyc_gene_to_reactions)
    
    reactions_files={}
    # Set up a command to run through each of the hits files
    commands=[]
    for hit in hits_files:
        
        hits_file=hits_files[hit]
        
        # Create a temp file for the reactions results
        file_out, reactions_file=tempfile.mkstemp()
        os.close(file_out)

        reactions_files[hit]=reactions_file

        commands.append([config.humann1_script_hits_to_reactions,
        [gene_to_reactions],[hits_file],[],
        reactions_file, hits_file])
    
    utilities.command_multiprocessing(threads, commands)

    # Remove the temp pickled hits files
    for file in hits_files:
        utilities.remove_file(hits_files[file])

    return reactions_files

def pathways(threads, reactions_files):
    """
    Compute the pathways for the reactions found
    """

    # Just run the all reactions for now
    reactions_file=reactions_files["all"]

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

    # Remove the temp reactions files
    for file in reactions_files:
        utilities.remove_file(reactions_files[file])

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
