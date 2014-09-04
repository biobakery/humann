""" 
Generate pathway coverage and abundance
"""
import os, shutil, tempfile, math
import utilities, config, store

def reactions_by_bug(args):
    """
    Identify the reactions from the hits found for a specific bug
    """
    reactions_database, hits, reactions_file, bug = args
    
    if config.verbose:
        print "Identify reactions for: " + bug
    
    file_handle=open(reactions_file,"w")
    
    alignments=store.alignments()
    
    query_index=alignments.find_index("query")
    reference_index=alignments.find_index("reference")
    evalue_index=alignments.find_index("evalue")
    
    # Group hits by query
    hits_by_query={}
    
    i=0
    for hit in hits:
        if hit[query_index] in hits_by_query:
            hits_by_query[hit[query_index]]+=[i]
        else:
            hits_by_query[hit[query_index]]=[i]
        i+=1
    
    if config.verbose:
        print "Total reads mapped to bug: " + str(len(hits_by_query))
    
    # Loop through hits by query
    # Identify all hits that have a reference which is part of the database
    total_genes={}
    for query in hits_by_query:
        genes=[]
        total_score=0
        for index in hits_by_query[query]:
            hit=hits[index]
            if reactions_database.gene_present(hit[reference_index]):
                score=math.exp(-hit[evalue_index])
                genes.append([hit[reference_index],score])
                total_score+=score
        # add these scores to the total gene scores
        for gene, score in genes:
            total_genes[gene]=score/total_score+total_genes.get(gene,0)
                
    # Merge the gene scores to reaction scores
    total_reactions={}
    for reaction in reactions_database.list_reactions():
        genes_list=reactions_database.find_genes(reaction)
        abundance=0
        # Add the scores for each gene to the total score for the reaction
        for gene in genes_list:
            abundance+=total_genes.get(gene,0)  
        
        # Only write out reactions where the abundance is greater than 0
        if abundance>0: 
            file_handle.write(reaction+config.output_file_column_delimiter
                +str(abundance)+"\n")
        
    file_handle.close()
    # Return the name of the temp file with the reactions
    return reactions_file
    
def reactions(threads, alignments):
    """
    Identify the reactions from the hits found
    """
    
    # load in the reactions database
    gene_to_reactions=os.path.join(config.data_folder,
        config.metacyc_gene_to_reactions)
    reactions_database=store.reactions_database(gene_to_reactions)
    
    reactions_files={}
    # Set up a command to run through each of the hits by bug
    # Also run on all of the hits at once
    args=[]
    for bug in alignments.bug_list()+["all"]:
        
        hits=alignments.hits_for_bug(bug)
        
        # Create a temp file for the reactions results
        file_out, reactions_file=tempfile.mkstemp()
        os.close(file_out)

        reactions_files[bug]=reactions_file

        args.append([reactions_database, hits, reactions_file, bug])
        
    utilities.command_multiprocessing(threads, args, function=reactions_by_bug)

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
