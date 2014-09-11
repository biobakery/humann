""" 
Generate pathway coverage and abundance
"""
import os, shutil, tempfile, math, re, sys, subprocess
import utilities, config, store, MinPath12hmp


def install_minpath():
    """
    Download and install the minpath software
    """
    
    # Download the minpath software v1.2
    # Check to see if already downloaded
    
    fullpath_scripts=os.path.dirname(os.path.realpath(__file__))
    minpath_exe=os.path.join(fullpath_scripts,config.minpath_folder,
        config.minpath_script)

    if not os.path.isfile(minpath_exe):
        utilities.download_tar_and_extract(config.minpath_url, 
            os.path.join(fullpath_scripts, config.minpath_file),fullpath_scripts)
        
    return minpath_exe


def run_minpath(reactions_file,metacyc_datafile):
    """
    Run minpath on the reactions file using the datafile of pathways
    """
    
    # Create temp files for the results
    file_out, tmpfile=tempfile.mkstemp()
    os.close(file_out)
    
    # Bypass minpath run if reactions file is empty
    if os.path.getsize(reactions_file):
    
        file_out2, tmpfile2=tempfile.mkstemp()
        os.close(file_out2)
    
        if config.verbose:
            print "\nRun MinPath .... "
    
        # Redirect stdout
        sys.stdout=open(os.devnull,"w")
    
        # Call minpath to identify pathways
        MinPath12hmp.Orth2Path(infile = reactions_file, reportfile = "/dev/null", 
            detailfile = tmpfile, whichdb = "ANY", mapfile=metacyc_datafile,
            mpsfile = tmpfile2)
    
        # Undo stdout redirect
        sys.stdout=sys.__stdout__
    
        utilities.remove_file(tmpfile2)
    
    return tmpfile


def pathways_by_bug(args):
    """
    Identify the reactions from the hits found for a specific bug
    """
    reactions_database, hits, reactions_file, bug = args
    
    if config.verbose:
        print "Identify reactions for: " + bug
    
    file_handle=open(reactions_file,"w")
    
    # Group hits by query
    hits_by_query={}
    
    for index, hit in enumerate(hits):
        if hit.get_query() in hits_by_query:
            hits_by_query[hit.get_query()]+=[index]
        else:
            hits_by_query[hit.get_query()]=[index]
    
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
            if reactions_database.gene_present(hit.get_reference()):
                score=math.exp(-hit.get_evalue())
                genes.append([hit.get_reference(),score])
                total_score+=score
        # add these scores to the total gene scores
        for gene, score in genes:
            total_genes[gene]=score/total_score+total_genes.get(gene,0)
                
    # Merge the gene scores to reaction scores
    reactions_store={}
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
            # Store the abundance data to compile with the minpath pathways
            reactions_store[reaction]=abundance
        
    file_handle.close()
    
    metacyc_datafile=os.path.join(config.data_folder,
        config.metacyc_reactions_to_pathways)

    # Run minpath to identify the pathways
    tmpfile=run_minpath(reactions_file, metacyc_datafile)
    
    # Overwrite the reactions file with the pathway data
    pathways_file=reactions_file
    
    # Process the minpath results
    file_handle_read=open(tmpfile, "r")
    file_handle_write=open(pathways_file, "w")
    
    line=file_handle_read.readline()
    
    pathways_store={}
    while line:
        data=line.strip().split(config.minpath_pathway_delimiter)
        if re.search("^path",line):
            current_pathway=data[config.minpath_pathway_index]
        else:
            current_reaction=data[config.minpath_reaction_index]
            # store the pathway and reaction
            pathways_store[current_reaction]=pathways_store.get(
                current_reaction,[]) + [current_pathway]      
        line=file_handle_read.readline()

    file_handle_read.close()
    
    # Write the pathway abundance for each reaction
    for current_reaction in reactions_store:
        # Find the pathways associated with reaction
        for current_pathway in pathways_store.get(current_reaction,[""]):
            new_line=config.output_file_column_delimiter.join([current_reaction,
                current_pathway, str(reactions_store[current_reaction]),"\n"])
            file_handle_write.write(new_line)
    file_handle_write.close()
   
    # Return the name of the temp file with the pathways
    return pathways_file
    
    
def pathways(threads, alignments):
    """
    Identify the reactions and then pathways from the hits found
    """
    
    # Install minpath
    minpath_exe=install_minpath()
    
    # load in the reactions database
    gene_to_reactions=os.path.join(config.data_folder,
        config.metacyc_gene_to_reactions)
    reactions_database=store.reactions_database(gene_to_reactions)
    
    pathways_files={}
    # Set up a command to run through each of the hits by bug
    # Also run on all of the hits at once
    args=[]
    for bug in alignments.bug_list()+["all"]:
        
        hits=alignments.hits_for_bug(bug)
        
        # Create a temp file for the reactions results
        file_out, pathways_file=tempfile.mkstemp()
        os.close(file_out)

        pathways_files[bug]=pathways_file

        args.append([reactions_database, hits, pathways_file, bug])
        
    results=utilities.command_multiprocessing(threads, args, function=pathways_by_bug)

    return pathways_files

def pathways_coverage_by_bug(args):
    """
    Compute the coverage of pathways for one bug
    """
    pathways_file, pathways_database, bug = args
    
    if config.verbose:
        print "Compute pathway coverage for bug: " + bug
    
    file_handle=open(pathways_file,"r")
    
    pathways_reactions_scores={}
    all_scores=[]
    for line in file_handle:
        data=line.strip().split(config.output_file_column_delimiter)
        if len(data) == 3:
            reaction, pathway, score = data
            
            # Bypass lines where the pathway is not listed
            if len(pathway) > 0:
                if pathway in pathways_reactions_scores:
                    pathways_reactions_scores[pathway][reaction]=float(score)
                else:
                    pathways_reactions_scores[pathway]={ reaction : float(score) }
                all_scores.append(float(score))        
    file_handle.close()

    # Find the median score value
    all_scores.sort()
    median_score_value=0
    if all_scores:
        median_score_value=all_scores[len(all_scores)/2]

    # Process through each pathway to compute coverage
    pathways_coverages={}
    xipe_input=[]
    for pathway, reaction_scores in pathways_reactions_scores.items():
        
        # Initialize any reactions in the pathway not found to 0
        for reaction in pathways_database.find_reactions(pathway):
            reaction_scores.setdefault(reaction, 0)
            
        # Count the reactions with scores greater than the median
        count_greater_than_median=0.0
        count_all=0.0
        for reaction, score in reaction_scores.items():
            if score > median_score_value:
               count_greater_than_median+=1.0
            count_all+=1.0
        
        # Compute coverage
        coverage=count_greater_than_median/count_all
        
        pathways_coverages[pathway]=str(coverage)
        xipe_input.append(config.xipe_delimiter.join([pathway,str(coverage)]))
    
    # Run xipe
    xipe_exe=os.path.join(os.path.dirname(os.path.realpath(__file__)),
        config.xipe_script)
    
    cmmd=[xipe_exe,"--file2",config.xipe_percent]
    xipe_subprocess = subprocess.Popen(cmmd, stdin = subprocess.PIPE,
        stdout = subprocess.PIPE, stderr = subprocess.PIPE )
    xipe_stdout, xipe_stderr = xipe_subprocess.communicate( "\n".join(xipe_input))
    
    # Record the pathways to remove based on the xipe error messages
    pathways_to_remove=[]
    for line in xipe_stderr.split("\n"):
        data=line.strip().split(config.xipe_delimiter)
        if len(data) == 2:
            pathways_to_remove.append(data[1])
    
    # Keep some of the pathways to remove based on their xipe scores
    for line in xipe_stdout.split("\n"):
        data=line.strip().split(config.xipe_delimiter)
        if len(data) == 2:
            pathway, pathway_data = data
            if pathway in pathways_to_remove:
                score, bin = pathway_data[1:-1].split(", ")
                if float(score) >= config.xipe_probability and int(bin) == config.xipe_bin:
                    pathways_to_remove.remove(pathway)
            
    # Remove the selected pathways
    for pathway in pathways_to_remove:
        del pathways_coverages[pathway]
    
    # Return a dictionary with a single key of the bug name
    return { bug: pathways_coverages }

def pathways_abundance_by_bug(args):
    """
    Compute the abundance of pathways for one bug
    """
    pathways_file, pathways_database, bug = args
    
    if config.verbose:
        print "Compute pathway abundance for bug: " + bug
    
    file_handle=open(pathways_file,"r")
    
    pathways_reactions_scores={}
    for line in file_handle:
        data=line.strip().split(config.output_file_column_delimiter)
        if len(data) == 3:
            reaction, pathway, score = data
            
            # Bypass lines where the pathway is not listed
            if len(pathway) > 0:
                if pathway in pathways_reactions_scores:
                    pathways_reactions_scores[pathway][reaction]=float(score)
                else:
                    pathways_reactions_scores[pathway]={ reaction : float(score) }        
    file_handle.close()

    # Process through each pathway to compute abundance
    pathways_abundances={}
    for pathway, reaction_scores in pathways_reactions_scores.items():
        
        # Initialize any reactions in the pathway not found to 0
        for reaction in pathways_database.find_reactions(pathway):
            reaction_scores.setdefault(reaction, 0)
            
        # Sort the scores for all of the reactions in the pathway from low to high
        sorted_reaction_scores=sorted(reaction_scores.values())
            
        # Select the second half of the list of reaction scores
        abundance_set=sorted_reaction_scores[(len(sorted_reaction_scores)/ 2):]
        
        # Compute abundance
        abundance=sum(abundance_set)/len(abundance_set)
        
        pathways_abundances[pathway]=str(abundance)
    
    # Return a dictionary with a single key of the bug name
    return { bug: pathways_abundances }
    
def print_pathways(pathways, file, header):
    """
    Print the pathways data to a file organized by pathway
    """
    
    delimiter=config.output_file_column_delimiter
    category_delimiter=config.output_file_category_delimiter
    
    file_handle=open(file,"w")
    
    # Write the header
    file_handle.write("Pathway"+ delimiter + header +"\n")
    
    # Unpack the list of dictionaries to a single dictionary
    # with the all pathways as a separate dictionary
    bug_pathways={}
    all_pathways={}
    for pathway in pathways:
        for bug, pathway_abundances in pathway.items():
            if bug == "all":
                all_pathways=pathway_abundances
            else:
                bug_pathways[bug]=pathway_abundances
    
    for pathway, score in all_pathways.items():
        
        # Write the pathway and score for all bugs
        # Only write pathways where the score is > 0
        if float(score) > 0:
            file_handle.write(delimiter.join([pathway,score])+"\n")
                                          
        # Identify if the pathway is present for each of the bugs
        # If present then print with bug identifier
        for bug in bug_pathways:
            if pathway in bug_pathways[bug]:
                # Write pathway if score is > 0
                if float(bug_pathways[bug][pathway]) > 0:
                    file_handle.write(pathway+category_delimiter+bug
                        +delimiter+bug_pathways[bug][pathway]+"\n")
                    
    file_handle.close()
    

def pathways_abundance_and_coverage(threads, pathways_files):
    """
    Compute the abundance and coverage of the pathways
    """

    # Load in the pathways database
    pathways_database=store.pathways_database(os.path.join(config.data_folder,
        config.metacyc_reactions_to_pathways))
    
    # Compute abundance for all pathways
    args=[]
    for bug, file in pathways_files.items():
         args.append([file, pathways_database, bug])
        
    pathways_abundance=utilities.command_multiprocessing(threads, args, 
        function=pathways_abundance_by_bug)

    # Print the pathways abundance data to file
    print_pathways(pathways_abundance, config.pathabundance_file, "Abundance")

    # Compute coverage 
    pathways_coverage=utilities.command_multiprocessing(threads, args, 
        function=pathways_coverage_by_bug)
    
    # Print the pathways abundance data to file
    print_pathways(pathways_coverage, config.pathcoverage_file, "Coverage")
    
    # Remove the temp pathways files
    for file in pathways_files:
        utilities.remove_file(pathways_files[file])

    return config.pathabundance_file, config.pathcoverage_file
