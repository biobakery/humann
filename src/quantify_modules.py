""" 
Generate pathway coverage and abundance
"""
import os, shutil, tempfile, math, re, sys, subprocess
import utilities, config, store, MinPath12hmp


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


def identify_reactions_and_pathways_by_bug(args):
    """
    Identify the reactions from the hits found for a specific bug
    """
    
    reactions_database, pathways_database, hits, bug = args
    
    if config.verbose:
        print "Compute gene scores by query ..."
    
    # Loop through hits to organize scores by query
    genes_by_query={}
    total_scores_by_query={}
    for hit in hits:
        bug, reference, query, evalue=hit
        score=math.exp(-evalue)
        if query in genes_by_query:
            genes_by_query[query].append([reference,score])
        else:
            genes_by_query[query]=[[reference,score]]
        total_scores_by_query[query]=total_scores_by_query.get(query,0)+score
    
    if config.verbose:
        print "Total reads mapped to " + bug + " : " + str(len(hits))            
    
    # add these scores by query to the total gene scores
    total_genes={}
    for query, genes in genes_by_query.items():
        for gene, score in genes:
            total_genes[gene]=total_genes.get(gene,0)+score/total_scores_by_query[query]
                
    # Merge the gene scores to reaction scores   
    reactions={}
    reactions_file_lines=[]
    for reaction in reactions_database.reaction_list():
        genes_list=reactions_database.find_genes(reaction)
        abundance=0
        # Add the scores for each gene to the total score for the reaction
        for gene in genes_list:
            abundance+=total_genes.get(gene,0)  
        
        # Only write out reactions where the abundance is greater than 0
        if abundance>0: 
            reactions_file_lines.append(reaction+config.output_file_column_delimiter
                +str(abundance)+"\n")
            # Store the abundance data to compile with the minpath pathways
            reactions[reaction]=abundance

    pathways={}
    if config.minpath_toggle == config.toggle_on:   

        # Create a temp file for the reactions results
        file_descriptor, reactions_file=tempfile.mkstemp()
        os.write(file_descriptor, "".join(reactions_file_lines))
        os.close(file_descriptor)
 
        metacyc_datafile=os.path.join(config.data_folder,
            config.metacyc_reactions_to_pathways)
    
        # Run minpath to identify the pathways
        tmpfile=run_minpath(reactions_file, metacyc_datafile)
        
        # Remove the temp reactions file
        utilities.remove_file(reactions_file)
        
        # Process the minpath results
        file_handle_read=open(tmpfile, "r")
        
        line=file_handle_read.readline()
        
        while line:
            data=line.strip().split(config.minpath_pathway_delimiter)
            if re.search(config.minpath_pathway_identifier,line):
                current_pathway=data[config.minpath_pathway_index]
            else:
                current_reaction=data[config.minpath_reaction_index]
                # store the pathway and reaction
                pathways[current_reaction]=pathways.get(
                    current_reaction,[]) + [current_pathway]      
            line=file_handle_read.readline()
    
        file_handle_read.close()
        
        # Remove the minpath results file
        utilities.remove_file(tmpfile)
    else:
        # Add all pathways associated with each reaction if not using minpath
        for current_reaction in reactions:
            pathways[current_reaction]=pathways.get(
                current_reaction, []) + pathways_database.find_pathways(current_reaction)
        
    
    pathways_and_reactions_store=store.pathways_and_reactions(bug)
    # Store the pathway abundance for each reaction
    for current_reaction in reactions:
        # Find the pathways associated with reaction
        for current_pathway in pathways.get(current_reaction,[""]):
            # Only store data for items with pathway names
            if current_pathway:
                pathways_and_reactions_store.add(current_reaction, current_pathway, 
                    reactions[current_reaction])
   
    # Return the name of the temp file with the pathways
    return pathways_and_reactions_store
    
    
def identify_reactions_and_pathways(threads, alignments, reactions_database, pathways_database):
    """
    Identify the reactions and then pathways from the hits found
    """
        
    # Remove the hits from the alignments that are not associated with a gene
    # in the gene to reactions database
    if config.verbose:
        print "Remove hits to genes not in reactions database ..."

    for gene in alignments.gene_list():
        if not reactions_database.gene_present(gene):
            alignments.delete_gene_and_hits(gene)
            
    # Update the bugs index list to remove any indexes that point to deleted hits        
    alignments.update_hits_for_bugs()
            
    print "Total bugs after filtering: " + str(alignments.count_bugs())
    alignments.print_bugs()            
            
    print "Total gene families after filtering: " + str(len(alignments.gene_list()))
            
    # Set up a command to run through each of the hits by bug
    args=[]
    for bug in alignments.bug_list():
        
        hits=alignments.hits_for_bug(bug)

        args.append([reactions_database, pathways_database, hits, bug])
        
    pathways_and_reactions_store=utilities.command_multiprocessing(threads, args, 
        function=identify_reactions_and_pathways_by_bug)

    return pathways_and_reactions_store

def pathways_coverage_by_bug(args):
    """
    Compute the coverage of pathways for one bug
    """
    
    pathways_and_reactions_store, pathways_database = args
    
    if config.verbose:
        print "Compute pathway coverage for bug: " + pathways_and_reactions_store.get_bug()
    
    # Process through each pathway to compute coverage
    pathways_coverages={}
    xipe_input=[]
    median_score_value=pathways_and_reactions_store.median_score()
    
    for pathway, reaction_scores in pathways_and_reactions_store.get_items():
        
        # Initialize any reactions in the pathway not found to 0
        for reaction in pathways_database.find_reactions(pathway):
            reaction_scores.setdefault(reaction, 0)
            
        # Count the reactions with scores greater than the median
        count_greater_than_median=0
        for reaction, score in reaction_scores.items():
            if score > median_score_value:
               count_greater_than_median+=1
        
        # Compute coverage
        coverage=count_greater_than_median/float(len(reaction_scores.keys()))
        
        pathways_coverages[pathway]=coverage
        xipe_input.append(config.xipe_delimiter.join([pathway,str(coverage)]))
    
    # Check config to determine if xipe should be run
    if config.xipe_toggle == config.toggle_on:
        # Run xipe
        xipe_exe=os.path.join(os.path.dirname(os.path.realpath(__file__)),
            config.xipe_script)
        
        cmmd=[xipe_exe,"--file2",config.xipe_percent]
        
        if config.verbose:
            print "Run xipe ...."
        xipe_subprocess = subprocess.Popen(cmmd, stdin = subprocess.PIPE,
            stdout = subprocess.PIPE, stderr = subprocess.PIPE )
        xipe_stdout, xipe_stderr = xipe_subprocess.communicate("\n".join(xipe_input))
        
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
    
    return store.pathways(pathways_and_reactions_store.get_bug(), pathways_coverages)

def pathways_abundance_by_bug(args):
    """
    Compute the abundance of pathways for one bug
    """
    
    pathways_and_reactions_store, pathways_database = args
    
    if config.verbose:
        print "Compute pathway abundance for bug: " + pathways_and_reactions_store.get_bug()

    # Process through each pathway to compute abundance
    pathways_abundances={}
    for pathway, reaction_scores in pathways_and_reactions_store.get_items():
        
        # Initialize any reactions in the pathway not found to 0
        for reaction in pathways_database.find_reactions(pathway):
            reaction_scores.setdefault(reaction, 0)
            
        # Sort the scores for all of the reactions in the pathway from low to high
        sorted_reaction_scores=sorted(reaction_scores.values())
            
        # Select the second half of the list of reaction scores
        abundance_set=sorted_reaction_scores[(len(sorted_reaction_scores)/ 2):]
        
        # Compute abundance
        abundance=sum(abundance_set)/len(abundance_set)
        
        pathways_abundances[pathway]=abundance
    
    # Return a dictionary with a single key of the bug name
    return store.pathways(pathways_and_reactions_store.get_bug(), pathways_abundances)
    
def print_pathways(pathways, file, header):
    """
    Print the pathways data to a file organized by pathway
    """
    
    delimiter=config.output_file_column_delimiter
    category_delimiter=config.output_file_category_delimiter
    

    
    # Compile data for all bugs by pathways
    all_pathways_scores={}
    all_pathways_scores_by_bug={}
    for bug_pathways in pathways:
        bug=bug_pathways.get_bug()
        # Add up all scores based on score for each bug for each pathway
        for pathway, score in bug_pathways.get_items():
            if score>0:
                all_pathways_scores[pathway]=all_pathways_scores.get(pathway,0)+score
                if not pathway in all_pathways_scores_by_bug:
                    all_pathways_scores_by_bug[pathway]={bug: score}
                else:
                    all_pathways_scores_by_bug[pathway][bug]=score             
 
    # Write the header
    file_handle=open(file,"w")
    file_handle.write("Pathway"+ delimiter + header +"\n")       
    
    # Print out the pathways with those with the highest scores first
    for pathway in sorted(all_pathways_scores, key=all_pathways_scores.get, reverse=True):
        # Print the sum of all bugs for pathway
        file_handle.write(pathway+delimiter+str(all_pathways_scores[pathway])+"\n")
        # Print scores per bug for pathway ordered with those with the highest values first
        for bug in sorted(all_pathways_scores_by_bug[pathway], key=all_pathways_scores_by_bug[pathway].get, reverse=True):
            file_handle.write(pathway+category_delimiter+bug+delimiter
                              +str(all_pathways_scores_by_bug[pathway][bug])+"\n")
                    
    file_handle.close()
    

def compute_pathways_abundance_and_coverage(threads, pathways_and_reactions_store, pathways_database):
    """
    Compute the abundance and coverage of the pathways
    """
    
    # Compute abundance for all pathways
    args=[]
    for bug_pathway_and_reactions_store in pathways_and_reactions_store:
         args.append([bug_pathway_and_reactions_store, pathways_database])
        
    pathways_abundance=utilities.command_multiprocessing(threads, args, 
        function=pathways_abundance_by_bug)

    # Print the pathways abundance data to file
    print_pathways(pathways_abundance, config.pathabundance_file, "Abundance")

    # Compute coverage 
    pathways_coverage=utilities.command_multiprocessing(threads, args, 
        function=pathways_coverage_by_bug)
    
    # Print the pathways abundance data to file
    print_pathways(pathways_coverage, config.pathcoverage_file, "Coverage")

    return config.pathabundance_file, config.pathcoverage_file
