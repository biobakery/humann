""" 
HUMAnN2: quantify_modules module
Generate pathway coverage and abundance

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
import shutil
import math
import re
import sys
import subprocess
import logging

import utilities
import config
import store
import quantify_families

# name global logging instance
logger=logging.getLogger(__name__)

def minpath_command(reactions_file,metacyc_datafile):
    """
    Return the minpath command for the reactions file using the datafile of pathways
    """
    
    # Create temp files for the results
    tmpfile=utilities.unnamed_temp_file()
    tmpfile2=utilities.unnamed_temp_file()
    tmpfile3=utilities.unnamed_temp_file()     

    minpath_script=os.path.join(os.path.dirname(os.path.realpath(__file__)),
        config.minpath_script)
    args=[minpath_script,"-any",reactions_file]
    args+=["-map",metacyc_datafile,"-report",tmpfile3]
    args+=["-details",tmpfile,"-mps",tmpfile2]
        
    command=[sys.executable,args,[reactions_file,metacyc_datafile],[],None,None,True,None]
    
    return tmpfile, command

def identify_reactions_and_pathways(gene_scores, reactions_database, pathways_database):
    """
    Identify the reactions and then pathways from the hits found
    """
            
    if config.minpath_toggle == "on":
        # Write a flat reactions to pathways file
        logger.debug("Write flat reactions to pathways file for Minpath")
        pathways_database_file=utilities.unnamed_temp_file()
        file_handle=open(pathways_database_file,"w")
        file_handle.write(pathways_database.get_database())
        file_handle.close()
    
    # Create a store for the pathways and reactions by bug
    pathways_and_reactions_store=store.PathwaysAndReactions()
    reactions={}
    
    minpath_results={}
    minpath_commands=[]
    # Run through each of the score sets by bug
    for bug in gene_scores.bug_list():
        gene_scores_for_bug=gene_scores.scores_for_bug(bug)
        
        # Merge the gene scores to reaction scores   
        message="Compute reaction scores for bug: " + bug
        logger.info(message)
        
        reactions[bug]={}
        reactions_file_lines=[]
        for reaction in reactions_database.reaction_list():
            genes_list=reactions_database.find_genes(reaction)
            abundance=0
            # Add the scores for each gene to the total score for the reaction
            for gene in genes_list:
                abundance+=gene_scores_for_bug.get(gene,0)  
            
            # Only write out reactions where the abundance is greater than 0
            if abundance>0: 
                reactions_file_lines.append(reaction+config.output_file_column_delimiter
                    +str(abundance)+"\n")
                # Store the abundance data to compile with the minpath pathways
                reactions[bug][reaction]=abundance
    
        pathways={}
        # Run minpath if toggle on and also if there is more than one reaction   
        if config.minpath_toggle == "on" and len(reactions_file_lines)>3:   
    
            # Create a temp file for the reactions results
            reactions_file=utilities.unnamed_temp_file()
            file_handle=open(reactions_file,"w")
            file_handle.write("".join(reactions_file_lines))
            file_handle.close()
        
            # Run minpath to identify the pathways
            logger.info("Run MinPath on " + bug)
                
            tmpfile, command=minpath_command(reactions_file, pathways_database_file)
            minpath_results[bug]=tmpfile
            minpath_commands.append(command)
            
    # Run through the minpath commands if minpath is to be run
    if minpath_commands:
        utilities.command_threading(config.threads,minpath_commands)
    
    # Link the pathways to reactions
    for bug in gene_scores.bug_list():
        pathways={}
        if bug in minpath_results:
            tmpfile=minpath_results[bug]
            # Process the minpath results
            if os.path.isfile(tmpfile):
        
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
            else:
                message="Empty results file from MinPath run for bug: " + bug
                print(message)
                logger.warning(message)
        else:
            # Add all pathways associated with each reaction if not using minpath
            for current_reaction in reactions.get(bug,{}):
                pathways[current_reaction]=pathways.get(
                    current_reaction, []) + pathways_database.find_pathways(current_reaction) 
         
        # Store the pathway abundance for each reaction
        for current_reaction in reactions.get(bug,{}):
            # Find the pathways associated with reaction
            for current_pathway in pathways.get(current_reaction,[""]):
                # Only store data for items with pathway names
                if current_pathway:
                    pathways_and_reactions_store.add(bug,current_reaction, current_pathway, 
                        reactions[bug][current_reaction])
   
    return pathways_and_reactions_store

def compute_pathways_coverage(pathways_and_reactions_store,pathways_database):
    """
    Compute the coverage of pathways for each bug
    """

    pathways_coverage_store=store.Pathways()
    for bug in pathways_and_reactions_store.bug_list():
    
        logger.debug("Compute pathway coverage for bug: " + bug)
        
        # Process through each pathway to compute coverage
        xipe_input=[]
        median_score_value=pathways_and_reactions_store.median_score(bug)
        
        for pathway in pathways_and_reactions_store.pathway_list(bug):
                
            reaction_scores=pathways_and_reactions_store.reaction_scores(bug,pathway)
            # Count the reactions with scores greater than the median
            count_greater_than_median=0
            for reaction, score in reaction_scores.items():
                if score > median_score_value:
                   count_greater_than_median+=1
            
            # Compute coverage
            coverage=0
            total_reactions_for_pathway=len(pathways_database.find_reactions(pathway))
            if total_reactions_for_pathway:
                coverage=count_greater_than_median/float(total_reactions_for_pathway)
            
            pathways_coverage_store.add(bug,pathway,coverage)
            xipe_input.append(config.xipe_delimiter.join([pathway,str(coverage)]))
        
        # Check config to determine if xipe should be run
        if config.xipe_toggle == "on":
            # Run xipe
            xipe_exe=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                config.xipe_script)
            
            cmmd=[xipe_exe,"--file2",config.xipe_percent]
            
            logger.info("Run xipe")
            
            # Create temp files for input and output
            infile=utilities.unnamed_temp_file()
            
            # Write the input to xipe
            file_handle=open(infile,"w")
            file_handle.write("\n".join(xipe_input))
            file_handle.close()
            
            stdout_file=utilities.unnamed_temp_file()
            stderr_file=utilities.unnamed_temp_file()
            
            utilities.execute_command(sys.executable,cmmd,[infile],[],stdout_file,
                infile,None,stderr_file)
            
            try:
                xipe_stderr=open(stderr_file,"r")
                xipe_stdout=open(stdout_file,"r")
            except EnvironmentError:
                logger.debug("Empty results file from Xipe")
                continue
            
            # Record the pathways to remove based on the xipe error messages
            pathways_to_remove=[]
            for line in xipe_stderr:
                data=line.strip().split(config.xipe_delimiter)
                if len(data) == 2:
                    pathways_to_remove.append(data[1])
            
            # Keep some of the pathways to remove based on their xipe scores
            for line in xipe_stdout:
                data=line.strip().split(config.xipe_delimiter)
                if len(data) == 2:
                    pathway, pathway_data = data
                    if pathway in pathways_to_remove:
                        score, bin = pathway_data[1:-1].split(", ")
                        if float(score) >= config.xipe_probability and int(bin) == config.xipe_bin:
                            pathways_to_remove.remove(pathway)
                    
            xipe_stderr.close()
            xipe_stdout.close()        
            
            # Remove the selected pathways
            for pathway in pathways_to_remove:
                pathways_coverage_store.delete(bug,pathway)

    return pathways_coverage_store

def compute_pathways_abundance(pathways_and_reactions_store, pathways_database):
    """
    Compute the abundance of pathways for each bug
    """
    
    # Process through each pathway for each bug to compute abundance
    pathways_abundance_store=store.Pathways()
    for bug in pathways_and_reactions_store.bug_list():
        
        logger.debug("Compute pathway abundance for bug: " + bug)
        
        for pathway in pathways_and_reactions_store.pathway_list(bug):
            
            reaction_scores=pathways_and_reactions_store.reaction_scores(bug,pathway)
            # Initialize any reactions in the pathway not found to 0
            for reaction in pathways_database.find_reactions(pathway):
                reaction_scores.setdefault(reaction, 0)
                
            # Sort the scores for all of the reactions in the pathway from low to high
            sorted_reaction_scores=sorted(reaction_scores.values())
                
            # Select the second half of the list of reaction scores
            abundance_set=sorted_reaction_scores[(len(sorted_reaction_scores)/ 2):]
            
            # Compute abundance
            abundance=sum(abundance_set)/len(abundance_set)
            
            # Store the abundance
            pathways_abundance_store.add(bug, pathway, abundance)
    
    return pathways_abundance_store
    
    
def print_pathways(pathways, file, header):
    """
    Print the pathways data to a file organized by pathway
    """
    
    logger.debug("Print pathways %s", header)
    
    delimiter=config.output_file_column_delimiter
    category_delimiter=config.output_file_category_delimiter            
 
    # Create the header
    tsv_output=["# Pathway"+ delimiter + config.file_basename]       
            
    # Print out the pathways with those with the highest scores first
    for pathway in pathways.get_pathways_double_sorted():
        all_score=pathways.get_score(pathway)
        if all_score>0:
            # Print the computation of all bugs for pathway
            tsv_output.append(pathway+delimiter+utilities.format_float_to_string(all_score))
            # Process and print per bug if selected
            if not config.remove_stratified_output:
                # Print scores per bug for pathway ordered with those with the highest values first
                for bug in pathways.get_bugs_double_sorted(pathway):
                    score=pathways.get_score_for_bug(bug,pathway)
                    if score>0:
                        tsv_output.append(pathway+category_delimiter+bug+delimiter
                            +utilities.format_float_to_string(score))
 
    if config.output_format == "biom":
        # Open a temp file if a conversion to biom is selected
        tmpfile=utilities.unnamed_temp_file()
        file_out=open(tmpfile,"w")
        file_out.write("\n".join(tsv_output))
        file_out.close()
        
        utilities.tsv_to_biom(tmpfile,file,"Pathway")
            
    else:
        # Write the final file as tsv
        file_handle=open(file,"w") 
        file_handle.write("\n".join(tsv_output))                  
        file_handle.close()
    

def compute_pathways_abundance_and_coverage(pathways_and_reactions_store, pathways_database):
    """
    Compute the abundance and coverage of the pathways
    """
    
    # Compute abundance for all pathways
    pathways_abundance=compute_pathways_abundance(pathways_and_reactions_store,
        pathways_database)

    # Print the pathways abundance data to file
    print_pathways(pathways_abundance, config.pathabundance_file, "Abundance (reads per kilobase)")

    # Compute coverage for all pathways
    pathways_coverage=compute_pathways_coverage(pathways_and_reactions_store,
        pathways_database)
    
    # Print the pathways abundance data to file
    print_pathways(pathways_coverage, config.pathcoverage_file, "Coverage")

    return config.pathabundance_file, config.pathcoverage_file
