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

from . import chi2cdf

from .. import utilities
from .. import config
from .. import store

# name global logging instance
logger=logging.getLogger(__name__)

def minpath_command(reactions_file,metacyc_datafile):
    """
    Return the minpath command and the name of the output file
    """
    
    # Create temp files for the results
    tmpfile=utilities.unnamed_temp_file()
    tmpfile2=utilities.unnamed_temp_file()
    tmpfile3=utilities.unnamed_temp_file()     
    tmpfile4=utilities.unnamed_temp_file()  

    minpath_script=os.path.join(os.path.dirname(os.path.abspath(__file__)),
        config.minpath_script)
    args=[minpath_script,"-any",reactions_file]
    args+=["-map",metacyc_datafile,"-report",tmpfile3]
    args+=["-details",tmpfile,"-mps",tmpfile2]
        
    command=[sys.executable,args,[reactions_file,metacyc_datafile],[],tmpfile4,None,True,None]
    
    return tmpfile, command

def xipe_command(infile):
    """
    Return the xipe command and the name of the output files
    """
    
    xipe_exe=os.path.join(os.path.dirname(os.path.abspath(__file__)),
    config.xipe_script)
            
    args=[xipe_exe,"--file1",infile,"--file2",config.xipe_percent]
            
    stdout_file=utilities.unnamed_temp_file()
    stderr_file=utilities.unnamed_temp_file()
    
    command=[sys.executable,args,[infile],[],stdout_file,None,True,stderr_file]
    
    return stdout_file, stderr_file, command
    

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
        if reactions_database:
            for reaction in sorted(reactions_database.reaction_list()):
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
        else:
            for gene in gene_scores_for_bug:
                score=gene_scores_for_bug[gene]
                
                if score>0:
                    reactions_file_lines.append(gene+config.output_file_column_delimiter
                        +str(score)+"\n")
                    # Store the abundance data to compile with the minpath pathways
                    reactions[bug][gene]=score
    
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
        
                file_handle_read=open(tmpfile, "rt")
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
    xipe_stdout_results={}
    xipe_stderr_results={}
    xipe_commands=[]
    for bug in pathways_and_reactions_store.bug_list():
    
        logger.debug("Compute pathway coverage for bug: " + bug)
        
        # Process through each pathway to compute coverage
        xipe_input=[]
        median_score_value=pathways_and_reactions_store.median_score(bug)
        
        for pathway in pathways_and_reactions_store.pathway_list(bug):
                
            reaction_scores=pathways_and_reactions_store.reaction_scores(bug,pathway)
            
            # Check if the pathways database is structured
            if pathways_database.is_structured():
                structure=pathways_database.get_structure_for_pathway(pathway)
                key_reactions=pathways_database.get_key_reactions_for_pathway(pathway)
                # Apply gap fill
                reaction_scores=gap_fill(key_reactions, reaction_scores)
                # Compute the structured pathway coverage
                coverage=compute_structured_pathway_abundance_or_coverage(structure,
                    key_reactions,reaction_scores,True,median_score_value)
            else:
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

            # Create temp file for input
            infile=utilities.unnamed_temp_file()
            
            # Write the input to xipe
            file_handle=open(infile,"w")
            file_handle.write("\n".join(xipe_input))
            file_handle.close()
            
            stdout_file, stderr_file, command = xipe_command(infile)
            
            xipe_commands.append(command)
            xipe_stdout_results[bug]=stdout_file
            xipe_stderr_results[bug]=stderr_file
            
    # Run xipe
    if xipe_commands:
        utilities.command_threading(config.threads,xipe_commands)
            
    # Process the xipe output
    for bug in xipe_stdout_results:
            
        try:
            xipe_stderr=open(xipe_stderr_results[bug],"rt")
            xipe_stdout=open(xipe_stdout_results[bug],"rt")
        except EnvironmentError:
            logger.debug("Empty results file from xipe")
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

def harmonic_mean(values):
    """
    Return the harmonic mean for the values
    """
    
    # If there are no values or if one of the values is zero, then the harmonic mean is zero
    mean=0
    if values and min(values) > 0:
        reciprocal_sum=sum((1.0/v) for v in values)
        mean=len(values)/reciprocal_sum
    
    return mean

def compute_structured_pathway_abundance_or_coverage(structure, key_reactions, reaction_scores, 
    coverage_computation, median_value):
    """
    Compute the abundance or coverage for a structured pathway
    """
    
    # Process through the structure to compute the abundance
    required_reaction_abundances=[]
    optional_reaction_abundances=[]
    # Select the join instead of removing from the list to not alter the list for
    # the calling function
    join=structure[0]
    for item in structure[1:]:
        if isinstance(item, list):
            required_reaction_abundances.append(compute_structured_pathway_abundance_or_coverage(item, 
                key_reactions, reaction_scores, coverage_computation, median_value))
        else:
            score=reaction_scores.get(item,0)
                
            # Update the score for the reaction if this is a coverage computation
            if coverage_computation:
                score=chi2cdf.chi2cdf(score,median_value)

            # Check if this is an optional reaction
            if item in key_reactions:
                required_reaction_abundances.append(score)
            else:
                optional_reaction_abundances.append(score)
    
    # If this is an OR join then use the max of all of the reaction abundances
    if join == config.pathway_OR:
        all_reaction_abundances=required_reaction_abundances + optional_reaction_abundances
        abundance=0
        if all_reaction_abundances:
            abundance = max(all_reaction_abundances)
    else:
        # If this is not an OR, then take the harmonic mean of the reactions
        abundance=harmonic_mean(required_reaction_abundances)
        # Add the optional reactions if they are present
        if optional_reaction_abundances:
            # Filter the optional abundances to only include those that are greater than the abundance
            # from the required reactions
            optional_reaction_abundances_filtered=[value for value in optional_reaction_abundances if value > abundance]
            abundance=harmonic_mean(required_reaction_abundances + optional_reaction_abundances_filtered)
        
    return abundance

def gap_fill(key_reactions, reaction_scores):
    """
    If all but one of the key reactions have abundance scores, then fill gap
    Boost the lowest abundance score
    """
    
    reaction_scores_gap_filled=reaction_scores.copy()
    
    # do not apply gap fill, if set to off
    if config.gap_fill_toggle == "off":
        return reaction_scores_gap_filled
    
    # get the scores for all of the key reactions
    key_reactions_nonzero_scores=[]
    for reaction in key_reactions:
        score=reaction_scores.get(reaction,0)
        if score > 0:
            key_reactions_nonzero_scores.append(score)

    if len(key_reactions)-len(key_reactions_nonzero_scores) == 1:
        # fill single zero gap with lowest key reaction score
        min_score=min(key_reactions_nonzero_scores)
        for reaction in key_reactions:
            score=reaction_scores.get(reaction,0)
            if score == 0:
                reaction_scores_gap_filled[reaction]=min_score
    elif len(key_reactions)-len(key_reactions_nonzero_scores) == 0:
        # boost lowest abundance score
        sorted_key_reactions_nonzero_scores=sorted(key_reactions_nonzero_scores)
        for reaction in key_reactions:
            score=reaction_scores.get(reaction,0)
            if score == sorted_key_reactions_nonzero_scores[0]:
                try:
                    reaction_scores_gap_filled[reaction]=sorted_key_reactions_nonzero_scores[1]
                except IndexError:
                    pass

    return reaction_scores_gap_filled
    

def compute_pathways_abundance(pathways_and_reactions_store, pathways_database):
    """
    Compute the abundance of pathways for each bug
    Also find the set of the reactions with abundance in all pathways present
    """
    
    # Store the reactions which have abundance in the pathways with abundance
    reactions_in_pathways_present={}
    
    # Process through each pathway for each bug to compute abundance
    pathways_abundance_store=store.Pathways()
    for bug in pathways_and_reactions_store.bug_list():
        
        logger.debug("Compute pathway abundance for bug: " + bug)
        
        reactions_in_pathways_present[bug]=set()
        for pathway in pathways_and_reactions_store.pathway_list(bug):
            
            reaction_scores=pathways_and_reactions_store.reaction_scores(bug,pathway)
            
            # Check if the pathways database is structured
            if pathways_database.is_structured():
                structure=pathways_database.get_structure_for_pathway(pathway)
                key_reactions=pathways_database.get_key_reactions_for_pathway(pathway)
                # Apply gap fill
                reaction_scores_gap_filled=gap_fill(key_reactions, reaction_scores)
                # Compute the structured pathway abundance
                abundance=compute_structured_pathway_abundance_or_coverage(structure,
                    key_reactions,reaction_scores_gap_filled,False,0)
            
            else:
                # Initialize any reactions in the pathway not found to 0
                for reaction in pathways_database.find_reactions(pathway):
                    reaction_scores.setdefault(reaction, 0)
                    
                # Sort the scores for all of the reactions in the pathway from low to high
                sorted_reaction_scores=sorted(reaction_scores.values())
                    
                # Select the second half of the list of reaction scores
                abundance_set=sorted_reaction_scores[int(len(sorted_reaction_scores)/ 2):]
                
                # Compute abundance
                abundance=sum(abundance_set)/len(abundance_set)
                
            # If this pathway is present, store those reactions with abundance
            if abundance > 0:
                for reaction,score in reaction_scores.items():
                    if score > 0:
                        reactions_in_pathways_present[bug].add(reaction)
            
            # Store the abundance
            pathways_abundance_store.add(bug, pathway, abundance)
    
    return pathways_abundance_store, reactions_in_pathways_present
    
    
def print_pathways(pathways, file, header, pathway_names, sorted_pathways_and_bugs,
                   unmapped_all, unintegrated_all, unintegrated_per_bug):
    """
    Print the pathways data to a file organized by pathway
    """
    
    logger.debug("Print pathways %s", header)
    
    delimiter=config.output_file_column_delimiter
    category_delimiter=config.output_file_category_delimiter            
 
    # Create the header
    column_name=config.file_basename + header
    if config.remove_column_description_output:
        column_name=config.file_basename
    tsv_output=["# Pathway"+ delimiter + column_name]
    
    # Add the unmapped and unintegrated values
    tsv_output.append(config.unmapped_pathway_name+delimiter+utilities.format_float_to_string(unmapped_all))
    tsv_output.append(config.unintegrated_pathway_name+delimiter+utilities.format_float_to_string(unintegrated_all))
    # Process and print per bug if selected
    if not config.remove_stratified_output:
        for bug in utilities.double_sort(unintegrated_per_bug):
            tsv_output.append(config.unintegrated_pathway_name+category_delimiter+bug
                              +delimiter+utilities.format_float_to_string(unintegrated_per_bug[bug]))
            
    # Print out all pathways sorted
    for pathway, bugs_list in sorted_pathways_and_bugs:
        all_score=pathways.get_score(pathway)
        pathway_name=pathway_names.get_name(pathway)
        # Print the computation of all bugs for pathway
        tsv_output.append(pathway_name+delimiter+utilities.format_float_to_string(all_score))
        # Process and print per bug if selected
        if not config.remove_stratified_output:
            # Print scores for all bugs sorted
            for bug in bugs_list:
                score=pathways.get_score_for_bug(bug,pathway)
                tsv_output.append(pathway_name+category_delimiter+bug+delimiter
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
    
def compute_gene_abundance_in_pathways(gene_scores, reactions_database, reactions_in_pathways_present):
    """
    Compute the abundance of genes present in pathways found
    Also compute the remaining gene abundance that did not contribute to any pathways present
    """
    
    # From each of the reactions present in pathways, find the list of all genes
    # that contributed to the pathway abundance (for each bug found in community)
    genes_in_pathways_present={}
    for bug in reactions_in_pathways_present:
        genes_in_pathways_present[bug]=set()
        for reaction in reactions_in_pathways_present[bug]:
            if reactions_database:
                # Find the list of genes for this reaction
                gene_list=reactions_database.find_genes(reaction)
                genes_in_pathways_present[bug].update(gene_list)
            else:
                # If the reactions database is not provided, then the pathway is defined
                # in terms of genes
                genes_in_pathways_present[bug].add(reaction)
        
    # Compute the abundance for the genes present
    gene_abundance_in_pathways={}
    remaining_gene_abundance={}
    for bug in genes_in_pathways_present:
        for gene in genes_in_pathways_present[bug]:
            gene_abundance_in_pathways[bug]=gene_abundance_in_pathways.get(bug,0)+gene_scores.get_score(bug,gene)
            
        # Compute the remaining gene abundance
        total_gene_abundance=sum(gene_scores.scores_for_bug(bug).values())
        remaining_gene_abundance[bug]=total_gene_abundance-gene_abundance_in_pathways.get(bug,0)
        
    return gene_abundance_in_pathways, remaining_gene_abundance

def compute_unmapped_and_unintegrated(gene_abundance_in_pathways, remaining_gene_abundance, unaligned_reads_count, pathways_abundance):
    """
    Compute the unmapped and unintegrated pathway values
    The compression constant is defined as the total abundance of all pathways divided by the total abundance of genes in pathways
    """
    
    # Compute the overall compression constant
    pathways_list=pathways_abundance.get_pathways_list()
    total_abundance_all_pathways=sum([pathways_abundance.get_score(pathway) for pathway in pathways_list])
    try:
        compression_all=total_abundance_all_pathways / ( gene_abundance_in_pathways.get("all",0) * 1.0 )
    except ZeroDivisionError:
        compression_all=0
    
    # Compute unmapped
    unmapped_all=compression_all * unaligned_reads_count
    
    # Compute unintegrated
    unintegrated_all=compression_all * remaining_gene_abundance.get("all",0)
    unintegrated_per_bug={}
    for bug in pathways_abundance.get_bugs_list():
        total_abundance_all_pathways=sum([pathways_abundance.get_score_for_bug(bug,pathway) for pathway in pathways_list])
        try:
            compression = total_abundance_all_pathways / ( gene_abundance_in_pathways.get(bug,0) * 1.0 )
        except ZeroDivisionError:
            compression=0 
        
        unintegrated_per_bug[bug]=compression * remaining_gene_abundance.get(bug,0)
    
    return unmapped_all, unintegrated_all, unintegrated_per_bug
    

def compute_pathways_abundance_and_coverage(gene_scores, reactions_database, 
                                            pathways_and_reactions_store, pathways_database, unaligned_reads_count):
    """
    Compute the abundance and coverage of the pathways
    """
    
    # Read in and store the pathway id to name mappings
    pathway_names=store.Names(config.pathway_name_mapping_file)
    
    # Compute abundance for all pathways
    pathways_abundance, reactions_in_pathways_present=compute_pathways_abundance(
        pathways_and_reactions_store, pathways_database)
    
    # Compute the abundance of genes in pathways and not in pathways
    gene_abundance_in_pathways, remaining_gene_abundance=compute_gene_abundance_in_pathways(
        gene_scores, reactions_database, reactions_in_pathways_present)
    
    # Compute the unmapped and unintegrated values
    unmapped_all, unintegrated_all, unintegrated_per_bug=compute_unmapped_and_unintegrated(
        gene_abundance_in_pathways, remaining_gene_abundance, unaligned_reads_count, pathways_abundance)
    
    # Compute coverage for all pathways
    pathways_coverage=compute_pathways_coverage(pathways_and_reactions_store,
        pathways_database)

    # Get the sorted list of pathways and bugs from the abundance values
    # This same sorting will be used for both the abundance and coverage output files
    sorted_pathways_and_bugs=pathways_abundance.get_pathways_and_bugs_nonzero_sorted()

    # Print the pathways abundance data to file
    print_pathways(pathways_abundance, config.pathabundance_file, "_Abundance", 
                   pathway_names, sorted_pathways_and_bugs, unmapped_all,
                   unintegrated_all, unintegrated_per_bug)
    
    # Set the unmapped and unintegrated to one for coverage output
    unmapped_all=1
    unintegrated_all=1
    unintegrated_per_bug={bug:1 for bug in unintegrated_per_bug.keys()}
    
    # Print the pathways coverage data to file
    print_pathways(pathways_coverage, config.pathcoverage_file, "_Coverage", 
                   pathway_names, sorted_pathways_and_bugs, unmapped_all,
                   unintegrated_all, unintegrated_per_bug)

    return config.pathabundance_file, config.pathcoverage_file
