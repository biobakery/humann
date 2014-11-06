"""
HUMAnN2: quantify_families module
Compute alignments by gene family

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
import logging
import math

import config
import utilities

# name global logging instance
logger=logging.getLogger(__name__)

def compute_gene_scores(hits):
    """
    Compute the scores for each gene for a set of hits
    """
 
    # Loop through hits to organize scores by query
    genes_by_query={}
    total_scores_by_query={}
    for hit in hits:
        try:
            bug, reference, reference_length, query, evalue=hit
            score=math.exp(-evalue)
            if query in genes_by_query:
                genes_by_query[query].append([reference,score,reference_length])
            else:
                genes_by_query[query]=[[reference,score,reference_length]]
            total_scores_by_query[query]=total_scores_by_query.get(query,0)+score
        except ValueError:
            logger.debug("Tried to use a hit that is not complete")
    
    # add these scores by query to the total gene scores
    total_gene_scores={}
    for query, genes in genes_by_query.items():
        for gene, score, gene_length in genes:
            total_gene_scores[gene]=total_gene_scores.get(gene,0)+(score/total_scores_by_query[query])/gene_length    
    
    return total_gene_scores

def compute_gene_scores_by_bug(args):
    """
    Compute the gene scores for a specific bug
    """
    
    hits, bug = args
    
    return [bug, compute_gene_scores(hits)]
    

def gene_families(alignments):
    """
    Compute the gene families from the alignments
    """
    
    logger.debug("Compute gene families")
    
    # Compute scores for each gene family per bug
    args=[]
    for bug in alignments.bug_list():
        hits=alignments.hits_for_bug(bug)
        args.append([hits, bug])
        
    # also run for all bugs
    hits=alignments.all_hits()
    args.append([hits, "all"])
    
    gene_scores_array=utilities.command_multiprocessing(config.threads, args, 
        function=compute_gene_scores_by_bug)
    
    # find the all scores and compile by gene
    all_scores={}
    all_scores_by_bug={}
    for bug, scores in gene_scores_array:
        if bug == "all":
            all_scores=scores
        else:
            for gene, score in scores.items():
                if score>0:
                    if not gene in all_scores_by_bug:
                        all_scores_by_bug[gene]={bug: score}
                    else:
                        all_scores_by_bug[gene][bug]=score 
        
    # Write the scores ordered with the top first
    tsv_output=["# Gene Family"+config.output_file_column_delimiter+"Abundance (reads per kilobase)"]
    
    delimiter=config.output_file_column_delimiter
    category_delimiter=config.output_file_category_delimiter

    # Print out the gene families with those with the highest scores first
    for gene in utilities.double_sort(all_scores):
        all_score=all_scores[gene]
        if all_score>0:
            # Print the computation of all bugs for gene family
            tsv_output.append(gene+delimiter+str(all_score))
            # Print scores per bug for family ordered with those with the highest values first
            if gene in all_scores_by_bug:
                for bug in utilities.double_sort(all_scores_by_bug[gene]):
                    tsv_output.append(gene+category_delimiter+bug+delimiter
                                      +str(all_scores_by_bug[gene][bug]))       
        
    if config.output_format=="biom":
        # Open a temp file if a conversion to biom is selected
        tmpfile=utilities.unnamed_temp_file()
        file_handle=open(tmpfile,'w')
        file_handle.write("\n".join(tsv_output))
        file_handle.close()
        
        utilities.tsv_to_biom(tmpfile,config.genefamilies_file)
        
    else:
        # Write output as tsv format
        file_handle = open(config.genefamilies_file, "w")
        file_handle.write("\n".join(tsv_output))
        file_handle.close()

    return config.genefamilies_file

