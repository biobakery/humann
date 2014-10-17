"""
Compute alignments by gene family
"""

import os
import logging

import config
import utilities

# name global logging instance
logger=logging.getLogger(__name__)

def gene_families(alignments):
    """
    Compute the gene families from the alignments
    """
    
    logger.debug("Compute gene families")
    
    # Compute hits by gene family
    gene_scores={}
    gene_output_lines={}
    for gene in alignments.gene_list():
        hit_list=alignments.hits_for_gene(gene)
        total_hits_for_gene=len(hit_list)
        
        if gene in gene_scores:
            gene_scores[gene]+=total_hits_for_gene
            gene_output_lines[gene].append(gene+config.output_file_column_delimiter+
            str(total_hits_for_gene))
        else:
            gene_scores[gene]=total_hits_for_gene
            gene_output_lines[gene]=[gene+config.output_file_column_delimiter+
            str(total_hits_for_gene)]
        
        # merge hits with the same bug
        hits_by_bug={}
        for hit in hit_list:
            bug, reference, query, evalue=hit
            hits_by_bug[bug]=hits_by_bug.get(bug,0)+1          
        
        # record the hits by bug for this gene
        for bug in sorted(hits_by_bug, key=hits_by_bug.get, reverse=True):
            gene_output_lines[gene].append(gene+config.output_file_category_delimiter+
                bug+config.output_file_column_delimiter+str(hits_by_bug[bug]))
        
    # Write the scores ordered with the top first
    tsv_output=["#Gene Family"+config.output_file_column_delimiter+"Hits"]
    for gene in sorted(gene_scores, key=gene_scores.get, reverse=True):
        tsv_output.append("\n".join(gene_output_lines[gene]))    
        
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

