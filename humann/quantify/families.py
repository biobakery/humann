"""
HUMAnN: quantify_families module
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

from .. import config
from .. import utilities
from .. import store

# name global logging instance
logger=logging.getLogger(__name__)

def gene_families(alignments,gene_scores,unaligned_reads_count):
    """
    Compute the gene families from the alignments
    """
    
    logger.debug("Compute gene families")
    
    # Compute scores for each gene family for each bug set
    alignments.convert_alignments_to_gene_scores(gene_scores)
        
    # Process the gene id to names mappings
    gene_names=store.Names(config.gene_family_name_mapping_file)
     
    delimiter=config.output_file_column_delimiter
    category_delimiter=config.output_file_category_delimiter     

    # Write the scores ordered with the top first
    column_name=config.file_basename+"_Abundance-RPKs"
    if config.remove_column_description_output:
        column_name=config.file_basename
    tsv_output=["# Gene Family"+delimiter+column_name]
    
    # Add the unaligned reads count
    tsv_output.append(config.unmapped_gene_name+delimiter+utilities.format_float_to_string(unaligned_reads_count))  

    # Print out the gene families with those with the highest scores first
    for gene in gene_scores.gene_list_sorted_by_score("all"):
        all_score=gene_scores.get_score("all",gene)
        if all_score>0:
            gene_name=gene_names.get_name(gene)
            # Print the computation of all bugs for gene family
            tsv_output.append(gene_name+delimiter+utilities.format_float_to_string(all_score))
            # Process and print per bug if selected
            if not config.remove_stratified_output:
                # Print scores per bug for family ordered with those with the highest values first
                scores_by_bug=gene_scores.get_scores_for_gene_by_bug(gene)
                for bug in utilities.double_sort(scores_by_bug):
                    if scores_by_bug[bug]>0:
                        tsv_output.append(gene_name+category_delimiter+bug+delimiter
                            +utilities.format_float_to_string(scores_by_bug[bug]))       
        
    if config.output_format=="biom":
        # Open a temp file if a conversion to biom is selected
        tmpfile=utilities.unnamed_temp_file()
        file_handle=open(tmpfile,'w')
        file_handle.write("\n".join(tsv_output))
        file_handle.close()
        
        utilities.tsv_to_biom(tmpfile,config.genefamilies_file,"Gene")
        
    else:
        # Write output as tsv format
        file_handle = open(config.genefamilies_file, "w")
        file_handle.write("\n".join(tsv_output))
        file_handle.close()

    return config.genefamilies_file

