"""
Compute alignments by gene family
"""

import config, store

def gene_families(alignments):
    """
    Compute the gene families from the alignments
    """
    
    file_handle = open(config.genefamilies_file, "w")

    file_handle.write("Gene Family"+config.output_file_column_delimiter+"Hits\n")
    
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
            bug=hit.get_bug()          
            if bug in hits_by_bug:
                hits_by_bug[bug]+=1
            else:
                hits_by_bug[bug]=1
        
        # record the hits by bug for this gene
        for bug in hits_by_bug.keys():
            gene_output_lines[gene].append(gene+config.output_file_category_delimiter+
                bug+config.output_file_column_delimiter+str(hits_by_bug[bug]))
        
    # Write the scores ordered with the top first
    for gene in sorted(gene_scores, key=gene_scores.get, reverse=True):
        file_handle.write("\n".join(gene_output_lines[gene])+"\n")    
        
    file_handle.close()

    return config.genefamilies_file

