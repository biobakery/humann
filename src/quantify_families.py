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
    for gene in alignments.gene_list():
        hit_list=alignments.hits_for_gene(gene)
        total_hits_for_gene=len(hit_list)
        
        file_handle.write(gene+config.output_file_column_delimiter+
            str(total_hits_for_gene)+"\n")
        
        # merge hits with the same bug
        hits_by_bug={}
        index=alignments.find_index("bug")
        for hit in hit_list:
            bug=hit[index]
            
            if bug in hits_by_bug:
                hits_by_bug[bug]+=1
            else:
                hits_by_bug[bug]=1
        
        # print the hits by bug for this gene
        for bug in hits_by_bug.keys():
            file_handle.write(gene+config.output_file_category_delimiter+
                bug+config.output_file_column_delimiter+str(hits_by_bug[bug])+"\n")
        
    file_handle.close()

    return config.genefamilies_file

