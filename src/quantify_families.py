"""
Identify hits from alignment steps
Compute gene family abundance
"""

import os, tempfile
import config, utilities, store

def hits(bam_output, blast_output):
    """
    Identify the hits from the alignment results
    """

    file_handle, hits_file=tempfile.mkstemp()
    os.close(file_handle)

    # find hits for the bam output
    if bam_output != "Empty":
        # call humann1 bam2hits.py
        utilities.execute_command(
            config.humann1_script_bam_to_hits,
            [], [bam_output], [], hits_file, bam_output)
   
        # call humann1 blast2hits.py with results
        if os.stat(hits_file).st_size:
            file_handle, hits_file2=tempfile.mkstemp()
            os.close(file_handle)
            utilities.execute_command(
                config.humann1_script_blast_to_hits,
                ["blastx",0,hits_file],[blast_output],[],
                    hits_file2, blast_output)
            utilities.remove_file(hits_file)
            hits_file=hits_file2
        else:
            if os.stat(blast_output).st_size:
                utilities.execute_command(
                    config.humann1_script_blast_to_hits,
                    [],[blast_output],[],hits_file, blast_output)

    else:
        # call humann1 blast2hits.py
        if os.stat(blast_output).st_size:
            utilities.execute_command(
                config.humann1_script_blast_to_hits,
                [],[blast_output],[],hits_file, blast_output)
        
    return hits_file

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

