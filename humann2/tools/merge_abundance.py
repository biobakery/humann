#!/usr/bin/env python

"""
Merge a gene families and pathway abundance file

This module will join gene and pathway tables output by HUMAnN2. 

To Run: 
$ ./merge_abundance.py --input-genes genefamilies.tsv --input-pathways pathabundance.tsv --output merged.tsv

"""

import argparse
import sys
import os
import gzip

import util

TABLE_DELIMITER="\t"
BUG_DELIMITER="|"
GENE_NAME_DELIMITER=": "
EC_DELIMITER=","

try:
    from humann2 import config
except ImportError:
    sys.exit("ERROR: Unable to find the HUMAnN2 install.")
    
def merge_abundances(gene_table,pathways_to_genes,input_pathways,output,additional_gene_info):
    """
    Read through the pathway abundances file
    Write the merged abundances
    """
    
    file_handle,header,line=util.process_gene_table_header(input_pathways, allow_for_missing_header=True)
    
    # open the output file
    try:
        write_file_handle=open(output,"w")
    except EnvironmentError:
        sys.exit("ERROR: Unable to write to output file: " + output)
        
    if header:
        write_file_handle.write(header)
        
    # read through the pathways mapping file, then merge in the gene families information
    # read thorugh the file first to allow for unordered input files
    pathways={}
    pathways_by_bug={}
    while line:
        data=line.rstrip().split(TABLE_DELIMITER)
        bug=""
        pathway=data[0]
        if BUG_DELIMITER in line:
            pathway,bug=data[0].split(BUG_DELIMITER)
        
        # remove pathway name if present
        if GENE_NAME_DELIMITER in pathway:
            pathway=pathway.split(GENE_NAME_DELIMITER)[0]
       
        try:
            data_point=float(data[-1])
        except ValueError:
            data_point=0
 
        if data_point > 0:
            if bug:
                if not pathway in pathways_by_bug:
                    pathways_by_bug[pathway]=[]
                pathways_by_bug[pathway].append([bug,data[-1]])
            else:
                pathways[pathway]=data[-1]

        line=file_handle.readline()
                        
    file_handle.close()

    for pathway in sorted(pathways):
        write_file_handle.write(TABLE_DELIMITER.join([pathway,pathways[pathway]])+"\n")
        for item in pathways_by_bug.get(pathway,[]):
            bug, data = item
            output_line=TABLE_DELIMITER.join([BUG_DELIMITER.join([pathway,bug]), data])
            write_file_handle.write(output_line+"\n")
            # check for gene families for this pathway for this bug
            for gene in pathways_to_genes.get(pathway,[]):
                if gene in gene_table:
                    if bug in gene_table[gene]:
                        # check if additional information should be added to the gene
                        gene_with_info=gene
                        if additional_gene_info:
                            if gene in additional_gene_info:
                                gene_with_info+=GENE_NAME_DELIMITER+additional_gene_info[gene]
                        # print out the gene for this bug
                        output_line=TABLE_DELIMITER.join([BUG_DELIMITER.join([pathway,bug,gene_with_info]),gene_table[gene][bug]])
                        write_file_handle.write(output_line+"\n")
    write_file_handle.close()
    
def read_mapping(gene_mapping_file,pathway_mapping_file):
    """
    Read the gene to reaction and then reaction to pathway mappings
    """
    
    # read the gene to reaction mappings
    reactions_to_genes={} 
    genes_to_ecs={}
    try:
        if gene_mapping_file.endswith(".gz"):
            file_handle=gzip.open(gene_mapping_file)
        else:
            file_handle=open(gene_mapping_file)
    except EnvironmentError:
        sys.exit("ERROR: Unable to open the mapping file: " + gene_mapping_file)
            
    # read through the mapping lines
    for line in file_handle:
        data=line.rstrip().split(TABLE_DELIMITER)
        reaction=data.pop(0)
        ecs=data.pop(0).split(EC_DELIMITER)
        reactions_to_genes[reaction]=data
        for gene in data:
            if not gene in genes_to_ecs:
                genes_to_ecs[gene]=set()
            genes_to_ecs[gene].update(ecs)
                
    file_handle.close()
    
    # convert the genes to ecs from sets to strings
    for gene in genes_to_ecs:
        genes_to_ecs[gene]=EC_DELIMITER.join(list(genes_to_ecs[gene]))
    
    # read the reaction to pathway mappings
    pathways_to_genes={}
    try:
        file_handle=open(pathway_mapping_file)
    except EnvironmentError:
        sys.exit("ERROR: Unable to open the mapping file: " + pathway_mapping_file)
        
    for line in file_handle:
        info=line.rstrip().split(TABLE_DELIMITER)
        if len(info) == 2:
            pathway, structure = info
            # get the reactions from the structure
            genes=set()
            for item in structure.split(" "):
                genes.update(reactions_to_genes.get(item,[]))
            pathways_to_genes[pathway]=genes
        
    file_handle.close()
  
    return pathways_to_genes,genes_to_ecs
        
def read_gene_table(gene_table):
    """
    Read the input gene table
    Store the abundances related to each bug
    """
    
    gene_table_data={}
    gene_names={}
    
    file_handle,header,line=util.process_gene_table_header(gene_table, allow_for_missing_header=True)
    
    while line:
        if BUG_DELIMITER in line:
            data=line.rstrip().split(TABLE_DELIMITER)
            try:
                gene,bug=data[0].split(BUG_DELIMITER)
                # remove the gene name if present and store
                if GENE_NAME_DELIMITER in gene:
                    info=gene.split(GENE_NAME_DELIMITER)
                    if len(info) > 1 :
                        gene = info[0]
                        gene_names[gene]=" ".join(info[1:])
                data_point=data[1]
            except IndexError:
                gene=""

            if gene:
                if not gene in gene_table_data:
                    gene_table_data[gene]={}
                gene_table_data[gene][bug]=data_point

        line=file_handle.readline()
            
    file_handle.close()
        
    return gene_table_data, gene_names

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Merge gene and pathway abundance tables\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--input-genes",
        help="the gene family abundance file\n",
        required=True)
    parser.add_argument(
        "--input-pathways",
        help="the pathway abundance file\n",
        required=True)
    parser.add_argument(
        "-a","--add",
        help="add gene information\n",
        choices=["ec","name"])
    parser.add_argument(
        "--gene-mapping",
        help="gene family to reaction mapping file\n")
    parser.add_argument(
        "--pathway-mapping",
        help="reaction to pathway mapping file\n")
    parser.add_argument(
        "-o","--output",
        help="the table to write\n",
        required=True)

    return parser.parse_args()


def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    # get the absolute paths
    input_genes=os.path.abspath(args.input_genes)
    input_pathways=os.path.abspath(args.input_pathways)
    
    output=os.path.abspath(args.output)
    
    # check input files
    for file in [input_genes,input_pathways]:
        if not os.path.isfile(file):
            sys.exit("The input provided is not a file. Please enter a new file. " + file)
        
        if not os.access(file, os.R_OK):
            sys.exit("The input provided is not readable. Please update the permissions. " + file)
        
    print("Reading the gene abundances.")
    gene_table,gene_names=read_gene_table(input_genes)
    
    # use the default mapping if not set
    if args.gene_mapping:
        args.gene_mapping=os.path.abspath(args.gene_mapping)
    else:
        args.gene_mapping=config.metacyc_gene_to_reactions
        
    if args.pathway_mapping:
        args.pathway_mapping=os.path.abspath(args.pathway_mapping)
    else:
        args.pathway_mapping=config.metacyc_reactions_to_pathways
    
    print("Reading the gene to pathway mapping.")
    pathways_to_genes,genes_to_ecs=read_mapping(args.gene_mapping,args.pathway_mapping)
    
    additional_info=None
    if args.add == "name":
        additional_info=gene_names
    elif args.add == "ec":
        additional_info=genes_to_ecs
    
    print("Merging abundances.")
    merge_abundances(gene_table,pathways_to_genes,input_pathways,output,additional_info)
    

if __name__ == "__main__":
    main()
