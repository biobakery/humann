#!/usr/bin/env python

"""
Merge a gene families and pathway abundance file

This module will join gene and pathway tables output by HUMAnN. 

To Run: 
$ ./merge_abundance.py --input-genes genefamilies.tsv --input-pathways pathabundance.tsv --output merged.tsv

"""

import argparse
import sys
import os
import gzip
import bz2
import re

from humann.tools import util

TABLE_DELIMITER="\t"
BUG_DELIMITER="|"
GENE_NAME_DELIMITER=": "
EC_DELIMITER=","

try:
    from humann import config
except ImportError:
    sys.exit("ERROR: Unable to find the HUMAnN install.")
    
def merge_abundances(gene_table,pathways_to_genes,input_pathways,output,additional_gene_info,remove_taxonomy):
    """
    Read through the pathway abundances file
    Write the merged abundances
    """
    
    lines=util.process_gene_table_with_header(input_pathways, allow_for_missing_header=True)
    header=next(lines)
    
    # open the output file
    try:
        write_file_handle=open(output,"w")
    except EnvironmentError:
        sys.exit("ERROR: Unable to write to output file: " + output)
        
    if header:
        write_file_handle.write(header+"\n")
        
    # read through the pathways mapping file, then merge in the gene families information
    # read thorugh the file first to allow for unordered input files
    pathways={}
    pathways_by_bug={}
    pathway_names={}
    for line in lines:
        data=line.split(TABLE_DELIMITER)
        bug=""
        pathway=data[0]
        if BUG_DELIMITER in line:
            pathway,bug=data[0].split(BUG_DELIMITER)
        
        # remove pathway name if present
        if GENE_NAME_DELIMITER in pathway:
            info=pathway.split(GENE_NAME_DELIMITER)
            pathway = info[0]
            pathway_names[pathway]=" ".join(info[1:])
       
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

    for pathway in sorted(pathways):
        # check if the name should be added to the pathway
        pathway_with_name=pathway
        if pathway in pathway_names:
            pathway_with_name+=GENE_NAME_DELIMITER+pathway_names[pathway]
        
        write_file_handle.write(TABLE_DELIMITER.join([pathway_with_name,pathways[pathway]])+"\n")
        
        # do not write the taxonomy if set
        if remove_taxonomy:
            for gene in pathways_to_genes.get(pathway,[]):
                if gene in gene_table:
                    # check if additional information should be added to the gene
                    gene_with_info=gene
                    if gene in additional_gene_info:
                        gene_with_info+=GENE_NAME_DELIMITER+additional_gene_info[gene]
                    
                     # look for the overall gene abundance
                    abundance=gene_table[gene].get("all",0)
                    if abundance > 0:
                        output_line=TABLE_DELIMITER.join([BUG_DELIMITER.join([pathway,gene_with_info]),str(abundance)])
                        write_file_handle.write(output_line+"\n")
        else:
            for item in pathways_by_bug.get(pathway,[]):
                bug, data = item
                output_line=TABLE_DELIMITER.join([BUG_DELIMITER.join([pathway,bug]), data])
                write_file_handle.write(output_line+"\n")
                # check for gene families for this pathway for this bug
                for gene in pathways_to_genes.get(pathway,[]):
                    if gene in gene_table:
                        # check if additional information should be added to the gene
                        gene_with_info=gene
                        if gene in additional_gene_info:
                            gene_with_info+=GENE_NAME_DELIMITER+additional_gene_info[gene]
                    
                        abundance=gene_table[gene].get(bug,0)
                        if abundance > 0:
                            output_line=TABLE_DELIMITER.join([BUG_DELIMITER.join([pathway,bug,gene_with_info]),str(abundance)])
                            write_file_handle.write(output_line+"\n")
                            
    write_file_handle.close()
    
def read_mapping(gene_mapping_file,pathway_mapping_file):
    """
    Read the gene to reaction and then reaction to pathway mappings
    """
    
    # read the gene to reaction mappings
    reactions_to_genes={} 
    reactions_to_ecs={}
    try:
        if gene_mapping_file.endswith(".gz"):
            file_handle=gzip.open(gene_mapping_file, "rt")
            readlines = file_handle.readlines()
        elif gene_mapping_file.endswith(".bz2"):
            file_handle = bz2.BZ2File(gene_mapping_file, "r")
            readlines = [line.decode('utf-8') for line in file_handle.readlines()]
        else:
            file_handle=open(gene_mapping_file, "rt")
            readlines = file_handle.readlines()
    except EnvironmentError:
        sys.exit("ERROR: Unable to open the mapping file: " + gene_mapping_file)
            
    file_handle.close()
    
    # read through the mapping lines
    for line in readlines:
        data=line.rstrip().split(TABLE_DELIMITER)
        reaction=data.pop(0)
        ecs=data.pop(0).split(EC_DELIMITER)
        reactions_to_genes[reaction]=data
        reactions_to_ecs[reaction]=ecs
                
    # read the reaction to pathway mappings
    pathways_to_genes={}
    pathways_to_ecs={}
    try:
        file_handle=open(pathway_mapping_file, "rt")
    except EnvironmentError:
        sys.exit("ERROR: Unable to open the mapping file: " + pathway_mapping_file)
        
    for line in file_handle:
        info=line.rstrip().split(TABLE_DELIMITER)
        if len(info) == 2:
            pathway, structure = info
            # get the reactions from the structure
            genes=set()
            ecs=set()
            for item in structure.split(" "):
                genes.update(reactions_to_genes.get(item,[]))
                ecs.update(reactions_to_ecs.get(item,[]))
            pathways_to_genes[pathway]=genes
            pathways_to_ecs[pathway]=ecs
        
    file_handle.close()
  
    return pathways_to_genes,pathways_to_ecs
        
def read_gene_table(gene_table):
    """
    Read the input gene table
    Store the abundances related to each bug
    """
    
    gene_table_data={}
    gene_names={}
    
    lines=util.process_gene_table_with_header(gene_table, allow_for_missing_header=True)
    header=next(lines)
    
    for line in lines:
        # process if not a comment
        if not re.match("#",line):
            data=line.split(TABLE_DELIMITER)
            
            try:
                gene,bug=data[0].split(BUG_DELIMITER)
            except ValueError:
                gene=data[0]
                bug="all"
                
            try:
                # remove the gene name if present and store
                if GENE_NAME_DELIMITER in gene:
                    info=gene.split(GENE_NAME_DELIMITER)
                    gene = info[0]
                    gene_names[gene]=" ".join(info[1:])
                data_point=float(data[1])
            except (IndexError,ValueError):
                gene=""

            if gene and data_point > 0:
                # check for an EC set
                genes=[gene]
                if EC_DELIMITER in gene:
                    genes=gene.split(EC_DELIMITER)
                    
                for gene in genes:
                    if not gene in gene_table_data:
                        gene_table_data[gene]={}
                    gene_table_data[gene][bug]=gene_table_data[gene].get(bug,0)+data_point
        
    return gene_table_data, gene_names

def determine_mapping_type(gene_table,pathways_to_genes,pathways_to_ecs):
    """
    Determine if the input file is of gene families or EC abundance type
    """
    
    all_genes=set()
    all_ecs=set()
    for pathway,genes in pathways_to_genes.items():
        all_genes.update(genes)
        all_ecs.update(pathways_to_ecs[pathway])
    
    for gene in gene_table:
        if gene in all_genes:
            return pathways_to_genes
        elif gene in all_ecs:
            return pathways_to_ecs
    

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Unpack pathway abundances to show genes included\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--input-genes",
        help="the gene family or EC abundance file\n",
        required=True)
    parser.add_argument(
        "--input-pathways",
        help="the pathway abundance file\n",
        required=True)
    parser.add_argument(
        "--gene-mapping",
        help="gene family to reaction mapping file\n")
    parser.add_argument(
        "--pathway-mapping",
        help="reaction to pathway mapping file\n")
    parser.add_argument(
        "-r","--remove-taxonomy",
        action="store_true",
        help="remove the taxonomy from the output file\n")
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
    pathways_to_genes,pathways_to_ecs=read_mapping(args.gene_mapping,args.pathway_mapping)
    
    pathways_to_input=determine_mapping_type(gene_table, pathways_to_genes, pathways_to_ecs)
    
    print("Merging abundances.")
    merge_abundances(gene_table,pathways_to_input,input_pathways,output,gene_names,args.remove_taxonomy)
    

if __name__ == "__main__":
    main()
