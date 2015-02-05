#!/usr/bin/env python

"""
Join a set of gene tables into a single table

This module will join gene tables output by HUMAnN2. 

Dependencies: Biom (only required if running with .biom files)

To Run: 
$ ./join_gene_tables.py -i <input_dir> -o <gene_table.{tsv,biom}>

"""

import argparse
import sys
import tempfile
import os
import shutil
import re

import util

GENE_TABLE_DELIMITER="\t"
GENE_TABLE_COMMENT_LINE="^#"
BIOM_FILE_EXTENSION=".biom"
TSV_FILE_EXTENSION=".tsv"
        
def join_gene_tables(gene_tables,output):
    """
    Join the gene tables to a single gene table
    """
    
    gene_table_data={}
    start_column_id=""
    genes={}
    for gene_table in gene_tables:
        
        file_handle,header,line=util.process_gene_table_header(gene_table)
    
        sample_names=header.rstrip().split(GENE_TABLE_DELIMITER)
        start_column_id=sample_names.pop(0)
        
        sample_data={}
        while line:
            data=line.rstrip().split(GENE_TABLE_DELIMITER)
            gene=data.pop(0)
            genes[gene]=1
            for i, data_point in enumerate(data):
                if sample_names[i] in sample_data:
                    sample_data[sample_names[i]][gene]=data_point
                else:
                    sample_data[sample_names[i]]={gene: data_point}
            line=file_handle.readline()
            
        file_handle.close()
    
        # check this sample id is unique
        for sample in sample_data:
            if sample in gene_table_data:
                sys.exit("Duplicate sample name: " + sample + ". Please remove " +
                    "duplicate sample names.")
            else:
                gene_table_data[sample]=sample_data[sample]
                
                
    # write the joined gene table
    sorted_sample_list=sorted(list(gene_table_data))
    sample_header=[start_column_id]+sorted_sample_list
    sorted_gene_list=sorted(list(genes))
    try:
        file_handle=open(output,"w")
        file_handle.write(GENE_TABLE_DELIMITER.join(sample_header)+"\n")
    except EnvironmentError:
        sys.exit("Unable to write file: " + file)  
        
    for gene in sorted_gene_list:
        write_data=[gene]
        for sample in sorted_sample_list:
            write_data.append(gene_table_data[sample].get(gene,"0"))
        file_handle.write(GENE_TABLE_DELIMITER.join(write_data)+"\n")
    
    file_handle.close()

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Join gene tables from output of HUMAnN2\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-v","--verbose", 
        help="additional output is printed\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "-i","--input",
        help="the directory of gene tables\n",
        required=True)
    parser.add_argument(
        "-o","--output",
        help="the gene table to write\n",
        required=True)
    parser.add_argument(
        "--file_name",
        help="only join gene tables with this string included in the file name")

    return parser.parse_args()


def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    # check for format of the gene tables
    input_dir=os.path.abspath(args.input)
    
    # check the directory exists
    if not os.path.isdir(input_dir):
        sys.exit("The input directory provided can not be found." + 
            "  Please enter a new directory.")
    
    biom_flag=False
    for file in os.listdir(input_dir):
        if file.endswith(BIOM_FILE_EXTENSION):
            biom_flag=True
            
    # Check for the biom software if running with a biom input file
    if biom_flag:
        if not util.find_exe_in_path("biom"):
            sys.exit("Could not find the location of the biom software."+
            " This software is required since the input file is a biom file.")       
    
    args.output=os.path.abspath(args.output)
    output_dir=os.path.dirname(args.output)
    
    # Create a temp folder for the biom conversions
    if biom_flag:
        temp_dir=tempfile.mkdtemp( 
            prefix='humann2_split_gene_tables_', dir=output_dir)
        if args.verbose:
            print("Temp folder created: " + temp_dir)
        
    gene_tables=[]
    file_list=os.listdir(input_dir)
    
    # filter out files which do not meet the name requirement if set
    if args.file_name:
        reduced_file_list=[]
        for file in file_list:
            if re.search(args.file_name,file):
                reduced_file_list.append(file)
        file_list=reduced_file_list
    
    for file in file_list:
        if file.endswith(BIOM_FILE_EXTENSION):
            # create a new temp file
            file_out, new_file=tempfile.mkstemp(dir=temp_dir)
            os.close(file_out)
        
            # convert biom file to tsv
            if args.verbose:
                print("Processing file: " + os.path.join(input_dir,file))
            util.biom_to_tsv(os.path.join(input_dir,file),new_file)
            gene_tables.append(new_file)
        else:
            gene_tables.append(os.path.join(input_dir,file))
        
    # split the gene table
    if args.verbose:
        print("Joining gene table")
        
    if biom_flag:
        # create a new temp file
        file_out, new_file=tempfile.mkstemp(dir=temp_dir)
        os.close(file_out)
        join_gene_tables(gene_tables,new_file)
        util.tsv_to_biom(new_file, args.output)
    else:
        join_gene_tables(gene_tables,args.output)
            
    # deleting temp folder with all files
    if biom_flag:
        if args.verbose:
            print("Deleting temp files in temp folder: " + temp_dir)
        shutil.rmtree(temp_dir)
    
    print("Gene table created: " + args.output)

if __name__ == "__main__":
    main()
