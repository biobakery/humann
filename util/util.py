"""
Utility functions
"""

import os
import sys
import re
import subprocess

def find_exe_in_path(exe):
    """
    Check that an executable exists in $PATH
    """
    
    paths = os.environ["PATH"].split(os.pathsep)
    for path in paths:
        fullexe = os.path.join(path,exe)
        if os.path.exists(fullexe):
            if os.access(fullexe,os.X_OK):
                return True
    return False

def biom_to_tsv(biom_file, new_tsv_file, taxonomy=None):
    """
    Convert from a biom to tsv file
    """

    cmd=["biom","convert","-i",biom_file,"-o",new_tsv_file,"--to-tsv"]

    # try to convert
    if taxonomy:
        cmd+=["--header-key","taxonomy"]
    
    try:
        if os.path.isfile(new_tsv_file):
            os.unlink(new_tsv_file)
        p_out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except (EnvironmentError, subprocess.CalledProcessError):
        sys.exit("Unable to convert biom file to tsv")
        
def tsv_to_biom(tsv_file, biom_file):
    """
    Convert from a biom to tsv file
    """
    
    cmd=["biom","convert","-i",tsv_file,"-o",biom_file,"--table-type","Gene table","--to-hdf5"]

    try:
        if os.path.isfile(biom_file):
            os.unlink(biom_file)
        p_out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except (EnvironmentError, subprocess.CalledProcessError):
        sys.exit("Unable to convert tsv file to biom")
        
def process_gene_table_header(gene_table):
    """
    Process through the header portion of the gene table file
    """
    
    # indicator of a comment line or the header
    # the last line in the file with this indicator is the header
    GENE_TABLE_COMMENT_LINE="^#"
    
    # try to open the file
    try:
        file_handle=open(gene_table,"r")
        line=file_handle.readline()
    except EnvironmentError:
        sys.exit("Unable to read file: " + gene_table)
            
    # find the headers
    header_flag=False
    header=line
    while re.match(GENE_TABLE_COMMENT_LINE, line):
        header_flag=True
        header = line
        line = file_handle.readline()
            
    if not header_flag:
        sys.exit("File does not have a required header: " + gene_table +
            " . Please add a header which includes the indicator: " +
            GENE_TABLE_COMMENT_LINE)
        
    return file_handle, header, line
    