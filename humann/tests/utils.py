
import os
import shutil
import subprocess
import tempfile

import cfg

def create_temp_folder(suffix):
    """ Create a temp folder """
    
    return tempfile.mkdtemp(prefix="humann_test_",suffix="_"+suffix)

def remove_temp_file(file):
    """ Remove a temp file """
    
    try:
        os.remove(file)
    except EnvironmentError:
        print("Warning: Unable to remove temp file: " +file)
        
def remove_temp_folder(folder):
    """ Remove a temp folder """
    
    try:
        shutil.rmtree(folder)
    except EnvironmentError:
        print("Warning: Unable to remove temp folder: " +folder)

def return_exe_path(exe):
    """
    Return the location of the exe in $PATH
    """
    paths = os.environ["PATH"].split(os.pathsep)
    full_path=""
    for path in paths:
        fullexe = os.path.join(path,exe)
        if os.path.exists(fullexe):
            if os.access(fullexe,os.X_OK):
                full_path=path
    return full_path

def run_humann(command):
    """
    Run the humann command 
    Use the demo chocophlan and uniref databases """
    
    command+=["--nucleotide-database",cfg.chocophlan_example_demo_folder,
              "--protein-database", cfg.uniref_example_demo_folder]
    run_command(command)

def run_command(command):
    """ Run the command """
    
    print("\nTesting command: ")
    print(" ".join(command))
    
    try:
        subprocess.check_call(command)
    except (EnvironmentError, subprocess.CalledProcessError):
        raise EnvironmentError("Warning: Unable to execute command in test.\n"+" ".join(command))    
        
def remove_temp_folder(tempdir):
    """ Remove the temp folder """
    
    try:
        shutil.rmtree(tempdir)
    except EnvironmentError:
        print("Warning: Unable to remove temp directory: " + tempdir)
        
def check_output(output_files_expected,output_folder=None):
    """ Check the output folder has the expected file and they are all non-zero """
    
    for file in output_files_expected:
        if output_folder:
            expected_file = os.path.join(output_folder,file)
        else:
            expected_file = file
        # check the file exists
        yield (os.path.isfile(os.path.join(expected_file)), "File does not exist: " + file)
        
        # check the file is not empty
        yield (os.stat(expected_file).st_size > 0, "File is empty: " + file)
        
def read_table_rows(file):
    """ Read in the table from a file, storing by rows """
    
    # The first line of the file is the header of column ids
    data=[]
    rows=[]
    with open(file) as file_handle:
        columns=file_handle.readline().rstrip()
        for line in file_handle:
            line_info=line.rstrip().split("\t")
            rows.append(line_info[0])
            data.append(line_info[1:])
    
    # reduce rows to string so it will be the same format as columns
    rows="\t".join(rows)   
      
    return columns, rows, data

def files_almost_equal(file1, file2, precision=None):
    """ Check that the files have data that is almost equal to allow for rounding """
    
    if not precision:
        precision=7
        
    # read the files
    columns1, rows1, data1 = read_table_rows(file1)
    columns2, rows2, data2 = read_table_rows(file2)    
    
    # check the headers are the same
    if not columns1 == columns2:
        return False, "Column names in files differ"
    
    # check the row names are the same
    if not rows1 == rows2:
        return False, "Row names in files differ"
    
    # check the data is almost the same to allow for rounding
    for row1_data, row2_data in zip(data1, data2):
        for data1_point, data2_point in zip(row1_data, row2_data):
            if not round(float(data1_point),precision) == round(float(data2_point),precision):
                return False, "Data points in files differ"
    
    return True, "NA"
        
def file_basename(file):
    """ Get the basename for the file without the extension """
    
    return os.path.splitext(os.path.basename(file))[0]

def read_biom_table( path ):
    """
    return the lines in the biom file
    """

    try:
        import biom
    except ImportError:
        raise ImportError("ERROR: Could not find the biom software."+
            " This software is required since the input file is a biom file.")
        
    try:
        tsv_table = biom.load_table( path ).to_tsv().split("\n")
    except (EnvironmentError, TypeError):
        raise EnvironmentError("ERROR: Unable to read biom input file.")
        
    return tsv_table
