
import os
import shutil
import subprocess
import tempfile

def create_temp_folder(suffix):
    """ Create a temp folder """
    
    return tempfile.mkdtemp(prefix="humann2_test_",suffix="_"+suffix)

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

def run_command(command):
    """ Run the humann2 command """
    
    print("\nTesting command: ")
    print(" ".join(command))
    
    try:
        subprocess.check_call(command)
    except (EnvironmentError, subprocess.CalledProcessError):
        print("Warning: Unable to execute humann2 in test.\n"+" ".join(command))    
        
def remove_temp_folder(tempdir):
    """ Remove the temp folder """
    
    try:
        shutil.rmtree(tempdir)
    except EnvironmentError:
        print("Warning: Unable to remove temp directory: " + tempdir)
        
def check_output(output_files_expected,output_folder):
    """ Check the output folder has the expected file and they are all non-zero """
    
    for file in output_files_expected:
        expected_file = os.path.join(output_folder,file)
        # check the file exists
        yield (os.path.isfile(os.path.join(expected_file)), "File does not exist: " + file)
        
        # check the file is not empty
        yield (os.stat(expected_file).st_size > 0, "File is empty: " + file)
        
def file_basename(file):
    """ Get the basename for the file without the extension """
    
    return os.path.splitext(os.path.basename(file))[0]

