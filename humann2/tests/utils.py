
import os
import shutil

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
