"""
HUMAnN2 setup

To run: python setup.py install

"""

import sys

try:
    import setuptools
except ImportError:
    sys.exit("Please install setuptools.")
    
    
import distutils
import os
import urllib
import tarfile
import subprocess
import shutil

VERSION = "0.1.4"

def install_tar(url,download_file,folder):
    """
    Download the url to the file and decompress into the folder
    """
    
    # Check for write permission to the target folder
    if not os.access(folder, os.W_OK):
        sys.exit("ERROR: The directory set to install is not writeable: "+
            folder + " . Please modify the permissions.")
    
    try:
        file, headers = urllib.urlretrieve(url,download_file)
        tarfile_handle=tarfile.open(download_file)
        tarfile_handle.extractall(path=folder)
    except EnvironmentError:
        sys.exit("ERROR: Unable to install minpath.")
        
    try:
        os.unlink(download_file)
    except EnvironmentError:
        print("warning: unable to remove the temp minpath download: " + downloaded_file)
        
def install_glpk(minpath_folder):
    """
    Download and install the most recent glpk for minpath
    """
    
    glpk_url="http://ftp.gnu.org/gnu/glpk/glpk-4.55.tar.gz"
    glpk_folder="glpk-4.55"
    
 
    # install the most recent gplk software
    glpk_download=os.path.join(minpath_folder,glpk_url.split('/')[-1])
    print("installing latest glpk")
    install_tar(glpk_url,glpk_download,minpath_folder)
        
    # move to the glpk directory
    current_working_directory=os.getcwd()
    glpk_install_folder=os.path.join(minpath_folder,glpk_folder)
    os.chdir(glpk_install_folder)
        
    glpk_install_error=False
    try:
        # run configure
        subprocess.call(["./configure"])
        # run make
        subprocess.call(["make"])
        
    except (EnvironmentError,subprocess.CalledProcessError):
        glpk_install_error=True
            
    # return to original working directory
    os.chdir(current_working_directory)
        
    # remove the latest install if needed
    if glpk_install_error:
        shutil.rmtree(glpk_install_folder,ignore_errors=True)    
        
def install_minpath(replace_install=None,update_glpk=None):
    """ 
    Download and install the minpath software if not already installed
    """
    
    # Download the minpath software v1.2
    # Check to see if already downloaded
    
    fullpath_scripts=os.path.join(os.path.dirname(os.path.abspath(__file__)),"humann2","quantify")

    minpath_file="minpath1.2.tar.gz"
    minpath_url="http://omics.informatics.indiana.edu/mg/get.php?" + \
    "justdoit=yes&software=" + minpath_file
    minpath_install_folder="MinPath"
    
    # install minpath if not already installed
    minpath_exe=os.path.join(fullpath_scripts,minpath_install_folder,"MinPath1.2.py")
    if not os.path.isfile(minpath_exe) or replace_install:
        download_file=os.path.join(fullpath_scripts, minpath_file)
        print("installing minpath")
        install_tar(minpath_url,download_file,fullpath_scripts)
        
        if update_glpk:
            minpath_folder=os.path.join(fullpath_scripts,minpath_install_folder)
            install_glpk(minpath_folder)
        
class InstallMinpath(distutils.cmd.Command):
    """
    Custom distutils command to install minpath
    """
    
    description = "install minpath"
    
    user_options = [('update-glpk', None, 'option to update glpk')]
    
    def initialize_options(self):
        self.update_glpk=False
        
    def finalize_options(self):
        pass
    
    def run(self):
        install_minpath(replace_install=True,update_glpk=self.update_glpk)
    
setuptools.setup(
    name="humann2",
    version=VERSION,
    license="MIT",
    description="HUMAnN2 is a pipeline for efficiently and accurately determining " + \
        "the coverage and abundance of microbial pathways in a community " + \
        "from metagenomic data. Sequencing a metagenome typically produces millions " + \
        "of short DNA/RNA reads.",
    maintainer="Lauren McIver",
    maintainer_email="lauren.j.mciver@gmail.com",
    url="http://huttenhower.sph.harvard.edu/humann2",
    keywords=["microbial","pathways","metabolic","analysis","metagenomic","metatranscriptomic"],
    classifiers=[
        "Programming Language :: Python",
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "License :: MIT License",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
    long_description=open('readme.md').read(),
    packages=setuptools.find_packages(),
    cmdclass={'minpath': InstallMinpath},
    package_data={
        'humann2' : [
            'humann2.cfg',
            'data/pathways/*',
            'data/misc/*',
            'data/uniref_DEMO/*',
            'data/chocophlan_DEMO/*',
            'tests/data/*',
            'quantify/MinPath/data/*',
            'quantify/MinPath/glpk-*/examples/glp*',
            'quantify/MinPath/glpk-*/examples/.libs/*'
        ]},
    entry_points={
        'console_scripts': [
            'humann2 = humann2.humann2:main',
            'humann2_databases = humann2.tools.humann2_databases:main',
            'humann2_config = humann2.tools.humann2_config:main',
            'humann2_join_tables = humann2.tools.join_tables:main',
            'humann2_split_table = humann2.tools.split_table:main',
            'humann2_rename_table = humann2.tools.rename_table:main',
            'humann2_renorm_table = humann2.tools.renorm_table:main',
            'humann2_regroup_table = humann2.tools.regroup_table:main',
            'humann2_humann1_kegg = humann2.tools.humann1_kegg:main',
            'humann2_rna_dna_norm = humann2.tools.rna_dna_norm',
        ]},
    test_suite= 'humann2.tests.humann2_test.main',
    zip_safe = False
 )
