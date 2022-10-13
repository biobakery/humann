"""
HUMAnN setup

To run: python setup.py install

"""

import sys

# required python versions (2.7+ or 3.0+)
required_python_version_major = [2,3]
required_python_version_minor = [7,0]
    
# check for either of the required versions
pass_check=False
try:
    for major, minor in zip(required_python_version_major, required_python_version_minor):
        if (sys.version_info[0] == major and sys.version_info[1] >= minor):
            pass_check=True
except (AttributeError,IndexError):
    sys.exit("CRITICAL ERROR: The python version found (version 1) " +
        "does not match the version required (version "+
        str(required_python_version_major)+"."+
        str(required_python_version_minor)+"+)")    

if not pass_check:
    sys.exit("CRITICAL ERROR: The python version found (version "+
        str(sys.version_info[0])+"."+str(sys.version_info[1])+") "+
        "does not match the version required (version "+
        str(required_python_version_major)+"."+
        str(required_python_version_minor)+"+)")


try:
    import setuptools
except ImportError:
    sys.exit("Please install setuptools.")
    
# check setuptools version    
required_setuptools_version_major = 1
try:
    setuptools_version = setuptools.__version__
    setuptools_version_major = int(setuptools_version.split(".")[0])
    if setuptools_version_major < required_setuptools_version_major:
        sys.exit("CRITICAL ERROR: The setuptools version found (version "+
                 setuptools_version+") does not match the version required "+
                 "(version "+str(required_setuptools_version_major) +"+)."
                 " Please upgrade your setuptools version.")
except (ValueError, IndexError, NameError):
    sys.exit("CRITICAL ERROR: Unable to call setuptools version. Please upgrade setuptools.")
    
from setuptools.command.install import install as _install
    
import distutils
import os

# try to import urllib.request.urlretrieve for python3
try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve

import tarfile
import subprocess
import shutil
import zipfile
import tempfile
import re
import time


VERSION = "3.6.1"

AUTHOR = "HUMAnN Development Team"
AUTHOR_EMAIL = "humann-users@googlegroups.com"
MAINTAINER = "Lauren McIver"
MAINTAINER_EMAIL = "lauren.j.mciver@gmail.com"

def byte_to_megabyte(byte):
    """
    Convert byte value to megabyte
    """
    
    return byte / (1024.0**2)

class ReportHook():
    def __init__(self):
        self.start_time=time.time()
        
    def report(self, blocknum, block_size, total_size):
        """
        Print download progress message
        """
        
        if blocknum == 0:
            self.start_time=time.time()
            if total_size > 0:
                print("Downloading file of size: " + "{:.2f}".format(byte_to_megabyte(total_size)) + " MB\n")
        else:
            total_downloaded=blocknum*block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))
                    
            if total_size > 0:
                percent_downloaded=total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stdout to overwrite stdout
                try:
                    download_rate=total_downloaded/(time.time()-self.start_time)
                    estimated_time=(total_size-total_downloaded)/download_rate
                except ZeroDivisionError:
                    download_rate=0
                    estimated_time=0
                estimated_minutes=int(estimated_time/60.0)
                estimated_seconds=estimated_time-estimated_minutes*60.0
                status +="{:3.2f}".format(percent_downloaded) + " %  " + \
                    "{:5.2f}".format(byte_to_megabyte(download_rate)) + " MB/sec " + \
                    "{:2.0f}".format(estimated_minutes) + " min " + \
                    "{:2.0f}".format(estimated_seconds) + " sec "
            status+="        \r"
            sys.stdout.write(status)

def download(url, download_file):
    """
    Download a file from a url
    """

    try:
        print("Downloading "+url)
        file, headers = urlretrieve(url,download_file,reporthook=ReportHook().report)
        # print final return to start new line of stdout
        print("\n")
    except EnvironmentError:
        print("WARNING: Unable to download "+url)
    

def download_unpack_tar(url,download_file_name,folder,software_name):
    """
    Download the url to the file and decompress into the folder
    """
    
    # Check for write permission to the target folder
    if not os.access(folder, os.W_OK):
        print("WARNING: The directory is not writeable: "+
            folder + " . Please modify the permissions.")
    
    download_file=os.path.join(folder, download_file_name)
    
    download(url, download_file)
    
    error_during_extract=False
    
    try:
        tarfile_handle=tarfile.open(download_file)
        tarfile_handle.extractall(path=folder)
        tarfile_handle.close()
    except (EnvironmentError,tarfile.ReadError):
        print("WARNING: Unable to extract "+software_name+".")
        error_during_extract=True
        
    if not error_during_extract:
        try:
            os.unlink(download_file)
        except EnvironmentError:
            print("WARNING: Unable to remove the temp download: " + download_file)
        
def download_unpack_zip(url,download_file_name,folder,software_name):
    """
    Download the url to the file and decompress into the folder
    """
    
    # Check for write permission to the target folder
    if not os.access(folder, os.W_OK):
        print("WARNING: The directory is not writeable: "+
            folder + " . Please modify the permissions.")
    
    download_file=os.path.join(folder, download_file_name)
    
    download(url, download_file)
    
    error_during_extract=False
    
    try:
        zipfile_handle=zipfile.ZipFile(download_file)
        zipfile_handle.extractall(path=folder)
        zipfile_handle.close()
    except EnvironmentError:
        print("WARNING: Unable to extract "+software_name+".")
        error_during_extract=True
        
    if not error_during_extract:
        try:
            os.unlink(download_file)
        except EnvironmentError:
            print("WARNING: Unable to remove the temp download: " + download_file)
        
def install_glpk(install_directory, replace_install=None):
    """
    Download and install the most recent glpk for minpath
    """
    
    glpk_url="http://ftp.gnu.org/gnu/glpk/glpk-4.55.tar.gz"
    glpk_folder="glpk-4.55"
    
    # Check if glpk is already installed
    glpk_installed=find_exe_in_path("glpsol")
    
    if not glpk_installed or replace_install:
        
        # get the prefix for the install directory
        install_directory=os.path.abspath(install_directory)
        prefix=os.path.abspath(os.path.join(install_directory,os.path.pardir))
        
        humann_source_folder=os.path.dirname(os.path.abspath(__file__))        
        tempfolder=tempfile.mkdtemp(prefix="glpk_download_",dir=humann_source_folder)
 
        # install the most recent gplk software
        glpk_download=os.path.join(tempfolder,glpk_url.split('/')[-1])
        print("Installing latest glpk.")
        download_unpack_tar(glpk_url,glpk_download,tempfolder,"glpk")
            
        # move to the glpk directory
        current_working_directory=os.getcwd()
        glpk_install_folder=os.path.join(tempfolder,glpk_folder)
        
        try:
            os.chdir(glpk_install_folder)
        except EnvironmentError:
            print("WARNING: glpk was not downloaded.")
        
        # test for gcc
        try:
            subprocess_output=subprocess.check_output(["gcc","--version"])
        except (EnvironmentError,subprocess.CalledProcessError):
            print("WARNING: Please install gcc.")
            
        # test for make
        try:
            subprocess_output=subprocess.check_output(["make","--version"])
        except (EnvironmentError,subprocess.CalledProcessError):
            print("WARNING: Please install make.")
            
        try:
            # run configure
            subprocess.call(["./configure","--prefix",prefix])
            # run make
            subprocess.call(["make"])
            # run make install
            subprocess.call(["make","install"])
        except (EnvironmentError,subprocess.CalledProcessError):
            print("WARNING: Errors installing new glpk version.")
                
        # return to original working directory
        os.chdir(current_working_directory)  
                
        # remove the temp download folder
        try:
            shutil.rmtree(tempfolder)
        except EnvironmentError:
            print("WARNING: Unable to remove temp install folder.")

    else:
        print("Found glpk install at "+glpk_installed)
        
        
def install_minpath(replace_install=None):
    """ 
    Download and install the minpath software if not already installed
    """
    
    # Download the minpath software v1.2
    # Check to see if already downloaded
    
    fullpath_scripts=os.path.join(os.path.dirname(os.path.abspath(__file__)),"humann","quantify")

    minpath_file="minpath1.2.tar.gz"
    minpath_url="http://omics.informatics.indiana.edu/mg/get.php?" + \
    "justdoit=yes&software=" + minpath_file
    minpath_install_folder="MinPath"
    
    # install minpath if not already installed
    minpath_exe=os.path.join(fullpath_scripts,minpath_install_folder,"MinPath1.2.py")
    if not os.path.isfile(minpath_exe) or replace_install:
        download_file=os.path.join(fullpath_scripts, minpath_file)
        print("Installing minpath.")
        download_unpack_tar(minpath_url,download_file,fullpath_scripts,"minpath")
    else:
        print("Found minpath install.")
        
def find_exe_in_path(exe):
    """
    Check that an executable exists in $PATH
    """
    
    paths = os.environ["PATH"].split(os.pathsep)
    for path in paths:
        fullexe = os.path.join(path,exe)
        if os.path.exists(fullexe):
            if os.access(fullexe,os.X_OK):
                return path
    return None

def install_boost(folder):
    """ Install boost locally in folder indicated
    """

    boost_file="boost_1_57_0.tar.gz"
    boost_url="http://sourceforge.net/projects/boost/files/boost/1.57.0/boost_1_57_0.tar.gz"

    print("Installing boost.")
    download_unpack_tar(boost_url, boost_file, folder, "boost")

    # get current working directory
    current_working_directory=os.getcwd()
    boost_build_dir=os.path.join(folder,"boost_1_57_0")

    try:
        os.chdir(boost_build_dir)
    except EnvironmentError:
        print("WARNING: boost directory does not exist")

    try:
        subprocess.call(["./bootstrap.sh","--with-libraries=timer,chrono,system,program_options,thread,iostreams","--prefix=../boost"])
        subprocess.call(["./b2","install"])
    except (EnvironmentError,subprocess.CalledProcessError):
        print("WARNING: Errors installing boost.")

    # return to original working directory
    os.chdir(current_working_directory)

        
def install_diamond(final_install_folder, build, replace_install=None):
    """ 
    Download and install the diamond software if not already installed
    """
    
    # Check if diamond is already installed
    diamond_installed=find_exe_in_path("diamond")
    
    if not diamond_installed or replace_install:
        diamond_version="2.0.15"
        diamond_exe="diamond"
        diamond_file="diamond-linux64.tar.gz"
        diamond_url="http://github.com/bbuchfink/diamond/releases/download/v{0}/diamond-linux64.tar.gz".format(diamond_version)
        
        # download source if build selected
        if build:
            diamond_file="v{0}.tar.gz".format(diamond_version)
            diamond_url="http://github.com/bbuchfink/diamond/archive/v{0}.tar.gz".format(diamond_version)

        humann_source_folder=os.path.dirname(os.path.abspath(__file__))        
        tempfolder=tempfile.mkdtemp(prefix="diamond_download_",dir=humann_source_folder)
        
        # install the diamond software
        print("Installing diamond.")
        error_during_install=False
        download_unpack_tar(diamond_url, diamond_file, tempfolder, diamond_exe)
        
        # compile diamond if build selected
        if build:
            # get the current directory
            current_working_directory=os.getcwd()
            diamond_build_dir=os.path.join(tempfolder,"diamond-{0}".format(diamond_version),"src")
            
            try:
                os.chdir(diamond_build_dir)
            except EnvironmentError:
                print("WARNING: diamond directory does not exist")
                
            # test for gcc
            try:
                subprocess_output=subprocess.check_output(["gcc","--version"])
            except (EnvironmentError,subprocess.CalledProcessError):
                print("WARNING: Please install gcc.")
                
            # test for make
            make_installed=True
            try:
                subprocess_output=subprocess.check_output(["make","--version"])
            except (EnvironmentError,subprocess.CalledProcessError):
                make_installed=False
                print("INFO: make is not installed.")
                
            # test for cmake
            cmake_installed=True
            try:
                subprocess_output=subprocess.check_output(["cmake","--version"])
            except (EnvironmentError,subprocess.CalledProcessError):
                cmake_installed=False
                print("INFO: cmake is not installed.")

            final_install_prefix=os.path.join(final_install_folder,os.pardir)
            if make_installed and cmake_installed:
                # if make and cmake are installed, run the standard build
                print("INFO: Installing diamond with cmake method.")
                diamond_install_folder=os.path.join(tempfolder,"diamond-{0}".format(diamond_version),"bin")
                try:
                    # make the install bin directory and change directories
                    os.mkdir(diamond_install_folder)
                    os.chdir(diamond_install_folder)
                    # run cmake and make
                    subprocess.call(["cmake","..","-DCMAKE_INSTALL_PREFIX="+final_install_prefix])
                    subprocess.call(["make","install"])
                except (EnvironmentError,subprocess.CalledProcessError):
                    print("WARNING: Errors installing diamond.")
                    error_during_install=True
            else:
                # if make or cmake are not installed, run the simple build
                print("INFO: Installing diamond with simple build method.")
                
                diamond_install_folder=os.path.join(tempfolder,"diamond-{0}".format(diamond_version))
                try:
                    # change directories to the diamond source folder
                    os.chdir(diamond_install_folder)
                
                    # remove the static flag from the build script
                    subprocess.call(["sed","-i","s/ -static//", "build_simple.sh"])
                
                    # run the build script
                    subprocess.call(["bash","-x", "build_simple.sh"])
                except (EnvironmentError,subprocess.CalledProcessError):
                    print("WARNING: Errors installing diamond.")
                    error_during_install=True
                # set the location of the diamond executable
                diamond_exe_full_path=os.path.join(diamond_install_folder, diamond_exe) 
               
                # copy the installed software to the final bin location
                try:
                    # copy to the install folder
                    shutil.copy(diamond_exe_full_path, final_install_folder)
                    # add executable permissions
                    os.chmod(os.path.join(final_install_folder,diamond_exe), 0o755)
                except (EnvironmentError, shutil.Error):
                    error_during_install=True 
                    
            # return to original working directory
            os.chdir(current_working_directory)
        else:
            # copy the installed software to the final bin location
            diamond_exe_full_path=os.path.join(tempfolder, diamond_exe)
            try:
                # copy to the install folder
                shutil.copy(diamond_exe_full_path, final_install_folder)
                # add executable permissions
                os.chmod(os.path.join(final_install_folder,diamond_exe), 0o755)
            except (EnvironmentError, shutil.Error) as e:
                print(e)
                error_during_install=True
            
        # remove the local diamond install
        try:
            shutil.rmtree(tempfolder)
        except EnvironmentError:
            print("WARNING: Unable to remove temp install folder.")
        
        if error_during_install:
            print("WARNING: Unable to install diamond. Please install diamond.")
        else:
            print("Installed diamond at "+final_install_folder)
    else:
        print("Found diamond install at "+diamond_installed)
        
        
def install_bowtie2(final_install_folder, mac_os, replace_install=None):
    """ 
    Download and install the bowtie2 software if not already installed
    """
    
    # Check if bowtie2 is already installed
    bowtie2_installed=find_exe_in_path("bowtie2")
    
    if not bowtie2_installed or replace_install:
        bowtie2_exe="bowtie2"
        bowtie2_file="bowtie2-2.2.3-linux-x86_64.zip"
        bowtie2_url="http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip"

        # if this is a MAC OS, select a different binary download
        if mac_os:
            bowtie2_file="bowtie2-2.2.3-macos-x86_64.zip"
            bowtie2_url="http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-macos-x86_64.zip"
            
        bowtie2_folder="bowtie2-2.2.3"
    
        humann_source_folder=os.path.dirname(os.path.abspath(__file__))
        tempfolder=tempfile.mkdtemp(prefix="bowtie2_download_",dir=humann_source_folder)

        # install the bowtie2 software
        print("Installing bowtie2.")
        error_during_install=False
        download_unpack_zip(bowtie2_url, bowtie2_file, tempfolder, bowtie2_exe)
        
        # copy the installed software to the final bin location
        # copy all bowtie2* executables
        fullpath_bowtie2_exe=os.path.join(tempfolder, bowtie2_folder)
        
        files=[]
        try:
            files=os.listdir(fullpath_bowtie2_exe)
        except EnvironmentError:
            print("WARNING: Bowtie2 files not found.")
            error_during_install=True
        
        for file in files:  
            # check if this file is one of the bowtie2* executables      
            if re.match(bowtie2_exe,file):  
                try:   
                    # copy to the install folder
                    shutil.copy(os.path.join(fullpath_bowtie2_exe,file), final_install_folder)
                    # add executable permissions
                    os.chmod(os.path.join(final_install_folder,file), 0o755)
                except (EnvironmentError, shutil.Error):
                    error_during_install=True
            
        # remove the local bowtie2 install
        try:
            shutil.rmtree(tempfolder)
        except EnvironmentError:
            print("WARNING: Unable to remove temp install folder.")
        
        if error_during_install:
            print("WARNING: Unable to install bowtie2. Please install bowtie2.")
        else:
            print("Installed bowtie2 at "+final_install_folder)
    else:
        print("Found bowtie2 install at "+bowtie2_installed)

        
class Install(_install):
    """
    Custom setuptools install command, set executable permissions for glpk
    """
    
    _install.user_options=_install.user_options+[('bypass-dependencies-install', 
        None, 'bypass install of dependencies'),('build-diamond',None,'build diamond'),
        ('replace-dependencies-install',None,'replace dependencies already installed')]
    
    def initialize_options(self):
        self.bypass_dependencies_install=False
        self.build_diamond=False
        self.replace_dependencies_install=False
        _install.initialize_options(self)
    
    def finalize_options(self):
        _install.finalize_options(self)
    
    def run(self):
        
        # install minpath if not already installed
        install_minpath(replace_install=self.replace_dependencies_install)
        
        _install.run(self)
        
        # find the current install folder
        current_install_folder=None
        for item in os.listdir(self.install_lib):
            full_path_item=os.path.join(self.install_lib, item)
            if os.path.isdir(full_path_item):
                if "humann" == item:
                    current_install_folder=full_path_item
        
        # find all glpsol executables
        if current_install_folder is None:
            print("Unable to find install folder at: " + self.install_lib)
        else:
            glpsols=[]
            minpath_folder=os.path.join(current_install_folder,"quantify","MinPath")
            for root, directories, files in os.walk(minpath_folder):
                for filename in files:
                    if filename == "glpsol":
                        glpsols.append(os.path.join(root,filename))
            # change the permissions of the glpk modules to make sure
            # they are executable
            for file in glpsols:
                try:
                    os.chmod(file,0o755)
                except EnvironmentError:
                    print("Unable to add execute permissions for file: " + file)
                    
        # find out the platform
        mac_os=False
        if sys.platform in ["darwin","os2","os2emx"]:
            mac_os=True
        
        # install dependencies if not already installed
        if not self.bypass_dependencies_install:
            # build diamond if set, or if on a Mac
            build_diamond=self.build_diamond
            if mac_os:
                build_diamond=True

            # get the final script install location for both source and wheel
            install_scripts=os.path.join(self.install_base, "bin")
            install_glpk(install_scripts,replace_install=self.replace_dependencies_install)
            install_diamond(install_scripts,build_diamond,replace_install=self.replace_dependencies_install)
            install_bowtie2(install_scripts,mac_os,replace_install=self.replace_dependencies_install)
        else:
            print("Bypass install of dependencies")
        
    
setuptools.setup(
    name="humann",
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    version=VERSION,
    license="MIT",
    description="HUMAnN: The HMP Unified Metabolic Analysis Network, version 3",
    long_description="HUMAnN v3 is a pipeline for efficiently and accurately determining " + \
        "the coverage and abundance of microbial pathways in a community " + \
        "from metagenomic data. Sequencing a metagenome typically produces millions " + \
        "of short DNA/RNA reads. This process, referred to as functional profiling, " + \
        "aims to describe the metabolic potential of a microbial community and its " + \
        "members. More generally, functional profiling answers the question: What " + \
        "are the microbes in my community-of-interest doing (or capable of doing)?",
    url="http://huttenhower.sph.harvard.edu/humann",
    keywords=['microbial','microbiome','bioinformatics','microbiology','metagenomic','metatranscriptomic','humann','humann'],
    platforms=['Linux','MacOS'],
    classifiers=[
        "Programming Language :: Python",
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.4",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
    packages=setuptools.find_packages(),
    cmdclass={'install': Install},
    package_data={
        'humann' : [
            'humann.cfg',
            'data/pathways/*',
            'data/misc/*',
            'data/uniref_DEMO/*',
            'data/chocophlan_DEMO/*',
            'tests/data/*.*',
            'tests/data/tooltest*/*',
            'quantify/MinPath/data/*',
            'quantify/MinPath/glpk-*/examples/glp*',
            'quantify/MinPath/glpk-*/examples/.libs/*'
        ]},
    entry_points={
        'console_scripts': [
            'humann = humann.humann:main',
            'humann3 = humann.humann:main',
            'humann_databases = humann.tools.humann_databases:main',
            'humann3_databases = humann.tools.humann_databases:main',
            'humann_config = humann.tools.humann_config:main',
            'humann_join_tables = humann.tools.join_tables:main',
            'humann_split_table = humann.tools.split_table:main',
            'humann_rename_table = humann.tools.rename_table:main',
            'humann_renorm_table = humann.tools.renorm_table:main',
            'humann_regroup_table = humann.tools.regroup_table:main',
            'humann_infer_taxonomy = humann.tools.infer_taxonomy:main',
            'humann_reduce_table = humann.tools.reduce_table:main',
            'humann_unpack_pathways = humann.tools.merge_abundance:main',
            'humann_test = humann.tests.humann_test:main',
            'humann_build_custom_database = humann.tools.build_custom_database:main',
            'humann_genefamilies_genus_level = humann.tools.genefamilies_genus_level:main',
            'humann_split_stratified_table = humann.tools.split_stratified_table:main',
            'humann_barplot = humann.tools.humann_barplot:main',
            'humann_benchmark = humann.tools.humann_benchmark:main',
            'humann_expand_cluster = humann.tools.expand_cluster:main'
        ]},
    test_suite= 'humann.tests.humann_test.unittests_suite_only',
    zip_safe = False
 )
