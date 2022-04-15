#!/usr/bin/env python

"""
HUMAnN: humann_databases module
Download databases an update config settings

Dependencies: None

To Run: humann_databases --download <database> <build> <install_location>

Copyright (c) 2014 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import sys
    
# Try to load one of the humann src modules to check the installation
try:
    from .. import check
except ImportError:
    sys.exit("CRITICAL ERROR: Unable to find the HUMAnN python package." +
        " Please check your install.") 

# Check the python version
check.python_version()

import argparse
import os

from .. import config
from .. import utilities

# the locations of the current databases to download
current_downloads={
    "chocophlan" : 
        {
            "full" : "http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v201901_v31.tar.gz",
            "DEMO" : "http://huttenhower.sph.harvard.edu/humann_data/chocophlan/DEMO_chocophlan.v201901_v31.tar.gz"
        },
    "uniref" : 
        {
            "uniref50_diamond" : "http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref50_annotated_v201901b_full.tar.gz",
            "uniref90_diamond" : "http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_annotated_v201901b_full.tar.gz",
            "uniref50_ec_filtered_diamond" : "http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_ec_filtered/uniref50_ec_filtered_201901b_subset.tar.gz",
            "uniref90_ec_filtered_diamond" : "http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_ec_filtered/uniref90_ec_filtered_201901b_subset.tar.gz",
            "DEMO_diamond" : "http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_DEMO_diamond_v201901b.tar.gz"
        },
    "utility_mapping" :
        {
            "full" : "http://huttenhower.sph.harvard.edu/humann_data/full_mapping_v201901b.tar.gz"
         }
}

database_type={
    "chocophlan" : "nucleotide",
    "uniref" : "protein",
    "utility_mapping" : "utility_mapping"
}

def check_user_database(original, user):
    """ Check that the user database is of the expected version """

    if original.split("/")[-1] == user.split("/")[-1]:
        return True
    else:
        sys.exit("The user database selected does not match that expected: "+original.split("/")[-1])

def download_database(database, build, location, database_location):
    """
    Download and decompress the selected database
    """
    
    install_location=""
    if database in current_downloads:
        if build in current_downloads[database]:
            # create a subfolder to hold the contents of the database
            install_location=os.path.join(location,database)
            if not os.path.isdir(install_location):
                try:
                    print("Creating subdirectory to install database: " + install_location)
                    os.mkdir(install_location)
                except EnvironmentError:
                    sys.exit("CRITICAL ERROR: Unable to create directory: " + install_location)

            # download the database
            downloaded_file=os.path.join(location,current_downloads[database][build].split('/')[-1])

            if database_location:
                check_user_database(current_downloads[database][build],database_location)
                utilities.download_tar_and_extract_with_progress_messages(database_location,
                    downloaded_file, install_location)
            else:
                utilities.download_tar_and_extract_with_progress_messages(current_downloads[database][build], 
                    downloaded_file, install_location)

            # remove the download (if not a local file provided by the user)
            if not database_location or (database_location and not os.path.isfile(database_location)):
                try:
                    os.unlink(downloaded_file)
                except EnvironmentError:
                    print("Unable to remove file: " + downloaded_file)
            
            print("\nDatabase installed: " + install_location + "\n")
        else:
            sys.exit("ERROR: Please select an available build.")
    else:
        sys.exit("ERROR: Please select an available database.")
        
    return install_location

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "HUMAnN Databases\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--available", 
        action="store_true",
        help="print the available databases\n")
    parser.add_argument(
        "--download", 
        nargs=3,
        metavar=("<database>","<build>","<install_location>"),
        help="download the selected database to the install location\n")
    parser.add_argument(
        "--update-config", 
        default="yes",
        choices=["yes","no"],
        help="update the config file to set the new database as the default [DEFAULT: yes]\n")
    parser.add_argument(
        "--database-location",
        help="location (local or remote) to pull the database")

    return parser.parse_args()

def main():
    # Parse arguments from the command line
    args=parse_arguments(sys.argv)
    
    if args.download:
        # download the database
        database=args.download[0]
        build=args.download[1]
        location=os.path.abspath(args.download[2])
        
        # create the install location if it does not already exist
        if not os.path.isdir(location):
            try:
                print("Creating directory to install database: " + location)
                os.mkdir(location)
            except EnvironmentError:
                sys.exit("CRITICAL ERROR: Unable to create directory: " + location)
        
        install_location=download_database(database,build,location,args.database_location)
        
        if args.update_config == "yes":
            # update the config file with the installed location
            config.update_user_edit_config_file_single_item("database_folders",
                database_type[database],install_location)
    
    if args.available or not args.download:
        # print the available databases
        current_config_items=config.read_user_edit_config_file()
        print("HUMAnN Databases ( database : build = location )")
        for database in current_downloads:
            for build, location in current_downloads[database].items():
                print(database+" : "+build+" = "+location)
                

