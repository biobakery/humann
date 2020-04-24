#!/usr/bin/env python

"""
HUMAnN: humann_config module
Configuration settings print and update

Dependencies: None

To Run: humann_config --update <section> <name> <value>

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

from .. import config

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "HUMAnN Configuration\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--print", 
        dest="print_config",
        action="store_true",
        help="print the configuration\n")
    parser.add_argument(
        "--update", 
        nargs=3,
        metavar=("<section>","<name>","<value>"),
        help="update the section : name to the value provided\n")
    
    return parser.parse_args()

def main():
    # Parse arguments from the command line
    args=parse_arguments(sys.argv)
    
    if args.update:
        # update the config file
        config.update_user_edit_config_file_single_item(args.update[0],args.update[1],args.update[2])
    
    if args.print_config or not args.update:
        # print the current configuration
        current_config_items=config.read_user_edit_config_file()
        print("HUMAnN Configuration ( Section : Name = Value )")
        for section in current_config_items:
            for name,value in current_config_items[section].items():
                print(section+" : "+name+" = "+str(value))
                

