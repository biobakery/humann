#!/usr/bin/python
"""
Create a mapping file of metacyc transporters and uniref50 ids.

To Run: 
$ python make_map_transporter_uniref50.py --input-transporters transporters.col --input-proteins proteins.dat 
--input-id-mapping map_uniprot_UniRef50.dat.gz --output map_uniref50_transporter.tar.gz

"""

import sys

try:
    import argparse
except ImportError:
    sys.exit("ERROR: Please upgrade to at least python v2.7")
    
import re
import os
import gzip

COLUMN_DELIMITER="\t"
ITEM_DELIMITER=","
OPTIONAL_TAG="-"
COMMENT_LINE="#"
METACYC_ID="UNIQUE-ID"
METACYC_ID_DELIMITER=" - "
METACYC_NAME_ID="COMMON-NAME"
METACYC_COMPONENT_ID="COMPONENT-OF"
METACYC_DATABASE_LINK_ID="DBLINKS"

def write_mapping_file(transporters,uniprot_to_transporters,uniref_to_uniprot,output_file):
    """
    Write the mapping file of uniref50 to transporters
    """

    try:
        file_handle=open(output_file,"w")
    except EnvironmentError:
        sys.exit("Unable to read file: " + file)
     
    transporter_to_uniref={}

    for uniref in uniref_to_uniprot:
        uniprots=uniref_to_uniprot[uniref]
        for uniprot in uniprots:
            for transporter in uniprot_to_transporters.get(uniprot,[]):
                name=transporters.get(transporter,None)
                if not name is None:
                    # if there is a name that is not "", add to the transporter id
                    if name:
                        transporter+=": " + name
                    if not transporter in transporter_to_uniref:
                        transporter_to_uniref[transporter]=set()
                    transporter_to_uniref[transporter].add(uniref)

    for transporter in transporter_to_uniref:
        file_handle.write(COLUMN_DELIMITER.join([transporter]+list(transporter_to_uniref[transporter]))+"\n")

    file_handle.close()

def read_uniprot(verbose,file,uniprot_to_transporters):
    """
    Read the uniprot AC to UniRef50 mapping
    """

    uniref_to_uniprot={}

    try:
        file_handle=gzip.open(file)
        line=file_handle.readline()
    except EnvironmentError:
        sys.exit("Unable to read file: " + file)

    count=0
    while line:
        data=line.rstrip().split(COLUMN_DELIMITER)
        count+=1
        if len(data) == 2:
            uniprot, uniref = data
            if uniprot in uniprot_to_transporters:
                if not uniref in uniref_to_uniprot:
                    uniref_to_uniprot[uniref]=set()
                uniref_to_uniprot[uniref].add(uniprot)

        if verbose and count % 100000 == 0:
            print("Read "+str(count)+" lines")
        
        line=file_handle.readline()

    file_handle.close()

    return uniref_to_uniprot

def read_proteins(file):
    """
    Read the MetaCyc proteins.dat file
    """

    uniprot_to_transporters={}
    
    try:
        file_handle=open(file,"r")
        line=file_handle.readline()
    except EnvironmentError:
        sys.exit("Unable to read file: " + file)
        
    id=""
    dbids=set()
    transporters=set()
    while line:
        if not re.match(COMMENT_LINE, line):
            # find the id
            if re.match(METACYC_ID, line):
                if id:
                    # record the last transporter
                    for uniprot_id in dbids:
                        if not uniprot_id in uniprot_to_transporters:
                            uniprot_to_transporters[uniprot_id]=set()
                        uniprot_to_transporters[uniprot_id].update(transporters)
                transporters=set()
                dbids=set()
                id=compound=line.rstrip().split(METACYC_ID_DELIMITER)[-1]
            # find the transporters
            elif re.match(METACYC_COMPONENT_ID, line):
                transporters.add(line.rstrip().split(METACYC_ID_DELIMITER)[-1])
            elif re.match(METACYC_DATABASE_LINK_ID, line) and "UNIPROT" in line:
                # example database line
                # DBLINKS - (HMDB "HMDB01319" NIL |kothari| 3608602403 NIL NIL)
                data=re.sub('[\(\)]',"",line.rstrip().split(METACYC_ID_DELIMITER)[-1]).split(" ")
                database=re.sub("\|","",data[0])
                id=re.sub("\"","",data[1])
                dbids.add(id)
            
        line=file_handle.readline()
        
    # store the last id
    if id:
        for uniprot_id in dbids:
            if not uniprot_id in uniprot_to_transporters:
                uniprot_to_transporters[uniprotid_id]=set()
            uniprot_to_transporters[uniprot_id].update(transporters)
        
    file_handle.close()
    
    return uniprot_to_transporters

def read_transporters(file):
    """
    Read the MetaCyc transporters.col file
    """

    transporters={}
    
    try:
        file_handle=open(file,"r")
        line=file_handle.readline()
    except EnvironmentError:
        sys.exit("Unable to read file: " + file)
        
    name=""
    id=""
    while line:
        if not re.match(COMMENT_LINE, line):
            # find the id and name
            data=line.rstrip().split(COLUMN_DELIMITER)
            if len(data)>2:
                transporters[data[0]]=re.sub("</*[a-z]*[A-Z]*>","",data[1])
            
        line=file_handle.readline()
        
    file_handle.close()
    
    return transporters

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Create mapping of MetaCyc transporters to UniRef50\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-v","--verbose",
        help="additional output is printed\n",
        action="store_true",
        default=False)
    parser.add_argument(
        "--input-transporters",
        help="the MetaCyc transporters.col file\n",
        required=True)
    parser.add_argument(
        "--input-proteins",
        help="the MetaCyc proteins.dat file\n",
        required=True)
    parser.add_argument(
        "--input-id-mapping",
        help="the UniProt file mapping AC to UniRef50\n",
        required=True)
    parser.add_argument(
        "-o","--output",
        help="the file to write the mapping\n",
        required=True)

    return parser.parse_args()


def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
     
    input_transporters=os.path.abspath(args.input_transporters)
    input_proteins=os.path.abspath(args.input_proteins)
    input_id_mapping=os.path.abspath(args.input_id_mapping)
    output_file=os.path.abspath(args.output)
    
    if args.verbose:
        print("Reading transporters file")
    transporters=read_transporters(input_transporters)    
    if args.verbose:
        print("Total transporters: "+str(len(transporters)))

    if args.verbose:
        print("Reading proteins file")
    uniprot_to_transporters=read_proteins(input_proteins)
    if args.verbose:
        print("Total uniprots mapping to transporters: " + str(len(uniprot_to_transporters)))   
 
    if args.verbose:
        print("Reading id mapping file")
    uniref_to_uniprot=read_uniprot(args.verbose,input_id_mapping,uniprot_to_transporters)        

    if args.verbose:
        print("Writing mapping file")
    write_mapping_file(transporters,uniprot_to_transporters,uniref_to_uniprot,output_file)
        
if __name__ == "__main__":
    main()
        
    
