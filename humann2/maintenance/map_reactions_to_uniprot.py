
"""
This script will create a reactions.dat formatted file that includes the direct AC mappings for genes.
It will also include all of the ECs for each reaction (only of level4). It can optionally
write out GO mappings. It can also only print out direct to AC mappings.

To run: python map_reactions_to_uniprot.py --input-reactions reactions.dat --input-enzrxn enzrxns.dat 
--input-proteins proteins.dat --input-gene-links gene-links.dat --output reactions_merged_eclevel4_only.dat

The reactions*.dat output file of this script can be used as input to Reaction_to_Uniref5090.py .

"""

import sys
import re
import argparse
import os

COMMENT_LINE="#"
METACYC_ID="UNIQUE-ID"
METACYC_DATABASE_LINK_ID="DBLINKS"
METACYC_ID_DELIMITER=" - "
METACYC_EC_ID="EC-NUMBER"
METACYC_ENZ_ID="ENZYMATIC-REACTION"
METACYC_EC_TAG="EC-"
METACYC_ENZYME_ID="ENZYME"
METACYC_GENE_ID="GENE"
METACYC_DATABASE_LINK_ID="DBLINKS"
METACYC_GO_ID="GO-TERMS"
METACYC_SUBREACTION_ID="REACTION-LIST"

def read_metacyc_gene_links(file):
    """
    Process the metacyc gene_links.dat file
    """
    
    genes={}
    
    try:
        file_handle=open(file,"r")
        line=file_handle.readline()
    except EnvironmentError:
        sys.exit("Unable to read file: " + file)
        
    while line:
        if not re.match(COMMENT_LINE, line):
            data=line.rstrip().split("\t")
            if len(data) == 3:
                gene, uniprot, name = data
                if uniprot:
                    genes[gene]=uniprot
        
        line=file_handle.readline()
        
    file_handle.close()
    
    return genes

def read_metacyc_reactions(file):
    """
    Process the metacyc reactions.dat file
    """
    
    reactions_ec={}
    reactions_enzymes={}
    reactions_uniprot={}
    reactions_subreactions={}
    
    try:
        file_handle=open(file,"r")
        line=file_handle.readline()
    except EnvironmentError:
        sys.exit("Unable to read file: " + file)
        
    reaction=""
    ec=[]
    enzyme=[]
    uniprot=[]
    subreactions=[]
    while line:
        if not re.match(COMMENT_LINE, line):
            # store the reaction id
            if re.match(METACYC_ID, line):
                # store the latest set of reaction information
                if reaction:
                    if ec:
                        reactions_ec[reaction]=ec
                    if enzyme:
                        reactions_enzymes[reaction]=enzyme
                    if uniprot:
                        reactions_uniprot[reaction]=uniprot
                    if subreactions:
                        reactions_subreactions[reaction]=subreactions
                                          
                ec=[]
                enzyme=[]
                uniprot=[]
                subreactions=[]
                reaction=line.rstrip().split(METACYC_ID_DELIMITER)[-1].strip()
            
            elif re.match(METACYC_EC_ID,line):
                current_ec=line.rstrip().split(METACYC_ID_DELIMITER)[-1].strip()
                current_ec=re.sub(METACYC_EC_TAG,"",re.sub("\|","",current_ec))
                # only output reactions of at least 4 levels
                if current_ec.count(".")>2:
                    ec.append(current_ec)
            elif re.match(METACYC_ENZ_ID,line):
                current_enzyme=line.rstrip().split(METACYC_ID_DELIMITER)[-1].strip()
                enzyme.append(current_enzyme)
            elif re.match(METACYC_SUBREACTION_ID,line):
                current_subreaction=line.rstrip().split(METACYC_ID_DELIMITER)[-1].strip()
                subreactions.append(current_subreaction)
            # find the database ids
            elif re.match(METACYC_DATABASE_LINK_ID, line) and "UNIPROT" in line:
                # example database line
                # DBLINKS - (HMDB "HMDB01319" NIL |kothari| 3608602403 NIL NIL)
                data=re.sub('[\(\)]',"",line.rstrip().split(METACYC_ID_DELIMITER)[-1]).split(" ")
                database=re.sub("\|","",data[0])
                id=re.sub("\"","",data[1])
                uniprot.append(id)
                
        line=file_handle.readline()
        
    file_handle.close()
    
    # store the latest set of reaction information
    if reaction:
        if ec:
            reactions_ec[reaction]=ec
        if enzyme:
            reactions_enzymes[reaction]=enzyme
        if uniprot:
            reactions_uniprot[reaction]=uniprot
        if subreactions:
            reactions_subreactions[reaction]=subreactions
    
    return reactions_ec, reactions_enzymes, reactions_uniprot, reactions_subreactions

def read_metacyc_enzrxns(file):
    """
    Process the metacyc enzrxns.dat file
    """
    
    enzrxns_enzymes={}
    
    try:
        file_handle=open(file,"r")
        line=file_handle.readline()
    except EnvironmentError:
        sys.exit("Unable to read file: " + file)
        
    id=""
    enzyme=[]
    while line:
        if not re.match(COMMENT_LINE, line):
            # store the reaction id
            if re.match(METACYC_ID, line):
                # store the latest set of information
                if id and enzyme:
                    enzrxns_enzymes[id]=enzyme
                                          
                enzyme=[]
                id=line.rstrip().split(METACYC_ID_DELIMITER)[-1].strip()
            
            elif re.match(METACYC_ENZYME_ID,line):
                current_enzyme=line.rstrip().split(METACYC_ID_DELIMITER)[-1].strip()
                enzyme.append(current_enzyme)
                
        line=file_handle.readline()
        
        # store the latest set of information
        if id and enzyme:
            enzrxns_enzymes[id]=enzyme
        
    file_handle.close()
    
    return enzrxns_enzymes

def read_metacyc_proteins(file):
    """
    Process the metacyc proteins.dat file
    """
    
    proteins_uniprot={}
    proteins_genes={}
    proteins_go={}
    
    try:
        file_handle=open(file,"r")
        line=file_handle.readline()
    except EnvironmentError:
        sys.exit("Unable to read file: " + file)
        
    id=""
    uniprot=[]
    gene=[]
    go=[]
    while line:
        if not re.match(COMMENT_LINE, line):
            # store the reaction id
            if re.match(METACYC_ID, line):
                # store the latest set of reaction information
                if id:
                    if uniprot:
                        proteins_uniprot[id]=uniprot
                    if gene:
                        proteins_genes[id]=gene
                    if go:
                        proteins_go[id]=go
                uniprot=[]
                gene=[]
                go=[]
                id=line.rstrip().split(METACYC_ID_DELIMITER)[-1].strip()
            
            elif re.match(METACYC_GENE_ID,line):
                current_gene=line.rstrip().split(METACYC_ID_DELIMITER)[-1].strip()
                gene.append(re.sub(METACYC_EC_TAG,"",re.sub("\|","",current_gene)))
            # find the database ids
            elif re.match(METACYC_DATABASE_LINK_ID, line) and "UNIPROT" in line:
                # example database line
                # DBLINKS - (HMDB "HMDB01319" NIL |kothari| 3608602403 NIL NIL)
                data=re.sub('[\(\)]',"",line.rstrip().split(METACYC_ID_DELIMITER)[-1]).split(" ")
                database=re.sub("\|","",data[0])
                db_id=re.sub("\"","",data[1])
                uniprot.append(db_id)
            elif re.match(METACYC_GO_ID,line):
                current_go=line.rstrip().split(METACYC_ID_DELIMITER)[-1].strip()
                go.append(re.sub("\|","",current_go))
                
        line=file_handle.readline()
        
    file_handle.close()
    
    if id:
        if uniprot:
            proteins_uniprot[id]=uniprot
        if gene:
            proteins_genes[id]=gene
        if go:
            proteins_go[id]=go
    
    return proteins_genes, proteins_uniprot, proteins_go

def write_reactions_mapping(reactions_subreactions, reactions_ec, reactions_enzymes, 
        reactions_uniprot, enzrxn_enzyme, proteins_genes, proteins_uniprot, proteins_go, 
        genes_uniprot, output_file, uniprot_only, add_go_mappings):
    """
    Merge all of the datasets and write output file
    """
    
    reactions_go={}
    
    for reaction, enzrxns in reactions_enzymes.items():
        # check for subreactions
        subreactions=reactions_subreactions.get(reaction,[])
        all_enzrxns=set()
        all_enzrxns.update(enzrxns)
        new_uniprots=set()
        for rxn in subreactions:
            all_enzrxns.update(reactions_enzymes.get(rxn,[]))
            new_uniprots.update(reactions_uniprot.get(rxn,[]))
        
        for enzrxn in enzrxns:
            for enzyme in enzrxn_enzyme.get(enzrxn,[]):
                uniprot=proteins_uniprot.get(enzyme,set())
                if uniprot:
                    new_uniprots.update(uniprot)
                if enzyme in proteins_go:
                    reactions_go[reaction]=proteins_go[enzyme]
                for gene in proteins_genes.get(enzyme,[]):
                    uniprot=genes_uniprot.get(gene,"")
                    if uniprot:
                        new_uniprots.add(uniprot)

        new_uniprots.update(reactions_uniprot.get(reaction,[]))
        if new_uniprots:           
            reactions_uniprot[reaction]=new_uniprots
        
    try:
        file_handle=open(output_file,"w")
    except EnvironmentError:
        sys.exit("Unable to write file: " + output_file)
       
    total_reactions_with_mappings=set(list(reactions_uniprot))
    total_reactions_with_mappings.update(list(reactions_ec))
    if add_go_mappings:
        total_reactions_with_mappings.update(list(reactions_go))
       
    # Write out file in format of reactions.dat 
    if uniprot_only:
        total_reactions_with_mappings=list(reactions_uniprot)
    
    for reaction in total_reactions_with_mappings:
        file_handle.write(METACYC_ID+METACYC_ID_DELIMITER+reaction+"\n")
        
        # write out the uniprots
        for uniprot in reactions_uniprot.get(reaction,[]):
            # example format
            # DBLINKS - (UNIPROT "P0A0E0" RELATED-TO |paley| 3439827039 NIL NIL)
            file_handle.write(METACYC_DATABASE_LINK_ID+METACYC_ID_DELIMITER+"(UNIPROT \""+uniprot+"\" )\n")
        
        if not uniprot_only:
            for ec in reactions_ec.get(reaction,[]):
                file_handle.write(METACYC_EC_ID+METACYC_ID_DELIMITER+METACYC_EC_TAG+ec+"\n")
                
            if add_go_mappings:
                for go in reactions_go.get(reaction,[]):
                    file_handle.write("GO"+METACYC_ID_DELIMITER+go+"\n")
            
        file_handle.write("//\n")
    
    return total_reactions_with_mappings
        
def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Map reactions to uniprot\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--input-reactions", required=True)
    parser.add_argument("--input-enzrxn", required=True)
    parser.add_argument("--input-proteins", required=True)
    parser.add_argument("--input-gene-links", required=True)
    parser.add_argument("--output",required=True)
    parser.add_argument("--uniprot-only",action="store_true",default=False)
    parser.add_argument("--add-go-mappings",action="store_true",default=False)

    return parser.parse_args()

def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    output_file=os.path.abspath(args.output)
    
    print("Processing reactions.dat")
    reactions_ec, reactions_enzymes, reactions_uniprot, reactions_subreactions = read_metacyc_reactions(args.input_reactions)
    print("Found "+str(len(reactions_ec))+" reactions with ECs.")
    print("Found "+str(len(reactions_enzymes))+" reactions with enzymes.")
    print("Found "+str(len(reactions_uniprot))+"  reactions with UniProt.")
    print("Found "+str(len(reactions_subreactions))+" reactions with subreactions.")
    all_reactions_with_mappings=set(list(reactions_ec))
    all_reactions_with_mappings.update(list(reactions_enzymes))
    all_reactions_with_mappings.update(list(reactions_uniprot))
    print("Found "+str(len(all_reactions_with_mappings))+" reactions with mappings to EC or enzymes or UniProt.")
    
    print("Processing enzrxns.dat")
    enzrxn_enzyme=read_metacyc_enzrxns(args.input_enzrxn)
    print("Found "+str(len(enzrxn_enzyme))+ " enzyme reactions with mappings to enzymes.")
    
    print("Processing proteins.dat")
    proteins_genes, proteins_uniprot, proteins_go=read_metacyc_proteins(args.input_proteins)
    print("Found "+str(len(proteins_genes))+ " proteins with mappings to genes.")
    print("Found "+str(len(proteins_go))+ " proteins with mappings to go.")
    print("Found "+str(len(proteins_uniprot))+ " proteins with mappings to UniProt.")
    
    print("Processing gene-links.dat")
    genes_uniprot = read_metacyc_gene_links(args.input_gene_links)
    print("Found "+str(len(genes_uniprot))+ " genes with mappings to UniProt.")
    
    reactions_with_mappings=write_reactions_mapping(reactions_subreactions, reactions_ec, reactions_enzymes, 
        reactions_uniprot, enzrxn_enzyme, proteins_genes, proteins_uniprot, proteins_go,
        genes_uniprot, output_file, args.uniprot_only, args.add_go_mappings)
    
    if args.uniprot_only:
        print("Found "+str(len(reactions_with_mappings))+" reactions with mappings to UniProt.")
    else:
        print("Found "+str(len(reactions_with_mappings))+" reactions with mappings to EC or UniProt.")
    
if __name__ == "__main__":
    main()
