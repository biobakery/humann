"""
Stores the alignments identified
Stores the reactions from the database selected
"""

import config

class alignments:
    def __init__(self):
        self.hits=[]
        self.bugs={}
        self.genes={}
        self.hit_index=["bug","reference","query","evalue","identity","coverage"]
        
    def add(self, reference, query, evalue, identity, coverage, bug): 
        """ 
        Add the hit to the list
        Add the index of the hit to the bugs list and gene list
        """
        self.hits.append([bug, reference, query, evalue, identity, coverage]) 
        
        index=len(self.hits)-1
        if bug in self.bugs:
            self.bugs[bug].append(index)
        else:
            self.bugs[bug]=[index]
            
        if reference in self.genes:
            self.genes[reference].append(index)
        else:
            self.genes[reference]=[index]
            
    def count_bugs(self):
        """ 
        Return total number of bugs
        """
        return len(self.bugs.keys())
    
    def count_genes(self):
        """ 
        Return total number of genes
        """
        return len(self.genes.keys())    
            
    def print_bugs(self):
        """
        Print out the bugs and total number of hits
        """
        for bug in self.bugs.keys():
            print bug + ": " + str(len(self.bugs[bug])) + " hits"
            
    def print_genes(self):
        """
        Print out the genes and the total number of hits
        """
        for gene in self.genes.keys():
            print gene + ": " + str(len(self.genes[gene])) + " hits"         
            
    def print_hits(self):
        """
        Print out the list of hits
        """
        for hit in self.hits:
            print ",".join(str(item) for item in hit)
            
    def gene_list(self):
        """
        Return a list of all of the gene families
        """
        
        return self.genes.keys()
    
    def bug_list(self):
        """
        Return a list of all of the bugs
        """
        
        return self.bugs.keys()
    
    def hits_for_gene(self,gene):
        """
        Return the alignments for the selected gene
        """
        
        hit_list=[]        
        for index in self.genes[gene]:
            hit_list.append(self.hits[index])
            
        return hit_list
    
    def hits_for_bug(self,bug):
        """
        Return the alignments for the selected bug
        """
        if bug=="all":
            hit_list=self.hits
        else:
            hit_list=[]        
            for index in self.bugs[bug]:
                hit_list.append(self.hits[index])
            
        return hit_list
       
    def find_index(self,item):
        """
        Return the index number to the item in a hit
        """
        
        return self.hit_index.index(item)
    
class reactions_database:
    def __init__(self, database):
        """
        Load in the reactions data from the database
        """
        self.reactions_to_genes={}
        self.genes_to_reactions={}
         
        file_handle=open(database,"r")
         
        line=file_handle.readline()
         
        # database is expected to contain a single line per reaction
        # this line begins with the reaction name and ec number and is followed 
        # by all genes associated with the reaction
         
        while line:
            data=line.rstrip().split(config.reactions_database_delimiter)
            if len(data)>2:
                reaction=data.pop(0)
                ec_number=data.pop(0)
             
                # store the data
                self.reactions_to_genes[reaction]=data
             
                for gene in data:
                    if gene in self.genes_to_reactions:
                        self.genes_to_reactions[gene]+=[reaction]
                    else:
                        self.genes_to_reactions[gene]=[reaction]    
             
            line=file_handle.readline()
        
        file_handle.close()
        
    def find_reactions(self,gene):
        """
        Return the list of reactions associated with the gene
        """
        
        list=[]
        if self.gene_present(gene):
            list=self.genes_to_reactions[gene]
            
        return list
    
    def find_genes(self,reaction):
        """
        Return the list of genes associated with the reaction
        """
        
        list=[]
        if self.reaction_present(reaction):
            list=self.reactions_to_genes[reaction]
    
        return list
    
    def gene_present(self, gene):
        """
        Return true if gene is part of database
        """
        
        found=False
        if gene in self.genes_to_reactions.keys():
            found=True
        
        return found
    
    def reaction_present(self, reaction):
        """
        Return true if reaction is part of database
        """
        
        found=False
        if reaction in self.reactions_to_genes.keys():
            found=True
        
        return found
         
    def list_reactions(self):
        """
        Return the list of all the reactions in the database
        """
           
        return self.reactions_to_genes.keys()
