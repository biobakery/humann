"""
Stores the alignments identified
"""

class alignments:
    def __init__(self):
        self.hits=[]
        self.bugs={}
        self.genes={}
        
    def add(self, reference, query, evalue, identity, coverage, bug): 
        """ 
        Add the hit to the list
        Add the index of the hit to the bugs list and gene list
        """
        self.hits.append([reference, query, evalue, identity, coverage]) 
        
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

