"""
HUMAnN2: store module
Stores the alignments identified
Stores the genes/reactions from the database selected
Stores the reactions/pathways from the database selected
Stores the pathways identified
Stores the unaligned reads

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
import os
import re
import logging
import copy

import config
import utilities

# name global logging instance
logger=logging.getLogger(__name__)

class Alignments:
    """
    Holds all of the alignments for all bugs
    """
    
    def __init__(self):
        self.__hits=[]
        self.__bugs={}
        self.__genes={}
        
    def add(self, reference, reference_length, query, evalue, bug): 
        """ 
        Add the hit to the list
        Add the index of the hit to the bugs list and gene list
        """
        
        if reference_length==0:
            reference_length=1000
            logger.debug("Reference length of 0 found for gene: " + reference)
        
        # store the reference length as per kilobase
        length=reference_length/1000.0
            
        self.__hits.append([bug, reference, length, query, evalue]) 
        
        index=len(self.__hits)-1
        if bug in self.__bugs:
            self.__bugs[bug].append(index)
        else:
            self.__bugs[bug]=[index]
            
        if reference in self.__genes:
            self.__genes[reference].append(index)
        else:
            self.__genes[reference]=[index]
            
    def count_bugs(self):
        """ 
        Return total number of bugs
        """
        return len(self.__bugs.keys())
    
    def count_genes(self):
        """ 
        Return total number of genes
        """
        return len(self.__genes.keys())    
            
    def counts_by_bug(self):
        """
        Return each bug and the total number of hits
        """
        lines=[]
        for bug in self.__bugs.keys():
            lines.append(bug + ": " + str(len(self.__bugs[bug])) + " hits")
            
        return "\n".join(lines)
            
    def gene_list(self):
        """
        Return a list of all of the gene families
        """
        
        return self.__genes.keys()
    
    def bug_list(self):
        """
        Return a list of all of the bugs
        """
        
        return self.__bugs.keys()
    
    def hits_for_gene(self,gene):
        """
        Return the alignments for the selected gene
        """
        
        hit_list=[]        
        for index in self.__genes[gene]:
            hit_list.append(self.__hits[index])
            
        return hit_list
    
    def hits_for_bug(self,bug):
        """
        Return the alignments for the selected bug
        """
        hit_list=[]        
        for index in self.__bugs[bug]:
            hit_list.append(self.__hits[index])
            
        return hit_list
    
    def all_hits(self):
        """
        Return all of the hits
        """
        hit_list=[]
        for hit in self.__hits:
            if hit:
                hit_list.append(hit)
                
        return hit_list
    
    def delete_gene_and_hits(self,gene):
        """
        Remove the gene and all data for all hits associated with the gene
        """
        
        if gene in self.__genes:
            for index in self.__genes[gene]:
                # Remove the data for the hit
                self.__hits[index]=[]
            # Remove the gene entry
            del self.__genes[gene]
            
    def update_hits_for_bugs(self):
        """
        Update the hit indexes associated with all of the bugs to remove
        any indexes that point to deleted hits
        """
        
        bugs_to_delete=[]
        for bug in self.__bugs:
            updated_indexes=[]
            for index in self.__bugs[bug]:
                # Check if the hit is empty
                if self.__hits[index]:
                    updated_indexes.append(index)
            if updated_indexes:
                self.__bugs[bug]=updated_indexes
            else:
                # if there are no hits remaining for the bug,
                # then indicate the bug should be deleted
                bugs_to_delete.append(bug)
        
        # delete the bugs that do not have any hits
        for bug in bugs_to_delete:
            logger.debug("Bug removed because no longer associated with any hits: %s", bug)
            del self.__bugs[bug]
            
    def filter_hits(self,reactions_database):
        """
        Remove hits from the alignments database that are not associated with
        a gene in the reactions database
        """
        
        # Remove the hits from the alignments that are not associated with a gene
        # in the gene to reactions database
        message="Remove hits to genes not in reactions database ..."
        logger.info(message)
        if config.verbose:
            print(message)
    
        logger.debug("Remove gene and hits")
        for gene in self.gene_list():
            if not reactions_database.gene_present(gene):
                self.delete_gene_and_hits(gene)
                
        # Update the bugs index list to remove any indexes that point to deleted hits
        logger.debug("Update hits for bugs to reflect removed hits")        
        self.update_hits_for_bugs()
                
        message="Total bugs after filtering: " + str(self.count_bugs())
        logger.info(message)        
        print(message)
        
        message=self.counts_by_bug()
        logger.info(message)
        print(message)            
                
        message="Total gene families after filtering: " + str(len(self.gene_list()))
        logger.info(message)
        print(message)
        
class GeneScores:
    """
    Holds scores for all of the genes
    """
    
    def __init__(self):
        self.__scores={}
        
    def add(self,gene_scores,bug):
        """ 
        Add gene scores for a specific bug
        """
        
        if bug in self.__scores:
            self.__scores[bug]=dict(self.__scores[bug].items() + gene_scores.items())
        else:
            self.__scores[bug]=gene_scores
        
    def count_genes_for_bug(self,bug):
        """
        Count the total number of genes stored for all bugs
        """
        
        return len(self.__scores[bug].keys())
        
    def get_score(self,bug,gene):
        """
        Get the score for a gene for bug
        """
        score=0
        if bug in self.__scores:
            score=self.__scores[bug].get(gene,0)
        
        return score

    def bug_list(self):
        """
        Return a list of the bugs including "all"
        """
        
        return self.__scores.keys()
    
    def scores_for_bug(self,bug):
        """
        Return the gene scores for a specific bug
        """
        
        scores={}
        if bug in self.__scores:
            scores=copy.copy(self.__scores[bug])
        else:
            logger.debug("Request for scores for bug that does not exist.")
        
        return scores
    
class PathwaysAndReactions:
    """
    Holds all of the pathways and reaction scores for one bug
    """
    
    def __init__(self, bug):
        self.__pathways={}
        self.__bug=bug
        
    def add(self, reaction, pathway, score): 
        """ 
        Add the pathway data to the dictionary
        """
        
        if pathway in self.__pathways:
            if reaction in self.__pathways[pathway]:
                logger.debug("Overwrite of pathway/reaction score: %s %s", pathway, reaction)
            self.__pathways[pathway][reaction]=score
        else:
            self.__pathways[pathway]={ reaction : score }
    
    def get_bug(self):
        """
        Return the bug associated with the pathways data
        """
        
        return str(self.__bug)
    
    def get_pathways(self):
        """
        Return the keys in the pathways dictionary
        """
        
        return self.__pathways.keys()
    
    def get_reactions(self, pathway):
        """
        Return the reactions in the pathways dictionary for a pathway
        """
        
        return copy.copy(self.__pathways.get(pathway,{}))
    
    def median_score(self):
        """
        Compute the median score for all scores in all pathways
        """
        
        # Create a list of all of the scores in all pathways
        all_scores=[]
        for item in self.__pathways.values():
            all_scores+=item.values()
        
        all_scores.sort()
        
        # Find the median score value
        median_score_value=0
        if all_scores:
            if len(all_scores) % 2 == 0:
                index1=len(all_scores)/2
                index2=index1-1
                median_score_value=(all_scores[index1]+all_scores[index2])/2.0
            else:
                median_score_value=all_scores[len(all_scores)/2]
            
        return median_score_value
    
    def count_pathways(self):
        """
        Return the total number of pathways
        """
        
        return len(self.__pathways.keys())
    
class Pathways:
    """
    Holds the pathways coverage or abundance data for a bug
    """
    
    def __init__(self, bug="None", pathways={}):
        self.__bug=bug
        self.__pathways=pathways

    def get_bug(self):
        """
        Return the bug associated with the pathways data
        """
        
        return str(self.__bug)
    
    def get_score(self, pathway):
        """
        Return the score for the pathway
        If the pathway does does not have a score, return 0
        """
        
        try:
            score=float(self.__pathways.get(pathway,0))
        except ValueError:
            score=0
            logger.debug("Non-float value found for pathway: " + pathway)
        return score
    
    def get_pathways(self):
        """
        Return the keys in the pathways dictionary
        """
        
        return self.__pathways.keys()
    
    def get_pathways_double_sorted(self):
        """
        Return the pathways sorted by value
        """
        
        return utilities.double_sort(self.__pathways)
        
    
class ReactionsDatabase:
    """
    Holds all of the genes/reactions data from the file provided
    """
    
    def __init__(self, database):
        """
        Load in the reactions data from the database
        """
        self.__reactions_to_genes={}
        self.__genes_to_reactions={}
        
        # Check the database file exists and is readable
        utilities.file_exists_readable(database)
         
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
                self.__reactions_to_genes[reaction]=data
             
                for gene in data:
                    if gene in self.__genes_to_reactions:
                        self.__genes_to_reactions[gene]+=[reaction]
                    else:
                        self.__genes_to_reactions[gene]=[reaction]    
             
            line=file_handle.readline()
        
        file_handle.close()
        
    def find_reactions(self,gene):
        """
        Return the list of reactions associated with the gene
        """
            
        return copy.copy(self.__genes_to_reactions.get(gene,[]))
    
    def find_genes(self,reaction):
        """
        Return the list of genes associated with the reaction
        """
        
        return copy.copy(self.__reactions_to_genes.get(reaction,[]))
         
    def reaction_list(self):
        """
        Return the list of all the reactions in the database
        """
           
        return self.__reactions_to_genes.keys()
    
    def gene_list(self):
        """
        Return the list of all the genes in the database
        """
           
        return self.__genes_to_reactions.keys()
    
    def gene_present(self, gene):
        """
        Check if the gene is included in the database
        """
        
        present=False
        if gene in self.__genes_to_reactions:
            present=True
            
        return present
    
class PathwaysDatabase:
    """
    Holds all of the reactions/pathways data from the file provided
    """
    
    def _is_pathway(self, item):
        """
        Determine if the item is a pathway or reaction
        """
        
        pathway=False
        # identifier can be at the beginning or end of the string
        if re.search("^"+config.pathway_identifier, 
            item) or re.search(config.pathway_identifier+"$", item):
            pathway=True
        
        return pathway    
    
    def _return_reactions(self, pathway, reactions):
        """
        Search recursively to find the reactions associated with the given pathway
        """
        
        reactions_for_pathway=[]
        for item in reactions.get(pathway,[]):
            # go through items to look for pathways to resolve
            if self._is_pathway(item):
                # find the reactions for the pathway
                reactions_for_pathway+=self._return_reactions(item, reactions)
            else:
                reactions_for_pathway+=[item]
                
        return reactions_for_pathway

    def __init__(self, database, recursion):
        """
        Load in the pathways data from the database
        """
        self.__pathways_to_reactions={}
        self.__reactions_to_pathways={}
        
        # Check the database file exists and is readable
        utilities.file_exists_readable(database)
        
        file_handle=open(database,"r")
         
        line=file_handle.readline()
         
        # database is expected to contain a single line per pathway
        # this line begins with the pathway name and is followed 
        # by all reactions and/or pathways associated with the pathway
         
        reactions={}
        while line:
            data=line.strip().split(config.pathways_database_delimiter)
            if len(data)>1:
                pathway=data.pop(0)
                reactions[pathway]=data
                
            line=file_handle.readline()
        
        file_handle.close()
        
        # process recursive pathways
        for pathway in reactions:
            for item in reactions[pathway]:
                # go through items to look for pathways to resolve
                reaction=[item]
                # identifier can be at the start or the end of the item name
                if self._is_pathway(item) and recursion:
                    # find the reactions for the pathway
                    reaction=self._return_reactions(item, reactions)
                
                self.__pathways_to_reactions[pathway]=self.__pathways_to_reactions.get(
                    pathway,[]) + reaction
                    
        # store all pathways associated with a reaction
        for pathway in self.__pathways_to_reactions:
            for reaction in self.__pathways_to_reactions[pathway]:
                self.__reactions_to_pathways[reaction]=self.__reactions_to_pathways.get(
                    reaction,[]) + [pathway]
    
    def find_reactions(self,pathway):
        """
        Return the list of reactions associated with the pathway
        """
         
        return copy.copy(self.__pathways_to_reactions.get(pathway, []))

    def find_pathways(self,reaction):
        """
        Return the list of pathways associated with the reaction
        """
         
        return copy.copy(self.__reactions_to_pathways.get(reaction, []))
    
    def reaction_list(self):
        """
        Return the list of reactions included in the database
        """
        
        return self.__reactions_to_pathways.keys()
    
    def pathway_list(self):
        """
        Return the list of pathways included in the database
        """
        
        return self.__pathways_to_reactions.keys()
    
    def get_database(self):
        """
        Return the database as a flat file with a single pathway per line
        """
        
        data=[]
        for pathway in self.__pathways_to_reactions:
            data.append(pathway+config.pathways_database_delimiter+
                config.pathways_database_delimiter.join(self.__pathways_to_reactions[pathway]))
            
        return "\n".join(data)
    
class Reads:
    """
    Holds all of the reads data to create a fasta file
    """
    
    def add(self, id, sequence):
        """
        Store the sequence and id which should correspond to the following:
        >id
        sequence
        """
        
        self.__reads[id]=sequence
    
    def __init__(self, file=None):
        """
        Create initial data structures and load if file name provided
        """
        self.__reads={}
              
        if file:
            
            # Check the file exists and is readable
            utilities.file_exists_readable(file)
            
            # Check that the file of reads is fasta
            # If it is fastq, then convert the file to fasta
            temp_file=""
            if utilities.fasta_or_fastq(file) == "fastq":
                input_fasta=utilities.fastq_to_fasta(file)
                temp_file=input_fasta
            else:
                input_fasta=file
                       
            file_handle=open(input_fasta,"r")
            
            sequence=""
            id=""
            for line in file_handle:
                if re.search("^>", line):
                    # store the prior sequence
                    if id:
                        self.add(id, sequence)
                    id=line.rstrip().replace(">","")
                    sequence=""
                else:
                    sequence+=line.rstrip()
            
            # add the last sequence
            self.add(id, sequence)
                
            file_handle.close()
            
            # Remove the temp fasta file if exists
            if temp_file:
                utilities.remove_file(temp_file)
    
    def remove_id(self, id):
        """
        Remove the id and sequence from the read structure
        """
        
        if id in self.__reads:
            del self.__reads[id]
                
    def get_fasta(self):
        """ 
        Return a string of the fasta file sequences stored
        """
        
        fasta=[]
        for id, sequence in self.__reads.items():
            fasta.append(">"+id+"\n"+sequence)
        
        return "\n".join(fasta)
    
    def id_list(self):
        """
        Return a list of all of the fasta ids
        """
        
        return self.__reads.keys()
    
    def count_reads(self):
        """
        Return the total number of reads stored
        """
            
        return len(self.__reads.keys())
    