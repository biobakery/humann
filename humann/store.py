"""
HUMAnN: store module
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
import math
import sys
import gzip
import bz2
import sqlite3

from . import config
from . import utilities

# name global logging instance
logger=logging.getLogger(__name__)

def store_id_mapping(file):
    """
    Store the id mapping data from the tab delimited file
    """

    id_mapping={}

    # Check the file exists and is readable
    utilities.file_exists_readable(file)
         
    file_handle=open(file,"rt")
         
    line=file_handle.readline()
    while line:
        # Ignore comment lines
        if not re.search(config.id_mapping_comment_indicator,line):
            data=line.rstrip().split(config.id_mapping_delimiter) 
            # set the default values for the mapping
            reference=""
            gene=""
            length=0
            bug="unclassified"
            try:
                reference=data[config.id_mapping_reference_index]
                gene=data[config.id_mapping_gene_index]
                length=data[config.id_mapping_gene_length_index]
                bug=data[config.id_mapping_bug_index]  
            except IndexError:
                logger.debug("Unable to read full mapping for id")
                    
            try:
                length=int(length)
            except ValueError:
                length=0
                    
            # if the reference and gene are found, store the mapping
            if reference and gene:
                id_mapping[reference]=[gene,length,bug]
                
        line=file_handle.readline()
        
    file_handle.close() 
    
    return id_mapping 

def normalized_gene_length(gene_length, read_length):
    """
    Compute the normalized gene length with the average read length if set
    Report in reads per kilobase
    """
    
    # if read length is not provided, default to 1
    if read_length < 1:
        read_length = 1

    return (abs(gene_length - read_length)+1)/1000.0

class SqliteStore:
    def __init__(self, minimize_memory_use = None):
        self.__minimize_memory_use=minimize_memory_use
        self.__dbpath = None
        self.__conn = None
        self.__is_within_transaction = False
        self.__stateful_ops_since_commit = None


    def connect(self):
        """
        Open the sqlite3 connection
        """
        if self.__conn:
            return

        if not self.__dbpath:
            store_name=type(self).__name__
            if self.__minimize_memory_use:
                self.__dbpath = utilities.unnamed_temp_file(store_name + ".sqlite")
                logger.debug("Initializing {0} store backed by a temporary file to minimize memory use".format(store_name))
            else:
                self.__dbpath = ":memory:"
                logger.debug("Initializing {0} store in-memory".format(store_name))
        
        self.__conn = sqlite3.connect(self.__dbpath, isolation_level=None)
        
    def do(self, *args):
        """
        Run a stateful statement like add or delete
        If within a transaction, commit and reopen every 100k operations
        """
        self.__conn.execute(*args)
        if self.__is_within_transaction:
            self.__stateful_ops_since_commit +=1
            if self.__stateful_ops_since_commit % 100000 == 0:
                self.__conn.execute("commit transaction")
                self.__conn.execute("begin transaction")

    def query(self, *args):
        """
        Use the sqlite3 connection
        """
        return self.__conn.execute(*args)

    def clear(self):
        """
        Clear all of the stored data
        """
        
        self.__conn.close()
        self.__conn = None

    def start_bulk_write(self):
        self.__is_within_transaction = True
        self.__stateful_ops_since_commit = 0
        self.__conn.execute("begin transaction")

    def end_bulk_write(self):
        self.__conn.execute("commit transaction")
        self.__is_within_transaction = False
        self.__stateful_ops_since_commit = None


class Alignments(SqliteStore):

    """
    Holds all of the alignments for all bugs
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.connect()
        self.do('''create table alignment (
              query text not null,
              bug text not null,
              reference text not null,
              score real not null,
              length real not null
            );''')
        self.__id_mapping={}
        
    def process_id_mapping(self,file):
        """
        Process the id mapping file
        """
        
        self.__id_mapping=store_id_mapping(file)
        
    def process_chocophlan_length(self,location,gene):
        """
        Return the length given the sequence location
        """
        
        try:
            if config.chocophlan_multiple_location_delimiter in location:
                locations=location.split(config.chocophlan_multiple_location_delimiter)
            else:
                locations=[location]
            length=0
            for location in locations:
                start, end = re.sub(config.chocophlan_location_extra_characters,
                    '',location).split(config.chocophlan_location_delimiter)
                length=length+abs(int(end)-int(start))+1
        except (ValueError, IndexError):
            length=0
            logger.debug("Unable to compute length for gene: " + gene)
        
        return length

    def process_reference_annotation(self,reference):
        """
        Process the reference string for information on gene, gene length, and bug
        Allow for chocophlan annotations, gene|gene_length, gene_length|gene, and gene
        Also use id mapping if provided
        """
        
        # if id mapping is provided first try to use it for the annotation data
        gene=""
        if self.__id_mapping:
            if reference in self.__id_mapping:
                [gene,length,bug]=self.__id_mapping[reference]
                
        # if id mapping is not provided or not found for the reference then
        # try to process the reference string
        if not gene:
            reference_info=reference.split(config.chocophlan_delimiter)
            
            # identify bug and gene families
            length=0
            gene=reference
            try:
                full_taxonomy=reference_info[config.chocophlan_bug_index]
                # Limit to species/genera
                bug_info=full_taxonomy.split(".")
                bug=".".join([bug_info[config.chocophlan_bug_genera_index],bug_info[config.chocophlan_bug_species_index]])

                # Join all genes selected
                gene_set=[]
                for index in config.chocophlan_gene_indexes:
                    gene_set.append(reference_info[index])
                if gene_set:
                    gene=config.chocophlan_delimiter.join(gene_set)
                length=int(reference_info[config.chocophlan_length_index])
            except (IndexError, ValueError):
                # try to find gene length if present
                bug="unclassified"
                # check for gene|gene_length|taxonomy
                if (len(reference_info)==3 and re.search("^[0-9]+$",reference_info[1])
                    and not re.search("^[0-9]+$",reference_info[2])):
                    bug=reference_info[2]
                    length=int(reference_info[1])
                    gene=reference_info[0]
                elif len(reference_info)==2:
                    if re.search("^[0-9]+$",reference_info[1]):
                        length=int(reference_info[1])
                        gene=reference_info[0]
                    elif re.search("^[0-9]+$",reference_info[0]):
                        length=int(reference_info[0])
                        gene=reference_info[1]
                    
        return [gene,length,bug]

    def add_annotated(self, query, matches, annotated_reference, read_length=None):
        """
        Add an alignment with an annotated reference
        """
        
        # Obtain the reference id length and bug
        [referenceid,length,bug]=self.process_reference_annotation(annotated_reference)
        
        self.add(referenceid, length, query, matches, bug, read_length)

    def add(self, reference, reference_length, query, matches, bug, read_length=None): 
        """ 
        Add the hit to the list
        Add the index of the hit to the bugs list and gene list
        """
        
        # set default read length
        if read_length is None:
            read_length = 1
        
        if reference_length==0:
            reference_length=config.default_reference_length
            logger.debug("Default gene length used for alignment to gene: " + reference)
        
        # store the score instead of the number of matches
        try:
            score=math.pow(matches,config.match_power)
        except ValueError:
            logger.debug("Could not convert the number of matches to score: " +  str(matches))
            score=0.0

        # Store the scores by bug and gene
        normalized_reference_length=normalized_gene_length(reference_length, read_length)
            
        # write the information for the hit
        self.do('insert into alignment (query, bug, reference, score, length) values (?,?,?,?,?)', [query, bug, reference, score, normalized_reference_length])

    def count_bugs(self):
        """ 
        Return total number of bugs
        """

        return self.query("select count (distinct bug) from alignment").fetchone()[0]
    
    def count_genes(self):
        """ 
        Return total number of genes
        """

        return self.query("select count (distinct reference) from alignment").fetchone()[0]
            
    def counts_by_bug(self):
        """
        Return each bug and the total number of hits
        """
 
        return "\n".join(["{0}: {1} hits".format(row[0], row[1]) for row in self.query('select bug, count(*) as c from alignment group by bug order by -c')])
            
    def gene_list(self):
        """
        Return a list of all of the gene families
        """
        
        return [row[0] for row in self.query('select distinct reference from alignment')]
    
    def bug_list(self):
        """
        Return a list of all of the bugs
        """

        return [row[0] for row in self.query('select distinct bug from alignment')]
    
    def get_hit_list(self):
        """
        Return a list of all of the hits
        """

        return [row for row in self.query('select query, bug, reference, score, length from alignment')]
    
    def hits_for_gene(self,gene):
        """
        Return a list of all of the hits for a specific gene
        """
        
        return [row for row in self.query('select query, bug, reference, score, length from alignment where reference=?', [gene])]
    
    def convert_alignments_to_gene_scores(self,gene_scores_store):
        """
        Computes the scores for all genes per bug
        Add to the gene_scores store
        """
        
        # calculate the score per bug and gene
        # for a single query result, it is 1/a.length
        # if a query matches multiple bugs and genes, its score is distributed by weighted average
        # scores from multiple queries are added up per bug and gene
        # see unit tests for examples
        result={}
        resultAll={}
        for bug, gene, score in self.query('''
              select bug, reference, sum(normalized_score_partial) as score from (
                  select
                    a.query,
                    a.bug,
                    a.reference,
                    sum(a.score / a.length ) / (total_score_for_query) as normalized_score_partial
                  from
                  alignment a join (
                    select query, sum(score) as total_score_for_query
                    from alignment group by query
                    ) as m
                  where a.query = m.query
                  group by a.bug, a.reference, a.query
               ) group by bug, reference
               order by bug, reference
            '''):
            if bug not in result:
                result[bug]={}
            result[bug][gene]=score
            if gene not in resultAll:
                resultAll[gene]=0
            resultAll[gene]+=score

        # Add to the store
        for bug in result:
            gene_scores_store.add(result[bug], bug)
        gene_scores_store.add(resultAll, "all")

        # Log a summary, and print if in verbose mode
        message="\n".join(["{0} : {1} gene families".format(bug, len(result[bug])) for bug in result])
        message="Total gene families  : " +str(sum([len(result[bug]) for bug in result]))+"\n"+message
        if config.verbose:
            print(message)
        logger.info("\n"+message)

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
            self.__scores[bug]=dict(list(self.__scores[bug].items()) + list(gene_scores.items()))
        else:
            self.__scores[bug]=gene_scores

    def add_single_score(self,bug,gene,score):
        """ 
        Add a score for a specific bug and gene
        """

        if bug in self.__scores:
            self.__scores[bug][gene]=score
        else:
            self.__scores[bug]={gene:score}
        
    def count_genes_for_bug(self,bug):
        """
        Count the total number of genes stored for all bugs
        """
        
        return len(self.__scores[bug])
        
    def get_score(self,bug,gene):
        """
        Get the score for a gene for bug
        """
        score=0
        if bug in self.__scores:
            score=self.__scores[bug].get(gene,0)
        
        return score
    
    def get_scores_for_gene_by_bug(self,gene):
        """
        Get the score for all bugs for a gene
        """
        scores={}
        for bug in self.__scores:
            score=self.__scores[bug].get(gene,0)
            if score>0 and bug != "all":
                scores[bug]=score
        
        return scores

    def bug_list(self):
        """
        Return a list of the bugs including "all"
        """
        
        return list(self.__scores.keys())

    def gene_list(self):
        """
        Return a list of the genes
        """
        genes={}
        for bug in self.__scores:
            for gene in self.__scores[bug]:
                genes[gene]=1
        
        return list(genes.keys())
    
    def gene_list_sorted_by_score(self,bug):
        """
        Return a list of the genes sorted by score for the bug
        """
        
        return utilities.double_sort(self.__scores.get(bug,{}))
    
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
    
    def add_from_file(self,file,id_mapping_file=None):
        """
        Add all of the gene scores from the file
        Use id mapping if provided
        """
        
        # Process the id mapping file if present
        id_mapping={}
        if id_mapping_file:
            id_mapping=store_id_mapping(id_mapping_file)
            
        unaligned_reads_count=0
        
        # Check the file exists and is readable
        utilities.file_exists_readable(file)
         
        file_handle=open(file,"rt")
         
        line=file_handle.readline()
        while line:
            # Ignore comment lines
            if not re.search(config.gene_table_comment_indicator,line):
                data=line.rstrip().split(config.gene_table_delimiter)
                gene=""
                bug="all"
                
                # Use id mapping if present
                if id_mapping:
                    if data[config.gene_table_gene_index] in id_mapping:
                         [gene,length,bug]=id_mapping[data[config.gene_table_gene_index]]
                
                # If gene not set with id mapping, then process
                if not gene:
                    if config.gene_table_category_delimiter in data[config.gene_table_gene_index]:
                        gene_data=data[config.gene_table_gene_index].split(
                            config.gene_table_category_delimiter)
                        gene=gene_data[0]
                        bug=gene_data[1]
                    else:
                        gene=data[config.gene_table_gene_index]
                    
                # remove the name of the gene if present
                if gene:
                    gene=gene.split(config.name_mapping_join)[0]
                    
                try:
                    value=float(data[config.gene_table_value_index])
                except (ValueError, IndexError):
                    value=0
                    if line.rstrip():
                        logger.debug("Unable to convert gene table value to float: %s",
                            line.rstrip())
                if gene == config.unmapped_gene_name:
                    unaligned_reads_count = value
                else:
                    self.add_single_score(bug,gene,value)
            line=file_handle.readline()

        file_handle.close()
        
        return unaligned_reads_count
    
class PathwaysAndReactions:
    """
    Holds all of the pathways and reaction scores for all bugs
    """
    
    def __init__(self):
        self.__pathways={}
        
    def add(self, bug, reaction, pathway, score): 
        """ 
        Add the pathway data to the dictionary
        """
        
        if not bug in self.__pathways:
            self.__pathways[bug]={}
        
        if pathway in self.__pathways[bug]:
            if reaction in self.__pathways[bug][pathway]:
                logger.debug("Overwrite of pathway/reaction score: %s %s", pathway, reaction)
            self.__pathways[bug][pathway][reaction]=score
        else:
            self.__pathways[bug][pathway]={ reaction : score }
            
    def bug_list(self):
        """
        Get a list of the bugs
        """
        
        return list(self.__pathways.keys())
    
    def pathway_list(self, bug):
        """
        Return the keys in the pathways dictionary for a bug
        """
        
        return list(self.__pathways.get(bug,{}).keys())
    
    def reaction_scores(self, bug, pathway):
        """
        Return the reactions in the pathways dictionary for a pathway and bug
        """
        
        return copy.copy(self.__pathways.get(bug,{}).get(pathway,{}))
    
    def median_score(self, bug):
        """
        Compute the median score for all scores in all pathways for one bug
        """
        
        # Create a list of all of the scores in all pathways
        all_scores=[]
        for item in self.__pathways.get(bug,{}).values():
            all_scores+=item.values()
        
        all_scores.sort()
        
        # Find the median score value
        median_score_value=0
        if all_scores:
            if len(all_scores) % 2 == 0:
                index1=int(len(all_scores)/2)
                index2=index1-1
                median_score_value=(all_scores[index1]+all_scores[index2])/2.0
            else:
                median_score_value=all_scores[int(len(all_scores)/2)]
            
        return median_score_value
    
    def max_median_score(self,bug):
        """
        Compute an alternative median score using the max values for the reactions for each pathway
        """
        
        # Create a list of the max scores for each of the pathways
        all_scores=[]
        for item in self.__pathways.get(bug,{}).values():
            all_scores.append(max(item.values()))
        
        all_scores.sort()
        
        # Find the median score value from the list of max scores per reaction
        median_score_value=0
        if all_scores:
            if len(all_scores) % 2 == 0:
                index1=int(len(all_scores)/2)
                index2=index1-1
                median_score_value=(all_scores[index1]+all_scores[index2])/2.0
            else:
                median_score_value=all_scores[int(len(all_scores)/2)]
            
        return median_score_value
    
    def count_pathways(self,bug):
        """
        Return the total number of pathways for a bug
        """
        
        return len(self.__pathways.get(bug,{}).keys())
    
class Pathways:
    """
    Holds the pathways coverage or abundance data for a bug
    """
    
    def __init__(self):
        self.__pathways={}
        self.__pathways_per_bug={}
        
    def add(self, bug, pathway, score):
        """
        Add the pathway score for the bug
        """
        
        # Try to convert the score to a float
        try:
            score=float(score)
        except ValueError:
            score=0
            logger.debug("Non-float value found for pathway: " + pathway)
        
        # Store all scores greater than 0
        if score>0:
            if bug == "all":
                self.__pathways[pathway]=score
            else:
                if pathway in self.__pathways_per_bug:
                    self.__pathways_per_bug[pathway][bug]=score
                else:
                    self.__pathways_per_bug[pathway]={bug:score}
            
    def delete(self, bug, pathway):
        """
        Delete the pathway for the bug
        """
        
        try:
            if bug == "all":
                del self.__pathways[pathway]
            else:
                del self.__pathways_per_bug[pathway][bug]
        except (KeyError,TypeError):
            pass
    
    def get_score(self, pathway):
        """
        Return the score for the pathway
        If the pathway does does not have a score, return 0
        """
        
        return self.__pathways.get(pathway,0)
    
    def get_score_for_bug(self, bug, pathway):
        """
        Return the score for the pathway for bug
        If the pathway does does not have a score, return 0
        """
        
        score=0
        if bug == "all":
            score=self.__pathways.get(pathway,0)
        elif pathway in self.__pathways_per_bug:
            score=self.__pathways_per_bug[pathway].get(bug,0)
         
        return score
    
    def get_pathways_double_sorted(self):
        """
        Return the pathways sorted by value
        """
        
        return utilities.double_sort(self.__pathways)
    
    def get_pathways_list(self):
        """
        Return the pathways list
        """
        
        return list(self.__pathways.keys())
    
    def get_bugs_double_sorted(self,pathway):
        """
        Return the bugs sorted by score for pathway
        """
        
        return utilities.double_sort(self.__pathways_per_bug.get(pathway,{}))

    def get_bugs_list(self):
        """
        Return the list of all of the bugs for all pathways
        """
        
        all_bugs=set()
        for pathway in self.__pathways_per_bug:
            all_bugs.update(self.__pathways_per_bug[pathway])
        
        return list(all_bugs)
        
    def get_pathways_and_bugs_nonzero_sorted(self):
        """
        Return lists that contains the pathways sorted with the
        largest abundance first, zero abundance pathways removed, with zero
        relative to the precision selected for the output file. Also include those
        bugs with each pathway, non-zero, sorted by abundance.
        """
        
        # Remove any pathways from the sorted list that are equivalent to zero based
        # on the precision selected for the output file
        # Double sort so that pathways with the same values are then sorted by name
        sorted_pathways_list=list(filter(lambda pathway: float(utilities.format_float_to_string(self.get_score(pathway))) > 0, self.get_pathways_double_sorted()))
        
        sorted_pathways_and_bugs=[]
        # Get the bugs for each of the pathways
        for pathway in sorted_pathways_list:
            # Remove any zero values based on precision selected and double sort
            bugs=list(filter(lambda bug: float(utilities.format_float_to_string(self.get_score_for_bug(bug, pathway))) > 0, self.get_bugs_double_sorted(pathway)))
            sorted_pathways_and_bugs.append([pathway,bugs])
            
        return sorted_pathways_and_bugs
    
class ReactionsDatabase:
    """
    Holds all of the genes/reactions data from the file provided
    """
    
    def __init__(self, database=None):
        """
        Load in the reactions data from the database
        """
        self.__reactions_to_genes={}
        self.__genes_to_reactions={}
        
        if not database is None:
            # Check the database file exists and is readable
            utilities.file_exists_readable(database)
            
            if database.endswith(".gz"):
                file_handle = gzip.open(database, "rt")
                readlines = file_handle.readlines()
            elif database.endswith(".bz2"):
                # read the bz2 file in binary and then decode for python2/3 compatibility
                file_handle = bz2.BZ2File(database, "r")
                readlines = [line.decode('utf-8') for line in file_handle.readlines()]
            else:
                file_handle=open(database,"rt")
                readlines = file_handle.readlines()
             
            # database is expected to contain a single line per reaction
            # this line begins with the reaction name and ec number and is followed 
            # by all genes associated with the reaction
             
            for line in readlines:
                data=line.rstrip().split(config.reactions_database_delimiter)
                if len(data)>2:
                    reaction=data.pop(0)
                    
                    if config.pathways_ec_column:
                        ec_number=data.pop(0)
                 
                    # store the data
                    self.__reactions_to_genes[reaction]=data
                 
                    for gene in data:
                        if gene in self.__genes_to_reactions:
                            self.__genes_to_reactions[gene]+=[reaction]
                        else:
                            self.__genes_to_reactions[gene]=[reaction]    
                 
            file_handle.close()
        
    def add_reactions(self, reactions):
        """
        Add these reactions and genes
        """
        
        for reaction in reactions:
            self.__reactions_to_genes[reaction]=self.__reactions_to_genes.get(reaction,[])+reactions[reaction]
            
            for gene in reactions[reaction]:
                self.__genes_to_reactions[gene]=self.__genes_to_reactions.get(gene,[])+[reaction]
        
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
    
    def _is_optional_reaction(self, item, reaction_names=None):
        """
        Check if this reaction is optional
        """
        
        if reaction_names is None:
            reaction_names=[]
        
        # count the number of indicators at the beginning of the name
        char=item[0]
        index=0
        optional_indicator_count=0
        while char == config.pathway_reaction_optional: 
            optional_indicator_count+=1
            index+=1
            try:
                char=item[index]
            except IndexError:
                break
            
        # this reaction has an optional indicator only if it is
        # not part of the original reaction name
        # original names can start with "--"  
        
        # first check for the original name in the reactions names if provided
        optional=False
        if item in reaction_names:
            optional=False
        elif optional_indicator_count in [1,3]:
            optional=True
            
        return optional
    
    def _find_reaction_list_and_key_reactions(self,items,reaction_names=None):
        """
        Find the reactions in the pathways items and also the key reactions
        """
        
        reaction_list=[]
        key_reactions=[]
        for item in items:
            # ignore items that are not reactions as they are part of pathway structure
            if not item in ["(",")","",config.pathway_AND,config.pathway_OR]:
                # check if the item name indicates an optional reaction
                # if so remove the optional reaction indicator at the beginning of the name                
                if self._is_optional_reaction(item, reaction_names):
                    item=item[1:]
                else:
                    # record that this is a key reaction
                    key_reactions.append(item)
                        
                reaction_list.append(item)      
        
        return reaction_list, key_reactions
        
    def _find_structure(self,items,reaction_names=None):
        """
        Find the structure of the pathway from the string
        """
        
        structure=[config.pathways_database_stucture_delimiter]
        levels={ 0: structure}
        current_level=0
            
        # Process through the list of strings
        for item in items:
            if item:
                # check if the item name indicates an optional reaction
                # if so remove the optional reaction indicator at the beginning of the name                
                if self._is_optional_reaction(item, reaction_names):
                    item=item[1:]
                        
                # Check if this is the start of a list
                if item == "(":
                    new_list=[config.pathways_database_stucture_delimiter]
                    levels[current_level+1]=new_list
                    # add the new list to the structure
                    levels[current_level].append(new_list)
                    # update the current level
                    current_level+=1
                # Check if this is the end of a list
                elif item == ")":
                    # Update the current level to close the list
                    current_level-=1
                # Check if this is a delimiter
                elif item in [config.pathway_AND,config.pathway_OR]:
                    # Update the delimiter at the beginning of the list
                    levels[current_level][0]=item
                else:
                    levels[current_level].append(item)

        return structure
        
    
    def _set_pathways_structure(self,reactions,reaction_names=None):
        """
        Determine the pathways structure from the input string
        """
        
        for pathway in reactions:
            # Check if the item is a list of items
            if isinstance(reactions[pathway], list):
                reactions[pathway]=config.pathways_database_stucture_delimiter.join(reactions[pathway])
            
            # Split the reactions information by the structured pathways delimiter
            reactions[pathway]=reactions[pathway].split(config.pathways_database_stucture_delimiter)
            
            # Find and store the structure for the pathway
            structure=self._find_structure(reactions[pathway],reaction_names)
            self.__pathways_structure[pathway]=structure
            
            # Find the list of reactions and the key reactions
            reaction_list, key_reactions = self._find_reaction_list_and_key_reactions(reactions[pathway],reaction_names)

            # Store the list of key reactions for the pathway
            self.__key_reactions[pathway]=key_reactions
            
            # Update the reactions dictionary to contain the list of reactions instead of the structure string
            reactions[pathway]=reaction_list
        
        return reactions
    
    def _store_pathways(self, reactions):
        """
        Create the dictionaries of reactions to pathways and pathways to reactions
        """
        
        for pathway in reactions:
            for reaction in reactions[pathway]:
                self.__pathways_to_reactions[pathway]=self.__pathways_to_reactions.get(
                    pathway,[]) + [reaction]
                self.__reactions_to_pathways[reaction]=self.__reactions_to_pathways.get(
                    reaction,[]) + [pathway]

    def __init__(self, database=None, reactions_database=None):
        """
        Load in the pathways data from the database
        """
        self.__pathways_to_reactions={}
        self.__reactions_to_pathways={}
        self.__pathways_structure={}
        self.__key_reactions={}
        
        reaction_names=None
        if not reactions_database is None:
            reaction_names=reactions_database.reaction_list()
        
        if not database is None:
            # Check the database file exists and is readable
            utilities.file_exists_readable(database)
            
            file_handle=open(database,"rt")
             
            line=file_handle.readline()
             
            # database is expected to contain a single line per pathway
            # this line begins with the pathway name and is followed 
            # by all reactions and/or pathways associated with the pathway
             
            reactions={}
            structured_pathway=False
            while line:
                data=line.strip().split(config.pathways_database_delimiter)
                if len(data)>1:
                    # replace any white spaces with underscores
                    pathway=data.pop(0).replace(" ","_")
                    reactions[pathway]=data
    
                    # check to see if the pathway has structure
                    if "(" in data[0]:
                        structured_pathway=True
                    
                line=file_handle.readline()
            
            file_handle.close()
            
            # if this is a structured pathways set, then store the structure
            if structured_pathway:
                reactions=self._set_pathways_structure(reactions, reaction_names)
            
            self._store_pathways(reactions)
            
    def is_structured(self):
        """
        Return True if this is a set of structured pathways
        """
        
        if self.__pathways_structure:
            return True
        else:
            return False
                    
    def add_pathway_structure(self, pathway, structure, reactions_database=None):
        """
        Add the string structure for a pathway
        """
        
        reaction_names=None
        if not reactions_database is None:
            reaction_names=reactions_database.reaction_list()
        
        reactions=self._set_pathways_structure({pathway: structure},reaction_names)
        self._store_pathways(reactions)
        
    def add_pathway(self, pathway, reactions):
        """
        Add the unstructured pathway
        """
        
        self._store_pathways({pathway: reactions})
        
    def get_structure_for_pathway(self,pathway):
        """ 
        Return the structure for a pathway
        """
        
        return copy.deepcopy(self.__pathways_structure.get(pathway, [])) 
        
    def get_key_reactions_for_pathway(self,pathway):
        """
        Return the key reactions for a pathway
        """
        
        return copy.copy(self.__key_reactions.get(pathway, []))
        
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
        
        return list(self.__reactions_to_pathways.keys())
    
    def pathway_list(self):
        """
        Return the list of pathways included in the database
        """
        
        return list(self.__pathways_to_reactions.keys())
    
    def get_database(self):
        """
        Return the database as a flat file with a single pathway per line
        """
        
        data=[]
        for pathway in sorted(self.__pathways_to_reactions):
            data.append(pathway+config.pathways_database_delimiter+
                config.pathways_database_delimiter.join(self.__pathways_to_reactions[pathway]))
        return "\n".join(data)
    
class Reads(SqliteStore):
    """
    Holds all of the reads data to create a fasta file
    """
    def __init__(self, file=None, **kwargs):
        super().__init__(**kwargs)
        self.connect()
        self.do('''create table read (
              id text not null primary key,
              sequence text not null
            );''')
    
    def add(self, id, sequence):
        """
        Store the sequence and id which should correspond to the following:
        >id
        sequence
        """
        self.do('insert or ignore into read (id, sequence) values (?,?)', [id, sequence])
            
    def process_file(self, file):
        """
        Process the file and yield ids and sequences
        """
        
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
                       
        file_handle=open(input_fasta,"rt")
            
        sequence=""
        id=""
        for line in file_handle:
            if re.search("^>", line):
                # store the prior sequence
                if id:
                    yield (id,sequence)
                id=line.rstrip().replace(">","")
                # only store the first id if multiple separated by spaces
                if " " in id:
                    ids=id.split(" ")
                    id=ids[0]
                sequence=""
            else:
                sequence+=line.rstrip()
            
        # add the last sequence
        yield (id, sequence)
                
        file_handle.close()
            
        # Remove the temp fasta file if exists
        if temp_file:
            utilities.remove_file(temp_file)
              
    def add_from_fasta(self, fasta_path):
        for (id,sequence) in self.process_file(fasta_path):
            self.add(id, sequence)

    def remove_id(self, id):
        """
        Remove the id and sequence from the read structure
        """
        self.do('delete from read where id = ? or id = ?', [id, utilities.remove_length_annotation(id)])

    def get_fasta(self):
        """ 
        Return a string of the fasta file sequences stored
        """
        for row in self.query('select id, sequence from read'):
            yield ">{0}\n{1}".format(*row)

    def id_list(self):
        """
        Return a list of all of the fasta ids
        """
        return [row[0] for row in self.query('select id from read')]
    
    def count_reads(self):
        """
        Return the total number of reads stored
        """
        
        return self.query("select count(*) from read").fetchone()[0]
         
class Names:
    """ 
    Holds all of the names that map to ids from a given file 
    """
    
    def __init__(self,file=None):
        """
        Read in the names and store
        """
        
        self.__names={}
        
        # Check the file exists and is readable
        unreadable_file=False
        if file is None:
            unreadable_file=True
        else:
            try:
                utilities.file_exists_readable(file, raise_IOError=True)
            except IOError:
                unreadable_file=True
                logger.debug("Unable to read Names file: " + file)
            
                       
        # Test if this is a gzipped file
        if not unreadable_file:
            if file.endswith(".gz"):
                file_handle = gzip.open(file, "rt")
                readlines = file_handle.readlines()
            elif file.endswith(".bz2"):
                # read the file in binary and then decode for python2/3 compatibility
                file_handle = bz2.BZ2File(file, "r")
                readlines = [line.decode('utf-8') for line in file_handle.readlines()]
            else:
                file_handle = open(file,"rt")
                readlines = file_handle.readlines()
                
            for line in readlines:
                try:
                    data = line.rstrip().split(config.name_mapping_file_delimiter)
                    id = data[0]
                    name = data[1]
                except IndexError:
                    id = ""
                
                if id:
                    self.__names[id]=name
                
            file_handle.close()
            
    def get_name(self,id):
        """
        Return the name for the given id
        """
        
        name = self.__names.get(id,"")
        
        if name:
            name = id + config.name_mapping_join + name
        else:
            # If there is not a name for the id, then just use the id
            name = id
        
        return name
    
