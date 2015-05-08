#!/usr/bin/python
"""
Create a structured pathways file from the MetaCyc flat pathways file (pathways.dat)

Dependencies: MetaCyc flat files

To Run: 
$ ./create_metacyc_structured_pathways_database.py --input pathways.dat --output metacyc_pathways.structured


There are a variety of different pathway structures to process
1) structures that do not have predecessors
2) super-pathways
3) recursive structures

"""

import sys

try:
    import argparse
except ImportError:
    sys.exit("ERROR: Please upgrade to at least python v2.7")

import os
import re
import copy

COMMENT_LINE="#"
METACYC_ID="UNIQUE-ID"
METACYC_PATHWAY_ID="PWY"
METACYC_ID_DELIMITER=" - "
METACYC_REACTION_LAYOUT_ID="^REACTION-LAYOUT"
METACYC_LEFT_PRIMARY=":LEFT-PRIMARIES"
METACYC_RIGHT_PRIMARY=":RIGHT-PRIMARIES"
METACYC_DIRECTION=":DIRECTION"
METACYC_RIGHT2LEFT=":R2L"
METACYC_LEFT2RIGHT=":L2R"
METACYC_KEYREACTION_ID="KEY-REACTIONS"
METACYC_PREDECESSORS="PREDECESSORS"
DELIMITER="\t"

OPTIONAL_REACTION_TAG="-"
OR_DEMILITER=" , "
AND_DEMILITER=" + "

def node_is_in_list(node,list_of_nodes):
    """
    Check if the node is present in the list of nodes
    """
    for node2 in list_of_nodes:
        if node is node2:
            return True
    return False

class PathwayStructure():
    def __init__(self):
        self.nodes=[]
        self.reactions=[]
        self.start_nodes=[]
        self.key_reactions=set()
        
    def add_node(self, data):
        if len(data) == 2:
            reaction, predecessor = data
            
            # check to see if predecessor exists
            predecessor_node = self.find_node(predecessor)
            
            if not predecessor_node:
                predecessor_node=Node(predecessor)
                self.nodes.append(predecessor_node)
            
            # check to see if the reaction exists
            reaction_node = self.find_node(reaction)
            
            if not reaction_node:
                reaction_node = Node(reaction)
                self.nodes.append(reaction_node)
                
            # set the predecessor and leaves
            predecessor_node.add_leaf(reaction_node)
            reaction_node.add_predecessor(predecessor_node)
            
        else:
            # add a new node that does not have a predecessor
            reaction = data[0]
            node=self.find_node(reaction)
            if not node:
                node=Node(reaction)
                self.nodes.append(node)
                
            self.start_nodes.append(node)
        
    def find_node(self, name):
        """ Find the node using the name """
        for node in self.nodes:
            if node.get_name() == name:
                return node
        return None
    
    def set_key_reactions(self, key_reactions):
        self.key_reactions=key_reactions
        
    def set_reactions(self, reactions):
        self.reactions=reactions
    
    def structure_to_string(self, item):
        """ Convert the structure of nested lists into a string """
        
        # test if item is a string
        if isinstance(item,basestring):
            return item
        # test if item is a list
        else:
            if item:
                if item[0] in [OR_DEMILITER,AND_DEMILITER]:
                    # then join the items with the delimiter
                    delimiter=item.pop(0)
                    string_start=" ( "
                    string_end=" ) "
                else:
                    delimiter=" "
                    string_start=" "
                    string_end=" "
            
                items_to_join=[]
                for i in item:
                    if isinstance(i, list):
                        items_to_join.append(self.structure_to_string(i))
                    else:
                        items_to_join.append(i) 
                        
                # remove any items that are empty
                items_to_join=filter(lambda x:x,items_to_join)
                if items_to_join:
                    if len(items_to_join) == 1:
                        string_start=" "
                        string_end=" "
                    return string_start+delimiter.join(items_to_join)+string_end
                else:
                    return ""
            else:
                return ""
        
    def create_structure(self):
        """ Get the structure for the pathway based on the levels of the nodes """
        
        # update the names of the nodes if key reactions are included
        if self.key_reactions:
            # find each of those nodes which are NOT key reactions and update their names
            for node in self.nodes:
                if not node.get_name() in self.key_reactions:
                    node.set_name(OPTIONAL_REACTION_TAG+node.get_name())
            #  also update the names of the reactions
            new_reactions={}
            for reaction in self.reactions:
                if not reaction in self.key_reactions:
                    # add the optional reaction identifier
                    new_reactions[OPTIONAL_REACTION_TAG+reaction]=self.reactions[reaction]
                else:
                    new_reactions[reaction]=self.reactions[reaction]
            self.reactions=new_reactions
        
        # set the levels for the nodes
        # all nodes that are specified by metacyc to not have predecessors are level 1
        for node in self.start_nodes:
            node.set_level(1)
            
        # identify nodes that do not currently have levels and set the levels
        for node in self.nodes:
            if node.get_level() is None and node.count_predecessors() == 0:
                node.set_level(1)
                
        # organize the nodes by levels
        nodes_by_levels={}
        for node in self.nodes:
            level=node.get_level()
            nodes_by_levels[level]=nodes_by_levels.get(level,[])+[node]
            
        # write out the nodes in order by level
        structure=[]
        prior_nodes=[]
        for level in sorted(nodes_by_levels):
            # start with the node with the most leaves
            for node in nodes_by_levels[level]:
                # only process if this node has not already been processed before
                if not node_is_in_list(node, prior_nodes):
                    structure_for_node, prior_nodes=node.create_structure(self.reactions,prior_nodes=prior_nodes)
                    structure+=structure_for_node
                    
        # add in nodes not accounted for already
        for node in self.nodes:
            if not node_is_in_list(node, prior_nodes):
                structure+=node.get_name()
                
        string_structure=self.structure_to_string(structure)
               
        # Check all reactions are included in structure
        for reaction in self.reactions:
            if not reaction in string_structure:
                string_structure+=" "+reaction
        
        # Check there are not any duplicate reactions
        check_string_structure=string_structure
        # Remove one of each of the reactions from the list
        removed_reactions=[]
        for node in self.nodes:
            name=node.get_name()
            check_string_structure=check_string_structure.replace(name,"",1)
            removed_reactions.append(name)
        for reaction in self.reactions:
            if not reaction in removed_reactions:
                check_string_structure=check_string_structure.replace(reaction,"",1)
                removed_reactions.append(reaction)
        
        # remove the join characters, parenthesis, and spaces
        check_string_structure_no_whitespace=check_string_structure.translate(None,
            "()"+OR_DEMILITER+AND_DEMILITER+OPTIONAL_REACTION_TAG).strip()
        
        # remove any duplicate reactions
        if not check_string_structure_no_whitespace == "":
            # search for the duplicated reactions
            for reaction in removed_reactions:
                occurances=check_string_structure.count(reaction)+1
                # remove the duplicated reactions starting with the last instances
                while occurances>1:
                    string_structure="".join(string_structure.rsplit(reaction,1))
                    # remove any extra AND/OR
                    string_structure=string_structure.replace(",    ,",",")
                    string_structure=string_structure.replace("+    +","+")
                    occurances-=1
                            
            # remove any unneeded AND/OR
            # check for at the end of parenthesis to remove extra joins
            string_structure=re.sub("\s*\,\s*\)"," )",string_structure)
            string_structure=re.sub("\s*\+\s*\)"," )",string_structure)
            string_structure=re.sub("\(\s*\)","",string_structure)
            
        # remove any extra whitespace
        string_structure=" ".join(string_structure.split())
    
        return string_structure

class Node():
    def __init__(self, name, level = None, predecessor = None, leaves = None):
        self.name = name
        self.predecessors = []
        if predecessor:
            self.predecessors = predecessors
        self.leaves = []
        if leaves:
            self.leaves = leaves
        self.level= None
        if level:
            self.level=level
            
    def get_name(self):
        return self.name
    
    def set_name(self, name):
        self.name=name
    
    def add_leaf(self,node):
        self.leaves.append(node)
        
    def add_predecessor(self,node):
        self.predecessors.append(node)
        
    def get_leaves(self):
        return copy.copy(self.leaves)
        
    def count_predecessors(self):
        return len(self.predecessors)
    
    def count_leaves(self):
        return len(self.leaves)
    
    def set_level(self,level,prior_nodes=None):
        """ Set the level of the node and the levels of the leaves and predecessors """
        if self.level is None:
            self.level=level
        
        if not prior_nodes:
            prior_nodes=[]
        
        # set the level of the leaves
        prior_nodes+=[self]
        for node in self.leaves:
            if not node_is_in_list(node, prior_nodes):
                node.set_level(level+1,prior_nodes)
                
        # set the level of the predecessors
        for node in self.predecessors:
            if not node_is_in_list(node, prior_nodes):
                node.set_level(level-1,prior_nodes)
        
    def get_level(self):
        return self.level
                
    def create_structure(self,reactions,prior_nodes=None):
        """
        Get the structure for the node and leaves
        """
        
        if not prior_nodes:
            prior_nodes=[]
        
        # search to see if self has already been visited in prior nodes
        self_visited_prior=False
        if node_is_in_list(self, prior_nodes):
            self_visited_prior=True
        else:
            prior_nodes+=[self]
            
        # remove recursive nodes from leaves if present
        non_recursive_leaves=[]
        for node in self.leaves:
            if not node_is_in_list(node, prior_nodes):
                non_recursive_leaves.append(node)
                
        # look ahead to see if this should be a contraction
        # see if any leaves have multiple precedessors
        multi_predecessor_leaves=[]
        leaves_reactants=set()
        for leaf in non_recursive_leaves:
            if leaf.count_predecessors() > 1:
                multi_predecessor_leaves+=leaf.predecessors
                # store the reactants from the leaves
                left,right=reactions.get(leaf.get_name(),[set(),set()])
                leaves_reactants.update(left)
                
        if multi_predecessor_leaves:
            # check if this a contraction or expansion
            products=set()
            for node in multi_predecessor_leaves:
                left,right=reactions.get(node.get_name(),[set(),set()])
                products.update(right)
                    
            new_list=[OR_DEMILITER]
            if len(products) > 1:
                # if the leaves only have a single reactant then this is an OR
                if len(leaves_reactants) == 1:
                    # check that this is not a space delimited list of reactants
                    if next(iter(leaves_reactants)).count(" ") > 0:
                        new_list=[AND_DEMILITER]
                    else:
                        new_list=[OR_DEMILITER]
                else:
                    new_list=[AND_DEMILITER]
                
            # add the name of this node if it was not visited prior
            if not self_visited_prior:
                new_list+=[self.name]
            
            for node in multi_predecessor_leaves:
                if not node_is_in_list(node, prior_nodes):
                    new_list.append(node.get_name())
                    # add the predecessors to the prior nodes
                    prior_nodes+=[node]
                
            structure=[new_list]    
        else:
            # add the name of this node if it was not visited prior
            if not self_visited_prior:
                structure=[self.name]
            else:
                structure=[]
        leaf_structures=[]
        # find the structures for the leaves
        for node in non_recursive_leaves:
            leaf_structure, prior_nodes=node.create_structure(reactions,prior_nodes)
            leaf_structures.append(leaf_structure)
            
        if len(non_recursive_leaves) == 1:
            structure+=leaf_structures[0]
        elif len(non_recursive_leaves) > 1:
            # OR expansion
            structure+=[[OR_DEMILITER]+leaf_structures]

        return structure, prior_nodes
        
def write_structured_pathways(metacyc_pathways, output_file):
    """
    Write the structure for each of the pathways
    """
    
    try:
        file_handle=open(output_file,"w")
    except EnvironmentError:
        sys.exit("Unable to write to file: "+output_file)
     
    for pathway in metacyc_pathways:
        structure=metacyc_pathways[pathway].create_structure()
        if not structure:
            sys.exit("MISSING structure for pathway: "+pathway)
        file_handle.write(pathway + DELIMITER + structure + "\n")
        
    file_handle.close()
            
def read_metacyc_pathways_structure(metacyc_pathways_file):
    """
    Process the metacyc pathways file to read the structure of the pathways
    """
    
    metacyc_pathways={}
    
    try:
        file_handle=open(metacyc_pathways_file,"r")
        line=file_handle.readline()
    except EnvironmentError:
        sys.exit("Unable to read file: " + metacyc_pathways_file)
        
    pathway=""
    key_reactions=set()
    reactions={}
    pathway_structure=PathwayStructure()
    while line:
        if not re.match(COMMENT_LINE, line):
            # find the pathway id
            if re.match(METACYC_ID, line):
                if pathway:
                    # record the last pathway
                    pathway_structure.set_key_reactions(key_reactions)
                    pathway_structure.set_reactions(reactions)
                    metacyc_pathways[pathway]=pathway_structure
                        
                pathway=""
                key_reactions=set()
                pathway_structure=PathwayStructure()
                reactions={}
                pathway=line.rstrip().split(METACYC_ID_DELIMITER)[-1]
                
            # find the ordering for the pathway
            elif re.match(METACYC_PREDECESSORS, line):
                pathway_structure.add_node(line.rstrip().split(METACYC_ID_DELIMITER)[-1].replace("(","").replace(")","").split(" "))
                
            # find the key reactions for the pathway
            elif re.match(METACYC_KEYREACTION_ID, line):
                key_reactions.add(line.rstrip().split(METACYC_ID_DELIMITER)[-1])
                
            # find the reaction information 
            elif re.match(METACYC_REACTION_LAYOUT_ID, line):
                # process the reaction layout
                # example format
                # REACTION-LAYOUT - (RXN-1401 (:LEFT-PRIMARIES TRYPTAMINE) (:DIRECTION :L2R) (:RIGHT-PRIMARIES INDOLE_ACETALDEHYDE))
                data=[re.sub("\)","",item).rstrip() for item in line.split("(")][1:]
                reaction=data.pop(0)
                
                left=set()
                right=set()
                direction=""
                for item in data:
                    # check for the type of identifier
                    if re.search(METACYC_LEFT_PRIMARY,item):
                        left.add(item.replace(METACYC_LEFT_PRIMARY,"").strip())
                    elif re.search(METACYC_RIGHT_PRIMARY,item):
                        right.add(item.replace(METACYC_RIGHT_PRIMARY,"").strip())                      
                    elif re.search(METACYC_DIRECTION,item):
                        direction=item.replace(METACYC_DIRECTION,"").strip()
                if re.search(METACYC_LEFT2RIGHT,direction):
                    reactions[reaction]=[left,right]
                elif re.search(METACYC_RIGHT2LEFT,direction):
                    reactions[reaction]=[right,left]
                    
        line=file_handle.readline()
        
    # store the last pathway
    if pathway:
        pathway_structure.set_key_reactions(key_reactions)
        pathway_structure.set_reactions(reactions)
        metacyc_pathways[pathway]=pathway_structure
        
    file_handle.close()
    
    return metacyc_pathways


def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Create structured MetaCyc pathways file\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-v","--verbose", 
        help="additional output is printed\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "-i","--input",
        help="the original MetaCyc pathways file\n",
        required=True)
    parser.add_argument(
        "-o","--output",
        help="the file to write the structured pathways\n",
        required=True)

    return parser.parse_args()


def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    input_file=os.path.abspath(args.input)
    output_file=os.path.abspath(args.output)

    # read the metacyc pathways file
    if args.verbose:
        print("Processing metacyc pathways data")
    metacyc_pathways=read_metacyc_pathways_structure(input_file)
    
    write_structured_pathways(metacyc_pathways, output_file)
    
    total_pathways=len(metacyc_pathways)
    if args.verbose:
        print("Total pathways identified from metacyc: " +str(total_pathways))
         

if __name__ == "__main__":
    main()
