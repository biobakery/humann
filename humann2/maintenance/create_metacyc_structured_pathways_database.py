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
OR_DELIMITER=" , "
AND_DELIMITER=" + "

def node_is_in_list(node,list_of_nodes):
    """
    Check if the node is present in the list of nodes
    """
    for node2 in list_of_nodes:
        if node is node2:
            return True
    return False

def find_list(input_list, structure):
    """ Find the list and index which contains the input list in the structure of lists """
            
    # Check through each of the lists and sub-lists to see if the input_list is included
    match=None,None
    for i in xrange(len(structure)):
        if isinstance(structure[i],list):
            if structure[i] is input_list:
                return structure, i
            else:
                if match[0] is None:
                    match=find_list(input_list,structure[i])
    return match

class PathwayStructure():
    def __init__(self):
        self.nodes=[]
        self.reactions=[]
        self.start_nodes=[]
        self.key_reactions=set()
        self.structure=[]
        self.subpathways=set()
        
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
        
    def get_structure(self):
        return self.structure
    
    def get_string_structure(self):
        return self.structure_to_string().replace("  "," ")
    
    def structure_to_string(self, item=None, subset=None):
        """ Convert the structure of nested lists into a string """
        
        # Check if the item is unset
        if item is None:
            item=self.structure
        
        # test if item is a string
        if isinstance(item,basestring):
            return item
        # test if item is a list
        else:
            if item:
                delimiter=" "
                string_start=" "
                string_end=" "
                next_subset=False
                if item[0] in [OR_DELIMITER,AND_DELIMITER]:
                    # if present, then join the items with the delimiter
                    delimiter=item.pop(0)
                    if len(item) > 1:
                        string_start=" ( "
                        string_end=" ) "
                        next_subset=True
                    else:
                        # do not use the delimiter if there is only one item
                        delimiter=" "
                elif subset:
                    string_start=" ( "
                    string_end=" ) "
            
                items_to_join=[]
                for i in item:
                    if isinstance(i, list):
                        items_to_join.append(self.structure_to_string(i,next_subset))
                    else:
                        items_to_join.append(i) 
                        
                # remove any items that are empty
                items_to_join=filter(lambda x:x,items_to_join)
                if items_to_join:
                    if len(items_to_join) == 1:
                        string_start=" "
                        string_end=" "
                        delimiter=" "
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
        structure_by_level={}
        prior_nodes={}
        insert=False
        for level in sorted(nodes_by_levels):
            # start with the node with the most leaves
            for node in nodes_by_levels[level]:
                # only process if this node has not already been processed before
                if not node_is_in_list(node, prior_nodes.keys()):
                    if insert:
                        structure_to_insert, insert_node, join_type=node.create_structure(self.reactions,prior_nodes)
                        # add this structure to the main structure
                        if not insert_node is None:
                            # check if the insert is into the main structure
                            if prior_nodes[insert_node] is structure:
                                # add this as a join to the top items in the main structure
                                # find those items in the beginning of the structure before any joins
                                # by looking for the first list
                                join_location=len(structure)
                                for i,item in enumerate(structure):
                                    if isinstance(item, list):
                                        join_location=i
                                        break
                                # this will always be a OR join as it is two ways to start the pathway
                                join_type=OR_DELIMITER
                                    
                                # if the join location is to use the full structure
                                if join_location == len(structure):
                                    new_structure=[join_type,structure,structure_to_insert]
                                    structure=new_structure
                                else:
                                    new_structure=[join_type,structure[0:join_location],structure_to_insert]
                                    del structure[0:join_location]
                                    structure=[new_structure]+structure
                            else:
                                # find the superlist that contains this node
                                insert_list, insert_location=find_list(prior_nodes[insert_node], structure)
                                
                                # check if an AND/OR needs to be added
                                if not insert_list[0] in [AND_DELIMITER,OR_DELIMITER]:
                                    new_list=[join_type]
                                    new_list.append(prior_nodes[insert_node])
                                    new_list.append(structure_to_insert)
                                    insert_list[insert_location]=new_list
                                else:
                                    insert_list.append(structure_to_insert)
                        else:
                            structure+=structure_to_insert
                    else:
                        # the first time this is run create the structure
                        # on other runs, insert into the structure
                        structure, insert_node, join_type=node.create_structure(self.reactions,prior_nodes,structure_by_level)
                        insert=True
                    
        # add in nodes not accounted for already
        for node in self.nodes:
            if not node_is_in_list(node, prior_nodes):
                structure+=node.get_name()
                
        self.structure=structure
        
    def replace_subpathway(self, subpathway, subpathway_structure, structure=None):
        """ Replace the subpathway id with the subpathway structure """
        
        if structure is None:
            structure=self.structure
        
        # For each instance of the subpathway in the structure, replace with the subpathway structure
        for i in xrange(len(structure)):
            if isinstance(structure[i],list):
                self.replace_subpathway(subpathway, subpathway_structure, structure[i])
            elif isinstance(structure[i],basestring):
                if structure[i] == subpathway:
                    structure[i]=subpathway_structure
                    
    def find_reaction(self, reaction, structure=None):
        """ Check if the reaction is included in the structure """
        
        if structure is None:
            structure=self.structure
        
        # Check through each of the lists and sub-lists to see if the reaction is included
        for i in xrange(len(structure)):
            if isinstance(structure[i],list):
                if self.find_reaction(reaction, structure[i]):
                    return True
            elif isinstance(structure[i],basestring):
                if structure[i] == reaction:
                    return True
                
        return False
    
    def count_reactions(self, reaction_count=None, structure=None):
        """ Count the number of times each reaction is included in the structure """
        
        if structure is None:
            structure=self.structure
            
        if reaction_count is None:
            reaction_count={}
        
        # Check through each of the lists and sub-lists to count the reactions
        for i in xrange(len(structure)):
            if isinstance(structure[i],list):
                reaction_count=self.count_reactions(reaction_count, structure[i])
            elif isinstance(structure[i],basestring):
                # do not count AND/OR
                if not structure[i] in [AND_DELIMITER,OR_DELIMITER]:
                    reaction_count[structure[i]]=reaction_count.get(structure[i],0)+1
                
        return reaction_count
    
    def resolve_subpathways(self,metacyc_pathway_structures):
        """ Resolve the subpathways in the metacyc structures """
        
        # Get a list of all of the subpathways for this pathway
        # Check both the reactions list and the nodes list
        # Pathways can be included in the predecessors in nodes but not in the
        # reactions layout sections
        self.subpathways=set()
        names_to_check=list(self.reactions)+[node.get_name() for node in self.nodes]
        for name in names_to_check:
            if ( re.match(OPTIONAL_REACTION_TAG+"*"+METACYC_PATHWAY_ID, name) 
                or re.search(METACYC_PATHWAY_ID+"$", name) ):
                self.subpathways.add(name)
                
        # Get the structure for the subpathway
        # Then replace the subpathway name with the structure for the subpathway
        for subpathway in self.subpathways:
            subpathway_id=subpathway
            # remove the optional reaction identifier if present
            if re.match(OPTIONAL_REACTION_TAG,subpathway_id):
                subpathway_id=subpathway_id.replace(OPTIONAL_REACTION_TAG,"",1)
            if subpathway_id in metacyc_pathway_structures:
                subpathway_structure=metacyc_pathway_structures[subpathway_id].get_structure()
            else:
                print("WARNING: Missing subpathway from database: " + subpathway_id)
                subpathway_structure=[""]
            # replace the pathway id with the structure for the pathway
            self.replace_subpathway(subpathway,subpathway_structure)
        
        # Check all reactions are included in structure (excluding subpathways)
        for reaction in self.reactions:
            if not reaction in self.subpathways:
                if not self.find_reaction(reaction):
                    self.structure+=[reaction]
                        
    def remove_duplicate_reaction(self, reaction, count, structure):
        """ Remove the duplicate reactions (keeping the last instance in the structure) """
        
        sublists=[]
        for i in xrange(len(structure)):
            if isinstance(structure[i],list):
                sublists.append(i)
            elif isinstance(structure[i],basestring):
                if structure[i] == reaction:
                    if count > 1:
                        structure[i]=""
                        count-=1
                        
        # process through any sublists in the main structure list
        for i in sublists:
            structure[i],count=self.remove_duplicate_reaction(reaction, count, structure[i])
                        
        return structure, count
                
    def resolve_duplicate_reactions(self):
        """ Check for duplicate reactions in the pathway and resolve """
        
        # Check there are not any duplicate reactions by counting the reactions
        reaction_count=self.count_reactions()
        
        # Check for duplicate reactions that are optional and key reactions in pathway
        duplicate_optional_reactions=set()
        for reaction,count in reaction_count.items():
            if OPTIONAL_REACTION_TAG+reaction in reaction_count:
                duplicate_optional_reactions.add(OPTIONAL_REACTION_TAG+reaction)
                
        # Update the structure to a deep copy of each object to
        # avoid recursive objects which reference themselves in the copy
        # This can happen with duplicate pathways included in the super-pathway

        # Also a copy of the structure is created so the subpathways that are included
        # do not have duplicate reactions removed from them that are only
        # considered duplicates if they are part of a super-pathway
        self.structure=deepcopy(self.structure)
                    
        # First remove the duplicate reactions that are optional that also have key reactions
        for reaction in duplicate_optional_reactions:
            # Remove all of these optional reactions which are duplicates to a key reaction
            self.structure,new_count=self.remove_duplicate_reaction(reaction,reaction_count[reaction]+1,self.structure)
            # Remove this reaction from the set of counts
            del reaction_count[reaction]
        
        # remove the duplicate reactions
        for reaction, count in reaction_count.items():
            if count > 1:
                self.structure,new_count=self.remove_duplicate_reaction(reaction,count,self.structure)


def deepcopy(structure):
    """ 
    Create a deepcopy of the structure. Each object referenced is created and included.
    This differs from copy.deepcopy in that there will not be duplicate references included in the copy returned.
    This is important for this structure as it could have duplicate references from duplicate subpathways.
    The duplicate references need to be resolved before the duplicate reactions are removed.
    This will prevent reaction removals from being propagated throughout duplicated references.
    """
        
    # Check if this is a list or string
    if isinstance(structure, basestring):
        return structure
    else:
        new_structure=[]
        for i in xrange(len(structure)):
            if isinstance(structure[i],list):
                new_structure.append(deepcopy(structure[i]))
            elif isinstance(structure[i],basestring):
                new_structure.append(structure[i])
                
    return new_structure
    

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
    
    def get_predecessors(self):
        return copy.copy(self.predecessors)
        
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
                
    def create_structure(self,reactions,prior_nodes,structure_by_level=None,current_level=None):
        """
        Get the structure for the node and leaves
        """
        
        if structure_by_level is None:
            structure_by_level={}
        
        if current_level is None:
            current_level=0
            structure_by_level[current_level]=[]
        
        # search to see if self has already been visited in prior nodes
        if not node_is_in_list(self, prior_nodes.keys()):
            # check if this node is from a contraction
            # check if the contraction is directly into an expansion
            if self.count_predecessors() > 1 and structure_by_level[current_level]:
                # if this is from a contraction then pop up one level before adding
                # if this is already at the starting level, then use the starting level
                if current_level > 0:
                    current_level-=1    
                    structure_by_level[current_level].append(self.name)
                    prior_nodes[self]=structure_by_level[current_level]
                else:
                    # create two lists
                    new_list=[structure_by_level[0],[self.name]]
                    structure_by_level[0]=new_list
                    current_level+=1
                    structure_by_level[current_level]=structure_by_level[0][1]
                    prior_nodes[self]=structure_by_level[0][1]
            else:
                structure_by_level[current_level].append(self.name)
                prior_nodes[self]=structure_by_level[current_level]
            
        # remove recursive nodes from leaves if present
        non_recursive_leaves=[]
        for node in self.leaves:
            if not node_is_in_list(node, prior_nodes.keys()):
                non_recursive_leaves.append(node)
            
        # find the structure for the leaves
        insert_node=None
        join_type=None
        if len(non_recursive_leaves) == 1:
            # Append this on to the same level as the predecessor
            leaf_structure, insert_node, join_type = non_recursive_leaves[0].create_structure(reactions,prior_nodes,structure_by_level,current_level)
        elif len(non_recursive_leaves) > 1:
            # This is an OR expansion
            # Add a set of lists each with a higher level value
            OR_list=[OR_DELIMITER]
            structure_by_level[current_level].append(OR_list)
            leaf_level=current_level+1
            for node in non_recursive_leaves:
                # Add a new list into the or list
                new_list=[]
                OR_list.append(new_list)
                structure_by_level[leaf_level]=new_list
                leaf_structure, insert_node, join_type = node.create_structure(reactions,prior_nodes,structure_by_level,leaf_level)
            
        # check if this is an insert and if there are not any other leaves remaining
        if len(non_recursive_leaves) == 0 and self.count_leaves() > 0:
            # return the node of the other predecessor to this leaf
            other_predecessors=[]
            for leaf in self.get_leaves():
                if insert_node:
                    break
                other_predecessors=leaf.get_predecessors()
                if len(other_predecessors) > 1:
                    for other_predecessor in other_predecessors:
                        if not other_predecessor is self and other_predecessor in prior_nodes:
                            insert_node=other_predecessor
                            break
                
            if insert_node:
                # determine the join type
                # if the insert node and this node have the same products and reactants then this is OR
                left,right=reactions.get(self.get_name(),[set(),set()])
                insert_node_left,insert_node_right=reactions.get(insert_node.get_name(),[set(),set()])
                
                join_type=AND_DELIMITER
                if left and left == insert_node_left and right and right == insert_node_right:
                    join_type=OR_DELIMITER

        return structure_by_level[0], insert_node, join_type
        
def write_structured_pathways(metacyc_pathways, output_file):
    """
    Write the structure for each of the pathways
    """
    
    try:
        file_handle=open(output_file,"w")
    except EnvironmentError:
        sys.exit("Unable to write to file: "+output_file)
        
    # print the structures
    for pathway in metacyc_pathways:
        structure=metacyc_pathways[pathway].get_string_structure()
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
    
    # create the initial structures for all of the pathways
    if args.verbose:
        print("Creating structures for pathways")
    for pathway in metacyc_pathways:
        metacyc_pathways[pathway].create_structure()
        
    # resolve any subpathways in the structures
    if args.verbose:
        print("Resolving subpathways in structures")
    for pathway in metacyc_pathways:
        metacyc_pathways[pathway].resolve_subpathways(metacyc_pathways)
        
    # resolve duplicate reactions in the structures
    if args.verbose:
        print("Resolving duplicate reactions in pathways")
    for pathway in metacyc_pathways:
        metacyc_pathways[pathway].resolve_duplicate_reactions()
    
    write_structured_pathways(metacyc_pathways, output_file)
    
    total_pathways=len(metacyc_pathways)
    if args.verbose:
        print("Total pathways identified from metacyc: " +str(total_pathways))
         

if __name__ == "__main__":
    main()
