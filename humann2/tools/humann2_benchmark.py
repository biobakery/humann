#!/usr/bin/env python

""" 
This script will capture the MaxRSS for a process and its children (and their children).
It transverses the pid tree to capture all processes spawned from the command. 

To run:
$ ./humann2_benchmark.py $COMMAND

Replacing $COMMAND with the command you would like to benchmark.

The MaxRSS of the processes are printed out in between the stdout from $COMMAND.

The final output has the system time in minutes and the MaxRSS sum of all processes in GB.
"""

import sys
import subprocess
import time

def process_ps_stdout(stdout):
    """ Process the stdout of the ps command """
    return [i.split()[0] for i in filter(lambda x: x, stdout.decode("utf-8").split("\n")[1:])]

def get_children(pid):
    """ Get the pids of the children of this pid """
    try:
        stdout=subprocess.check_output(["ps","--ppid",pid,"-o","pid"])
    except subprocess.CalledProcessError:
        stdout=[]

    pids=[]
    if stdout:
        pids=process_ps_stdout(stdout)

    return pids

def traverse_tree(pid,nodes):
    """ Get all of the children of the pid by walking the tree """

    for child in get_children(pid):
        nodes.update(traverse_tree(child,nodes))
    nodes.add(pid)

    return nodes

def get_pids(pid):
    """ Get the children and their children from parent pid """

    pids=set([pid])
    for child in get_children(pid):
        pids.update(traverse_tree(child,pids))
    
    return list(pids)

def main():
    """ Capture time and memory from command run """
    process = subprocess.Popen(sys.argv[1:],shell=False)
    pid=str(process.pid)
    start=time.time()
    max_memory=0
    while process.poll() is None:
        time.sleep(1)
        # while the process is running check on the memory use
        # get the pids of the main process and all children (and their children)
        pids=get_pids(pid)
        stdout=subprocess.check_output(["ps","--pid",",".join(pids),"-o","pid,rss,command"]).decode("utf-8")
        print("\n"+stdout+"\n")
        # remove the header from the process output
        status=[i.split() for i in filter(lambda x: x, stdout.split("\n")[1:])]
        # memory is the sum of all rss
        memory=sum(int(i[1]) for i in status)
        if memory > max_memory:
            max_memory=memory
    
    end=time.time()
    print("Time: {:.0f} minutes".format((end-start)/60))
    print("Max Memory (RSS): {:.1f} GB".format(max_memory*1.0/1024**2))

if __name__ == "__main__":
    main()

