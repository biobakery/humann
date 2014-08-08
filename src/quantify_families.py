"""
Identify hits from alignment steps
Compute gene family abundance
"""

import os, tempfile
import config, utilities

def find_hits(bam_output, blast_output):
    """
    Identify the hits from the alignment results
    """

    file_handle, hits_file=tempfile.mkstemp()
    os.close(file_handle)

    print "Running to find hits ........"

    # find hits for the bam output
    if bam_output != "Empty":
        # call humann1 bam2hits.py
        utilities.execute_command("bam2hits.py", [], [bam_output],
            [], hits_file, bam_output)
    
        # call humann1 blast2hits.py with results
        utilities.execute_command("blast2hits.py",
            ["blastx",0,hits_file],[blast_output],[],
                hits_file, blast_output)
    else:
        # call humann1 blast2hits.py
        utilities.execute_command("blast2hits.py",
            [],[blast_output],[],hits_file, blast_output)
        
    return hits_file
