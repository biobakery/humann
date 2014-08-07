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

    # find hits for the bam output
    if bam_output != "Empty":
        # call humann1 bam2hits.py
        exe="bam2hits.py"
    
        print "Running " + exe + " ........"
        utilities.execute_command(exe, [], [bam_output],
            [], hits_file, bam_output)
    

    # find this for the output in blast format
    exe="blast2hits.py"
    
    return hits_file
