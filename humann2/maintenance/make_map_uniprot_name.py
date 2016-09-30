#! /usr/bin/env python

"""
Extract information from UniRef50/90 fasta
to make a mapping from UniRef ID to english name
Exclude mappings to "Uncharacterized Protein" to save space
====
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

import os, sys, re, argparse, collections

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

# uninformative names to skip
c_disallow = """
MULTISPECIES: hypothetical protein
hypothetical protein, partial
Putative uncharacterized protein (Fragment)
Predicted protein
Uncharacterized protein (Fragment)
Putative uncharacterized protein
hypothetical protein
Uncharacterized protein
"""
c_disallow = {n:1 for n in c_disallow.split( "\n" ) if n != ""}

# characters not allowed to appear in names (used elsewhere in humann2)
c_bad_strings = ["|", ": ", ";"]

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument( 'fasta', help='uniprot/uniref fasta file' )
args = parser.parse_args()

outfile = "map_{}_name.txt".format( os.path.split( args.fasta )[1].split( "." )[0] )
print >>sys.stderr, "writing to:", outfile

report = collections.Counter()
with open( args.fasta ) as fhin, open( outfile, "w" ) as fhout:
    for line in fhin:
        if line[0] == ">":
            match = re.search( ">(.*?) (.*) n=", line )
            if match:
                code, name = match.groups()
                if name not in c_disallow:
                    for s in c_bad_strings:
                        name = name.replace( s, "_" )
                    print >>fhout, "\t".join( [code, name] )
                    report["Saved"] += 1
                else:
                    report[name] += 1                   
            else:
                print >>sys.stderr, "BAD UNI-FASTA LINE:", line.strip()
                report["Unparsed"] += 1

for k, v in report.items():
    print >>sys.stderr, v, "\t", k, "<- skipped" if k in c_disallow else ""
