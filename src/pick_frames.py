#! /usr/bin/env python

"""
Standalone pipe for converting DNA reads in FASTA
format to FASTA formatted peptides.

Examines all 6 reading frames and prints those that
(1) do not contain a stop codon and 
(2) do not contain a non-amino acid character (here, "X")

Usage: cat reads.fasta | python <this>
-----------------------------------------------
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

import sys

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

# used when we find an untranslatable codon
bad_aa_char = "X"

# for making reverse complement
switch = {"A":"T", "T":"A", "G":"C", "C":"G"}

# for translating dna to protein
genetic_code = """
Ala/A GCT GCC GCA GCG
Arg/R CGT CGC CGA CGG AGA AGG
Asn/N AAT AAC
Asp/D GAT GAC
Cys/C TGT TGC
Gln/Q CAA CAG
Glu/E GAA GAG
Gly/G GGT GGC GGA GGG
His/H CAT CAC
Ile/I ATT ATC ATA
Leu/L TTA TTG CTT CTC CTA CTG
Lys/K AAA AAG
Met/M ATG
Phe/F TTT TTC
Pro/P CCT CCC CCA CCG
Ser/S TCT TCC TCA TCG AGT AGC
Thr/T ACT ACC ACA ACG
Trp/W TGG
Tyr/Y TAT TAC
Val/V GTT GTC GTA GTG
Stp/* TAA TGA TAG
"""

# convert genetic code to dict
decode = {}
for line in genetic_code.split( "\n" ):
    items = line.split( )
    if len( items ) > 0:
        long, short = items[0].split( "/" )
        for codon in items[1:]:
            decode[codon] = short

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------
        
def reverse_complement ( dna ):
    return "".join( [switch[nuc] for nuc in dna[::-1]] )
    
def translate ( dna, frame=0 ):
    peptide = ""
    i = frame
    while i + 3 <= len( dna ):
        peptide += decode.get( dna[i:i+3], bad_aa_char )
        i += 3
    return peptide

def pick_frames ( header, sequence ):
    sequence_rc = reverse_complement( sequence )
    peptides = []
    # forward translations
    for i in range( 3 ):
        peptides.append( translate( sequence, frame=i ) )
    # reverse translations
    for i in range( 3 ):
        peptides.append( translate( sequence_rc, frame=i ) )
    for peptide in peptides:
        if not ( "*" in peptide or bad_aa_char in peptide ):
            print header
            print peptide

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main ( ):
    header, sequence = None, ""
    for line in sys.stdin:
        line = line.strip()
        if line[0] == ">":
            # new record, process active record and reset
            if header is not None:
                pick_frames( header, sequence )
            header, sequence = line, ""
        else:
            sequence += line
    # process last record
    pick_frames( header, sequence )

if __name__ == "__main__":
    main()
