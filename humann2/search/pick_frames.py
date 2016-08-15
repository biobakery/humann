"""
HUMAnN2: pick_frames module
Screen nucleotide reads for protein-like reading frames

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
    """ convert a dna strand to its reverse complement """
    return "".join( [switch.get(nuc,"N") for nuc in dna[::-1]] )
    
def translate ( dna, frame=0 ):
    """ translate a dna sequence in the desired frame (0,1,2) """
    peptide = ""
    i = frame
    while i + 3 <= len( dna ):
        peptide += decode.get( dna[i:i+3], bad_aa_char )
        i += 3
    return peptide

def pick_frames ( sequence ):
    """ identify +/- frames with no bad chars (including "*"=stop) """
    sequence = sequence.upper()
    sequence_rc = reverse_complement( sequence )
    valid_peptides = []
    # forward translations
    for i in range( 3 ):
        p = translate( sequence, frame=i )
        if "*" not in p and bad_aa_char not in p:
            valid_peptides.append( p )
    # reverse translations
    for i in range( 3 ):
        p = translate( sequence_rc, frame=i )
        if "*" not in p and bad_aa_char not in p:
            valid_peptides.append( p )
    # return / output
    return valid_peptides

def write ( header, sequence ):
    """ find a print valid frames from sequence in fasta format """
    valid_peptides = pick_frames( sequence )
    for p in valid_peptides:
        print(header)
        print(p)

# ---------------------------------------------------------------
# main (script mode)
# ---------------------------------------------------------------

def main ( ):
    """ called to stream dna -> protein in script mode """
    header, sequence = None, ""
    for line in sys.stdin:
        line = line.strip()
        if line[0] == ">":
            # new record, process active record and reset
            if header is not None:
                write( header, sequence )
            header, sequence = line, ""
        else:
            sequence += line
    # process last record
    write( header, sequence )

if __name__ == "__main__":
    main()
