#! /usr/bin/env python

"""
HUMAnN2 utility for renaming and renormalizing TSV files
Run ./rename_renorm.py -h for usage help
"""

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import csv
import re
import gzip
import argparse
import sys

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_strat_delim = "|"
c_name_delim = ": "
c_multiname_delim = ";"
c_str_unknown = "NO_NAME"

# ---------------------------------------------------------------
# classes
# ---------------------------------------------------------------

class CTable ( ):

    """ Very basic table class; would be more efficient using numpy 2D array """

    def __init__ ( self, path ):
        self.anchor = None
        self.colheads = []
        self.rowheads = []
        self.data = []
        self.is_stratified = False
        with try_zip_open( path ) as fh:
            for row in csv.reader( fh, dialect='excel-tab' ):
                if self.anchor is None:
                    self.anchor = row[0]
                    self.colheads = row[1:]
                else:
                    self.rowheads.append( row[0] )
                    self.data.append( row[1:] )
        for rowhead in self.rowheads:
            if c_strat_delim in rowhead:
                print( "Treating", path, "as stratified output, e.g.", 
                       rowhead.split( c_strat_delim ), file=sys.stderr )
                self.is_stratified = True
                break

    def normalize ( self, cpm=False ):
        divisor = 1 
        if self.is_stratified:
            divisor /= 2.0
        if cpm:
            divisor /= 1e6
        totals = [0 for k in range( len( self.colheads ) )]
        for i, row in enumerate( self.data ):
            self.data[i] = [float( k ) for k in row]
            totals = [k1 + k2 for k1, k2 in zip( totals, self.data[i] )]
        for i, row in enumerate( self.data ):
            self.data[i] = ["%.6g" % ( row[j] / totals[j] / divisor ) for j in range( len( totals ) )] 

    def rename ( self, renaming_dict ):
        seen = {}
        for i, rowhead in enumerate( self.rowheads ):
            items = rowhead.split( c_strat_delim )
            old_name = items[0]
            seen[old_name] = False
            new_name = renaming_dict.get( old_name, c_str_unknown )
            if new_name != c_str_unknown:
                seen[old_name] = True
            items[0] = c_name_delim.join( [old_name, new_name] )
            self.rowheads[i] = c_strat_delim.join( items )
        tcount = seen.values().count( True )
        print( "Renamed %d of %d entries (%.2f%%)" \
                   % ( tcount, len( seen ), 100 * tcount / float( len( seen ) ) ), 
               file=sys.stderr )

    def write ( self, fh ):
        writer = csv.writer( fh, dialect='excel-tab' )
        writer.writerow( [self.anchor] + self.colheads )
        for i in range( len( self.rowheads ) ):
            writer.writerow( [self.rowheads[i]] + self.data[i] )

# ---------------------------------------------------------------
# utility functions 
# ---------------------------------------------------------------

def get_args ():
    """ Get args from Argparse """
    parser = argparse.ArgumentParser()
    parser.add_argument( 
        "-i", "--input", 
        help="Original HUMAnN2 output table (.tsv format)",
        )
    parser.add_argument( 
        "-n", "--names", 
        help="Mapping of gene/pathway IDs to full names (.tsv or .tsv.gz)",
        )
    parser.add_argument( 
        "-r", "--norm", 
        choices=["off", "cpm", "relab"], 
        default="off", 
        help="Normalization scheme: [off], copies per million [cpm], relative abundance [relab]; default=[off]",
        )
    parser.add_argument( 
        "-o", "--output", 
        default=None,
        help="Path for modified output table; default=[STDOUT]",
        )
    args = parser.parse_args()
    return args

def try_zip_open( path ):
    """ open an uncompressed or gzipped file; fail gracefully """
    fh = None
    try:
        fh = open( path ) if not re.search( r".gz$", path ) else gzip.GzipFile( path )
    except:
        print( "Problem loading", path, file=sys.stderr )
    return fh

def load_renaming ( path, table ):
    """ load a tsv file mapping one name to another (e.g. uniref50 id to english name) """
    # get possible "old names" from table (avoids storing full mapping)
    table_names = {}
    for rowhead in table.rowheads:
        table_names[rowhead.split( c_strat_delim )[0]] = 1
    # load renaming dict
    renaming_dict = {}
    with try_zip_open( path ) as fh:
        for row in csv.reader( fh, dialect="excel-tab" ):
            old_name = row[0]
            if old_name in table_names:
                for new_name in row[1:]:
                    renaming_dict.setdefault( old_name, [] ).append( new_name.strip() )
    for old_name, new_names in renaming_dict.items():
        renaming_dict[old_name] = c_multiname_delim.join( new_names )
    print( "Loaded renaming dict from", path, file=sys.stderr )
    return renaming_dict

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main ( ):
    args = get_args()
    table = CTable( args.input )
    # renaming
    if args.names is not None:
        table.rename( load_renaming( args.names, table ) )
    # normalizing
    if args.norm != "off":
        table.normalize( cpm=True if args.norm == "cpm" else False )
    # output
    fh = open( args.output, "w" ) if args.output is not None else sys.stdout
    table.write( fh )

if __name__ == "__main__":
    main()
