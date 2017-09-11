#!/usr/bin/env python

"""
===============================================
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import os
import sys
import csv

# ---------------------------------------------------------------
# import checks
# ---------------------------------------------------------------

# humann2 checks
try:
    from humann2 import config
    from humann2.tools import util
except ImportError:
    sys.exit("CRITICAL ERROR: Unable to find the HUMAnN2 python package." +
        " Please check your install.")

# numpy checks
try:
    import numpy as np
except:
    sys.exit( "This module requires the Python scientific stack: numpy, scipy, and matplotlib." )

# ---------------------------------------------------------------
# main class
# ---------------------------------------------------------------

class Table:

    def __init__( self, source=None, headers=None, metadata=None ):

        self.data = {}
        self.anchor = "# SAMPLES"
        self.headers = []
        self.metadata = {}
        self.source = ""
        self.n = None
        self.verbose = True

        if type( source ) is dict:
            self.data = source
            self.headers = headers
            self.metadata = metadata
            self.verbose = False
        else:
            if source is None:
                self.source = "STDIN"
                handle = sys.stdin
            else:
                self.source = source
                handle = util.try_zip_open( source )
            print( "Reading table data from <{}>".format( self.source ), file=sys.stderr )
            IN_METADATA = False if metadata is None else True
            for row in csv.reader( handle, csv.excel_tab ):
                if self.headers == []:
                    self.anchor = row[0]
                    self.headers = row[1:]
                elif IN_METADATA:
                    self.metadata[row[0]] == row[1:]
                    if row[0] == metadata:
                        IN_METADATA = False
                else:
                    self.data[row[0]] = np.array( row[1:], dtype=np.float )
        
        self.check_lengths( )
        self.n = len( self.headers )
        if self.source is not None:
            self.report( )

    def zeros( self ):
        return np.zeros( self.n )

    def report( self ):
        if self.verbose:
            print( "  # of samples: {}".format( len( self.headers ) ), file=sys.stderr )
            print( "  # of metadata rows: {}".format( len( self.metadata ) ), file=sys.stderr )
            print( "  # of feature rows: {}".format( len( self.data ) ), file=sys.stderr )
            totals = {util.fsplit( f )[0] for f in util.fsort( self.data )}
            print( "  # of feature totals: {}".format( len( totals ) ), file=sys.stderr )

    def check_lengths( self ):
        row_lens = {len( row ) for name, row in self.data.items( )}
        if len( row_lens ) == 0:
            sys.exit( "No data rows loaded." )
        elif len( row_lens ) > 1:
            sys.exit( "Data rows have unequal lengths." )
        elif row_lens != set( [len( self.headers )] ):
            print( row_lens, set( [len( self.headers )] ) ) 
            sys.exit( "Data row length does not match number of headers." )
        elif self.metadata != {}:
            meta_lens = set( len( row ) for name, row in self.metadata.items( ) )
            if meta_lens != row_lens:
                sys.exit( "Metadata row lengths not consistent with data row lengths." )

    def iter_rows( self, unfloat=False ):
        # headers
        yield [self.anchor] + self.headers
        # metadata if any
        for name in sorted( self.metadata ):
            yield [name] + metadata[name]
        # feature rows (data)
        for name in util.fsort( self.data ):
            values = list( self.data[name] )
            if unfloat:
                values = ["0" if x == 0 else "%.6g" % ( x ) for x in values]
            yield [name] + values

    def write( self, path=None, unfloat=False ):
        if path is None:
            path = "STDOUT"
            fh = sys.stdout
        else:
            fh = util.try_zip_open( path, "w" )
        writer = csv.writer( fh, csv.excel_tab )
        for row in self.iter_rows( unfloat=unfloat ):
            writer.writerow( row )
        print( "Wrote table data to <{}>".format( path ), file=sys.stderr )
        self.verbose = True
        self.report( )
