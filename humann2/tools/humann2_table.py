#!/usr/bin/env python

"""
Updated HUMAnN2 table object
===============================================
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

from __future__ import print_function # PYTHON 2.7+ REQUIRED
import os
import sys
import csv
from collections import OrderedDict

# ---------------------------------------------------------------
# import checks
# ---------------------------------------------------------------

try:
    from humann2 import config
    from humann2.tools import util
except ImportError:
    sys.exit("CRITICAL ERROR: Unable to find the HUMAnN2 python package." +
        " Please check your install.")

try:
    import numpy as np
except:
    sys.exit( "This module requires the Python scientific stack: numpy, scipy, and matplotlib." )

# ---------------------------------------------------------------
# the table class
# ---------------------------------------------------------------

class Table:

    def __init__( self, source=None, headers=None, metadata=None, 
                  last_metadata=None ):

        # the col1, row1 entry
        self.anchor = "# SAMPLES"
        # list of the sample names (excludes anchor)
        self.headers = []
        # metadata rows (dict: header->list)
        self.metadata = OrderedDict( )
        # data rows (dict: feature->array)
        self.data = {}
        # properties of the table
        self.n = None
        self.source = ""
        self.verbose = True

        if type( source ) is dict:
            self.data = source
            self.headers = headers
            self.metadata = self.metadata if metadata is None else metadata
            self.verbose = False
        else:
            if source is None:
                self.source = "STDIN"
                handle = sys.stdin
            else:
                self.source = source
                handle = util.try_zip_open( source )
            print( "Reading table data from <{}>".format( self.source ), file=sys.stderr )
            IN_METADATA = False if last_metadata is None else True
            for row in csv.reader( handle, csv.excel_tab ):
                if self.headers == []:
                    self.anchor = row[0]
                    self.headers = row[1:]
                elif IN_METADATA:
                    self.metadata[row[0]] = row[1:]
                    if row[0] == last_metadata:
                        IN_METADATA = False
                else:
                    try:
                        self.data[row[0]] = np.array( row[1:], dtype=np.float )
                    except:
                        sys.exit( ("Died while trying to turn row <{}> into numbers;\n"
                                  "Use --last-metadata <row> to specify end of metadata.".format( row[0] )) )

        self.n = len( self.headers )
        self.check_lengths( )
        if self.source is not None:
            self.report( )

    def zeros( self ):
        # update in case table size changed?
        self.n = len( self.headers )
        return np.zeros( self.n )

    def report( self ):
        if self.verbose:
            print( "  # of samples: {:,}".format( len( self.headers ) ), file=sys.stderr )
            if len( self.metadata ) > 0:
                print( "  # of metadata rows: {:,}".format( len( self.metadata ) ), file=sys.stderr )
            totals = {util.fsplit( f )[0] for f in util.fsort( self.data )}
            print( "  # of feature totals: {:,}".format( len( totals ) ), file=sys.stderr )
            print( "  # of feature rows: {:,}".format( len( self.data ) ), file=sys.stderr )

    def check_lengths( self ):
        row_lens = {len( row ) for name, row in self.data.items( )}
        if len( row_lens ) == 0:
            # try this as a warning and not a failure
            print( "WARNING: No data rows loaded.", file=sys.stderr )
        elif len( row_lens ) > 1:
            sys.exit( "CRITICAL ERROR: Data rows have unequal lengths." )
        elif row_lens != set( [len( self.headers )] ):
            print( row_lens, set( [len( self.headers )] ) ) 
            sys.exit( "CRITICAL ERROR: Data row length does not match number of headers." )
        elif len( self.metadata ) > 0:
            meta_lens = set( len( row ) for name, row in self.metadata.items( ) )
            if meta_lens != row_lens:
                sys.exit( "CRITICAL ERROR: Metadata row lengths not consistent with data row lengths." )

    def resample( self, samples ):
        """ slice columns from table as new table """
        samples = set( samples )
        original = set( self.headers )
        for s in samples:
            if s not in original:
                print( "WARNING: Can't find sample <{}>".format( s ), file=sys.stderr )
        order = []
        for i, h in enumerate( self.headers ):
            if h in samples:
                order.append( i )
        new_headers = [self.headers[i] for i in order]
        new_data = {}
        for f, row in self.data.items( ):
            new_data[f] = row[order]
        new_metadata = OrderedDict( )
        for f, row in self.metadata.items( ):
            new_metadata[f] = [row[i] for i in order]
        return Table( new_data, headers=new_headers, metadata=new_metadata )

    def iter_rows( self, unfloat=False ):
        # headers
        yield [self.anchor] + self.headers
        # metadata if any
        for name, metarow in self.metadata.items( ):
            yield [name] + metarow
        # feature rows (data)
        for name in util.fsort( self.data ):
            values = list( self.data[name] )
            if unfloat:
                values = ["0" if x == 0 else "%.6g" % ( x ) for x in values]
            yield [name] + values

    def write( self, path=None, unfloat=False, report=True ):
        if path is None:
            path = "STDOUT"
            fh = sys.stdout
        else:
            fh = util.try_zip_open( path, "w" )
        writer = csv.writer( fh, csv.excel_tab )
        for row in self.iter_rows( unfloat=unfloat ):
            writer.writerow( row )
        print( "Wrote table data to <{}>".format( path ), file=sys.stderr )
        if report:
            self.verbose = True
            self.report( )
