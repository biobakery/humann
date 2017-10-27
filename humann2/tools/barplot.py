#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import csv
import argparse
import math

from collections import Counter

try:
    from humann2 import config
    from humann2.tools import util
    from humann2.tools.humann2_table import Table
    from humann2.tools.legendizer import Legendizer
except ImportError:
    sys.exit( "CRITICAL ERROR: Unable to find the HUMAnN2 python package." +
              " Please check your install." )

try:
    import matplotlib
    matplotlib.use( "Agg" )
    matplotlib.rcParams["pdf.fonttype"] = 42
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import numpy as np
    import scipy.cluster.hierarchy as sch
except ImportError:
    sys.exit( "This script requires the Python scientific stack: numpy, scipy, and matplotlib." )   

# ---------------------------------------------------------------
# global constants
# ---------------------------------------------------------------

c_hunits = 9
c_other = "Other"

# ---------------------------------------------------------------
# argument parsing
# ---------------------------------------------------------------

c_sort_help = """Sample sorting methods (can use more than one; will evaluate in order)

none        : Default
sum         : Sum of stratified values
dominant    : Value of the most dominant stratification
similarity  : Bray-Curtis agreement of relative stratifications
usimilarity : Bray-Curtis agreement of raw stratifications
metadata    : Given metadata label

"""

c_choices_help = """Scaling options for total bar heights (strata are always proportional to height)

none        : Default
pseudolog   : Total bar heights log10-scaled (strata are NOT log10-scaled)
pseudosqrt  : Total bar heights sqrt-scaled (strata are NOT sqrt-scaled)
normalize   : Bars all have height=1 (highlighting relative attribution)

"""

# ---------------------------------------------------------------
# command-line interface
# ---------------------------------------------------------------

description = util.wrap( """
HUMAnN2 utility for plotting a single stratified feature

Plots the stratified contributions of a specified function. Can optionally
sort samples to reveal ecological- or metadata-linked trends. Can perform custom
scaling to highlight stratifications even when community totals have high dynamic range.
""" )

def get_args( ):
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    util.attach_common_arguments( parser, no_output=True )
    parser.add_argument( "-f", "--focal-feature",
                         metavar = "<feature>",
                         help="Feature ID of interest (give ID not full name)", )
    parser.add_argument( "-t", "--top-strata",
                         metavar = "<int>",
                         type=int,
                         default=7,
                         help="Number of top stratifications to highlight (top = highest grand means)\n[Default=7]", )
    parser.add_argument( "-s", "--sort",
                         metavar = "<sorting methods>",
                         nargs="+",
                         default=["none"],
                         choices=["none", "sum", "dominant", "similarity", "usimilarity", "metadata"],
                         help=c_sort_help, )             
    parser.add_argument( "-m", "--focal-metadata",
                         metavar="<feature>",
                         default=None,
                         help="Indicate metadata to highlight / group by", )
    parser.add_argument( "-l", "--max-metalevels",
                         metavar="<int>",
                         type=int,
                         default=7,
                         help="Color the N most populous metadata levels; collapse others as 'Other'\n[Default=7]" )
    parser.add_argument( "-c", "--colormap",
                         metavar="<colormap>",
                         default="jet",
                         help="Color space for stratifications", )
    parser.add_argument( "-k", "--meta-colormap",
                         metavar="<colormap>",
                         default="Dark2",
                         help="Color space for metadata levels", )
    parser.add_argument( "-x", "--exclude-unclassified",
                         action="store_true",
                         help="Do not include the 'unclassified' stratum", )
    parser.add_argument( "-o", "--output",
                         metavar="<file.ext>",
                         default=None,
                         help="Where to save the figure", )
    parser.add_argument( "-a", "--scaling",
                         metavar="<choice>",
                         choices=["none", "normalize", "pseudolog", "pseudosqrt"],
                         default="none",
                         help=c_choices_help, )
    parser.add_argument( "-g", "--as-genera",
                         action="store_true",
                         help="Collapse species to genera", )
    parser.add_argument( "-r", "--grid",
                         action="store_true",
                         help="Add y-axis grid", )
    parser.add_argument( "-z", "--remove-zeroes",
                         action="store_true",
                         help="Do not plot samples with zero sum for this feature", )
    parser.add_argument( "-w", "--width",
                         metavar = "<int>",
                         default = 3,
                         type=int,
                         help="Relative width of the plot vs. legend (default: 3)", )
    parser.add_argument( "-d", "--dimensions",
                         metavar = "<size>",
                         nargs=2,
                         type=float,
                         default=[10.0,4.0],
                         help="Image height and width in inches (default: 8 4)", )
    parser.add_argument( "-y", "--ylims",
                         metavar = "<limit>",
                         nargs=2,
                         default=[None, None],
                         help="Fix limits for y-axis", )
    parser.add_argument( "-e", "--legend-stretch",
                         metavar = "<float>",
                         type=float,
                         default=1.0,
                         help="Stretch/compress legend elements", )
    return parser.parse_args( )

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

def subseq( seq, index ):
    """numpy-style slicing and indexing for lists"""
    return [seq[i] for i in index]

def bugname( rowhead ):
    """make bug names look nicer"""
    if "s__" in rowhead:
        return rowhead.split( "." )[1].replace( "s__", "" ).replace( "_", " " )
    elif "g__" in rowhead:
        return rowhead.replace( "g__", "" ).replace( "_", " " )
    else:
        return rowhead

def get_colors( colormap, n ):
    """utility for defining N evenly spaced colors across a color map"""
    cmap = plt.get_cmap( colormap )
    cmap_max = cmap.N
    if n == 1:
        return [cmap( int( 0.5 * cmap_max ) )]
    else:
        return [cmap( int( k * cmap_max / (n - 1) ) ) for k in range( n )]

def dummy( ax, border=False ):
    """dummy an axis"""
    ax.set_xticks( [] )
    ax.set_yticks( [] )
    if not border:
        for k, v in ax.spines.items():
            v.set_visible( False )

def tsv_reader( path ):
    """quick read"""
    with open( path ) as fh:
        for row in csv.reader( fh, csv.excel_tab ):
            yield row

def simplify_metadata( values, n=None ):
    counts = Counter( values )
    if len( counts ) > n:
        index = [[-c, k] for k, c in counts.items( )]
        index.sort( )
        allowed = {x[1] for x in index[0:n]}
        values = [k if k in allowed else c_other for k in values]
    return values
            
# ---------------------------------------------------------------
# Single-feature table class with numerical methods
# ---------------------------------------------------------------
    
class FeatureTable:

    def __init__( self, Table, focus, metaname=None, exclude_unclassified=False ):

        self.colheads = Table.headers[:]
        self.rowheads = []
        self.data = []
        self.metaname = metaname
        self.metarow = None
        self.title = None

        # isolate optional metadata row for plotting
        if metaname is not None:
            self.metarow = Table.metadata.get( metaname, None )

        # pull out relevant feature rows
        for f in util.fsort( Table.data ):
            fcode, fname, fstratum = util.fsplit( f )
            if fcode == focus and fstratum is not None:
                self.title = util.fjoin( fcode, fname )
                if fstratum != util.c_unclassified or not exclude_unclassified:
                    self.rowheads.append( fstratum )
                    self.data.append( Table.data[f] )
        self.data = np.array( self.data )
        self.update( )
        
    def update( self ):
        self.nrows, self.ncols = self.data.shape
        self.colsums = sum( self.data )
        self.colsums = np.array( [k if k > 0 else -1.0 for k in self.colsums] )
        self.rowmap = {}
        self.colmap = {}
        for i, h in zip( range( self.nrows ), self.rowheads ):
            self.rowmap[h] = i
        for i, h in zip( range( self.ncols ), self.colheads ):
            self.colmap[h] = i
        assert self.nrows == len( self.rowheads ) == len( self.rowmap ), "row dim failure"
        assert self.ncols == len( self.colheads ) == len( self.colmap ), "col dim failure"

    def as_genera( self ):
        print( "Regrouping to genera (before selecting/sorting strata)", file=sys.stderr )
        temp = {}
        for rowhead, row in zip( self.rowheads, self.data ):
            rowhead = rowhead.split( "." )[0]
            temp[rowhead] = temp.get( rowhead, np.zeros( self.ncols ) ) + row
        self.rowheads = list( temp.keys( ) )
        self.data = np.array( [temp[k] for k in self.rowheads] )
        self.update( )

    def remove_zeroes( self ):
        order = [i for i, k in enumerate( self.colsums ) if k > 0]
        self.colheads = subseq( self.colheads, order )
        self.data = self.data[:, order]
        if self.metarow is not None:
            self.metarow = subseq( self.metarow, order )
        self.update( )
            
    def sort( self, method=None ):
        if method == "none":
            order = range( self.ncols )
        elif method == "similarity":
            norm = self.data / ( self.colsums + util.c_eps * np.ones( self.ncols ) )
            # note: linkage assumes things to cluster = rows; we want cols
            order = sch.leaves_list( sch.linkage( norm.transpose(), metric="braycurtis" ) )
        elif method == "usimilarity":
            # note: linkage assumes things to cluster = rows; we want cols
            order = sch.leaves_list( sch.linkage( self.data.transpose(), metric="braycurtis" ) )
        elif method == "dominant":
            maxes = [(max( self.data[i, :] ), i) for i in range( self.nrows )]
            maxes.sort()
            ranks = [None for k in maxes]
            for i, (s, i2) in enumerate( maxes ):
                ranks[i2] = i
            argmax = []
            for c in range( self.ncols ):
                argmax.append( sorted( range( self.nrows ), key=lambda r: self.data[r][c] )[-1] )
            order = sorted( range( self.ncols ), key=lambda c: (ranks[argmax[c]], self.data[argmax[c], c]), reverse=True )
        elif method == "sum":
            order = sorted( range( self.ncols ), key=lambda c: self.colsums[c], reverse=True )
        elif method == "metadata":
            order = sorted( range( self.ncols ), key=lambda c: self.metarow[c] )
        else:
            sys.exit( "Can't sort with method: %s" % ( method ) )
        self.data = self.data[:, order]
        self.colheads = subseq( self.colheads, order )
        if self.metarow is not None:
            self.metarow = subseq( self.metarow, order )
        self.update( )
        
    def filter_top_strata( self, top_strata ):
        # save unclassified if present
        unclassified = None
        if "unclassified" in self.rowmap:
            unclassified = self.data[self.rowmap["unclassified"]]   
        # select features
        best = []
        for i, f in enumerate( self.rowheads ):
            frow = self.data[i]
            if f != "unclassified":
                best.append( [np.mean( frow ), f] )
        best.sort( reverse=True )
        top_strata = min( top_strata, len( best ) )
        select = [pair[1] for pair in best[0:top_strata]]
        # make a new row for "other"
        other = np.zeros( self.ncols )
        for i, f in enumerate( self.rowheads ):
            frow = self.data[i]
            if f != "unclassified" and f not in select:
                other = other + frow
        # subset and reorder the table
        index = [self.rowmap[k] for k in select]
        self.rowheads = subseq( self.rowheads, index )
        self.data = self.data[index, :]
        if sum( other ) > 0:
            self.rowheads += ["Other"]
            self.data = np.vstack( [self.data, other] )
        if unclassified is not None:
            self.rowheads += ["Unclassified"]
            self.data = np.vstack( [self.data, unclassified] )
        # reverse row order to simplify plotting
        self.rowheads.reverse( )
        self.data = self.data[::-1]
        self.update( )

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):

    args = get_args( )
    
    # treat "-" as none in ylims
    a, b = args.ylims
    args.ylims[0] = None if a in ["-", None] else float( a )
    args.ylims[1] = None if b in ["-", None] else float( b )

    # load table manipulation 
    table = Table( 
        args.input, 
        last_metadata=args.last_metadata,
        )
    table = FeatureTable( 
        table,
        focus=args.focal_feature,
        metaname=args.focal_metadata,
        exclude_unclassified=args.exclude_unclassified,
    )

    # collapse to genera?
    if args.as_genera:
        table.as_genera( )
    
    # remove zeroes?
    if args.remove_zeroes:
        table.remove_zeroes( )
        
    # apply one or more sorting methods
    for method in args.sort:
        table.sort( method )

    # filter/collapse features (moved to take place AFTER sorting)
    table.filter_top_strata( args.top_strata )
        
    # simplify metadata
    if table.metarow is not None:
        table.metarow = simplify_metadata( table.metarow, args.max_metalevels )

    # set up axis system
    wunits = args.width + 1
    fig = plt.figure( )
    fig.set_size_inches( *args.dimensions )
    if table.metarow is not None:
        main_ax = plt.subplot2grid( ( c_hunits, wunits ), ( 0, 0 ), rowspan=c_hunits-1, colspan=wunits-1 )
        anno_ax = plt.subplot2grid( ( c_hunits, wunits ), ( 0, wunits-1 ), rowspan=c_hunits, colspan=1 )
        meta_ax = plt.subplot2grid( ( c_hunits, wunits ), ( c_hunits-1, 0 ), rowspan=1, colspan=wunits-1 )
        dummy( meta_ax, border=True )
    else:
        main_ax = plt.subplot2grid( ( 1, wunits ), ( 0, 0 ), rowspan=1, colspan=wunits-1 )
        anno_ax = plt.subplot2grid( ( 1, wunits ), ( 0, wunits-1 ), rowspan=1, colspan=1 )
    dummy( anno_ax )
    anno_ax.set_xlim( 0, 1 )
    anno_ax.set_ylim( 0, 1 )
    
    # setup: strata colors
    cdict = {"Other":"0.5", "Unclassified":"0.8"}
    if os.path.exists( args.colormap ):
        print( "Reading strata colors from file:", args.colormap, file=sys.stderr )
        for item, color in tsv_reader( args.colormap ):
            if item not in cdict:
                cdict[item] = color
        for f in table.rowheads:
            cdict[f] = cdict.get( f, "black" )
    else:
        colors = get_colors( args.colormap, args.top_strata )
        for f in table.rowheads:
            if f not in cdict:
                cdict[f] = colors.pop( )

    # scaling options
    if args.scaling == "none":
        ylabel      = "Relative abundance"
        bottoms     = np.zeros( table.ncols )
        ymin        = 0 if args.ylims[0] is None else args.ylims[0]
        ymax        = max( sum( table.data ) ) if args.ylims[1] is None else args.ylims[1]
        main_ax.set_ylim( ymin, ymax )
    elif args.scaling == "normalize":
        ylabel      = "Relative contributions"
        table.data  = table.data / table.colsums
        bottoms     = np.zeros( table.ncols )
        main_ax.set_ylim( 0, 1 )
    elif args.scaling == "pseudolog":
        ylabel      = "log10( Relative abundance )"
        ymin        = min( [k for k in table.colsums if k > 0] ) if args.ylims[0] is None else args.ylims[0]
        floor       = math.floor( np.log10( ymin ) )
        #floor       = floor if abs( np.log10( ymin ) - floor ) >= util.c_eps else (floor - 1)
        floors      = floor * np.ones( table.ncols )
        crests      = np.array( [np.log10( k ) if k > 10**floor else floor for k in table.colsums] )
        heights     = crests - floors
        table.data  = table.data / table.colsums * heights
        ymax        = max( table.colsums ) if args.ylims[1] is None else args.ylims[1]
        ceil        = math.ceil( np.log10( ymax ) )
        bottoms     = floors
        main_ax.set_ylim( floor, ceil )
    elif args.scaling == "pseudosqrt":
        ylabel      = "sqrt( Relative abundance )"
        ymin        = 0 if args.ylims[0] is None else args.ylims[0]
        floor       = np.sqrt( ymin )
        floors      = floor * np.ones( table.ncols )
        crests      = np.array( [np.sqrt( k ) if k > floor**2 else floor for k in table.colsums] )
        heights     = crests - floors
        table.data  = table.data / table.colsums * heights
        ymax        = max( table.colsums ) if args.ylims[1] is None else args.ylims[1]
        ceil        = math.ceil( np.sqrt( ymax ) )
        bottoms     = floors
        main_ax.set_ylim( floor, ceil )
        
    # add bars
    series = []
    for i, f in enumerate( table.rowheads ):
        frow = table.data[table.rowmap[f]]
        series.append( main_ax.bar(
            range( table.ncols ),
            frow,
            width=1,
            bottom=bottoms,
            color=cdict[f],
            edgecolor="none", ) )
        bottoms += frow

    # setup: meta colors
    if table.metarow is not None:
        if os.path.exists( args.meta_colormap ):
            sys.stderr.write( "Reading meta colors from file: {}\n".format( args.meta_colormap ) )
            mcdict = {"Other":"0.5"}
            for item, color in tsv_reader( args.meta_colormap ):
                mcdict[item] = color
            for m in table.metarow:
                mcdict[m] = mcdict.get( m, "black" )
        else:
            unique = sorted( set( table.metarow ) )
            mcdict = {v:c for v, c in zip( unique, get_colors( args.meta_colormap, len( unique ) ) )}
        
    # plot metadata?
    if table.metarow is not None:
        meta_ax.set_xlim( 0, len( table.metarow ) )
        for i, v in enumerate( table.metarow ):
            meta_ax.bar(
                i,
                1,
                width=1,
                color=mcdict[v],
                edgecolor="none",
                )

    # attach metadata separators
    if table.metarow is not None and args.sort[-1] == "metadata":
        xcoords = []
        for i, value in enumerate( table.metarow ):
            if i > 0 and value != table.metarow[i-1]:
                main_ax.axvline( x=i, color="black", zorder=2 )
                meta_ax.axvline( x=i, color="black", zorder=2 )
                
    # axis limits
    main_ax.set_xlim( 0, table.ncols )

    # labels
    samp_ax = main_ax if table.metarow is None else meta_ax
    samp_ax.set_xlabel( "Samples (N=%d)" % ( table.ncols ) )
    main_ax.set_title( table.title, weight="bold", size=12 )
    main_ax.set_ylabel( ylabel, size=12 )
    # tick params
    main_ax.tick_params( axis="x", which="major", direction="out", bottom="on", top="off" )
    main_ax.tick_params( axis="y", which="major", direction="out", left="on", right="off" )
    main_ax.set_xticks( [] )

    # pseudoscaling note
    if args.scaling in ["pseudolog", "pseudosqrt"]:
        xmin, xmax = main_ax.get_xlim( )
        x = xmin + 0.01 * abs( xmax - xmin )
        ymin, ymax = main_ax.get_ylim( )
        y = ymax - 0.04 * abs( ymax - ymin )
        main_ax.text( x, y,
                      "*Stratifications are proportional",
                      va="top", size=11,
                      backgroundcolor="white", )

    # optional yaxis grid
    if args.grid:
        for ycoord in main_ax.yaxis.get_majorticklocs( ):
            main_ax.axhline( y=ycoord, color="0.75", ls=":", zorder=0 )

    # legend
    L = Legendizer( anno_ax, markersize=50 )

    # add legend for stratifications
    L.subhead( "Stratifications:" )
    for i in range( len( table.rowheads ) ):
        i = -(i+1)
        value = table.rowheads[i]
        name = bugname( value )
        color = cdict[value]
        L.element( marker="s", color=color, label=name,
                   label_style="italic" if "s__" in value else "normal",
                   )

    # add legend for metadata
    if table.metarow is not None:
        L.subhead( "Sample type:" )
        levels = sorted( set( table.metarow ), key=lambda x: (0, x) if x != c_other else (1, x) )
        for i in range( len( levels ) ):
            L.element( marker="s", color=mcdict[levels[i]], label=levels[i] )
                         
    # wrapup
    L.draw( )
    plt.tight_layout( )
    fig.subplots_adjust( hspace=0.2, wspace=0.02 )
    if args.output is not None:
        plt.savefig( args.output )
    else:
        plt.savefig( os.path.split( args.input )[1]+".pdf" )
        
if __name__ == "__main__":
    main( )
