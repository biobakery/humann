#!/usr/bin/env python

"""
HUMAnN2 Plotting Tool
===============================================
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

import os
import sys
import csv
import argparse
import math

try:
    import matplotlib
    matplotlib.use( "Agg" )
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import numpy as np
    import scipy.cluster.hierarchy as sch
except:
    sys.exit( "This script requires the Python scientific stack: numpy, scipy, and matplotlib." )
    
# ---------------------------------------------------------------
# global constants
# ---------------------------------------------------------------

c_epsilon = 1e-10
c_hunits = 9

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
pseudolog   : Total bar heights log-scaled (strata are NOT log scaled)
normalize   : Bars all have height=1 (highlighting relative attribution)

"""

def get_args( ):
    parser = argparse.ArgumentParser(
        description="HUMAnN2 plotting tool",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument( "-i", "--input",
                         metavar = "<input table>",
                         required = True, 
                         help="HUMAnN2 table with optional metadata", )
    parser.add_argument( "-f", "--focal-feature",
                         metavar = "<feature id>",
                         help="Feature ID of interest (give ID not full name)", )
    parser.add_argument( "-t", "--top-strata",
                         metavar = "<int>",
                         type=int,
                         default=7,
                         help="Number of top stratifications to highlight (top = highest grand means)", )
    parser.add_argument( "-s", "--sort",
                         metavar = "<sorting methods>",
                         nargs="+",
                         default=["none"],
                         choices=["none", "sum", "dominant", "similarity", "usimilarity", "metadata"],
                         help=c_sort_help, )             
    parser.add_argument( "-l", "--last-metadatum",
                         metavar = "<feature>",
                         default=None,
                         help="Indicate end of metadata rows", )
    parser.add_argument( "-m", "--focal-metadatum",
                         metavar="<feature>",
                         default=None,
                         help="Indicate metadatum to highlight / group by", )
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
                         choices=["none", "normalize", "pseudolog"],
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
                         help="Relative width of the plot vs. legend (default: 5)", )
    parser.add_argument( "-d", "--dimensions",
                         metavar = "<size>",
                         nargs=2,
                         type=float,
                         default=[10.0,4.0],
                         help="Image height and width in inches (default: 8 4)", )
    parser.add_argument( "-y", "--ylims",
                         metavar = "<limit>",
                         nargs=2,
                         type=float,
                         default=[None, None],
                         help="Fix limits for y-axis", )
    parser.add_argument( "-e", "--legend-stretch",
                         metavar = "",
                         type=float,
                         default=1.0,
                         help="Stretch/compress legend elements", )
    return parser.parse_args()

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
            
# ---------------------------------------------------------------
# table class
# ---------------------------------------------------------------
    
class FeatureTable:

    """ """
    
    def __init__( self, path, focus=None, last=None, metaheader=None, exclude_unclassified=False ):
        self.colheads = None
        self.rowheads = []
        self.data = []
        self.metarow = None
        self.fname = None
        in_features = False
        for row in tsv_reader( path ):
            rowhead, values = row[0], row[1:]
            if self.colheads is None:
                self.colheads = values
            else:
                if (in_features or last is None) and "|" in rowhead:
                    feature, stratum = rowhead.split( "|" )
                    feature_id = feature.split( ": " )[0]
                    if focus is None:
                        focus = feature_id
                        sys.stderr.write( "No feature selected; defaulting to 1st feature: " + feature + "\n" )
                    if focus != feature_id:
                        continue
                    else:
                        self.fname = feature
                        if stratum != "unclassified" or not exclude_unclassified:
                            self.rowheads.append( stratum )
                            self.data.append( list( map( float, values ) ) )
                if metaheader is not None and rowhead == metaheader:
                    self.metarow = values
                if last is not None and rowhead == last:
                    in_features = True
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
        sys.stderr.write( "Regrouping to genera (before selecting/sorting strata)\n" )
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
            norm = self.data / ( self.colsums + c_epsilon * np.ones( self.ncols ) )
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
        self.update()
        
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
        # reverse row order to simplfy plotting
        self.rowheads.reverse( )
        self.data = self.data[::-1]
        self.update()

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):

    args = get_args()
    
    # load table manipulation 
    table = FeatureTable(
        args.input,
        focus=args.focal_feature,
        last=args.last_metadatum,
        metaheader=args.focal_metadatum,
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
        
    # set up axis system
    wunits = args.width + 1
    fig = plt.figure()
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
        sys.stderr.write( "Reading strata colors from file: {}\n".format( args.colormap ) )
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
        ylabel      = "log10(Relative abundance)"
        ymin        = min( [k for k in table.colsums if k > 0] ) if args.ylims[0] is None else args.ylims[0]
        floor       = math.floor( np.log10( ymin ) )
        floor       = floor if abs( np.log10( ymin ) - floor ) >= c_epsilon else floor - 1
        floors      = floor * np.ones( table.ncols )
        crests      = np.array( [np.log10( k ) if k > 10**floor else floor for k in table.colsums] )
        heights     = crests - floors
        table.data  = table.data / table.colsums * heights
        ymax        = max( table.colsums ) if args.ylims[1] is None else args.ylims[1]
        ceil        = math.ceil( np.log10( ymax ) )
        bottoms = floors
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
            mcdict = {}
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
    #main_ax.set_title( "%s\nFrom: %s" % ( table.fname, os.path.split( args.input )[1] ) )
    main_ax.set_title( table.fname, weight="bold" )
    main_ax.set_ylabel( ylabel, size=12 )
    # tick params
    main_ax.tick_params( axis="x", which="major", direction="out", bottom="on", top="off" )
    main_ax.tick_params( axis="y", which="major", direction="out", left="on", right="off" )
    main_ax.set_xticks( [] )

    # pseudolog note
    if args.scaling == "pseudolog":
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
    scale = args.legend_stretch
    xmar = 0.01
    xsep = 0.03
    ybuf = 0.06 * scale
    ysep = 0.02 * scale
    yinc = 0.03 * scale
    rech = 0.03 * scale
    recw = 0.1
    big_font = 10 * scale
    sml_font = 8 * scale
    ydex = 1.0

    def add_items( title, labels, colors, ydex, bugmode=False ):
        ydex -= ybuf
        anno_ax.text( xmar, ydex, title, weight="bold", va="center", size=big_font )
        ydex -= ybuf
        for l, c in zip( labels, colors ):
            anno_ax.add_patch(
                patches.Rectangle(
                    (xmar, ydex-rech),
                    recw,
                    rech,
                    facecolor=c,
                    edgecolor="k",
                )
            )
            anno_ax.text(
                xsep + recw,
                ydex - 0.5 * rech,
                l,
                size=sml_font,
                va="center",
                style="italic" if (bugmode and l not in ["Unclassified", "Other"]) else "normal",
            )
            ydex -= rech + ysep
        ydex += ysep
        return ydex

    # add legend for stratifications
    ydex = add_items(
        "Stratifications:",
        map( bugname, table.rowheads[::-1] ),
        [cdict[r] for r in table.rowheads[::-1]],
        ydex,
        bugmode=True,
        )

    # add legend for metadata
    if table.metarow is not None:
        levels = sorted( set( table.metarow ) )
        add_items(
            "Sample label:",
            levels,
            [mcdict[k] for k in levels],
            ydex,
            )
                         
    # wrapup
    plt.tight_layout()
    fig.subplots_adjust( hspace=0.2, wspace=0.03 )
    if args.output is not None:
        plt.savefig( args.output )
    else:
        plt.savefig( os.path.split( args.input )[1]+".pdf" )

# ---------------------------------------------------------------
# boilerplate
# ---------------------------------------------------------------
        
if __name__ == "__main__":
    main()
