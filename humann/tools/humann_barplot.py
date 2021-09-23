#!/usr/bin/env python

import os
import sys
import re
import csv
import argparse
import math

from collections import Counter

#-------------------------------------------------------------------------------
# complex imports
#-------------------------------------------------------------------------------

try:
    from humann.tools import util
except ImportError:
    sys.exit( "CRITICAL ERROR: Unable to find the HUMAnN python package." +
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

#-------------------------------------------------------------------------------
# global constants
#-------------------------------------------------------------------------------

c_main_h             = 6
c_other_str          = "other"
c_other_color        = "0.50"
c_unclassified_str   = "unclassified"
c_unclassified_color = "0.75"
c_font1              = 10
c_font2              = 12
c_font3              = 14

#-------------------------------------------------------------------------------
# long argument help strings
#-------------------------------------------------------------------------------

c_sort_help = """Sample sorting methods (can use more than one; will evaluate in order)

none         : Maintains sample order from input file [default]
sum          : Sort on decreasing sum of stratified values
dominant     : Sort on samples' greatest taxon stratification
braycurtis   : Sort on Bray-Curtis agreement of taxon stratifications, unweighted
braycurtis_w : Sort on Bray-Curtis agreement of taxon stratifications, abundance-weighted
metadata     : Sort on specified metadata label
file         : Apply sorting order read in from a file

"""

c_choices_help = """Scaling options for total bar heights (taxa are always proportional to height)

original : Plot original units [default]
logstack : Community totals (stacked bar peaks) are log10-scaled
totalsum : Community totals (stacked bar peaks) are fixed at 1.0

"""

#-------------------------------------------------------------------------------
# command-line interface
#-------------------------------------------------------------------------------

description = util.wrap( """
HUMAnN utility for plotting a single stratified feature

Plots the taxon-stratified contributions of a specified function. Can optionally
sort samples to reveal ecological- or metadata-linked trends. Can perform custom
scaling to highlight stratifications even when community totals have high dynamic range.
""" )

def get_args( ):

    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # main arguments
    util.attach_common_arguments( parser, no_output=True )
    parser.add_argument( "-f", "--focal-feature",
                         metavar="<id>",
                         required=True,
                         help="Feature ID of interest (give ID not full name) [required]", )
    parser.add_argument( "-o", "--output",
                         metavar="<path.ext>",
                         default="humann_barplot.png",
                         help="Where and how to save the figure [humann_barplot.png]", )

    # manipulate stratified taxa data
    group1 = parser.add_argument_group( "manipulate species contributions" )
    group1.add_argument( "--top-taxa",
                         metavar="<int>",
                         type=int,
                         default=18,
                         help="Max taxon stratifications (by grand mean) to highlight [18]", )
    group1.add_argument( "--as-genera",
                         action="store_true",
                         help="Collapse species to genera [off]", )
    group1.add_argument( "--exclude-unclassified",
                         action="store_true",
                         help="Do not include the 'unclassified' taxon [off]", )
    group1.add_argument( "--remove-zeros",
                         action="store_true",
                         help="Do not analyze samples with zero sum for this feature [off]", )
    group1.add_argument( "--sort",
                         metavar = "<sorting methods>",
                         nargs="+",
                         default=["none"],
                         choices=["none", "sum", "dominant", "braycurtis", "braycurtis_w", "metadata", "file"],
                         help=c_sort_help, )
    group1.add_argument( "--taxa-colormap",
                         metavar="<named colormap OR colormap file>",
                         default=None,
                         help="Color space for taxa [automatic]", )
    group1.add_argument( "--write-taxa-colors",
                         metavar="<path>",
                         default=None,
                         help="Write taxa colors to a file for cross-plot consistency [off]", )

    # manipulate metadata
    group2 = parser.add_argument_group( "plot sample metadata" )
    group2.add_argument( "-m", "--focal-metadata",
                         metavar="<name>",
                         default=None,
                         help="Indicate metadata to highlight / group by [none]", )
    group2.add_argument( "--max-metalevels",
                         metavar="<int>",
                         type=int,
                         default=7,
                         help="Keep the most frequent metadata levels and collapse others [7]" )
    group2.add_argument( "--meta-colormap",
                         metavar="<named colormap OR colormap file>",
                         default=None,
                         help="Color space for metadata levels [automatic]", )

    # non-legend graphical tweaks
    group3 = parser.add_argument_group( "graphical tweaks" )
    group3.add_argument( "--scaling",
                         metavar="<choice>",
                         choices=["original", "logstack", "totalsum"],
                         default="original",
                         help=c_choices_help, )
    group3.add_argument( "--ylims",
                         metavar = "<limit>",
                         nargs=2,
                         default=[None, None],
                         help="Fix limits for y-axis [automatic]", )
    group3.add_argument( "--no-grid",
                         action="store_true",
                         help="Don't plot y-axis grid lines [on]", )
    group3.add_argument( "--dimensions",
                         metavar="<float>",
                         nargs=2,
                         type=float,
                         default=[11.0,6.0],
                         help="Image width and height in inches [11 6]", )
    group3.add_argument( "--units",
                         metavar="<text>",
                         default="Abundance (unspecified units)",
                         help="Name for y-axis abundance units [generic]", )

    # legend tweaks
    group4 = parser.add_argument_group( "legend layout" )
    group4.add_argument( "--legend-cols",
                         metavar="<int>",
                         type=int,
                         default=3,
                         help="Number of legend columns [3]", )
    group4.add_argument( "--legend-rows",
                         metavar="<int>",
                         type=int,
                         default=10,
                         help="Number of legend rows [10]", )
    group4.add_argument( "--legend-height",
                         metavar="<float>",
                         type=float,
                         default=1.0,
                         help="Ratio of legend to data axis height [1.0]", )

    # special
    group5 = parser.add_argument_group( "read or write sample order" )
    group5.add_argument( "--sample-order",
                         metavar="<path>",
                         default=None,
                         help="Read sample order from this file [none]", )
    group5.add_argument( "--write-sample-order",
                         metavar="<path>",
                         default=None,
                         help="Write sample order to this file [none]", )
    
    return parser.parse_args( )

#-------------------------------------------------------------------------------
# utilities
#-------------------------------------------------------------------------------

def subseq( seq, index ):
    """numpy-style slicing and indexing for lists"""
    return [seq[i] for i in index]

def taxname( rowhead ):
    """make g.s taxa names look nicer"""
    if "s__" in rowhead:
        return rowhead.split( "." )[1].replace( "s__", "" ).replace( "_", " " )
    elif "g__" in rowhead:
        return rowhead.replace( "g__", "" ).replace( "_", " " )
    else:
        return rowhead

def tsv_reader( path ):
    """quick read"""
    with util.try_zip_open( path ) as fh:
        for row in csv.reader( fh, csv.excel_tab ):
            yield row

def simplify_metadata( values, n=None ):
    """collapse rare metadata levels"""
    counts = Counter( values )
    if len( counts ) > n:
        index = [[-c, k] for k, c in counts.items( )]
        index.sort( )
        allowed = {x[1] for x in index[0:n]}
        values = [k if k in allowed else c_other_str for k in values]
    return values

def empty_axis( ax, border=False ):
    """turn an axis into a blank canvas"""
    ax.set_xticks( [] )
    ax.set_yticks( [] )
    ax.set_facecolor( "none" )
    if not border:
        for k, v in ax.spines.items( ):
            v.set_visible( False )

def ncolorlist( colormap, n ):
    """utility for defining N evenly spaced colors across a color map"""
    cmap = plt.get_cmap( colormap )
    cmap_max = cmap.N
    if n == 1:
        return [cmap( int( 0.5 * cmap_max ) )]
    else:
        return [cmap( int( k * cmap_max / (n - 1) ) ) for k in range( n )]

def generate_colormap( name, reference, items, overrides=None ):
    """generic colormap loader for taxa and metadata"""
    color_dict = {} if overrides is None else overrides
    novel = [k for k in items if k not in color_dict]
    # load colors from a file (default to black for missing)
    if reference is not None and os.path.exists( reference ):
        util.say( "Reading {} colors from file: {}".format( name, reference ) )
        for item, color in tsv_reader( reference ):
            if item not in color_dict:
                color_dict[item] = color
        for item in items:
            color_dict[item] = color_dict.get( item, "black" )
    # generate evenly spaced colors
    elif reference is not None:
        colors = ncolorlist( reference, len( novel ) )
        for item in novel:
            color_dict[item] = colors.pop( )
    # (default) assign "good enough" colors from a fixed set
    elif len( novel ) <= 18:
        colors = ncolorlist( "tab20", 20 )
        # ****skip two gray colors****
        colors = colors[0:14] + colors[16:]
        index = 0
        slide = 2 if len( novel ) <= 9 else 1
        for item in novel:
            color_dict[item] = colors[index]
            index += slide
    # too many colors to pick them automatically
    else:
        util.die( "Can't auto-color >18 {} taxa".format( name ) )
    return color_dict

#-------------------------------------------------------------------------------
# custom table class for this script
#-------------------------------------------------------------------------------
    
class BarplotTable:

    def __init__( self, path, focal_feature=None, 
        last_metadata=None, focal_metadata=None, 
        exclude_unclassified=False ):

        # table features
        self.colheads = None
        self.rowheads = []
        self.data = []
        self.metarow = None
        self.focus_name = None
        IN_FEATURES = False
        
        # pull relevant rows from input table
        for row in tsv_reader( path ):
            rowhead, values = row[0], row[1:]
            if self.colheads is None:
                self.colheads = values
                continue
            # ****focal meta and last meta can be the same thing****
            if focal_metadata is not None and rowhead == focal_metadata:
                self.metarow = values
            if last_metadata is not None and rowhead == last_metadata:
                IN_FEATURES = True
            if last_metadata is None or IN_FEATURES:
                code, name, stratum = util.fsplit( rowhead )
                if code == focal_feature and stratum is not None:
                    if stratum != c_unclassified_str or not exclude_unclassified:
                        self.focus_name = util.fjoin( code, name )
                        self.rowheads.append( stratum )
                        self.data.append( [float( k ) for k in values] )
                        
        # check that we found something
        if self.focus_name is None:
            util.die( "Requested feature <{}> was missing or not stratified".format( focal_feature ) )
    
        # update the table
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
        if not (self.nrows == len( self.rowheads ) == len( self.rowmap )):
            util.die( "Row dimension mismatch" )
        if not (self.ncols == len( self.colheads ) == len( self.colmap )):
            util.die( "Col dimension mismatch" )

    def as_genera( self ):
        temp = {}
        for rowhead, row in zip( self.rowheads, self.data ):
            rowhead = rowhead.split( "." )[0]
            temp[rowhead] = temp.get( rowhead, np.zeros( self.ncols ) ) + row
        self.rowheads = list( temp.keys( ) )
        self.data = np.array( [temp[k] for k in self.rowheads] )
        self.update( )

    def remove_zeros( self ):
        order = [i for i, k in enumerate( self.colsums ) if k > 0]
        self.colheads = subseq( self.colheads, order )
        self.data = self.data[:, order]
        if self.metarow is not None:
            self.metarow = subseq( self.metarow, order )
        self.update( )
            
    def sort( self, method=None, args=None ):
        # original order
        if method == "none":
            order = range( self.ncols )
        # bc similarity of the relative contributions
        elif method == "braycurtis":
            mat = self.data / self.colsums
            # ****linkage assumes things to cluster = rows; we want cols****
            mat = mat.transpose( )
            order = sch.leaves_list( sch.linkage( mat, method="average", metric="braycurtis" ) )
        # bc similarity of the actual values (weighted)
        elif method == "braycurtis_w":
            mat = self.data.transpose( )
            order = sch.leaves_list( sch.linkage( mat, method="average", metric="braycurtis" ) )
        # sort by dominant bug
        elif method == "dominant":
            # 1) rank bugs by row index according to their max value
            maxes = [(max( self.data[i, :] ), i) for i in range( self.nrows )]
            maxes.sort( )
            ranks = [None for k in maxes]
            for i, (s, i2) in enumerate( maxes ):
                ranks[i2] = i
            # 2) determine which bug is maximized in each sample
            argmax = []
            for c in range( self.ncols ):
                argmax.append( sorted( range( self.nrows ), key=lambda r: self.data[r][c] )[-1] )
            # 3) order samples by the rank of their max bug, then by their abundance for that bug
            order = sorted( range( self.ncols ), key=lambda c: (ranks[argmax[c]], self.data[argmax[c], c]), reverse=True )
        # sort by community totals
        elif method == "sum":
            order = sorted( range( self.ncols ), key=lambda c: self.colsums[c], reverse=True )
        # sort by metadata value
        elif method == "metadata":
            order = sorted( range( self.ncols ), key=lambda c: self.metarow[c] )
        # load sorting order from a file
        elif method == "file":
            index = []
            with open( args[0] ) as fh:
                for line in fh:
                    index.append( line.strip( ) )
            index = {line:i for i, line in enumerate( index )}
            order = sorted( range( self.ncols ), key=lambda c: index[self.colheads[c]] )
        else:
            util.die( "Can't sort with:", method )
        # reorder key data, then update
        self.data = self.data[:, order]
        self.colheads = subseq( self.colheads, order )
        if self.metarow is not None:
            self.metarow = subseq( self.metarow, order )
        self.update( )
        
    def filter_top_taxa( self, top_taxa ):
        # save unclassified if present
        unclassified = None
        if c_unclassified_str in self.rowmap:
            unclassified = self.data[self.rowmap[c_unclassified_str]]
        # select features
        best = []
        for i, f in enumerate( self.rowheads ):
            frow = self.data[i]
            if f != c_unclassified_str:
                best.append( [np.mean( frow ), f] )
        best.sort( reverse=True )
        top_taxa = min( top_taxa, len( best ) )
        select = [pair[1] for pair in best[0:top_taxa]]
        # make a new row for "other"
        other = np.zeros( self.ncols )
        for i, f in enumerate( self.rowheads ):
            frow = self.data[i]
            if f != c_unclassified_str and f not in select:
                other = other + frow
        # subset and reorder the table
        index = [self.rowmap[k] for k in select]
        self.rowheads = subseq( self.rowheads, index )
        self.data = self.data[index, :]
        if sum( other ) > 0:
            self.rowheads += [c_other_str]
            self.data = np.vstack( [self.data, other] )
        if unclassified is not None:
            self.rowheads += [c_unclassified_str]
            self.data = np.vstack( [self.data, unclassified] )
        # reverse row order to simplify plotting
        self.rowheads.reverse( )
        self.data = self.data[::-1]
        self.update( )        

#-------------------------------------------------------------------------------
# custom legend
#-------------------------------------------------------------------------------

class BarplotLegend( ):

    def __init__( self, ax, cols=None, rows=None ):
        self.ax    = ax
        self.cols  = cols
        # add in one row for the titles
        self.rows  = rows + 1
        self.cdex  = -1
        self.rdex  = 0
        self.vstep = 1 / (self.rows + 1)
        self.hstep = 0.01        
        self.htext = (1 - 3 * self.cols * self.hstep) / self.cols
        self.cinit = [k * (3 * self.hstep + self.htext) for k in range( self.cols )]

    def check( self ):
        if self.cdex + 1 > self.cols:
            util.die( "Ran out of legend space. Increase rows/cols or label fewer things." )

    def group( self, name ):
        self.rdex  = 0
        self.cdex += 1
        self.check( )
        self.ax.text(             
            self.cinit[self.cdex],
            1 - self.vstep,
            name,
            va = "center",
            ha = "left",
            size = c_font2,
            weight = "bold",
            clip_on = False,
        )

    def member( self, color="white", label="label", label_style="normal" ):
        # increment index
        if self.rdex + 1 < self.rows:
            self.rdex += 1
        else:
            self.rdex  = 1
            self.cdex += 1
        self.check( )
        # draw a rectangle for the color
        self.ax.bar(
            self.cinit[self.cdex],
            self.vstep / 2,
            bottom = 1 - (self.rdex + 1) * self.vstep - self.vstep / 4,
            align = "edge",
            width = 2 * self.hstep,
            color = color,
            edgecolor = "black",
            clip_on = False,
        )
        # draw the label
        self.ax.text(          
            self.cinit[self.cdex] + 3 * self.hstep,
            1 - (self.rdex + 1) * self.vstep,
            label,
            color = "black",
            va = "center",
            ha = "left",
            size = c_font1,
            style = label_style,
            clip_on = False,
        )

#-------------------------------------------------------------------------------
# main
#-------------------------------------------------------------------------------

def main( ):

    args = get_args( )
    
    # treat "-" as none in ylims
    a, b = args.ylims
    args.ylims[0] = None if a in ["-", None] else float( a )
    args.ylims[1] = None if b in ["-", None] else float( b )

    # load feature table 
    T = BarplotTable(
        args.input,
        focal_feature=args.focal_feature,
        focal_metadata=args.focal_metadata,
        last_metadata=args.last_metadata,
        exclude_unclassified=args.exclude_unclassified,
    )

    # collapse species to genera?
    if args.as_genera:
        T.as_genera( )
    
    # remove zero-valued samples?
    lost_samples = 0
    if args.remove_zeros:
        old = T.ncols
        T.remove_zeros( )
        new = T.ncols
        lost_samples = old - new

    # apply one or more sorting methods?
    for method in args.sort:
        if "braycurtis" in method and not args.remove_zeros:
            util.die( "Can't sort by <{}> without invoking <--remove-zeros>".format( method ) )
        T.sort( method, args=[args.sample_order if method == "file" else None] )

    # filter/collapse features (moved to take place AFTER sorting)
    T.filter_top_taxa( args.top_taxa )
        
    # simplify metadata?
    if T.metarow is not None:
        T.metarow = simplify_metadata( T.metarow, args.max_metalevels )

    # set up axis system
    main_h = c_main_h
    full_w = 1
    anno_h = math.ceil( main_h * args.legend_height )
    fig = plt.figure( )
    fig.set_size_inches( *args.dimensions )
    if T.metarow is not None:
        full_h = main_h + anno_h + 1
        main_ax = plt.subplot2grid( ( full_h, full_w ), ( 0,          0 ), rowspan=main_h, colspan=1 )
        meta_ax = plt.subplot2grid( ( full_h, full_w ), ( main_h,     0 ), rowspan=1,      colspan=1 )
        anno_ax = plt.subplot2grid( ( full_h, full_w ), ( main_h + 1, 0 ), rowspan=anno_h, colspan=1 )
        empty_axis( meta_ax, border=True )
        meta_ax.set_xlim( 0, T.ncols )
        meta_ax.set_ylim( 0, 1 )
    else:
        full_h = main_h + anno_h
        main_ax = plt.subplot2grid( ( full_h, full_w ), ( 0,          0 ), rowspan=main_h, colspan=1 )
        anno_ax = plt.subplot2grid( ( full_h, full_w ), ( main_h,     0 ), rowspan=anno_h, colspan=1 )
    main_ax.set_xlim( 0, T.ncols )
    empty_axis( anno_ax )
    anno_ax.set_xlim( 0, 1 )
    anno_ax.set_ylim( 0, 1 )

    # design taxa colors
    taxa_colors = generate_colormap(
        "taxa",
        args.taxa_colormap, 
        T.rowheads[::-1],
        {c_other_str: c_other_color, c_unclassified_str: c_unclassified_color},
    )

    # write taxa colors?
    if args.write_taxa_colors is not None:
        with open( args.write_taxa_colors, "w" ) as fh:
            for stratum in T.rowheads[::-1]:
                color = matplotlib.colors.to_hex( taxa_colors[stratum] )
                print( "{}\t{}".format( stratum, color ), file=fh )

    # scale abundance values
    if args.scaling == "original":
        ylabel  = args.units
        bottoms = np.zeros( T.ncols )
        ymin    = 0 if args.ylims[0] is None else args.ylims[0]
        ymax    = max( sum( T.data ) ) if args.ylims[1] is None else args.ylims[1]
        main_ax.set_ylim( ymin, ymax )
    elif args.scaling == "totalsum":
        ylabel  = "Relative contributions"
        T.data  = T.data / T.colsums
        bottoms = np.zeros( T.ncols )
        main_ax.set_ylim( 0, 1 )
    elif args.scaling == "logstack":
        # while plotting stacked bars...
        # 1) the top of the stacks ("crests") are log scaled
        # 2) the bot is ~arbitrary (smallest non-zero order of magnitude)
        # 3) taxa are fractions of the top-bot distance (not log'ed)
        ylabel  = args.units
        ymin    = min( [k for k in T.colsums if k > 0] ) if args.ylims[0] is None else args.ylims[0]
        floor   = math.floor( np.log10( ymin ) )
        floors  = floor * np.ones( T.ncols )
        crests  = np.array( [np.log10( k ) if k > 10**floor else floor for k in T.colsums] )
        heights = crests - floors
        T.data  = T.data / T.colsums * heights
        ymax    = max( T.colsums ) if args.ylims[1] is None else args.ylims[1]
        ceil    = math.ceil( np.log10( ymax ) )
        bottoms = floors
        main_ax.set_ylim( floor, ceil )
        ticks = list( range( floor, ceil + 1 ) )
        main_ax.set_yticks( ticks )
        main_ax.set_yticklabels( ["$10^{" + str( k ) + "}$" for k in ticks], fontsize=c_font2 )
        
    # add contribution bars
    series = []
    for i, f in enumerate( T.rowheads ):
        frow = T.data[T.rowmap[f]]
        series.append( main_ax.bar(
            range( T.ncols ),
            frow,
            align="edge",
            width=1,
            bottom=bottoms,
            color=taxa_colors[f],
            edgecolor="none", ) )
        bottoms += frow

    # plot metadata?
    if T.metarow is not None:
        # design colors
        meta_colors = generate_colormap(
                "metadata",
                args.meta_colormap,
                sorted( set( T.metarow ) ),
                {c_other_str: c_other_color},
        )
        # add bars
        for i, v in enumerate( T.metarow ):
            meta_ax.bar(
                i,
                1,
                align="edge",
                width=1,
                color=meta_colors[v],
                edgecolor="none",
                )
        # add level separators if samples grouped on metadata (via last sort)
        if args.sort[-1] == "metadata":
            xcoords = []
            for i, value in enumerate( T.metarow ):
                if i > 0 and value != T.metarow[i-1]:
                    main_ax.axvline( x=i, color="black", lw=1.0, zorder=2 )
                    meta_ax.axvline( x=i, color="black", lw=1.0, zorder=2 )

    # add plot title
    main_ax.set_title( T.focus_name, weight="bold", size=c_font3 )

    # label x-axis (indicate possible sample loss)
    samp_ax = main_ax if T.metarow is None else meta_ax
    xlabel = "{:,} samples".format( T.ncols )
    if lost_samples > 0:
        xlabel += " of {:,} total (removed {:,} with zero stratified abundance)".format(
            T.ncols + lost_samples,
            lost_samples,
        )
    samp_ax.set_xlabel( xlabel, size=c_font2 )

    # label y-axis (defined during scaling)
    main_ax.set_ylabel( ylabel, size=c_font2 )
    
    # modify tick params
    main_ax.tick_params( axis="x", which="major", direction="out", bottom=True, top=False )
    main_ax.tick_params( axis="y", which="major", direction="out", left=True, right=False )
    main_ax.set_xticks( [] )

    # add optional yaxis grid
    if not args.no_grid:
        for ycoord in main_ax.yaxis.get_majorticklocs( ):
            main_ax.axhline( y=ycoord, color="0.75", ls="--", lw=1.0, zorder=0 )

    # define the legend
    L = BarplotLegend( anno_ax, cols=args.legend_cols, rows=args.legend_rows )
    L.group( "Stratified contributions:" )
    for i in range( len( T.rowheads ) ):
        i = -(i+1)
        value = T.rowheads[i]
        name = taxname( value )
        color = taxa_colors[value]
        L.member( 
            color=color, 
            label=name,
            label_style="italic" if re.search( "^[gs]__", value ) else "normal",
        )
    if T.metarow is not None:
        L.group( "Sample label (metadata):" )
        levels = sorted( set( T.metarow ), key=lambda x: (0, x) if x != c_other_str else (1, x) )
        for i in range( len( levels ) ):
            L.member( color=meta_colors[levels[i]], label=levels[i] )

    # write sample order?
    if args.write_sample_order is not None:
        with open( args.write_sample_order, "w" ) as fh:
            for sampleid in T.colheads:
                print( sampleid, file=fh )
    # wrapup
    plt.tight_layout( )
    fig.subplots_adjust( hspace=1.0, wspace=0.0 )
    plt.savefig( args.output, dpi=300 )
        
if __name__ == "__main__":
    main( )
