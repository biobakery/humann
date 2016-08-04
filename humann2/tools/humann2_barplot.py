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
c_unit_width = 0.1
c_min_width = 6
c_height = 6

# ---------------------------------------------------------------
# argument parsing
# ---------------------------------------------------------------

def get_args( ):
    parser = argparse.ArgumentParser()
    parser.add_argument( "-i", "--input",
                         help="HUMAnN2 table with optional metadata", )
    parser.add_argument( "-f", "--focal-feature",
                         help="Feature ID of interest (give ID not full name)", )
    parser.add_argument( "-t", "--top-strata",
                         type=int,
                         default=7,
                         help="Number of top stratifications (species) to highlight by highest grand means", )
    parser.add_argument( "-s", "--sort",
                         nargs="+",
                         default=["none"],
                         choices=["none", "sum", "dominant", "similarity", "usimilarity", "metadata"],
                         help="Sorting methods (can use more than one; will evaluate in order).", )
    parser.add_argument( "-n", "--norm",
                         action="store_true",
                         help="Renormalize to highlight within-sample composition", )
    parser.add_argument( "-l", "--last-metadatum",
                         default=None,
                         help="Indicate end of metadata rows", )
    parser.add_argument( "-m", "--focal-metadatum",
                         default=None,
                         help="Indicate metadatum to highlight / group by", )
    parser.add_argument( "-c", "--colormap",
                         default="jet",
                         help="Color space for top stratifications", )
    parser.add_argument( "-k", "--meta-colormap",
                         default="jet",
                         help="Color space for metadata levels" )
    parser.add_argument( "-x", "--exclude-unclassified",
                         action="store_true",
                         help="Do not include the 'unclassified' stratum", )
    parser.add_argument( "-o", "--output",
                         default=None,
                         help="Where to save the figure", )
    parser.add_argument( "-y", "--ymax",
                         default=None,
                         type=float,
                         help="Force y-max", )
    return parser.parse_args()

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

def subseq( seq, index ):
    """numpy-style slicing and indexing for lists"""
    return [seq[i] for i in index]

def bugname( rowhead ):
    """make bug names look nicer"""
    if "." in rowhead:
        return rowhead.split( "." )[1].replace( "s__", "" ).replace( "_", " " )
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
        adding = False
        with open( path ) as fh:
            for row in csv.reader( fh, dialect="excel-tab" ):
                rowhead, values = row[0], row[1:]
                if self.colheads is None:
                    self.colheads = values
                else:
                    if (adding or last is None) and "|" in rowhead:
                        feature, stratum = rowhead.split( "|" )
                        feature_id = feature.split( ": " )[0]
                        if focus is None:
                            focus = feature_id
                            print >>sys.stderr, "No feature selected, so defaulting to first feature:", feature
                        if focus != feature_id:
                            continue
                        else:
                            self.fname = feature
                            if stratum != "unclassified" or not exclude_unclassified:
                                self.rowheads.append( stratum )
                                self.data.append( map( float, values ) )
                    if metaheader is not None and rowhead == metaheader:
                        self.metarow = values
                    if last is not None and rowhead == last:
                        adding = True
        self.rowheads = map( bugname, self.rowheads )
        self.data = np.array( self.data )
        self.update()
        
    def update( self ):
        self.nrows, self.ncols = self.data.shape
        self.colsums = sum( self.data )
        self.rowmap = {}
        self.colmap = {}
        for i, h in zip( range( self.nrows ), self.rowheads ):
            self.rowmap[h] = i
        for i, h in zip( range( self.ncols ), self.colheads ):
            self.colmap[h] = i
        assert self.nrows == len( self.rowheads ) == len( self.rowmap ), "row dim failure"
        assert self.ncols == len( self.colheads ) == len( self.colmap ), "col dim failure"

    def norm( self ):
        self.data = self.data / ( self.colsums + c_epsilon * np.ones( self.ncols ) )

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
        self.rowheads.reverse()
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
    # normalize  
    if args.norm:
        table.norm()
    # apply one or more sorting methods
    for method in args.sort:
        table.sort( method )
    # filter/collapse features (moved to take place AFTER sorting)
    table.filter_top_strata( args.top_strata )
        
    # set up axis system
    fig = plt.figure()
    fig.set_size_inches( 10, 5 )

    def dumb( ax ):
        ax.set_xticks( [] )
        ax.set_yticks( [] )
        for k, v in ax.spines.items():
            v.set_visible( False )

    if table.metarow is not None:
        main_ax = plt.subplot2grid( ( 9, 3 ), ( 0, 0 ), rowspan=8, colspan=2 )
        anno_ax = plt.subplot2grid( ( 9, 3 ), ( 0, 2 ), rowspan=9, colspan=1 )
        meta_ax = plt.subplot2grid( ( 9, 3 ), ( 8, 0 ), rowspan=1, colspan=2 )
        dumb( anno_ax )
        dumb( meta_ax )
    else:
        main_ax = plt.subplot2grid( ( 1, 3 ), ( 0, 0 ), rowspan=1, colspan=2 )
        anno_ax = plt.subplot2grid( ( 1, 3 ), ( 0, 2 ), rowspan=1, colspan=1 )
        dumb( anno_ax )
        
    # setup: colors
    cdict = {"Other":"0.5", "Unclassified":"0.8"}
    colors = get_colors( args.colormap, args.top_strata )
    for f in table.rowheads:
        if f not in cdict:
            cdict[f] = colors.pop()
    # add bars
    series = []
    bottoms = np.zeros( table.ncols )
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

    # plot metadata?
    if table.metarow is not None:
        unique = sorted( set( table.metarow ) )
        mc = {v:c for v, c in zip( unique, get_colors( args.meta_colormap, len( unique ) ) )}
        meta_ax.set_xlim( 0, len( table.metarow ) )
        for i, v in enumerate( table.metarow ):
            meta_ax.bar(
                i,
                1,
                width=1,
                color=mc[v],
                edgecolor="none",
                )
        
    # axis limits
    main_ax.set_xlim( 0, table.ncols )
    main_ax.set_ylim( 0, max( bottoms ) )  
    # labels
    samp_ax = main_ax if table.metarow is None else meta_ax
    samp_ax.set_xlabel( "Samples (N=%d)" % ( table.ncols ) )
    main_ax.set_title( "%s\nFrom: %s" % ( table.fname, os.path.split( args.input )[1] ) )
    main_ax.set_ylabel( "Relative Abundance" )
    # tick params
    main_ax.tick_params( axis="x", which="major", direction="out", bottom="on", top="off" )
    main_ax.tick_params( axis="y", which="major", direction="out", left="on", right="off" )
    main_ax.set_xticks( [] )

    # norm override   
    if args.norm:
        main_ax.set_yticks( [] )
        main_ax.set_ylabel( "" )

    if args.ymax is not None:
        main_ax.set_ylim( [0, args.ymax] )                        
        main_ax.text( 0, args.ymax, "clipped at y="+str( args.ymax ), va="top" )
        
    # legend
    xmar = 0.01
    xsep = 0.03
    ydex = 0.99
    ysep = 0.01
    yinc = 0.01
    rech = 0.025
    recw = 0.1

    def add_items( title, labels, colors, ydex, bugmode=False ):
        anno_ax.text( xmar, ydex, title, weight="bold", va="top" )
        ydex -= 3.5 * rech
        for l, c in zip( labels, colors ):
            anno_ax.add_patch(
                patches.Rectangle(
                    (xmar, ydex),
                    recw,
                    rech,
                    facecolor=c,
                    edgecolor="k",
                )
            )
            anno_ax.text(
                xsep + recw,
                ydex + 0.5 * rech,
                l,
                size=8,
                va="center",
                style="italic" if (bugmode and l not in ["Unclassified", "Other"]) else "normal",
            )
            ydex -= ysep + rech
        return ydex
        return 
        
    ydex = add_items(
        "Functional strata",
        table.rowheads[::-1],
        [cdict[r] for r in table.rowheads[::-1]],
        ydex,
        bugmode=True,
        )

    if table.metarow is not None:
        add_items(
            args.focal_metadatum,
            sorted( mc.keys() ),
            [mc[k] for k in sorted( mc.keys() )],
            ydex,
            )
    
    # output
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
