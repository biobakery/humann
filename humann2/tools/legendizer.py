#!/usr/bin/env python

"""
Adapted from zopy.mplutils2.legendizer
"""

class Legendizer( ):

    def __init__( self, ax, hscale=1.0, vscale=1.0, markersize=10, labelsize=8 ):
        self.ax = ax
        self.spacing = 0.05 * vscale
        self.start = 0.5 #1 - self.spacing / 2.0
        self.margin = 0.1 * hscale
        self.textbuffer = 0.05 * hscale
        self.labelsize = labelsize
        self.markersize = markersize
        self.height = self.start
        self.delta = 0
        self.ax.set_ylim( 0, 1 )
        self.ax.set_xlim( 0, 1 )
        self.colortext = False
        self.commands = []

    def element( self, *args, **kwargs ):
        self.delta += 0.5 * self.spacing
        self.commands.append( ["element", args, kwargs] )
        
    def draw_element( self, marker="o", color="white", edgecolor="black", 
                      label="NO_LABEL", label_style="normal", ):
        if "text:" in marker:
            self.ax.text( self.margin, self.height,
                          marker.replace( "text:", "" ),
                          va="center", ha="center",
                          size=self.labelsize,
                          color=color, weight="bold", )
        else:
            self.ax.scatter( [self.margin], [self.height],
                          marker=marker, s=self.markersize,
                          color=color, edgecolor=edgecolor, )
            # background?
            if color in ["white", "none"] and edgecolor in ["white", "none"]:
                self.ax.scatter( [self.margin], [self.height], marker="s",
                                 color="black", edgecolor="none",
                                 s=1.5 * self.markersize, zorder=0, )
        self.ax.text( self.margin + self.textbuffer,
                      self.height, label,
                      color="black" if not self.colortext else color,
                      va="center", ha="left",
                      size=self.labelsize, clip_on=False, style=label_style )
        self.height -= self.spacing

    def subhead( self, *args, **kwargs ):
        self.skip( )
        self.delta += 0.5 * self.spacing
        self.commands.append( ["subhead", args, kwargs] )                          

    def draw_subhead( self, text ):
        self.ax.text( 0.65 * self.margin, self.height,
                      text,
                      va="center", ha="left",
                      size=self.labelsize * 1.2,
                      weight="bold", )
        self.height -= self.spacing

    def skip( self, *args, **kwargs ):
        self.delta += 0.5 * 0.5 * self.spacing
        self.commands.append( ["skip", args, kwargs] )
                                  
    def draw_skip( self ):
        self.height -= self.spacing / 2.0

    def draw( self ):
        self.height += self.delta
        for name, args, kwargs in self.commands:
            if name == "element":
                self.draw_element( *args, **kwargs )
            elif name == "subhead":
                self.draw_subhead( *args, **kwargs )
            elif name == "skip":
                self.draw_skip( *args, **kwargs )
