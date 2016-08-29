#!/usr/bin/env python

from __future__ import print_function
import urllib2
import gzip

"""
Example "stanza" from the go-basic obo file:

[Term]
id: GO:0000011
name: vacuole inheritance
namespace: biological_process
def: "The distribution of vacuoles ... " [GOC:mcc, PMID:10873824, PMID:14616069]
is_a: GO:0007033 ! vacuole organization
is_a: GO:0048308 ! organelle inheritance
"""

url = "http://geneontology.org/ontology/go-basic.obo"

opener = urllib2.build_opener( )
opener.addheaders = [('User-agent', 'Mozilla/5.0')]

wh = opener.open( url )

data = {}
for line in wh:
    line = line.strip( )
    items = line.split( ":" )
    key = items[0]
    text = ":".join( items[1:] ).strip( )
    if key == "id":
        inner = data.setdefault( text, {} )
    elif len( data ) > 0:
        inner[key] = text

with gzip.GzipFile( "map_go_name.txt.gz", "w" ) as fh:
    for goid in sorted( data ):
        namespace = data[goid]["namespace"]
        namespace = [k[0] for k in namespace.split( "_" )]
        namespace = "".join( namespace ).upper( )
        name = data[goid]["name"]
        name = "[{}] {}".format( namespace, name )
        print( goid+"\t"+name, file=fh )
