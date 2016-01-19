#! /usr/bin/env python

"""
Make the mapping from EC # to EC Name
"""

import gzip
import subprocess

c_url = "ftp://ftp.expasy.org/databases/enzyme/enzyme.dat"
c_id_prefix = "ID   "
c_de_prefix = "DE   "
c_no_name = "NO_NAME"

mapping = {}
cmd = subprocess.Popen( "curl {}".format( c_url ), shell=True, stdout=subprocess.PIPE )
for line in cmd.stdout:
    line = line.strip()
    if c_id_prefix in line:
        other, name = line.split( c_id_prefix )
        mapping[name] = c_no_name
    if c_de_prefix in line:
        other, desc = line.split( c_de_prefix )
        # slice to remove trailing period
        mapping[name] = desc[:-1] 

with gzip.GzipFile( "map_ec_name.txt.gz", "w" ) as fh:
    for name in sorted( mapping ):
        print >>fh, "\t".join( [name, mapping[name]] )
