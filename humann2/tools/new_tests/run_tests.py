#!/usr/bin/env python

import os
import sys
import glob

os.chdir( sys.argv[1] )

outputs = glob.glob( "*.output" )
for o in outputs:
    os.remove( o )
outputs = glob.glob( "*.stderr" )
for o in outputs:
    os.remove( o )

files = set( glob.glob( "*.command" ) )
for f in files:
    print "--> Executing:", f
    os.system( "bash {} 2> {}.stderr".format( f, f ) )

files = set( glob.glob( "*.expected" ) )
for f in files:
    print "--> expected output", f
    f2 = f.replace( ".expected", "" )
    if os.path.exists( f2 ):
        print "--> found matching output", f2, "(comparing)"
        os.system( "diff {} {}".format( f, f2 ) )
    else:
        sys.exit( "LETHAL ERROR: missing expected output {}".format( f2 ) )
