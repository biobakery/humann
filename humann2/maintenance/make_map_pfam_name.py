#!/usr/bin/env python

import os

url = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz"

os.system( "curl {} | zcat | cut -f1,5 | gzip > map_pfam_name.txt.gz".format( url ) )
