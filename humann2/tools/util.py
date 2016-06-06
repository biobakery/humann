from __future__ import print_function # PYTHON 2.7+ REQUIRED
import os
import sys
import re
import subprocess
import csv
import gzip
import bz2

# ---------------------------------------------------------------
# utilities used by the split and join tables scripts
# ---------------------------------------------------------------

def find_exe_in_path(exe):
    """
    Check that an executable exists in $PATH
    """
    
    paths = os.environ["PATH"].split(os.pathsep)
    for path in paths:
        fullexe = os.path.join(path,exe)
        if os.path.exists(fullexe):
            if os.access(fullexe,os.X_OK):
                return True
    return False

def biom_to_tsv(biom_file, new_tsv_file, taxonomy=None):
    """
    Convert from a biom to tsv file
    """

    cmd=["biom","convert","-i",biom_file,"-o",new_tsv_file,"--to-tsv"]

    # check if taxonomy is set (can be set to zero)
    if taxonomy != None:
        cmd+=["--header-key","taxonomy"]
    
    try:
        if os.path.isfile(new_tsv_file):
            os.unlink(new_tsv_file)
        p_out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except (EnvironmentError, subprocess.CalledProcessError):
        command=" ".join(cmd)
        sys.exit("Unable to convert biom file to tsv"+"\n"+command)
        
def tsv_to_biom(tsv_file, biom_file):
    """
    Convert from a biom to tsv file
    """
    
    cmd=["biom","convert","-i",tsv_file,"-o",biom_file,"--table-type","Gene table","--to-hdf5"]

    try:
        if os.path.isfile(biom_file):
            os.unlink(biom_file)
        p_out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except (EnvironmentError, subprocess.CalledProcessError):
        command=" ".join(cmd)
        sys.exit("Unable to convert tsv file to biom"+"\n"+command)
        
def process_gene_table_header(gene_table, allow_for_missing_header=None):
    """
    Process through the header portion of the gene table file
    """
    
    # indicator of a comment line or the header
    # the last line in the file with this indicator is the header
    GENE_TABLE_COMMENT_LINE="^#"
    
    # try to open the file
    try:
        file_handle=open(gene_table,"r")
        line=file_handle.readline()
    except EnvironmentError:
        sys.exit("Unable to read file: " + gene_table)
            
    # find the headers
    header_flag=False
    header=line
    while re.match(GENE_TABLE_COMMENT_LINE, line):
        header_flag=True
        header = line
        line = file_handle.readline()
            
    if not header_flag and not allow_for_missing_header:
        sys.exit("File does not have a required header: " + gene_table +
            " . Please add a header which includes the indicator: " +
            GENE_TABLE_COMMENT_LINE)
        
    if not header_flag and allow_for_missing_header:
        header = ""
        
    return file_handle, header, line

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# utilities used by the rename, renorm, etc. scripts
# ---------------------------------------------------------------
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_strat_delim     = "|"
c_taxon_delim     = "."
c_name_delim      = ": "
c_multiname_delim = ";"
c_str_unknown     = "NO_NAME"
c_ungrouped       = "UNGROUPED"
c_unmapped        = "UNMAPPED"
c_unintegrated    = "UNINTEGRATED"
c_many_bytes      = 1e8
c_zip_multiplier  = 10

c_topsort = {
    c_unmapped:0,
    c_ungrouped:1,
    c_unintegrated:2,
    "UniRef50_unknown":3,
    "UniRef90_unknown":4,
}

# ---------------------------------------------------------------
# helper classes
# ---------------------------------------------------------------

class Table ( ):

    """ 
    Very basic table class; would be more efficient using numpy 2D array
    """

    def __init__ ( self, path ):
        self.anchor = None
        self.colheads = []
        self.rowheads = []
        self.data = []
        self.is_stratified = False
        if path is None:
            fh = sys.stdin
            path = "STDIN"
            print( "Loading table from: <STDIN>", file=sys.stderr )
        else:
            fh = try_zip_open( path )
            print( "Loading table from:", path, file=sys.stderr )
            size_warn( path )
        for row in csv.reader( fh, dialect='excel-tab' ):
            try:
                if self.anchor is None:
                    self.anchor = row[0]
                    self.colheads = row[1:]
                else:
                    self.rowheads.append( row[0] )
                    self.data.append( row[1:] )
            except IndexError:
                # ignore empty lines in input file
                pass
        for rowhead in self.rowheads:
            if c_strat_delim in rowhead:
                print( "  Treating", path, "as stratified output, e.g.", 
                       rowhead.split( c_strat_delim ), file=sys.stderr )
                self.is_stratified = True
                break

    def write ( self, path=None, unfloat=False ):
        fh = try_zip_open( path, "w" ) if path is not None else sys.stdout
        writer = csv.writer( fh, delimiter="\t", lineterminator="\n")
        writer.writerow( [self.anchor] + self.colheads )
        for i in range( len( self.rowheads ) ):
            values = self.data[i][:]
            if unfloat:
                values = map( lambda x: "%.6g" % ( x ), values )
            writer.writerow( [self.rowheads[i]] + values )

class Ticker( ):
    def __init__( self, iterable, step=100, pad="  " ):
        self.count = 0        
        self.total = len( iterable )
        self.step = 100
        self.pad = pad
    def tick( self ):
        self.count += 1
        if self.count % self.step == 0:
            self.report( )
    def report( self ):
        frac = self.count / float( self.total )
        print( self.pad+"{:.1f}%".format( 100 * frac ), file=sys.stderr, end="\r" )

# ---------------------------------------------------------------
# helper functions
# ---------------------------------------------------------------

def size_warn( path ):
    m = 1 if ".gz" not in path else c_zip_multiplier
    if m * os.path.getsize( path ) > c_many_bytes:
        print( "  This is a large file, one moment please...", file=sys.stderr )

def try_zip_open( path, *args ):
    """ 
    open an uncompressed or gzipped file; fail gracefully 
    """
    fh = None
    try:
        if path.endswith(".gz"):
            fh = gzip.GzipFile( path, *args )
        elif path.endswith(".bz2"):
            fh = bz2.BZ2File( path, *args )
        else:
            fh = open( path, *args )
    except EnvironmentError:
        sys.exit( "Problem opening file: " + path)
    return fh

def load_polymap ( path, start=0, skip=None, allowed_keys=None, allowed_values=None ):
    """
    Load a file like:
    A 1 2
    B 1
    B 3
    C 1 2 4
    To a nested dict structure:
    {A:{1:1, 2:1}, B:{1:1, 3:1}, C:{1:1, 2:2, 4:1}
    Inner values are not important (set to 1)
    """
    polymap = {}
    with try_zip_open( path ) as fh:
        print( "Loading mapping file from:", path, file=sys.stderr )
        size_warn( path )
        for row in csv.reader( fh, dialect="excel-tab" ):
            key = row[start]
            if allowed_keys is None or key in allowed_keys:
                for i, value in enumerate( row ):
                    if i != start and (skip is None or i not in skip):
                        if allowed_values is None or value in allowed_values:
                            polymap.setdefault( key, {} )[value] = 1
    return polymap

def fsplit( feature ):
    items = feature.split( c_strat_delim )
    stratum = None if len( items ) == 1 else items[1]
    items = items[0].split( c_name_delim )
    name = None if len( items ) == 1 else items[1]
    feature = items[0]
    return feature, name, stratum

def fjoin( feature, name=None, stratum=None ):
    if name is not None:
        feature = c_name_delim.join( [feature, name] )
    if stratum is not None:
        feature = c_strat_delim.join( [feature, stratum] )
    return feature

def fsort( features ):
    # force 1|A to come before 11
    features = sorted( features, key=lambda f: f.split( c_strat_delim ) )
    # force special features to the top (defined above)
    default = 1 + max( c_topsort.values() )
    features = sorted( features, key=lambda f: c_topsort.get( fsplit( f )[0], default ) )
    return features
