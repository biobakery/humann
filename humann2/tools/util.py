from __future__ import print_function # PYTHON 2.7+ REQUIRED
import os
import sys
import re
import subprocess
import csv
import gzip

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
# utilities used by the rename, renorm scripts
# ---------------------------------------------------------------

# constants
c_strat_delim = "|"
c_name_delim = ": "
c_multiname_delim = ";"
c_str_unknown = "NO_NAME"

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
            print( "Loading table from STDIN", file=sys.stderr )
        else:
            fh = try_zip_open( path )
        for row in csv.reader( fh, dialect='excel-tab' ):
            if self.anchor is None:
                self.anchor = row[0]
                self.colheads = row[1:]
            else:
                self.rowheads.append( row[0] )
                self.data.append( row[1:] )
        for rowhead in self.rowheads:
            if c_strat_delim in rowhead:
                print( "Treating", path, "as stratified output, e.g.", 
                       rowhead.split( c_strat_delim ), file=sys.stderr )
                self.is_stratified = True
                break

    def write ( self, path=None, unfloat=False ):
        fh = try_zip_open( path, "w" ) if path is not None else sys.stdout
        writer = csv.writer( fh, dialect='excel-tab' )
        writer.writerow( [self.anchor] + self.colheads )
        for i in range( len( self.rowheads ) ):
            values = self.data[i][:]
            if unfloat:
                values = map( lambda x: "%.6g" % ( x ), values )
            writer.writerow( [self.rowheads[i]] + values )

def try_zip_open( path, *args ):
    """ 
    open an uncompressed or gzipped file; fail gracefully 
    """
    fh = None
    try:
        fh = open( path, *args ) if not re.search( r".gz$", path ) else gzip.GzipFile( path, *args )
    except:
        print( "Problem opening", path, file=sys.stderr )
    return fh

def load_polymap ( path ):
    """ 
    load a tsv file mapping one name to another (e.g. uniref50 id to english name)
    """
    polymap = {}
    with try_zip_open( path ) as fh:
        for row in csv.reader( fh, dialect="excel-tab" ):
            old_name = row[0]
            for new_name in row[1:]:
                polymap.setdefault( old_name, {} )[new_name] = 1
    print( "Loaded polymap from", path, file=sys.stderr )
    return polymap
