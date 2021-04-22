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

c_unclassified="unclassified"

# indicator of a comment line or the header
# the last line in the file with this indicator is the header
GENE_TABLE_COMMENT_LINE="#"

# the extension used for biom files
BIOM_FILE_EXTENSION=".biom"

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
        
def process_gene_table_with_header(gene_table, allow_for_missing_header=None):
    """
    Process through the header portion of the gene table file
    """
    
    # try to open the file
    try:
        lines=gzip_bzip2_biom_open_readlines( gene_table )
    except EnvironmentError:
        sys.exit("Unable to read file: " + gene_table)
            
    # find the headers
    header=""
    first_data_line=""
    for line in lines:
        if line[0] == GENE_TABLE_COMMENT_LINE:
            header = line
        else:
            first_data_line = line
            break
            
    if not header and not allow_for_missing_header:
        sys.exit("File does not have a required header: " + gene_table +
            " . Please add a header which includes the indicator: " +
            GENE_TABLE_COMMENT_LINE)
        
    # provide the header, if one was found
    if header:
        yield header
    
    # provide the first data line
    yield first_data_line
    
    # now provide the remaining lines
    for line in lines:
        yield line

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
            rows = csv.reader( sys.stdin, dialect='excel-tab' )
            path = "STDIN"
            print( "Loading table from: <STDIN>", file=sys.stderr )
        else:
            rows = [line.split("\t") for line in process_gene_table_with_header( path, True )]
            print( "Loading table from:", path, file=sys.stderr )
            size_warn( path )
        for row in rows:
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
        """ If the output file has the biom extension, then write a biom output file """
        
        # get the rows of output to write
        rows = self.write_rows(unfloat)
        
        # check if the output file specified is a biom file based on the extension
        try:
            biom_file = path.endswith(BIOM_FILE_EXTENSION)
        except (AttributeError, ValueError):
            biom_file = False
        
        if biom_file:
            write_biom(path, rows)
        else:
            write_tsv(path, rows)
        
    def write_rows(self, unfloat=False):
        """ Yield the values to write to the output file """
        yield [self.anchor] + self.colheads
        for i in range( len( self.rowheads ) ):
            values = self.data[i][:]
            if unfloat:
                values = list( map( lambda x: "%.6g" % ( x ), values ))
            yield [self.rowheads[i]] + values

def write_tsv( path, rows ):
    """ Write the output in tsv (possibly compressed) format to a file or stdout """
    fh = try_zip_open( path, write=True) if path is not None else sys.stdout
    writer = csv.writer( fh, delimiter="\t", lineterminator="\n")

    for row in rows:
        writer.writerow(row)
            
def write_biom( path, rows ):
    """ Write the file in biom format """
        
    try:
        import biom
    except ImportError:
        sys.exit("Could not find the biom software."+
            " This software is required since the input file is a biom file.")
            
    try:
        import numpy
    except ImportError:
        sys.exit("Could not find the numpy software."+
            " This software is required since the input file is a biom file.")
            
    try:
        import h5py
    except ImportError:
        sys.exit("Could not find the h5py software."+
            " This software is required since the input file is a biom file.")
            
    # reformat the rows into a biom table
    samples=next(rows)[1:]
    ids=[]
    data=[]
    for row in rows:
        ids.append(row[0])
        data.append(row[1:])
            
    table=biom.Table(numpy.array(data), ids, samples)
        
    # write a h5py biom table
    with h5py.File(path, 'w') as file_handle:
        table.to_hdf5(file_handle, "humann utility script")
        

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

def try_zip_open( path, write=None ):
    """ 
    open an uncompressed or gzipped file; fail gracefully 
    """
    fh = None

    # set the open mode
    if write:
        open_mode = "w"
    elif path.endswith(".bz2"):
        open_mode = "rt"
    else:
        open_mode = "rt"

    try:
        if path.endswith(".gz"):
            fh = gzip.open( path, open_mode )
        elif path.endswith(".bz2"):
            fh = bz2.open( path, open_mode )
        else:
            fh = open( path, open_mode )
    except EnvironmentError:
        sys.exit( "Problem opening file: " + path)
    return fh

def read_biom_table( path ):
    """
    return the lines in the biom file
    """

    try:
        import biom
    except ImportError:
        sys.exit("Could not find the biom software."+
            " This software is required since the input file is a biom file.")
        
    try:
        tsv_table = biom.load_table( path ).to_tsv().split("\n")
    except (EnvironmentError, TypeError):
        sys.exit("ERROR: Unable to read biom input file.")
        
    return tsv_table

def gzip_bzip2_biom_open_readlines( path ):
    """
    return the lines in the opened file for tab delimited text, gzip, bzip2 and biom files
    """

    # if the file is biom, convert to text and return lines
    if path.endswith(BIOM_FILE_EXTENSION):
        for line in read_biom_table(path):
            yield line
    else:
        with try_zip_open( path ) as file_handle:
            for line in file_handle:
                yield line.rstrip()

def load_polymap( path, start=0, skip=None, allowed_keys=None, allowed_values=None ):
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
    print( "Loading mapping file from:", path, file=sys.stderr )
    size_warn( path )
    for line in gzip_bzip2_biom_open_readlines( path ):
        row = line.split("\t")
        key = row[start]
        if allowed_keys is None or key in allowed_keys:
            for i, value in enumerate( row ):
                if i != start and (skip is None or i not in skip):
                    if allowed_values is None or value in allowed_values:
                        polymap.setdefault( key, {} )[value] = 1
    return polymap

def fsplit( feature ):
    """
    Expected format is:
      CODE: NAME|STRATUM
    As in:
      GENE123: tRNA transferase|s__Bacteroides_dorei
    Robust to names containing ": ":
      GENE123: tRNA transferase: proline|s__Bacteroides_dorei
    Not robust to names containing "|" (lethal error):
      GENE123: tRNA transferase|proline|s__Bacteroides_dorei
    """
    code = None
    name = None
    stratum = None
    # extract stratum
    items = feature.split( c_strat_delim )
    if len( items ) == 1:
        # looks like an unstratified feature
        code = items[0]
    elif len( items ) == 2:
        # stratified feature
        code = items[0]
        stratum = items[1]
    else:
        # possible names with "|"
        sys.exit( "LETHAL ERROR: bad feature name: {}".format( f ) )
    # extract name
    items = code.split( c_name_delim )
    if len( items ) == 1:
        code = items[0]
    elif len( items ) == 2:
        code = items[0]
        name = items[1]
    elif len( items ) > 2:
        # possible name containing ": "
        code = items[0]
        name = c_name_delim.join( items[1:] )
    return code, name, stratum

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

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# new code from the tools branch
# ---------------------------------------------------------------
# ---------------------------------------------------------------

def say( *args ):
    print( " ".join( [str( k ) for k in args] ), file=sys.stderr )
    
def die( *args ):
    args = ["LETHAL ERROR:"] + list( args )
    say( *args )
    sys.exit( "EXITED" )

def attach_common_arguments( parser, no_output=False ):
    parser.add_argument(
        "-i", "--input",
        default=None,
        metavar="<path>",
        help="HUMAnN table (.tsv or .biom format) [stdin]",
        )
    parser.add_argument(
        "-l", "--last-metadata",
        default=None,
        metavar="<row>",
        help="The name (header) of the last row containing metadata, if any [none]",
        )
    # some scripts have custom output
    if not no_output:
        parser.add_argument(
            "-o", "--output",
            default=None,
            metavar="<path>",
            help="Path for output table [stdout]",
        )
    return None

def pretty_grid( items, cols=3, desc="Please select one of these options:" ):
    padding = 2 + max( [len( k ) for k in items] )
    counter = 0
    desc += "\n"
    desc += "-" * (len( desc ) - 1) + "\n"
    for k in items:
        desc += k.ljust( padding )
        counter += 1
        if counter == cols:
            desc += "\n"
            counter = 0
    return desc

def wrap( text, width=85 ):
    # remove flanking whitespace
    text = text.strip( )
    lines = text.split( "\n" )
    # title
    rule = "=" * width
    newlines = [rule, lines[0], rule, "\n"]
    newline = ""
    # description lines
    for line in lines[2:]:
        line = line.strip( )
        if line == "":
            newlines.append( newline )
            newlines.append( "\n" )
            newline = ""
            continue
        words = line.split( )
        for word in words:
            if len( word ) > width:
                newlines.append( newline )
                newlines.append( word )
                newline = ""
            elif len( newline + " " + word ) > width:
                newlines.append( newline )
                newline = word
            else:
                newline += (" " if newline != "" else "") + word
    if len( newline ) > 0:
        newlines.append( newline )
    newlines += ["\n", rule]
    return "\n".join( [k for k in newlines] )