#!/usr/bin/env python
"""
Utilities relating to executing third party software, file permissions,
and file formats
"""

import os, sys, subprocess, re

def file_exists_readable(file):
	"""
	Exit with error if file does not exist or is not readable
	"""
	if not os.path.isfile(file):
		sys.exit("ERROR: Can not find file "+ file)
		
	if not os.access(file, os.R_OK):
		sys.exit("ERROR: Not able to read file " + file)


def find_exe_in_path(exe, user_path):
	"""
	Check that an executable exists in the provided path
	
	Path provided should be of the same format as $PATH
	"""
	paths = user_path.split(os.pathsep)
	for path in paths:
		fullexe = os.path.join(path,exe)
		if os.path.exists(fullexe):
			if os.access(fullexe,os.X_OK):
				return True
	return False	

def execute_software(exe, args, infiles, outfiles):
	"""
	Execute third party software
	"""
	
	# determine which path to check for the software
	# if the exe is a .py, check the python path
	path = os.environ["PATH"]
	if os.path.splitext(exe)[1] in [".py"]:
		path = os.pathsep.join(sys.path)
	
	# check that the executable can be found
	if not find_exe_in_path(exe, path):
		sys.exit("ERROR: Can not find executable " + exe)
	
	# check that the files exist and are readable
	for file in infiles:
		file_exists_readable(file)
		
	if args:
		cmd=[exe,args]
	else:
		cmd=[exe]
	
	print "Running: " + exe + " " + args + " ........"
	
	try:
		p = subprocess.call(cmd)
	
	except OSError as e:
		sys.exit("Error: Problem executing " + exe + " " + args + "\n" +
			e.strerror)
		
	# check that the output files exist and are readable
	for file in outfiles:
		file_exists_readable(file)
	
def fasta_or_fastq(file):
	"""
	Check to see if a file is of fasta or fastq format
	
	Fastq format short example:
	@SEQ_ID
	GATCTGG
	+
	!''****
	
	Fasta format short example:
	>SEQ_INFO
	GATCTGG
	
	Returns error if not of fasta or fastq format
	"""
	
	format="error"
	
	# check file exists
	file_exists_readable(file)
	
	# read in first 4 lines of file to check format
	file_handle = open(file, "r")
	
	first_line = file_handle.readline()
	second_line = file_handle.readline()
	
	# check that second line is only nucleotides
	if re.search("^[aATtcCgGN]+$", second_line):
		# check first line to determine fasta or fastq format
		if re.search("^@",first_line):
			format="fastq"		
		if re.search("^>",first_line):
			format="fasta"
			
	return format
	
