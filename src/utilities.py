#!/usr/bin/env python
"""
Utilities relating to executing third party software
"""

import os

def find_exe_in_path(exe):
	"""
	Check that an executable exists in the users path
	"""
	paths = os.environ["PATH"].split(os.pathsep)
	for path in paths:
		fullexe = os.path.join(path,exe)
		if os.path.exists(fullexe):
			if os.access(fullexe,os.X_OK):
				return True
	return False	
