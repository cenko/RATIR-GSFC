"""
Translated from ratlist.sh by John Capone (jicapone@astro.umd.edu).

Notes:
	- added option to specify working directory
"""

import os
import fnmatch

"""
Purpose:	To create lists of detected fits files in each subdir of the specified working dir. Each camera has it's own list file.
Input:
	working_dir:	if walk_dirs=True, this is the directory with subdirectories containing data (ex. visit numbers).  if walk_dirs=False, this is the directory containing the data you would like to list by camera.
"""
def ratlist( working_dir='.', walk_dirs=True ):
	
	# move to working directory
	start_dir = os.getcwd()
	os.chdir( working_dir ) # move into working dir
	
	# get list of directories
	dirs = os.walk( '.' ).next()[1]
	
	# go through FITs files for each directory
	for d in dirs:
		os.chdir( d ) # move into subdir
		fns = []
		for i in range(4):
			fns.append( open( 'C{}_{}.list'.format( i, d ), 'w' ) ) # create list file for each camera
		for f in os.listdir( '.' ):
			for i in range(4):
				if fnmatch.fnmatch( f, '*C{}*.fits'.format(i) ):
					fns[i].write( f + '\n' ) # add detected file to correct list file
		os.chdir( '../' ) # move out of subdir
	os.chdir( start_dir ) # move back to starting dir
