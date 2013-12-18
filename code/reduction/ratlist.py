"""
Translated from ratlist.sh by John Capone (jicapone@astro.umd.edu).

Notes:
	- added option to specify working directory
"""

import os
import fnmatch

"""
Purpose:			To create lists of detected fits files in each subdir of the specified working dir. Each camera has it's own list file.
Input:
	workdir:		directory with subdirectories containing data (ex. visit numbers)
	cams:			lists created for specified camera numbers
"""
def ratlist( workdir='.', cams=[0,1,2,3] ):

	# move to working directory
	start_dir = os.getcwd()
	os.chdir( workdir ) # move into working dir
	
	d = os.getcwd().split('/')[-1] # name of current directory
	print "* * * creating lists of frames in {} * * *".format(d)
	fns = []
	for i in cams:
		fns.append( open( 'C{}_{}.list'.format( i, d ), 'w' ) ) # create list file for each camera
	for f in os.listdir( '.' ):
		for i in range(len(cams)):
			if fnmatch.fnmatch( f, '*C{}*.fits'.format(cams[i]) ):
				fns[i].write( f + '\n' ) # add detected file to correct list file
	os.chdir( start_dir ) # move back to starting dir
