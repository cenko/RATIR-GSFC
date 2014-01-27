"""
	Purpose:	repository for general purpose functions used by the RATIR pipeline.  also contains constants (it may be better to put these in a configuration file rather than in this python script).

	Usage:		
		*)	most of these functions and constants are usually called by other scripts, not by the user
		1)	enter python or ipython environment
		2)	can load all functions using:
			- "from rat_preproc import *" if you want to call functions using just the function's name
			- "import rat_preproc as rp" if you want to call functions using rp.function_name( args )

	Notes:
		- is numpy mean() buggy? found cases where mean(array) < min(array)
		- added show_list() to aid users in reviewing results

	Future Improvements:
		- 

"""

import os
import astropy.io.fits as pf
import matplotlib.pylab as pl
import numpy as np
import time
from scipy.ndimage.interpolation import zoom

# CONSTANTS
PROPOSALS = { 'NIRstandard': '0000', 'OPTstandard': '00001', 'cluster': '0002', 'galaxy': '0003', 'blank': '0004', 'pointing': '0005', 'bias': '0006', 'dark': '0007', 'flat': '0008', 'focus': '0009', 'misc': '0010', 'GRB': '1000' } # 2012 proposal names and id numbers. source: rsync://ratir.astroscu.unam.mx/public/proposalidentifiers.txt
CAM_NAMES = [ 'C0', 'C1', 'C2', 'C3' ] # RATIR camera names
CAM_ROTAT = [ 0, 0, 1, 1 ] # frames are rotated by value * 90 degrees
H2RG_FILTERS = [ 'Z', 'J', 'Y', 'H' ] # RATIR NIR bands.  0+2 are C2, 1+3 are C3
# RATIR H2RG filter slices
Z_SLICE = np.s_[1:1700,1:900]
Y_SLICE = np.s_[1:1700,1144:2043]
J_SLICE = np.s_[1:1700,1:900]
H_SLICE = np.s_[1:1700,1144:2043]
H2RG_SLICES = [ Z_SLICE, J_SLICE, Y_SLICE, H_SLICE ] # same order as H2RG_FILTERS.  0+2 are C2, 1+3 are C3
OBJ_NAME = 'img' # designator for object frames
SKY_NAME = 'sky' # designator for sky frames
FLAT_NAME = 'flat' # designator for flat frames
BIAS_NAME = 'bias' # designator for bias frames
CONFIG_LOCATION = 'astro_functs.py' # name of file containing configuration information, currently this file.
CAM_GAIN = [ lambda SOFTGAIN: 16.80/SOFTGAIN, lambda SOFTGAIN: 18.64/SOFTGAIN, lambda SOFTGAIN: 2.2/SOFTGAIN, lambda SOFTGAIN: 2.4/SOFTGAIN ] # gain of each camera as a function of the SOFTGAIN keyword extracted from a frame's header
CAM_SATUR = [ lambda SOFTGAIN: (2.**16/SOFTGAIN)-1, lambda SOFTGAIN: (2.**16/SOFTGAIN)-1, lambda SOFTGAIN: (32000./SOFTGAIN)-1, lambda SOFTGAIN: (32000./SOFTGAIN)-1 ] # saturation levels for each detector in DNs as a function of the SOFTGAIN keyword extracted from a frame's header
CENTER_KEY = 'STRRQAP' # RATIR header keyword specifying which H2RG filters the target is focused on

"""
	Written by John Capone (jicapone@astro.umd.edu).

	Purpose:		combine stack of frames

	Input:
		indata:		stack of frames to be combined
		type:		function used to combine the stack.  currently only mean or median
		ret_std:	if set to true, return the standard deviation of each pixel.  default is false
	
	Output:
		combined:	combined stack
		sigma:		standard deviation of each pixel (optional)

	Notes:
		- 

	Future Improvements:
		- add outlier rejection
"""
def imcombine( indata, type='median', ret_std=False ):
	if indata.ndim != 3:
		print "Warning: data should be 3D stack of frames."
	if type is 'mean':
		combined = np.mean( indata, axis=0 )
	else:
		combined = np.median( indata, axis=0 )
	if ret_std:
		sigma = np.std( indata, axis=0 )
		return combined, sigma
	else:
		return combined

"""
	Converted from IDL ROBUST_SIGMA function by John Capone (jicapone@astro.umd.edu).

	Purpose:		Calculate a resistant estimate of the dispersion of a distribution. For an uncontaminated distribution, this is identical to the standard deviation.

	Input:
		y:			Vector of quantity for which the dispersion is to be calculated
		zero:		if set, the dispersion is calculated w.r.t. 0.0 rather than the central value of the vector. If Y is a vector of residuals, this should be set.
	
	Output:			robust_sigma returns the dispersion. In case of failure, returns value of -1.0

	Notes:
		- 
"""
def robust_sigma( y, zero=False ):
	eps = 1.0e-20
	if zero:
		y0 = 0.
	else:
		y0 = np.median(y)
	# first, the median absolute deviation about the median:
	mad = np.median( np.abs(y-y0) )/.6745
	# if the mad=0, try the mean absolute deviation:
	if mad < eps:
		mad = np.average( np.abs(y-y0) )/.80
	if mad < eps:
		sigma = 0.
		return sigma
	# now the biweighted value:
	u   = (y-y0)/(6.*mad)
	uu  = u*u
	q = uu <= 1.0
	count = np.sum(q)
	if count < 3:
		print 'robust_sigma: this distribution is too weird! returning -1'
		sigma = -1.
		return sigma
	numerator = np.sum( (y[q] - y0)**2. * (1 - uu[q])**4. )
	n = y.size
	den1 = np.sum( (1. - uu[q]) * (1. - 5.*uu[q]) )
	sigma = n * numerator / (den1 * (den1 - 1.))
	if sigma > 0.:
		sigma = np.sqrt( sigma )
	else:
		sigma = 0.
	return sigma

"""
	Written by John Capone (jicapone@astro.umd.edu).

	Purpose:		displays images in specified list file.

	Input:
		list_fn:	file name of list file
		nx:			number of images to display simultaneously in x
		ny:			number of images to display simultaneously in y
		size_mult:	multiple to determine image sizes
		zoom_lvl:	amount to decrease image resolution by (to save memory)
		fontsize:	fontsize for plots

	Usage:
		1)	enter python or ipython environment
		2)	load function -> 'from rat_preproc import show_list'
		3)	run function -> 'show_list( list_fn='path/to/list/file.list' )'
			- decreasing zoom_lvl (i.e. from 0.5 to 0.1) decreases the size of the displayed image, thus decreasing the amount of memory required
		4)	function will display arrays of images in list file for inspection by user

	Notes:
		- added escape character
		- user can now change font size

	Future Improvements:
		* get data directory from list_fn
"""
def show_list( list_fn, nx=5, ny=3, size_mult=3.2, zoom_lvl=0.5, fontsize=8 ):

	nx = int(nx); ny = int(ny) # force parameter types to int

	try:
		fin = open( list_fn, 'r' ) # open list of FITs files
	except IOError:
		print "Error: {} not found.  Exiting...".format( list_fn )
		return
	
	# move to working directory
	start_dir = os.getcwd()
	workdir = os.path.split( list_fn )[0]
	if workdir == '':
		workdir = '.'
	os.chdir( workdir )

	fits_fns = fin.readlines() # read file names from list
	fin.close() # close files
	nfits = len(fits_fns) # number of fits files listed

	pl.ion() # pylab in interactive mode

	# create figures of subplots for review
	nfigs = int( np.ceil( nfits / float( nx * ny ) ) )
	for i in range( nfigs ):

		start_fits = i*nx*ny
		if (i + 1)*nx*ny <= nfits:
			stop_fits = (i + 1)*nx*ny - 1
			nsubplts = nx*ny
		else:
			stop_fits = nfits
			nsubplts = nfits - start_fits
		
		pl.figure( "FITS {} - {}".format( start_fits, stop_fits ), figsize=(nx*size_mult,ny*size_mult), tight_layout=True ) # create new figure

		# display image in each subplot
		for j in range( nsubplts ):

			print "{:<4}:{:<4}".format( i, start_fits + j )

			pl.subplot( ny, nx, j+1 ) # new subplot

			fits_fn = fits_fns[start_fits + j].rstrip() # current fits file name with return removed
			temp = fits_fn.split('.')
			if len( temp ) == 1:
				fits_fn += '.fits'
			elif len( temp ) == 2:
				if temp[1].lower() != 'fits':
					print "Error: invalid \"{}\" file type detected.  Should be \"fits\" file type. Exiting...".format( temp[1] )
					os.chdir( start_dir ) # move back to starting directory
					return
			else:
				print "Error: file names should not include \".\" except for file type extention.  Exiting..."
				os.chdir( start_dir ) # move back to starting directory
				return
			fits_id = temp[0] # fits file name with extention removed
			
			# open data
			hdulist = pf.open( fits_fn )
			im = hdulist[0].data
			h = hdulist[0].header
			hdulist.close() # close FITs file

			# image statistics
			m = np.median( im )
			s = robust_sigma( im )

			# display
			imdisp = zoom( im, zoom_lvl )
			axim = pl.imshow( imdisp, vmin=m-5*s, vmax=m+5*s, origin='lower', cmap=pl.cm.gray )
			axim.get_axes().get_xaxis().set_ticks([]) # remove numbers from x axis
			axim.get_axes().get_yaxis().set_ticks([]) # remove numbers from y axis
			pl.title( "{} - {} filter".format( fits_id, h['FILTER'] ), fontsize=fontsize ) # title with identifier
			pl.xlabel( "Median: {}".format( m ), fontsize=fontsize )

		usr_select = raw_input( "Press any key to continue or \"q\" to quit: " ) # prompt user to continue
		pl.close('all') # close image to free memory
		if usr_select.lower() == 'q':
			print "Exiting..."
			os.chdir( start_dir ) # move back to starting directory
			return

	os.chdir( start_dir ) # move back to starting directory
