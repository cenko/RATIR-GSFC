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
#CAM_ROTAT = [ 0, 0, 1, 1 ] # frames are rotated by value * 90 degrees
#CAM_XFLIP = [ False, False, False, False ] # frames are flipped in x
#CAM_YFLIP = [ False, False, False, True ] # frames are flipped in y
CAM_PXSCALE = [0.32, 0.32, 0.3, 0.3] # C0, C1, C2, C3 in arcsec/px
SPLIT_FILTERS = [ 'Z', 'J', 'Y', 'H' ] # RATIR NIR bands.  0+2 are C2, 1+3 are C3
CAM_BIAS   = [True, True, False, False]

FILTERS_CAM = [1, 1, 2, 2]
# RATIR H2RG filter slices
Z_SLICE = np.s_[4:975,100:2000]
Y_SLICE = np.s_[1135:2043,240:2043]
J_SLICE = np.s_[50:1000,4:2000]
H_SLICE = np.s_[1200:2043,4:1940]
#Z_SLICE = np.s_[1:1700,1:900]
#Y_SLICE = np.s_[1:1700,1144:2043]
#J_SLICE = np.s_[1:1700,1:900]
#H_SLICE = np.s_[1:1700,1144:2043]
H2RG_SLICES = [ Z_SLICE, J_SLICE, Y_SLICE, H_SLICE ] # same order as H2RG_FILTERS.  0+2 are C2, 1+3 are C3
OBJ_NAME = 'img' # designator for object frames
#SKY_NAME = 'sky' # designator for sky frames
FLAT_NAME = 'flat' # designator for flat frames
BIAS_NAME = 'bias' # designator for bias frames
CONFIG_LOCATION = 'astro_functs.py' # name of file containing configuration information, currently this file.

CAM_WAVE  = ['OPT', 'OPT', 'IR', 'IR']
CAM_SPLIT = [False, False, True, True]
#WCS relevant parameters (RATIR H2RGs have barrel distortions)
a = -19.60381671
b = -4128.15179797
CAM_SECPIX1  = [0.3168, 0.3171, 0.2988, 0.2983]
CAM_SECPIX2  = [0.3171, 0.3191, -0.2955, 0.2945]
CAM_THETA    = [0.60, 2.40, 92.55, -89.25]
CAM_X0       = [513, 513, 1177, 924]
CAM_Y0		 = [513, 513, 1031, 982]

H2RG_ASTR = {'PV1_1': 1.0, 'PV2_1': 1.0, 'PV1_17':a, 'PV2_17':a, 'PV1_19': 2.0*a, 'PV2_19':2.0*a, 'PV1_21':a, 'PV2_21':a, 'PV1_31':b, 'PV2_31': b, 'PV1_33':3.0*b, 'PV2_33':3.0*b, 'PV1_35':3.0*b, 'PV2_35':3.0*b, 'PV1_37':b, 'PV2_37': b}
RA_KEY = 'ETRRQRA'
DEC_KEY = 'ETRRQDE'


CAM_GAIN = [ lambda SOFTGAIN: 16.80/SOFTGAIN, lambda SOFTGAIN: 18.64/SOFTGAIN, lambda SOFTGAIN: 2.2/SOFTGAIN, lambda SOFTGAIN: 2.4/SOFTGAIN ] # gain of each camera as a function of the SOFTGAIN keyword extracted from a frame's header
CAM_SATUR = [ lambda SOFTGAIN: (2.**16/SOFTGAIN)-1, lambda SOFTGAIN: (2.**16/SOFTGAIN)-1, lambda SOFTGAIN: (36000./SOFTGAIN)-1, lambda SOFTGAIN: (36000./SOFTGAIN)-1 ] # saturation levels for each detector in DNs as a function of the SOFTGAIN keyword extracted from a frame's header
CENTER_KEY = 'STRRQAP' # RATIR header keyword specifying which H2RG filters the target is focused on
# frame corners in arcmin offset from center.  top-left, bottom-left, bottom-right, top-right.  top==north, left==east
CAMOFFS = np.array([[[2.785,2.632], [2.604,-2.775], [-2.800,-2.615], [-2.635,2.789]],	# C0 corner offsets in arcmin
					[[2.817,2.624], [2.607,-2.818], [-2.807,-2.624], [-2.616,2.818]],	# C1 corner offsets in arcmin
					[[5.012,6.229], [4.569,-3.905], [-0.227,-3.678], [0.228,6.453]],	# C2-Z corner offsets in arcmin
					[[-0.556,6.488],[-1.013,-3.642],[-5.488,-3.445], [-5.018,6.683]],	# C2-Y corner offsets in arcmin
					[[4.834,5.430], [4.720,-4.701], [-0.324,-4.647], [-0.185,5.494]],	# C3-J corner offsets in arcmin
					[[-0.916,5.503],[-1.059,-4.639],[-5.318,-4.594], [-5.154,5.557]]])	# C3-H corner offsets in arcmin
FRAMECENTER = CAMOFFS.mean( axis=1 ) # field centers in arcmin (E,N) offset from center
# Offsets of the pointing apertures east and north in arcmin
APOFFS	= {	"rcenter":		np.array([0,0]),
			"icenter":		np.array([0,0]),
			"ricenter":		np.array([0,0]),
			"riZJcenter":	np.array([1.2,0]),
			"riYHcenter":	np.array([-1.8,0]),
			"ZJcenter":		np.array([2.2,0.7]),
			"YHcenter":		np.array([-3.2,0.7])}
# return aperture center in pixels from frame center
#def APCENTER( apstr, fnum):
#	ao = APOFFS[apstr] # offset of pointing in arcmin
#	fc = FRAMECENTER[fnum] # offset of frame center in arcmin
#	if fnum < 2: ps = CAM_PXSCALE[fnum]
#	elif fnum < 4: ps = CAM_PXSCALE[2]
#	else: ps = CAM_PXSCALE[3]
#	return ((ao - fc)*60./ps).astype(np.int)

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
		- added normalization before combining

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
