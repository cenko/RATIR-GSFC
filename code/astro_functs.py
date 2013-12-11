"""
Written by John Capone (jicapone@astro.umd.edu).

Purpose:
	- this code is primarily used for the constants found below.  it may be better to put these in a configuration file rather than in this python script.

Notes:
	- the meadian_clip and image combining functions are not working well.  numpy masked arrays are slow.  there is definitely a better way to do this.
	- is numpy mean() buggy? found cases where mean(array) < min(array)
"""

import astropy.io.fits as pf
import numpy as np
import time

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

"""
Purpose:	load a ratir fits file
Input:
	fn:		filename
	clean:	return only those header keywords specified in the function (optional)
Output:
	header:	FITS header
	data:	FITS data
"""
def get_fits( fn, clean=False ):
	hdulist = pf.open( fn )
	if len(hdulist) != 1:
		print "Warning: HDU list contains more than one HDU.  Only the first HDU will be returned."
	header = hdulist[0].header
	data = hdulist[0].data
	if clean:
		ik_C0 = [	'JD', # exposure JD
					'EXPTIME', # exposure time in s
					'SCRC1', # ccd temperature in K
					'PRPSLID', # proposal id number
					'VSTID', # id number of the visit
					'FILTER', # camera filter
					'SOFTGAIN', # software gain of image (value image is divided by to reduce size)
					'AVERAGE', # pre-software gain average
					'STDEV'] # pre-software gain std
		ik_C1 = [	'JD', # exposure JD
					'EXPTIME', # exposure time in s
					'SCRC2', # ccd temperature in K
					'PRPSLID', # proposal id number
					'VSTID', # id number of the visit
					'FILTER', # camera filter
					'SOFTGAIN', # software gain of image (value image is divided by to reduce size)
					'AVERAGE', # pre-software gain average
					'STDEV'] # pre-software gain std
		ik_C2 = [	'ACQTIME', # exposure JD
					'EXPTIME', # exposure time in s
					'SCRC3', # H2RG temperature in K
					'PRPSLID', # proposal id number
					'VSTID', # id number of the visit
					'SOFTGAIN', # software gain of image (value image is divided by to reduce size)
					'AVERAGE', # pre-software gain average
					'STDEV'] # pre-software gain std
		ik_C3 = [	'ACQTIME', # exposure JD
					'EXPTIME', # exposure time in s
					'SCRC4', # H2RG temperature in K
					'PRPSLID', # proposal id number
					'VSTID', # id number of the visit
					'SOFTGAIN', # software gain of image (value image is divided by to reduce size)
					'AVERAGE', # pre-software gain average
					'STDEV'] # pre-software gain std
		ik = ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2','EXTEND'] # basic keys which are required
		if fn.find('C0') != -1:
			ik += ik_C0
		elif fn.find('C1') != -1:
			ik += ik_C1
		elif fn.find('C2') != -1:
			ik += ik_C2
		elif fn.find('C3') != -1:
			ik += ik_C3
		else:
			print "Warning: detector type not found in filename."
		htemp = pf.Header()
		cards = header.cards
		for c in cards:
			if ik.count( c[0] ) != 0:
				htemp.append(c)
		header = htemp
	return header, data

"""
Purpose:	iterative median sigma filter.  if axis=0 for a stack of images, then the clipping will be done pixel by pixel.
Input:
	indata:	data set to be filtered
******************* Not finished *******************
"""
def median_clip( indata, clipsig=5., maxiter=10, converge_num=0.01, axis=None ):
	data = np.copy( indata )
	i = 0; c1 = 1.0 ; c2 = 0.0
	clipped = np.zeros( data.shape, dtype=bool )
	ct = data.size
	while ( ct != 0 ) and ( c1 >= c2 ) and ( i < maxiter ):
		lastct = ct
		medianval = np.median( data, axis=axis )
		sigs = np.std( data, axis=axis )
		clipped |= np.abs( data - medianval ) > clipsig * sigs
		ct = np.sum( clipped )
		if ct > 0:
			data[clipped] = medianval # set outliers to median value
		c1 = np.abs(ct - lastct)
		c2 = converge_num * lastct
		i += 1
	return data, clipped

"""
Purpose:		combine stack of frames using mean
Input:
	indata:		stack of frames to be combined
	type:		function used to combine the stack.  currently only mean or median
Output:
	combined:	combined stack
	sigma:		standard deviation of each pixel
"""
def imcombine( indata, type='mean' ):
	if indata.ndim != 3:
		print "Warning: data should be 3D stack of frames."
	if type(indata) is np.ma.core.MaskedArray:
		if type is 'mean':
			combined = np.ma.mean( indata, axis=0 )
		else:
			combined = np.ma.median( indata, axis=0 )
		sigma = np.ma.std( indata, axis=0 )
	else:
		if type is 'mean':
			combined = np.mean( indata, axis=0 )
		else:
			combined = np.median( indata, axis=0 )
		sigma = np.std( indata, axis=0 )
	return combined, sigma

"""
Purpose:		mean add stack of frames with median filter
******************* Not finished *******************
"""
def clip_combine( indata, clipsig=5.0, maxiter=10, converge_num=0.01 ):
	if indata.ndim != 3:
		print "Warning: data should be 3D stack of frames."
	data = np.ma.zeros( indata.shape )
	for i in range(indata.shape[0]):
		data[i] = median_clip( indata[i], clipsig, maxiter, converge_num )
	return imcombine( data )

"""
* converted from IDL ROBUST_SIGMA function
Purpose:		Calculate a resistant estimate of the dispersion of a distribution. For an uncontaminated distribution, this is identical to the standard deviation.
Input:
	y:			Vector of quantity for which the dispersion is to be calculated
	zero:		if set, the dispersion is calculated w.r.t. 0.0 rather than the central value of the vector. If Y is a vector of residuals, this should be set.
Output:		robust_sigma returns the dispersion. In case of failure, returns value of -1.0
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




