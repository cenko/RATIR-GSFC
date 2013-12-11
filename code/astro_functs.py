"""
Written by John Capone (jicapone@astro.umd.edu).

Purpose:
	- this code is primarily used for the constants found below.  it may be better to put these in a configuration file rather than in this python script.

Notes:
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




