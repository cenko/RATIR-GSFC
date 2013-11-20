#! /usr/bin/env python

import astropy.io.fits as f
import numpy as np

##########

C1SAT = 16384.0
C1MIN = 2000.0
C1MAX = 13000.0

###########

def ratir_bias(images, camera="C1", outpre="Bias"):

	"""Construct a bias frame from a series of RATIR 0 s images."""
	
	barray = np.array([f.open(im)[0].data for im in images])

	# Median combination (no outlier rejection for the time being
	bmed = np.median(barray, axis=0)
	
	# Create new frame (using header from old one)
	im = f.open(images[0])
	im[0].data = bmed
	im[0].header["BITPIX"] = -32
	im[0].header["BZERO"] = 0
	im.writeto("%s-%s.fits" % (outpre, camera))
	
	return
	
##########

def ratir_flat(images, camera="C1", outpre="Flat", bias="Bias-C1.fits",
			   maximum=C1MAX, minimum=C1MIN, fkey="FILTER"):

	"""Construct a flat from a series of sky flats."""
	
	goodflats = []
	
	# Bias subtraction
	bdata = f.open(bias)[0].data

	for im in images:
		new = f.open(im)
		new[0].data -= bdata
			
		# Check for good values	
		newmed = np.median(new[0].data)
		if (newmed < maximum) and (newmed > minimum):
			goodflats.append(new[0].data / newmed)
			
	# Median combine flats
	gfarray = np.array(goodflats)
	gfmed = np.median(gfarray, axis=0)
	
	# Write results
	im = f.open(images[0])
	im[0].data = gfmed
	im[0].header["BITPIX"] = -32
	im[0].header["BZERO"] = 0
	im.writeto("%s-%s-%s.fits" % (outpre, camera, im[0].header[fkey]))
	
	return
	
##########

def ratir_detrend(images, bias="Bias-C1.fits", flat="Flat-C1-i.fits", 
                  outpre="d"):
                  
	"""Detrend RATIR data (bias subtract and flat field."""
	
	bdata = f.open(bias)[0].data
	fdata = f.open(flat)[0].data
	
	for im in images:
	
		new = f.open(im)
		new[0].data -= bdata
		new[0].data /= fdata
		new[0].header["BITPIX"] = -32
		new[0].header["BZERO"] = 0
		new.writeto("%s%s" % (outpre, im))
		
	return
	

	
	
	