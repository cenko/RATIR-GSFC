"""
Translated from icoords.pro by John Capone (jicapone@astro.umd.edu).

Notes:
	selection of x and y limits for differenct filters is a mess.  should automate
"""

import numpy as np
import os
import fnmatch
import astropy.io.fits as pf
import pylab as pl
import scipy as sp
from astropy import wcs

# find index of xarr and yarr with minimum RSS distance from x and y
#   Note: needs updated implementation optimized for numpy
def nearest( x, y, xarr, yarr, mindist ):
	index = -1
	imindist = mindist
	if np.shape(xarr) != np.shape(yarr):
		raise plotratirError( "xarr and yarr must have equal dimensions" )	
	for i in range( np.size(xarr) ):
		dist = np.sqrt( (x - xarr[i])**2. + (y - yarr[i])**2. )
		if dist < imindist:
			imindist = dist
			index = i
	return index

# make a circle for identifying sources in images
def circle( xcenter, ycenter, radius ):
	points = np.linspace( 0., 2.*np.pi, 100 )
	x = xcenter + radius * np.cos(points)
	y = ycenter + radius * np.sin(points)
	return np.transpose([x,y])

# Define a number of global variables
filters = ['r','i','z','y','J','H']

# returns files in directory "loc" which start with prefix and end with postfix
def get_files( selection, loc='.' ):
	matches = []
	for files in os.listdir(loc):
		if fnmatch.fnmatch( files, selection ):
			matches.append(files)
	return matches

# Identify files
prefchar = 'coadd'
wildcharimg = '?????-????-?????_?'
zffiles = get_files( prefchar + wildcharimg + '.fits' )
weightfiles = get_files( prefchar + wildcharimg + '.weight.fits' )

# hextract IDL function - Extract a subimage from an array and update astrometry in FITS header
#	Note: fits xy convention opposite of python
def hextract( oldim, oldhd, x0, x1, y0, y1 ):
	newim = np.copy( oldim[y0:y1,x0:x1] )
	newhd = oldhd.copy()
	# new reference pixels are center pixels of sub-image
	refpx1 = np.float(y0 + y1)/2.
	refpx2 = np.float(x0 + x1)/2.
	# Parse the WCS keywords in the primary HDU
	w = wcs.WCS(newhd)
	pixcrd = np.array( [[refpx1, refpx2]] )
	world = w.wcs_pix2world( pixcrd, 1 )
	# update header astrometry
	newhd.update( 'CRPIX1', pixcrd[0,0] )
	newhd.update( 'CRPIX2', pixcrd[0,1] )
	newhd.update( 'CRVAL1', world[0,0] )
	newhd.update( 'CRVAL2', world[0,1] )
	# update header array size
	newhd.update( 'NAXIS1', newim.shape[1] )
	newhd.update( 'NAXIS2', newim.shape[0] )
	return newim, newhd

def writefits( ofile, img, header ):
	hdu = pf.PrimaryHDU(img)
	hdu.header = header
	hdulist = pf.HDUList([hdu])
	if os.path.exists(ofile):
		os.remove(ofile)
	hdulist.writeto(ofile)

def icoords():
	
	# perform initial crop to remove noisy edges
	for i in range(np.size(zffiles)):
		hdulist = pf.open(zffiles[i])
		h = hdulist[0].header
		img = hdulist[0].data
		cfilter = zffiles[i].split('_')[1].split('.')[0] # extract filter label from file name
		x1 = x2 = y1 = y2 = 0
		if cfilter == 'H' and False:	# case is commented out in IDL code --> always False here
			x1 = 550.
			x2 = 1250.
			y1 = 1500.
			y2 = 3100.
		elif cfilter == 'i' and True:
			x1 = 1103.
			x2 = 2067.
			y1 = 238.
			y2 = 1160.
		elif cfilter == 'r' and True:
			x1 = 1090.
			x2 = 2070.
			y1 = 75.
			y2 = 1020.
		elif cfilter == 'J' and False:	# case is commented out in IDL code --> always False here
			x1 = 700.
			x2 = 1500.
			y1 = 200.
			y2 = 1800.
		elif cfilter == 'Y' and False:	# case is commented out in IDL code --> always False here
			x1 = 170.
			x2 = 900.
			y1 = 200.
			y2 = 1800.
		elif cfilter == 'Z' and False:	# case is commented out in IDL code --> always False here
			x1 = 125.
			x2 = 975.
			y1 = 150.
			y2 = 1950.
		img, h = hextract( img, h, x1, x2, y1, y2 )
		ofile = zffiles[i].split('.')[0] + '.crop.fits'
		writefits( ofile, img, h )
		hdulist = pf.open(weightfiles[i])
		wh = hdulist[0].header
		wimg = hdulist[0].data
		wimg, wh = hextract( wimg, wh, x1, x2, y1, y2 )
		ofile = weightfiles[i].split('.')[0] + '.crop.weight.fits'
		writefits( ofile, wimg, wh )

	# Resample all images using SWarp to a reference image called multicolor
	for i in range(np.size(zffiles)):
		ifile = zffiles[i].split('.')[0] + '.crop.fits'
		swarpstr = ''
		if i == 0:
			swarpstr = '"' + ifile + '"'
		else:
			swarpstr = swarpstr + ' "' + ifile + '"'

	stackcmd = 'swarp ' + swarpstr + ' -DELETE_TMPFILES N'
	stackcmd += ' -IMAGEOUT_NAME multicolor.fits -WEIGHTOUT_NAME multicolor.weight.fits'
	os.system( stackcmd )

	# Rename all the resampled files
	for i in range(np.size(zffiles)):
		tmp = zffiles[i].split('.')[0]
		ifile = tmp +'.crop.resamp.fits'
		ofile = tmp + '.crop.fits'
		mvcmd = 'mv -f ' + ifile + ' ' + ofile
		os.system( mvcmd )
		ifile = tmp + '.crop.resamp.weight.fits'
		ofile = tmp + '.crop.weight.fits'
		mvcmd = 'mv -f ' + ifile + ' ' + ofile
		os.system( mvcmd )

	# run sextractor on pipeline reduced files to identify point sources
	for i in range(np.size(zffiles)):
		cfilter = zffiles[i].split('_')[1].split('.')[0] # extract filter label from file name
		if cfilter == 'Z' or cfilter == 'Y':
			cfilter = cfilter.lower()
		# run SExtractor on each image
		ifile = zffiles[i].split('.')[0] + '.crop.fits'
		cmd = 'sex ' + ifile + ' -c "ratir_nir.sex"'
		print cmd
		os.system( cmd )
		# read in results from sextractor and produce IRAF coordinate file
		x,y,ra,dec,mag,magerr,e,fwhm = np.loadtxt( 'temp.cat', unpack=True )
		f = open( 'coords' + cfilter, 'w' )
		for j in range(np.size(ra)):
			f.write( '{:15.6f}{:15.6f}'.format( ra[j], dec[j] ) )
		f.close()
		cmd = 'mv -f temp.cat fluxes1_' + cfilter + '.txt'
		print cmd
		os.system( cmd )
		fwhm = np.median(fwhm)
		f = open( cfilter + '.aannulus', 'w' )
		f.write( '{:15.6f}'.format( fwhm ) )
		f.close()
