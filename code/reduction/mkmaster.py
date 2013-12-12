"""
Written by John Capone (jicapone@astro.umd.edu).

Notes:
	- 
"""

import os
import shutil
import astropy.io.fits as pf
import numpy as np
import matplotlib.pylab as pl
import astro_functs as af # contains basic functions and RATIR constants

# function to make FITS keywords to store file names of combined frames
FITS_IN_KEY = lambda n: 'IMCMB{:03}'.format(int(n))

"""
Purpose:		bias list comes directly from ratlist, so the list name must be changed to work with mkmaster. this function is called from mkmaster and should not be called by the user.
Input:
	workdir:	directory where function is to be executed
"""
def _prep_bias_list( cams, workdir='.' ):

	# move to working directory
	start_dir = os.getcwd()
	os.chdir( workdir )

	# rename bias list
	for cam in cams:
		shutil.copy( 'C{}_{}.list'.format( cam, os.getcwd().split('/')[-1] ), '{}_{}.list'.format( af.BIAS_NAME, cam ) )

"""
Purpose:		make master bias and master flat frames
				* currently no outlier rejection
Input:
	mtype:		type of master frame. should be either af.FLAT_NAME or af.BIAS_NAME
	workdir:	directory where function is to be executed
	bands:		camera number or H2RG filter name.  all by default
"""
def mkmaster( mtype, workdir='.', bands=['0','1']+af.H2RG_FILTERS ):

	# check for valid mtype
	if mtype not in [af.FLAT_NAME, af.BIAS_NAME]:
		print "Error: valid arguments for mtype are", af.FLAT_NAME, "and", af.BIAS_NAME, ". Exiting..."
		return

	if mtype is af.BIAS_NAME:
		_prep_bias_list( cams=bands, workdir=workdir )

	print "* * * making master {} frame * * *".format( mtype )

	# move to working directory
	start_dir = os.getcwd()
	os.chdir( workdir )

	# work on FITs files for specified cameras
	for band in bands:
		flist = open( '{}_{}.list'.format( mtype, band ), 'r' )
		hdu = pf.PrimaryHDU()
		data_arr = []
		i = 0
		for f in flist.readlines():
			fits_fn = f.rstrip() # current fits file name with return removed
			fits_id = fits_fn.split('.')[0] # fits file name with extention removed
			fn = fits_id + '.fits'
			hdu.header[FITS_IN_KEY(i)] = fn # add flat fn to master flat header
			print fn
			hdulist = pf.open( fn )
			data_arr.append(hdulist[0].data)
			i += 1
		data_arr = np.array( data_arr )
		master = af.imcombine( data_arr, type='median' )
		dtemp = master.astype(np.float)
		hdu.data = dtemp/np.median(dtemp) # scale is median
		hdulist = pf.HDUList( [hdu] )

		hdulist.writeto( '{}_{}.fits'.format( mtype, band ), clobber=True )

	# move back to starting directory
	os.chdir( start_dir )