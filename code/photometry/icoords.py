"""
NAME:
	icoords

PURPOSE:
	Creates a stacked image with all filters (saved to multicolor.fits and multicolor.weight.fits).  
	Then runs sextractor on each image (with weight file) and saves the RA and DEC values to 
	coords(FILTER) and the sextractor output to fluxes1_(FILTER).txt


OUTPUTS:
	multicolor.fits, multicolor.weight.fits - files with all filter images stacked
	coords(FILTER) - RA and DEC coordinates from sextractor from cropped images
	fluxes1_(FILTER).txt - sextractor output from cropped images
	imagelist(FILTER) - name of cropped file 
	
Translated from icoords.pro by John Capone (jicapone@astro.umd.edu).
Modified by Vicki Toy (vtoy@astro.umd.edu) removes need for crop using weight file
"""

import numpy as np
import os
import astropy.io.fits as pf
import photprocesslibrary as pplib

def icoords():

	#Identify files (must have same number of images files as weight files)
	prefchar    = 'coadd'
	zffiles     = pplib.choosefiles(prefchar + '*_?.fits')
	weightfiles = pplib.choosefiles(prefchar + '*_?.weight.fits')
	
	if len(zffiles) > len(weightfiles):
		print 'Must have matching weight file to each image file to run automatic crop.'
		print 'To use manual crop user manualcrop keyword and change crop values by hand'
		return -1
		
	numfiles = len(zffiles)
		
	#Resample all images using SWarp to a reference image called multicolor using weight files
	swarpstr = ''
	for i in range(numfiles):
		swarpstr = swarpstr + zffiles[i] + ' '

	stackcmd = 'swarp ' + swarpstr + '-DELETE_TMPFILES N -WEIGHT_TYPE MAP_WEIGHT -IMAGEOUT_NAME multicolor.fits -WEIGHTOUT_NAME multicolor.weight.fits'
		
	os.system( stackcmd )

	#Rename all the resampled files to crop files
	for i in range(numfiles):
		tmp = zffiles[i].split('.')[0]
		ifile = tmp + '.resamp.fits' #tmp +'.crop.resamp.fits'
		ofile = tmp + '.ref.fits'	 #tmp + '.crop.fits'
		mvcmd = 'mv -f ' + ifile + ' ' + ofile
		os.system(mvcmd)
		
		ifile = tmp + '.resamp.weight.fits' #tmp + '.crop.resamp.weight.fits'
		ofile = tmp + '.ref.weight.fits' #tmp + '.crop.weight.fits'
		mvcmd = 'mv -f ' + ifile + ' ' + ofile
		os.system(mvcmd)
	
	#Find directory where this python code is located
	propath = os.path.dirname(os.path.realpath(__file__))
	
	#Make sure configuration file is in current working directory, if not copy it from
	#location where this code is stored
	if not os.path.exists('ratir_weighted.sex'): 
		os.system('cp '+propath+'/defaults/ratir_weighted.sex .')
	
	if not os.path.exists('temp.param'): 
		os.system('cp '+propath+'/defaults/temp.param .')
		
	if not os.path.exists('ratir.conv'): 
		os.system('cp '+propath+'/defaults/ratir.conv .')
		
	if not os.path.exists('ratir_nir.nnw'): 
		os.system('cp '+propath+'/defaults/ratir_nir.nnw .')

	#Run sextractor on pipeline reduced files to identify point sources
	for i in range(numfiles):
		cfilter = zffiles[i].split('_')[1].split('.')[0] # extract filter label from file name
		
		#Create imagelist with cropped filename in it for each filter
		ofile = 'imagelist'+cfilter
		f = open(ofile, 'w')
		f.write(zffiles[i].split('.')[0] + '.ref.fits')
		f.close
		
		if cfilter == 'Z' or cfilter == 'Y':
			cfilter = cfilter.lower()
		
		ifile = zffiles[i].split('.')[0] + '.ref.fits ' + ' -WEIGHT_IMAGE ' + zffiles[i].split('.')[0] + '.ref.weight.fits'
		cmd = 'sex ' + ifile + ' -c "ratir_weighted.sex"'
		print cmd
		os.system(cmd)
		
		#Read in results from sextractor and produce IRAF coordinate file save as coords(FILTER)
		#Rename sextractor results to fluxes1_(FILTER).txt
		x,y,ra,dec,mag,magerr,e,fwhm,flag = np.loadtxt('temp.cat', unpack=True)
		f = open('coords' + cfilter, 'w')
		
		for j in range(len(ra)):
			f.write('{:15.6f}{:15.6f}\n'.format(ra[j], dec[j]))
		f.close()
		
		cmd = 'mv -f temp.cat fluxes1_' + cfilter + '.txt'
		print cmd
		os.system(cmd)
