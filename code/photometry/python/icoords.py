"""
NAME:
	icoords

PURPOSE:
	Automatically crops coadded fits files using weighted files.  Can run crop manually, but
	will have to alter script (default for manual crop is 10% off each side).
	Then creates a stacked image with all filters (saved to multicolor.fits and multicolor.weight.fits).  
	Then runs sextractor on each cropped image and saves the RA and DEC values to coords(FILTER) and 
	the sextractor output to fluxes1_(FILTER).txt

OPTIONAL KEYWORD INPUTS:
	manualcrop - set to anything but None if don't want to use weight images to crop, 
			     will have to go manually into program to change values for each filter

OUTPUTS:
	multicolor.fits, multicolor.weight.fits - files with all filter images stacked
	coords(FILTER) - RA and DEC coordinates from sextractor from cropped images
	fluxes1_(FILTER).txt - sextractor output from cropped images
	imagelist(FILTER) - name of cropped file 
	
Translated from icoords.pro by John Capone (jicapone@astro.umd.edu).
Modified by Vicki Toy (vtoy@astro.umd.edu) to include automatic cropping using weight file
"""

import numpy as np
import os
import sys
import astropy.io.fits as pf
import photprocesslibrary as pplib

def icoords(manualcrop=None):

	#Identify files (must have same number of images files as weight files)
	prefchar    = 'coadd'
	zffiles     = pplib.choosefiles(prefchar + '*_?.fits')
	weightfiles = pplib.choosefiles(prefchar + '*_?.weight.fits')
	
	if len(zffiles) > len(weightfiles):
		print 'Must have matching weight file to each image file to run automatic crop.'
		print 'To use manual crop user manualcrop keyword and change crop values by hand'
		return -1
		
	numfiles = len(zffiles)

	#Perform initial crop to remove noisy edges
	for i in range(numfiles):
	
		whdulist = pf.open(weightfiles[i])
		wh       = whdulist[0].header
		wimg     = whdulist[0].data

		xaxis = wh['NAXIS1']
		yaxis = wh['NAXIS2']
		
		#Finds points to crop to using weight files 
		#(99% level chosen after visual inspection, can be modified if not producing good crops)
		scale = 0.99
		x1 = pplib.weightedge(wimg, range(xaxis), scale=scale, column=1)
		x2 = pplib.weightedge(wimg, reversed(range(xaxis)), scale=scale, column=1)
		y1 = pplib.weightedge(wimg, range(yaxis), scale=scale, row=1)
		y2 = pplib.weightedge(wimg, reversed(range(yaxis)), scale=scale, row=1)
		
		if manualcrop is not None:
			hdulist = pf.open(zffiles[i])
			h       = hdulist[0].header
			img     = hdulist[0].data
	
			xaxis = h['NAXIS1']
			yaxis = h['NAXIS2']
			
			x1 = 0.1*xaxis
			x2 = 0.9*xaxis
			y1 = 0.1*yaxis
			y2 = 0.9*yaxis			

		#Crops both the weight and image fits files and changes headers so
		#center pixels are adjusted and adds information of crop
		wofile = weightfiles[i].split('.')[0] + '.crop.weight.fits'
		pplib.hextractlite(wofile, wimg, wh, x1, x2, y1, y2)

		hdulist = pf.open(zffiles[i])
		h       = hdulist[0].header
		img     = hdulist[0].data
		
		ofile = zffiles[i].split('.')[0] + '.crop.fits'
		pplib.hextractlite(ofile, img, h, x1, x2, y1, y2)
		
	#Resample all cropped images using SWarp to a reference image called multicolor
	swarpstr = ''
	for i in range(numfiles):
		ifile = zffiles[i].split('.')[0] + '.crop.fits'	
		swarpstr = swarpstr + ifile + ' '

	stackcmd = 'swarp ' + swarpstr + '-DELETE_TMPFILES N -IMAGEOUT_NAME multicolor.fits -WEIGHTOUT_NAME multicolor.weight.fits'
	os.system( stackcmd )
	
	#Rename all the resampled files to overwrite crop files
	for i in range(numfiles):
		tmp = zffiles[i].split('.')[0]
		ifile = tmp +'.crop.resamp.fits'
		ofile = tmp + '.crop.fits'
		mvcmd = 'mv -f ' + ifile + ' ' + ofile
		os.system(mvcmd)
		
		ifile = tmp + '.crop.resamp.weight.fits'
		ofile = tmp + '.crop.weight.fits'
		mvcmd = 'mv -f ' + ifile + ' ' + ofile
		os.system(mvcmd)
	
	#Find directory where this python code is located
	propath = os.path.dirname(os.path.realpath(__file__))
	
	#Make sure configuration file is in current working directory, if not copy it from
	#location where this code is stored
	if not os.path.exists('ratir_nir.sex'): 
		os.system('cp '+propath+'/defaults/ratir_nir.sex .')
	
	if not os.path.exists('temp.param'): 
		os.system('cp '+propath+'/defaults/temp.param .')
		
	if not os.path.exists('sex.conv'): 
		os.system('cp '+propath+'/defaults/sex.conv .')

	#Run sextractor on pipeline reduced files to identify point sources
	for i in range(numfiles):
		cfilter = zffiles[i].split('_')[1].split('.')[0] # extract filter label from file name
		
		#Create imagelist with cropped filename in it for each filter
		ofile = 'imagelist'+cfilter
		f = open(ofile, 'w')
		f.write(zffiles[i].split('.')[0] + '.crop')
		f.close
		
		if cfilter == 'Z' or cfilter == 'Y':
			cfilter = cfilter.lower()
			
		#Run SExtractor on each image
		ifile = zffiles[i].split('.')[0] + '.crop.fits'
		cmd = 'sex ' + ifile + ' -c "ratir_nir.sex"'
		print cmd
		os.system(cmd)
		
		#Read in results from sextractor and produce IRAF coordinate file save as coords(FILTER)
		#Rename sextractor results to fluxes1_(FILTER).txt
		#Removed FWHM aannulus files (all were 0...???)
		x,y,ra,dec,mag,magerr,e,fwhm = np.loadtxt('temp.cat', unpack=True)
		f = open('coords' + cfilter, 'w')
		
		for j in range(len(ra)):
			f.write('{:15.6f}{:15.6f}\n'.format(ra[j], dec[j]))
		f.close()
		
		cmd = 'mv -f temp.cat fluxes1_' + cfilter + '.txt'
		print cmd
		os.system(cmd)
