'''
NAME:
	assoc
	
PURPOSE:
	Crops all of the filter images and combined image to the same size and coordinates.
	Then finds sources using the master combined image, and calculates the magnitude based on
	just the filtered images using sextractor.  Saves new sextractor values to "fluxes2_(FILTER).txt'

OUTPUTS:
	coadd*.crop.multi.fits - files cropped so that they include all filters
	fluxes2_(FILTER).txt   - files with sextractor output using multicolor to find sources 
							 and filter fits files to find magnitudes
							 
Translated from assoc.pro and modified by Vicki Toy (vtoy@astro.umd.edu)
'''

import astropy.io.fits as pf
from string import index
from numpy import shape
from astropy import wcs
import os
import photprocesslibrary as pplib

def assoc():

	coaddfiles = pplib.choosefiles('coadd*.crop.fits')

	ra1arr  = []
	dec1arr = []
	ra2arr  = []
	dec2arr = []

	#Finds the RA and DEC of the first and the last pixel of each cropped coadded file
	for files in coaddfiles:

		fitsfile = pf.open(files)
		fitsheader = fitsfile[0].header
		data = fitsfile[0].data
	
		imSize = shape(data)
	
		#Converts pixel value to RA and DEC using header information (AstroPy function)
		w = wcs.WCS(fitsheader)
		pixcrd = [[0.,0.], [imSize[1]-1.0, imSize[0]-1.0]]
		[[ra1,dec1],[ra2,dec2]] = w.wcs_pix2world(pixcrd, 0)
	
		#Stores information into arrays
		ra1arr.append(ra1)
		dec1arr.append(dec1)
		ra2arr.append(ra2)
		dec2arr.append(dec2)

	#Finds the coordinates that fit all of the data
	raleft  = min(ra1arr)
	raright = max(ra2arr)
	decbot  = max(dec1arr)
	dectop  = min(dec2arr)

	#Crops data so the size of every filter image matches and saves to file 'coadd*.crop.multi.fits'
	for files in coaddfiles:

		newfile = files[:-4]+'multi.fits'
		fitsfile = pf.open(files)
		fitsheader = fitsfile[0].header
		data = fitsfile[0].data
	
		w = wcs.WCS(fitsheader)
		[[x1,y1],[x2,y2]] = w.wcs_world2pix([[raleft, decbot], [raright, dectop]], 0)
		
		pplib.hextractlite(newfile, data, fitsheader, x1+1, x2, y1+1, y2)

	#Crops the multicolor fits file so it matches the same coordinates
	mixfile = 'multicolor.fits'
	mixfitsfile = pf.open(mixfile)
	mixfitsheader = mixfitsfile[0].header
	mixdata = mixfitsfile[0].data

	mixw = wcs.WCS(mixfitsheader)
	[[mx1,my1],[mx2,my2]] = mixw.wcs_world2pix([[raleft, decbot], [raright, dectop]], 0)
	pplib.hextractlite(mixfile, mixdata, mixfitsheader, mx1+1, mx2, my1+1, my2)

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

	#Uses sextractor to find the magnitude and location of sources for each file
	#Saves this information into 'fluxes2_*.txt' files
	for files in coaddfiles:

		#Finds filter name and makes sure it is capitalized correctly		
		filter = files.split('_')[1].split('.')[0]
	
		if filter.lower() in ('j','h','k'):
			filter = filter.upper()
		else:
			filter = filter.lower()
		
		compfile = files[:-4]+'multi.fits'

		#Call to sextractor in double image mode (image1 used for detection of sources, image2 only for measurements - must be same size) 
		#"sex image1, image2 -c configuration file" uses multicolor file for source detection and filter file for magnitude measurements
	
		os.system('sex ' + mixfile + ', ' + compfile + ' -c "ratir_nir.sex"')
		os.system('mv -f temp.cat fluxes2_'+filter+'.txt')