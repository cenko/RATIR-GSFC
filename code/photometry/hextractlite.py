'''
Pared down version of hextract.pro from the IDL Astro Library
Translated to python by Vicki Toy (vtoy@astro.umd.edu

NAME:
	hextractlite
PURPOSE:
	Extracts subimage from an array and updates astrometry in FITS file.  Saves to newfile name
INPUTS:
	newfile     - file name to save altered fits file to
	data        - array to be altered
	fitsheader  - header to be altered
	x1,x2,y1,y2 - respectively: first and last x pixel, first and last y pixel to be
				  extracted from data array.  Need to be integer scalars
OUTPUTS:
	Saves altered data to newfile
EXAMPLE:
	hextractlite(newfile, data, fitsheader, 100,500,20,1900)
'''

import astropy.io.fits as pf

def hextractlite(newfile, data, fitsheader, x1, x2, y1, y2):

	#Truncates the coordinates
	x1 = int(x1)
	y1 = int(y1)
	x2 = int(x2)
	y2 = int(y2)

	fitsheader.add_history('HEXTRACT: Original image size was '+str(fitsheader['NAXIS2'])+' by '+str(fitsheader['NAXIS1']))
	fitsheader.add_history('Extracted Image: ['+ str(y1)+':'+str(y2+1)+','+str(x1)+':'+ str(x2+1)+']')

	#Naxis altered to reflect the new subarray
	fitsheader.update('naxis1', x2-x1+1)
	fitsheader.update('naxis2', y2-y1+1)
	
	#Crpix shifted by new origin
	oldcrpix1 = fitsheader['crpix1']
	oldcrpix2 = fitsheader['crpix2']
	fitsheader.update('crpix1', oldcrpix1-x1)
	fitsheader.update('crpix2', oldcrpix2-y1)

	#Extract subarray
	#Python has opposite coordinate convention
	newdata = data[y1:y2+1, x1:x2+1]	
	
	#Saves new fits file and header to specified fits file name
	pf.writeto(newfile, newdata, fitsheader, clobber=True)
	