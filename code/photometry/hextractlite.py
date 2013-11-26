'''
Pared down version of hextract.pro from the IDL Astro Library
Takes a subarray of fits file and alters the naxis and crpix header
Saves the new fits file and new header under the specified file name
'''

import astropy.io.fits as pf

def hextractlite(newfile, data, fitsheader, x1, x2, y1, y2):

	#Truncates the coordinates
	x1 = int(x1)
	y1 = int(y1)
	x2 = int(x2)
	y2 = int(y2)
	
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
	