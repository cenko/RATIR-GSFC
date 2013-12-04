#Photometry processing library (misc function that can be reused)

import astropy.io.fits as pf
import fnmatch
import os


"""
NAME:
	choosefiles
PURPOSE:
	Returns files in directory "loc" that contain "selection" text
INPUTS:
	selection - text required to be in filename
OPTIONAL KEYWORD INPUTS:
	loc - directory where to look, default is current working directory
OUTPUT:
	Matches to search
EXAMPLE:
	matches = choosefiles('coadd*.fits', loc='.')

Written by John Capone (jicapone@astro.umd.edu)
"""
#
def choosefiles( selection, loc='.' ):
	matches = []
	for files in os.listdir(loc):
		if fnmatch.fnmatch( files, selection ):
			matches.append(files)
	return matches

"""
NAME:
	weightedge
PURPOSE:
	To find the innermost corner of a weighted image for a good crop
	Looks for point where values are roughly constant between rows or columns
	Can change this with scale keyword	
INPUTS:
	array - array to search
	itarray - array of indices to search through	
OPTIONAL KEYWORD INPUTS:
	column - set True if iterating through x-axis
	row    - set True if iterating through y-axis
	scale  - set value to look for scaling to previous row or column
OUTPUT:
	Returns the column or row of "nonzero" component on weighted file
EXAMPLE:
	leftedge = weightedge(data, range(xaxis), column=1, scale=0.99)
	
Written by Vicki Toy (vtoy@astro.umd.edu)
"""
def weightedge(array, itarray, scale=1, column=None, row=None):
	oldsum = 0
	for i in itarray:
		if column is not None:
			newsum = sum(array[:,i])
		elif row is not None:
			newsum = sum(array[i,:])
			
		if scale*newsum >= oldsum:
			oldsum = newsum
		else:
			return i-1		


'''
Pared down version of hextract.pro from the IDL Astro Library
Translated to python by Vicki Toy (vtoy@astro.umd.edu)

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
	