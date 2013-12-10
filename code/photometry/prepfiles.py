'''
Creates files with the desired reference files within

OUTPUTS: 'imagelist*'
'''

from choosefiles import choosefiles
from string import index

coaddfiles = choosefiles('coadd*.crop.fits')

for file in coaddfiles:

	#Finds filter name from character after underscore
	underscore_loc = index(file, '_')
	filter = file[underscore_loc + 1]
	
	#Creates reference file name (remove .fits from end)
	#Saves reference file ('coadd*.crop') in imagelist file ('imagelist*')
	fits_loc = index(file, '.fits')
	reffile = file[:fits_loc]
	
	ofile = 'imagelist'+filter
	f = open(ofile, 'w')
	f.write(reffile)
	f.close()	