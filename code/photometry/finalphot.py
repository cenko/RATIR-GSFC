'''
NAME:
	finalphot

PURPOSE:
	Calculates absolute magnitude and errors based on 'fluxes2_(FILTER).txt' (magnitude from sextractor), 
	'irafoffset*' (correction for instrumentation magnitude).  Saves final magnitudes to 'finalphot(FILTER).am'

OUTPUTS: 
	finalphot(FILTER).am - file containing pixel and coordinate location and 
						   corrected magnitude of sextractor sources found (fluxes2_(FILTER).txt)
						   
Translated from finalphot.pro and modified by Vicki Toy (vtoy@astro.umd.edu)
'''

import numpy as np
import photprocesslibrary as pplib

def finalphot():
	#Finds matching flux files
	txtmatch  = 'fluxes2_*.txt'
	fluxfiles = pplib.choosefiles(txtmatch)

	for file in fluxfiles:
	
		#Finds the filter of each file from naming convention and
		#makes sure the filters have the right capitalization 
		filter = file.split('_')[1].split('.')[0]

		if filter.lower() in ('j','h','k'):
			filter = filter.upper()
		else:
			filter = filter.lower()

		#Grabs data from irafoffset and the current fluxfile
		ifile  = 'irafoffset'+filter
	
		[ABoffset, err] = np.loadtxt(ifile, unpack=True)
		[x,y,ra,dec,mag,magerr,e,fwhm] = np.loadtxt(file, unpack=True)
	
		#Converts to absolute magnitude (calibrates with catalog values from offset file)
		amfile = 'finalphot'+filter+'.am'
	
		aperid = np.arange(len(mag)) + 1
		aper_err = np.sqrt(magerr**2. + err**2.)
		aper_mag = mag - ABoffset
	
		#Creates Absolute Magnitude file with coordinates
		np.savetxt(amfile, np.transpose([aperid, x, y, ra, dec, aper_mag, aper_err]), fmt='%15.6f', 
			header='ID\t X\t Y\t RA\t DEC\t CAL_MAG\t CAL_MAG_ERR\t')

	
	