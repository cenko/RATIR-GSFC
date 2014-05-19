'''
NAME:
	calcoff
	
PURPOSE:
	Calculates the magnitude offset from image (using values from SExtractor saved in 
	'fluxes1_(FILTER).txt' and compares with a star catalog.  Saves the sigma clipped
	median magnitude offset to 'irafoffset(FILTER)' 
	
OUTPUTS:
	irafoffset(FILTER)  - file with the median magnitude offset and error 
						  (sigma of median/number of sextractor sources used in fitting)
	coadd*crop.fits.im  - file with RA, DEC, and magnitude of sextractor sources, with RA
						  and DEC recalculated with x and y coordinates using WCS conversions
	coadd*crop.fits.cat - file with RA, DEC, and magnitude of catalog sources

Requires cdsclient package installed on computer

Translated from calcoff.pro and modified by Vicki Toy (vtoy@astro.umd.edu)
'''

import get_SEDs
import astropy.io.fits as pf
from astropy import wcs
import numpy as np
import photprocesslibrary as pplib


def calcoff():
	#Finds all 'fluxes1_*.txt' (sextractor first results from cropped images) in current directory
	fluxfiles = pplib.choosefiles('fluxes1_*.txt')

	for files in fluxfiles:

		#Finds filter from filename convention
		filter = files.split('_')[1].split('.')[0]
	
		#Filter names must match conventional filter capitalization to use catalogs (r,i,z,y,J,H,K)
		if filter.lower() in ('j','h','k'):
			filter = filter.upper()
		else:
			filter = filter.lower()
	
		#Reads file and saves columns into variables
		[x,y,ra,dec,mag_aper,magerr_aper,e,fwhm,flag] = np.loadtxt(files, unpack=True)
	
		#Put limits on the magnitude aperture, magnitude aperture error, and fwhm
		tmp = (mag_aper > 10) & (mag_aper < 20) & (magerr_aper < 0.4) & (fwhm < 9)
	
  		x    = x[tmp]
   		y    = y[tmp]
   		mag  = mag_aper[tmp]
   		magerr = magerr_aper[tmp]
   	
   		#Finds correct cropped fits file to read
   		filename = open('imagelist'+filter, 'r').read().strip()
		fitsfile = pf.open(filename)
		fitsheader = fitsfile[0].header
		
		#Calculate coordinates from pixel position and header (slightly different than ra and dec) and saves to '*.im'
		w = wcs.WCS(fitsheader)
		pixcrd = []
		for ind in range(len(x)):
			pixcrd.append([x[ind],y[ind]])
		coord_out = w.wcs_pix2world(pixcrd, 0)
	
		ra_out  = coord_out[:,0]
		dec_out = coord_out[:,1]
	
		imfile = filename+'.im'
		np.savetxt(imfile, np.transpose([ra_out, dec_out, mag]), fmt='%15.6f', header='RA\t DEC\t INST_MAG\t')

		#Use the more precisely calculated RA and DEC from converting header information than from sextractor
		ra  = ra_out
		dec = dec_out

		#Calculates SEDs (for objects with matching coordinates with star catalog) and outputs '*.cat' 
		catfile = filename+'.cat'
		zp, mad = get_SEDs.zeropoint(imfile, filter, output_file=catfile)
	
		#Unpacks magnitude values
		[refra,refdec,u_mag,g_mag,r_mag,i_mag,
			z_mag,y_mag,bigB_mag,bigV_mag,bigR_mag,
			bigI_mag,J_mag,H_mag,K_mag,u_err,g_err,
			r_err,i_err,z_err,y_err,bigB_err,bigV_err,
			bigR_err,bigI_err,J_err,H_err,K_err,mode] = np.loadtxt(catfile, unpack=True)
	
		#Can expand to include other filters
		magdict = {'r_mag': r_mag, 'i_mag': i_mag, 'z_mag': z_mag, 'y_mag': y_mag, 'J_mag': J_mag, 'H_mag': H_mag, 'K_mag': K_mag} 
		errdict = {'r_err': r_err, 'i_err': i_err, 'z_err': z_err, 'y_err': y_err, 'J_err': J_err, 'H_err': H_err, 'K_err': K_err} 	
	
		#Only use values that have mode set to 0 (modeled from SDSS and 2-MASS)	and save the magnitude (of file's same filter
		use = (mode == 0)
	
		refra  = refra[use]
		refdec = refdec[use]
		refmag = magdict[filter+'_mag'][use]
		referr = errdict[filter+'_err'][use]

		#Find sources that are within one arcsecond of a catalog source
		difarr = []
	
		for ind in range(len(ra)):
			smatch = pplib.nearest( ra[ind]*np.cos(dec[ind]*np.pi/180.),dec[ind],refra*np.cos(refdec*np.pi/180.),refdec, 1./3600.)
			
			#If any of the sources match a non-
			if any(smatch) & (refmag[smatch] < 99):
				dif = mag[ind] - refmag[smatch]
				difarr.extend(dif)

		difarr = np.array(difarr)
	
		#Calculates statistics of magnitude difference based on iterative sigma clipping and save to 'irafoffset_*'
		[difMean, difSig, difMedian, difMask, difSaveMask] = pplib.djs_iterstat(difarr)
		err = difSig/np.sqrt(len(ra))
	
		s_difMedian = '%.6f'%(difMedian)
		s_err       = '%.6f'%(err)
	
		ofile = 'irafoffset'+filter
		f = open(ofile, 'w')
		f.write(s_difMedian+'\t'+s_err)
		f.close()