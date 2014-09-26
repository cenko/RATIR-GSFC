#######################################
# Routines for processing of RATIR data
# SBC - Started 2 May 2014
#######################################

import pyraf
from pyraf import iraf
import astropy.io.fits as pyfits
import astropy.coordinates
from scipy import optimize
import numpy as np
import glob, os, shutil, re
import iqpkg
from iqutils import *

# Necessary packages
iraf.images()
iraf.immatch()
#iraf.imfilter()
iraf.noao()
iraf.imred()
iraf.ccdred()
iraf.stsdas()
iraf.hst_calib()
iraf.nicmos()
iraf.imutil()

yes=iraf.yes
no=iraf.no
INDEF=iraf.INDEF
globclob=yes
globver=yes

RATIRFILTS = ["i"]
RATIRPIXSCALE =  {"C0": 0.32, "C1": 0.32, "C2": 0.30, "C3": 0.30}
RATIRFLAT = "evening flats (bright)"
RATIRBIAS = "biases"
RATIRSAT = {"C1": 12000.0}
C0FILTS = []
C1FILTS = ["i"]
C2FILTS = []
C3FILTS = []
PYASTROM = "/Users/scenko/python/RATIR-GSFC/code/reduction/astrom/vlt_autoastrometry.py"


##############################################################################

def iq_ratir(chip="C1", filt="i", reflist="sdss.reg"):

	'''Process RATIR data for given chip and filter'''
	
	# Check and see if bias frame exists.  If not, create one
	if not os.path.exists("Bias-%s.fits" % chip):
	
		imlist = glob.glob("[0-9]*T[0-9]*%sb.fits" % chip)
		hdr = pyfits.getheader(imlist[0])
		
		# Setup zerocombine
		zerocombine = iraf.ccdred.zerocombine
		zerocombine.combine = 'median'
		zerocombine.reject = 'avsigclip'
		zerocombine.ccdtype = ''
		zerocombine.process = no
		zerocombine.delete = no
		zerocombine.clobber = no
		zerocombine.scale = 'none'
		zerocombine.statsec = '*'
		zerocombine.nlow = 0
		zerocombine.nhigh = 1
		zerocombine.nkeep = 1
		zerocombine.mclip = yes
		zerocombine.lsigma = 3.0
		zerocombine.hsigma = 3.0
		zerocombine.rdnoise = 0.0
		zerocombine.gain = hdr['SOFTGAIN']
		zerocombine.snoise = 0
		zerocombine.pclip = -0.5
		zerocombine.blank = 0.0
		
		# Run zerocombine
		bstr = ",".join(imlist)
		zerocombine(input=bstr, output="Bias-%s.fits" % chip)
		
	# Bias subtract flat frames
	imlist = glob.glob("[0-9]*T[0-9]*%sf.fits" % chip)
	flis = []
	for im in imlist:
		hdr = pyfits.getheader(im)
		if (hdr["FILTER"] == filt) and (not os.path.exists("b%s" % im)):
			
			# Subtract bias frame
			ccdproc = iraf.ccdred.ccdproc	
			ccdproc.ccdtype = ""
			ccdproc.noproc = no
			ccdproc.fixpix = no
			ccdproc.overscan = no
			ccdproc.trim = no
			ccdproc.zerocor = yes
			ccdproc.darkcor = no
			ccdproc.flatcor = no
			ccdproc.illumcor = no
			ccdproc.fringecor = no
			ccdproc.readcor = no
			ccdproc.scancor = no
			ccdproc.readaxis = 'line'
			ccdproc.fixfile = ""
			ccdproc.zero = "Bias-%s.fits" % chip
			ccdproc.interactive = no
			ccdproc.function = "legendre"
			ccdproc.order = 1
			ccdproc.sample = "*"
			ccdproc.naverage = 1
			ccdproc.niterate = 1
			ccdproc.low_reject = 3.0
			ccdproc.high_reject = 3.0
			ccdproc.grow = 0
			ccdproc(images=im, output="b%s" % im)
				
			flis.append("b%s" % im)
		elif (hdr["FILTER"] == filt):
			flis.append("b%s" % im)

	# Create flat field
	if not os.path.exists("Flat-%s-%s.fits" % (chip, filt)):
	
		glis = []
		for im in flis:
			iraf.iterstat.nsigrej = 5.0
			iraf.iterstat.maxiter = 10
			iraf.iterstat.verbose = globver
			iraf.iterstat.lower = INDEF
			iraf.iterstat.upper = INDEF
			iraf.iterstat(im)
			if (iraf.iterstat.median < RATIRSAT[chip]) and (iraf.iterstat.median > 1000.0):
				glis.append(im)
		     	
		# Set up flatcombine
		flatcombine = iraf.ccdred.flatcombine
		flatcombine.combine = 'median'
		flatcombine.reject = 'avsigclip'
		flatcombine.ccdtype = ''
		flatcombine.scale = 'median'
		flatcombine.statsec = ''
		flatcombine.nlow = 1
		flatcombine.nhigh = 1
		flatcombine.nkeep = 1
		flatcombine.mclip = yes
		flatcombine.lsigma = 3.0
		flatcombine.hsigma = 3.0
		flatcombine.rdnoise = 0.0
		flatcombine.gain = hdr["SOFTGAIN"]
		flatcombine.snoise = 0.0
		flatcombine.pclip = -0.5
			
		# Run flatcombine
		fstr = ",".join(glis[:10])
		flatcombine(fstr, output="Flat-%s-%s.fits" % (chip, filt))
			
		# Normalize
		iraf.iterstat.nsigrej = 5.0
		iraf.iterstat.maxiter = 10
		iraf.iterstat.verbose = globver
		iraf.iterstat.lower = INDEF
		iraf.iterstat.upper = INDEF
		iraf.iterstat("Flat-%s-%s.fits" % (chip, filt))
		iraf.imarith("Flat-%s-%s.fits" % (chip, filt), "/", 
			         iraf.iterstat.median, "Flat-%s-%s.fits" % (chip, filt))
			         	
	# Bias subtract and flat-field science images
	imlist = glob.glob("[0-9]*T[0-9]*%so_img_1.fits" % chip)
	for im in imlist:
		hdr = pyfits.getheader(im)
		if (hdr["FILTER"] == filt) and (not os.path.exists("fb%s" % im)):
			
			hdr = pyfits.getheader(im)
			
			# Subtract bias frame
			ccdproc = iraf.ccdred.ccdproc	
			ccdproc.ccdtype = ""
			ccdproc.noproc = no
			ccdproc.fixpix = no
			ccdproc.overscan = no
			ccdproc.trim = no
			ccdproc.zerocor = yes
			ccdproc.darkcor = no
			ccdproc.flatcor = yes
			ccdproc.illumcor = no
			ccdproc.fringecor = no
			ccdproc.readcor = no
			ccdproc.scancor = no
			ccdproc.readaxis = 'line'
			ccdproc.fixfile = ""
			ccdproc.zero = "Bias-%s.fits" % chip
			ccdproc.flat = "Flat-%s-%s.fits" % (chip, filt)
			ccdproc.interactive = no
			ccdproc.function = "legendre"
			ccdproc.order = 1
			ccdproc.sample = "*"
			ccdproc.naverage = 1
			ccdproc.niterate = 1
			ccdproc.low_reject = 3.0
			ccdproc.high_reject = 3.0
			ccdproc.grow = 0
			ccdproc(images=im, output="fb%s" % im)
			
	# Create sky frame (including dark current!)
	if not os.path.exists("Sky-%s-%s.fits" % (chip, filt)):
		imlist = glob.glob("fb[0-9]*T[0-9]*%so_img_1.fits" % chip)
		hdr = pyfits.getheader(imlist[0])
		slis = []
		for im in imlist:
			if (hdr["FILTER"] == filt):
				slis.append(im)
				
		# Set up combine
		imcombine = iraf.immatch.imcombine
		imcombine.headers = ""
		imcombine.bpmasks = ""
		imcombine.rejmasks = ""
		imcombine.nrejmasks = ""
		imcombine.expmasks = ""
		imcombine.sigmas = ""
		imcombine.combine = "median"
		imcombine.reject = "avsigclip"
		imcombine.project = no
		imcombine.outtype = "real"
		imcombine.offsets = "none"
		imcombine.masktype = "none"
		imcombine.scale = "median"
		imcombine.zero = "none"
		imcombine.weight = "none"
		imcombine.statsec = ""
		imcombine.lsigma = 3.0
		imcombine.hsigma = 3.0
		imcombine.rdnoise = 0.0
		imcombine.gain = hdr["SOFTGAIN"]
		
		sstr = 	",".join(slis[:10])
		imcombine(sstr, "Sky-%s-%s.fits" % (chip, filt))
		
		# Normalize
		iraf.iterstat("Sky-%s-%s.fits" % (chip, filt))
		iraf.imarith("Sky-%s-%s.fits" % (chip, filt), "/", 
			         iraf.iterstat.median, "Sky-%s-%s.fits" % (chip, filt))
			         
	# Subtract sky frame (and dark!)
	imlist = glob.glob("fb[0-9]*T[0-9]*%so_img_1.fits" % chip)
	for im in imlist:
		if not os.path.exists("s%s" % im):
		
			iraf.iterstat(im)
			iraf.imarith("Sky-%s-%s.fits" % (chip, filt), "*", iraf.iterstat.median,
			             "temp.fits")
			iraf.imarith(im, "-", "temp.fits", "s%s" % im)
			os.remove("temp.fits")
			fimg = pyfits.open("s%s" % im, "update")
			fimg[0].header["SKYMED"] = iraf.iterstat.median
			fimg[0].header["SKYSIG"] = iraf.iterstat.sigma
			fimg[0].header["SKYSUB"] = 1
			fimg.flush()
								
	# Remove cosmic rays
	imlist = glob.glob("sfb[0-9]*T[0-9]*%so_img_1.fits" % chip)
	for im in imlist:
		if not os.path.exists("c%s" % im):
			fimg = pyfits.open(im)
			iraf.lacos_im(im, "c%s" % im, "c%s.mask.fits" % im[:-5],
			              gain=hdr['SOFTGAIN'], readn=0.0, statsec="*,*",
			              skyval=fimg[0].header["SKYMED"], sigclip=4.5, 
			              sigfrac=0.5, objlim=1.0, niter=3, verbose=no)

	# Add WCS keywords
	clis = glob.glob("csfb[0-9]*T[0-9]*%so_img_1.fits" % chip)
	for im in clis:
	
		fimg = pyfits.open(im, mode='update')
		if fimg[0].header.get("RA")==None:
		
			ra0 = fimg[0].header["STRCURA"]; dec0 = fimg[0].header["STRCUDE"]
			nax = [fimg[0].header["NAXIS1"], fimg[0].header["NAXIS2"]]
			fimg[0].header["RA"] = ra0
			fimg[0].header["DEC"] = dec0		
			fimg[0].header["PIXSCALE"] = RATIRPIXSCALE[chip]
			fimg[0].header["PIXSCAL1"] = RATIRPIXSCALE[chip]
			fimg[0].header["PIXSCAL2"] = RATIRPIXSCALE[chip]
			fimg[0].header["CTYPE1"] = "RA---TAN"
			fimg[0].header["CTYPE2"] = "DEC--TAN"
			fimg[0].header["WCSDIM"] = 2
			fimg[0].header["WAT0_001"] = "system=image"
			fimg[0].header["WAT1_001"] = "wtype=tan axtype=ra"
			fimg[0].header["WAT2_001"] = "wtype=tan axtype=dec"
			fimg[0].header["LTM1_1"] = 1.0
			fimg[0].header["LTM2_2"] = 1.0	
			fimg[0].header["CRPIX1"] = nax[0] / 2.0
			fimg[0].header["CRPIX2"] = nax[1] / 2.0
			fimg[0].header["CRVAL1"] = ra0
			fimg[0].header["CRVAL2"] = dec0
			fimg[0].header["CD1_1"] = -RATIRPIXSCALE[chip] / 3600.0
			fimg[0].header["CD1_2"] = 0.0
			fimg[0].header["CD2_1"] = 0.0
			fimg[0].header["CD2_2"] = RATIRPIXSCALE[chip] / 3600.0
			
		fimg.flush()
		
	# Crude astrometry
	for im in clis:
	
		if not os.path.exists("a%s" % im):
			os.system("python %s %s" % (PYASTROM, im))
				
	# Refine Astrometry
	o1 = open("daofind.param", "w")
	o1.write("NUMBER\nXWIN_IMAGE\nYWIN_IMAGE\nMAG_AUTO\nFLAGS\nA_IMAGE\nB_IMAGE\n")
	o1.write("ELONGATION\nFWHM_IMAGE\nXWIN_WORLD\nYWIN_WORLD\n")
	o1.write("ERRAWIN_IMAGE\nERRBWIN_IMAGE\nERRTHETAWIN_IMAGE\nERRAWIN_WORLD\n")
	o1.write("ERRBWIN_WORLD\nERRTHETAWIN_WORLD\nFLUX_AUTO\nFLUX_RADIUS\nFLUXERR_AUTO")
	o1.close()
	
	o2 = open("default.conv", "w")
	o2.write("CONV NORM\n# 5x5 convolution mask of a gaussian PSF with FWHM = 3.0 pixels.\n0.092163 0.221178 0.296069 0.221178 0.092163\n0.221178 0.530797 0.710525 0.530797 0.221178\n0.296069 0.710525 0.951108 0.710525 0.296069\n0.221178 0.530797 0.710525 0.530797 0.221178\n0.092163 0.221178 0.296069 0.221178 0.092163")
	o2.close()
	
	alis = glob.glob("acsfb[0-9]*T[0-9]*%so_img_1.fits" % chip)
	for im in alis:
	
		# Detect sources
		os.system("sex -CATALOG_NAME %s.cat -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME daofind.param -DETECT_THRESH 2.0 -ANALYSIS_THRESH 2.0 -GAIN_KEY SOFTGAIN -PIXEL_SCALE 0 %s" % (im[:-5], im))
		
	# Run scamp for alignment
	imlist = " ".join(alis)
	os.system("scamp -ASTREF_CATALOG SDSS-R7 -DISTORTDEG 1 -SOLVE_PHOTOM N -SN_THRESHOLDS 3.0,10.0 -CHECKPLOT_DEV NULL %s" % imlist.replace(".fits", ".cat"))
	
	# Update header
	for im in alis:
		os.system("missfits %s" % im)
		os.system("rm %s.head %s.back" % (im[:-5], im))
		
	# Find good images for coadd
	slis = []
	for im in alis:
		hdr = pyfits.getheader(im)
		rms1 = hdr["ASTRRMS1"]; rms2 = hdr["ASTRRMS2"]
		if (rms1 < 1.0e-4) and (rms1 > 5.0e-6) and (rms2 < 1.0e-4) and (rms2 > 5.0e-6):
			slis.append(im)
			
	# Initial (unweighted) coadd
	sstr = " ".join(slis)
	os.system("swarp -GAIN_KEYWORD SOFTGAIN %s" % sstr)
	
	# Identify sources in the coadded frame
	hdr = pyfits.getheader("coadd.fits")
	iqpkg.iqobjs("coadd.fits", 10.0, RATIRSAT[chip], skyval="0.0",
	             pix=RATIRPIXSCALE[chip], gain=hdr["GAIN"], aperture=20.0,
	             wtimage="coadd.weight.fits")
	hdr = pyfits.getheader("coadd.fits")
	cpsfdiam = 1.34 * float(hdr["SEEPIX"])
	iqpkg.iqobjs("coadd.fits", 10.0, RATIRSAT[chip], skyval="0.0",
	             pix=RATIRPIXSCALE[chip], gain=hdr["GAIN"], aperture=cpsfdiam,
	             wtimage="coadd.weight.fits")
	             
	refstars1 = Starlist("coadd.fits.stars")
	refstars1.pix2wcs("coadd.fits")
	truemags = np.array(refstars1.mags())
	maglist = []; errlist = []
	
	# (Relative) Zeropoints for individual images
	for im in slis:
		hdr = pyfits.getheader(im)
		iqpkg.iqobjs(im, 3.0, RATIRSAT[chip], skyval="!SKYMED", 
		       pix=RATIRPIXSCALE[chip], gain=hdr["SOFTGAIN"], aperture=20.0)
		#hdr = pyfits.getheader(im)
		#psfdiam = 1.34 * float(hdr["SEEPIX"])
		#iqpkg.iqobjs(im, 3.0, RATIRSAT[chip], skyval="!SKYMED", 
		       #pix=RATIRPIXSCALE[chip], gain=hdr["SOFTGAIN"], aperture=psfdiam)
		stars = Starlist("%s.stars" % im)
		refstars1.wcs2pix(im)
		a,b = stars.match(refstars1,tol=10.0,maxnum=1000)
		newmags = np.zeros(len(refstars1)); newwts = np.zeros(len(refstars1))
		for i in range(len(refstars1)):
			for j in range(len(a)):
				if b[j].mag == truemags[i]:
					newmags[i] = a[j].mag
					newwts[i] = 1.0 / np.power(np.maximum(a[j].magu, 0.01),2)
					continue
		maglist.append(newmags); errlist.append(newwts)
		
	obsmags = np.array(maglist)
	wts = np.array(errlist)
	
	[zpts, scatt, rms] = calc_zpts(truemags, obsmags, wts, sigma=3.0)
	
	# Update headers
	medzp = np.median(zpts)
	for i in range(len(slis)):
		im = slis[i]
		fimg = pyfits.open(im, mode="update")
		fimg[0].header["RELZPT"] = zpts[i]
		fimg[0].header["RELZPTSC"] = scatt[i]
		fimg[0].header["RELZPRMS"] = rms[i]
		fimg[0].header["FLXSCALE"] = 1.0 / np.power(10, (zpts[i] - medzp) / 2.5)
		fimg.flush()
		
	# Final coadd
	os.system("swarp -GAIN_KEYWORD SOFTGAIN %s" % sstr)
	
	# Calibration
	iqpkg.iqobjs("coadd.fits", 10.0, RATIRSAT[chip], skyval="0.0",
	             pix=RATIRPIXSCALE[chip], gain=hdr["GAIN"], aperture=cpsfdiam,
	             wtimage="coadd.weight.fits")
	stars = Starlist("coadd.fits.stars")
	refstars = Starlist(reflist)
	refstars.wcs2pix("coadd.fits")
	refstars.set_mag("%sMAG" % filt.upper())
	truemags = np.array(refstars.mags())
	maglist = []; errlist = []
	
	a,b = stars.match(refstars, tol=10.0, maxnum=1000)
	newmags = np.zeros(len(refstars)); newwts = np.zeros(len(refstars))
	for i in range(len(refstars)):
		for j in range(len(a)):
			if b[j].mag == truemags[i]:
				newmags[i] = a[j].mag
				newwts[i] = 1.0 / np.power(np.maximum(a[j].magu, 0.02),2)
				continue
	maglist.append(newmags); errlist.append(newwts)
	
	obsmags = np.array(maglist)
	wts = np.array(errlist)
	
	[zpts, scatt, rms] = calc_zpts(truemags, obsmags, wts, sigma=3.0)
	fimg = pyfits.open("coadd.fits", mode="update")
	fimg[0].header["ABSZPT"] = zpts[0] + 25
	fimg[0].header["ABSZPTSC"] = scatt[0]
	fimg[0].header["ABZPRMS"] = rms[0]
	
	print "Zeropoint for coadded image: %.3f" % (25.0+zpts[0])
	print "Robust scatter for coadded image: %.3f" % scatt[0]
	print "RMS for coadded image: %.3f" % rms[0]
	
    
#####

def calc_zpts(truemags, obsmags, wts, sigma=3.0):

	nobs = obsmags.shape[0]
	nstars = obsmags.shape[1]
	zpt0 = np.ones(nobs)
	
	residual = lambda z, M, m, w: np.ravel(np.power(np.tile(M, (m.shape[0],1)) - m - np.ones((m.shape[0],m.shape[1])) * np.reshape(z, (m.shape[0],1)), 2) * w)
	
	result = optimize.leastsq(residual, zpt0, args=(truemags, obsmags, wts),
	                          full_output=1)

	modmags = obsmags + result[0].reshape(nobs,1).repeat(nstars,axis=1)
	adiff1 = np.tile(truemags, (nobs,1)) - modmags
	
	for i in range(nobs):
		newd1 = adiff1[i,:][np.where(wts[i,:]!=0)]
		scat1 = sMAD(newd1)
		for j in range(nstars):
			if np.abs(adiff1[i,j]-np.median(newd1)) > sigma * scat1:
				wts[i,j] = 0
	
	result2 = optimize.leastsq(residual, result[0], args=(truemags, obsmags, wts),
	                          full_output=1)
	modmags2 = obsmags + result2[0].reshape(nobs,1).repeat(nstars,axis=1)
	adiff2 = np.tile(truemags, (nobs,1)) - modmags2
	
	scats = np.zeros(nobs)
	rmss = np.zeros(nobs)
	for i in range(nobs):
		newd2 = adiff2[i,:][np.where(wts[i,:]!=0)]
		scat2 = sMAD(newd2)
		scats[i] = scat2
		rms2 = np.std(newd2)
		rmss[i] = rms2
		
	return [result2[0], scats, rmss]
	
#####

def sMAD(vec):
    '''Calculates a robust scatter using the median of the absolute deviaiton'''
    return 1.48*np.median(np.abs(vec-np.median(vec)))
	
	
	
	                          
	
    
    
