#######################################
# Routines for processing of LMI data
# SBC - Started 8 February 2014
#######################################

import pyraf
from pyraf import iraf
import astropy.io.fits as pyfits
import astropy.coordinates
import numpy as np
import glob, os, shutil

from vlt_autoastrometry import autoastrometry

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

LMIFILTS = ["U", "B", "V", "R", "I", "SDSS-G", "SDSS-R", "SDSS-I", "SDSS-Z"]
LMIPIXSCALE = 0.240
ASTROMCMD = "/Users/scenko/python/astrometry/vlt_autoastrometry.py"

##############################################################################

def preproc(image, fkey="FILTER", ppre="p", bkey="BIASSEC", tkey="TRIMSEC",
            clobber=globclob, verbose=globver):

	'''Update header keywords, subtract overscan, and trim image.'''
	
	# Create 'FILTER' keyword
	pyim = pyfits.open(image)
	f1 = pyim[0].header["FILTER1"]; f2 = pyim[0].header["FILTER2"]
	if f1 == "OPEN":
		pyim[0].header[fkey] = f2
	elif f2 == "OPEN":
		pyim[0].header[fkey] = f1
	else:
		pyim[0].header[fkey] = "%s-%s" % (f1, f2)
		
	# Grab bias and trim sections for ccdproc
	bsec = pyim[0].header[bkey]
	tsec = pyim[0].header[tkey]
	
	# Write out results
	pyim.writeto("%s%s" % (ppre, image), clobber=clobber)
	
	# Configure ccdproc
	ccdproc = iraf.ccdred.ccdproc
	ccdproc.ccdtype = ""
	ccdproc.noproc = no
	ccdproc.fixpix = no
	ccdproc.overscan = yes
	ccdproc.trim = yes
	ccdproc.zerocor = no
	ccdproc.darkcor = no
	ccdproc.flatcor = no
	ccdproc.illumcor = no
	ccdproc.fringecor = no
	ccdproc.readcor = no
	ccdproc.scancor = no
	ccdproc.readaxis = 'line'
	ccdproc.fixfile = ""
	ccdproc.biassec = bsec
	ccdproc.trimsec = tsec
	ccdproc.interactive = no
	ccdproc.function = "legendre"
	ccdproc.order = 1
	ccdproc.sample = "*"
	ccdproc.naverage = 1
	ccdproc.niterate = 1
	ccdproc.low_reject = 3.0
	ccdproc.high_reject = 3.0
	ccdproc.grow = 0
	ccdproc(images="%s%s" % (ppre,image), output="")
	
	return
	
#############################################################################

def lmi_cals(imlist, dobias=yes, dobpm=yes, doflats=yes, btype="BIAS", 
             ftype="SKY FLAT", fkey="FILTER", ppre="p", bkey="BIASSEC", 
             tkey="TRIMSEC", bfile="Bias.fits", bpmfile="BPM.pl", 
             bpre="b", flatpre="Flat", clobber=globclob, verbose=globver):

	'''Process bias and twilight flats, creating relevant calibration files for
	   nightly processing of LMI data.'''
	   
	images = glob.glob(imlist)   
	blist = []; flist = []
	
	# Preprocess all the relevant images
	for im in images:
		hdr = pyfits.getheader(im)
		if hdr['OBSTYPE'] == btype:
			blist.append("%s%s" % (ppre, im))
			preproc(im, fkey=fkey, ppre=ppre, bkey=bkey, tkey=tkey, clobber=clobber,
		    	    verbose=verbose)
		elif hdr['OBSTYPE'] == ftype:
			flist.append("%s%s" % (ppre, im))
			preproc(im, fkey=fkey, ppre=ppre, bkey=bkey, tkey=tkey, clobber=clobber,
		    	    verbose=verbose)
		    	    
	if dobias:
		
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
		zerocombine.rdnoise = hdr['RDNOISE']
		zerocombine.gain = hdr['GAIN']
		zerocombine.snoise = 0
		zerocombine.pclip = -0.5
		zerocombine.blank = 0.0
		
		# Run zerocombine
		bstr = ",".join(blist)
		if os.path.exists(bfile):
			os.remove(bfile)
		zerocombine(input=bstr, output=bfile)

	if dobpm:
	
		# Setup imcombine
		imcombine = iraf.immatch.imcombine
		imcombine.sigmas = "bsigma.fits"
		imcombine.combine = "median"
		imcombine.reject = "none"
		imcombine.project = no
		imcombine.outtype = "real"
		imcombine.outlimits = ""
		imcombine.offsets = "none"
		imcombine.masktype = "none"
		imcombine.scale = "none"
		imcombine.zero = "none"
		imcombine.weight = "none"
		
		# Run imcombine
		bstr = ",".join(blist)
		if os.path.exists("junk.fits"):
			os.remove("junk.fits")
		if os.path.exists("bsigma.fits"):
			os.remove("bsigma.fits")
		imcombine(input=bstr, output="junk.fits")

		# Run ccdmask to create BPM file
		if os.path.exists(bpmfile):
			os.remove(bpmfile)		
		iraf.ccdred.ccdmask("bsigma.fits", bpmfile, ncmed=7, nlmed=7, ncsig=15,
		                    nlsig=15, lsigma=10, hsigma=10, ngood=1, linterp=1, 
		                    cinterp=1, eqinterp=1)

		os.remove("bsigma.fits")
		os.remove("junk.fits")
		
	if doflats:
	
		# Subtract bias images from all flats
		for im in flist:
			lmi_debias(im, bfile=bfile)
			
		# Loop over filters
		for filt in LMIFILTS:
				
			# ID files in appropriate filter
			flis = []
			for im in flist:
				hdr = pyfits.getheader(im)
				if hdr[fkey] == filt:
					flis.append("%s%s" % (bpre, im))
					
			if flis == []:
				continue
					
			# Set up flatcombine
			flatcombine = iraf.ccdred.flatcombine
			flatcombine.combine = 'median'
			flatcombine.reject = 'avsigclip'
			flatcombine.ccdtype = ''
			flatcombine.process = no
			flatcombine.scale = 'median'
			flatcombine.statsec = ''
			flatcombine.nlow = 1
			flatcombine.nhigh = 1
			flatcombine.nkeep = 1
			flatcombine.mclip = yes
			flatcombine.lsigma = 3.0
			flatcombine.hsigma = 3.0
			flatcombine.rdnoise = hdr["RDNOISE"]
			flatcombine.gain = hdr["GAIN"]
			flatcombine.snoise = 0.0
			flatcombine.pclip = -0.5
			flatcombine.blank = 1.0
			
			# Run flatcombine
			fstr = ",".join(flis)
			if os.path.exists("%s-%s.fits" % (flatpre, filt)):
				os.remove("%s-%s.fits" % (flatpre, filt))
			flatcombine(input=fstr, output="%s-%s.fits" % (flatpre, filt))
			
			# Normalize
			iraf.iterstat.nsigrej = 5.0
			iraf.iterstat.maxiter = 10
			iraf.iterstat.verbose = globver
			iraf.iterstat.lower = INDEF
			iraf.iterstat.upper = INDEF
			iraf.iterstat("%s-%s.fits" % (flatpre, filt))
			iraf.imarith("%s-%s.fits" % (flatpre, filt), "/", 
			             iraf.iterstat.median, "%s-%s.fits" % (flatpre, filt))
			
	return
	
###########################################################################

def lmi_detrend(imlist, otype="OBJECT", ppre="p", bkey="BIASSEC", 
		        tkey="TRIMSEC", bpre="b", bfile="Bias.fits", fkey="FILTER", 
		        fpre="f", flatpre="Flat", clobber=globclob, verbose=globver):

	'''Identify science frames, pre-process, subtract bias, divide by flat, 
	   and correct for non-linearity.'''
	   
	images = glob.glob(imlist)   
	
	# Loop through all images
	for im in images:
		hdr = pyfits.getheader(im)
		if hdr['OBSTYPE'] == otype:
		
			# Preprocess
			preproc(im, fkey=fkey, ppre=ppre, bkey=bkey, tkey=tkey, clobber=clobber,
		    	    verbose=verbose)
		
			# Bias subtraction
			lmi_debias("%s%s" % (ppre, im), bpre=bpre, bfile=bfile)
			
			# Flat field
			lmi_flat("%s%s%s" % (bpre, ppre, im), fpre=fpre, fkey=fkey,
			         flatpre="Flat") 
			         
			# Linearity correction
			lmi_lincor("%s%s%s%s" % (fpre, bpre, ppre, im))
	
	return
	
###########################################################################

def lmi_lincor(image):

	'''Apply (pre-calculated) linearity correction for LMI data.'''
	
	fimg = pyfits.open(image, mode='update')
	fimg[0].data = 1.00305 * fimg[0].data - 1.20752e-6 * np.power(fimg[0].data, 2) + 1.267053e-11 * np.power(fimg[0].data, 3)
	fimg.flush()
	return
	
###########################################################################

def lmi_astrom(imlist, wpre="w"):

	'''WCS fits for images.'''
	
	images = glob.glob(imlist)   
	
	# Loop through all images
	for image in images:
	
		# First need to update a bunch of keywords
		fimg = pyfits.open(image)
		fimg[0].header["PIXSCALE"] = LMIPIXSCALE
		fimg[0].header["PIXSCAL1"] = LMIPIXSCALE
		fimg[0].header["PIXSCAL2"] = LMIPIXSCALE
		fimg[0].header["CTYPE1"] = "RA---TAN"
		fimg[0].header["CTYPE2"] = "DEC--TAN"
		fimg[0].header["WCSDIM"] = 2
		fimg[0].header["WAT0_001"] = "system=image"
		fimg[0].header["WAT1_001"] = "wtype=tan axtype=ra"
		fimg[0].header["WAT2_001"] = "wtype=tan axtype=dec"
		fimg[0].header["LTM1_1"] = 1.0
		fimg[0].header["LTM2_2"] = 1.0
	
		nax1 = fimg[0].header["NAXIS1"]; nax2 = fimg[0].header["NAXIS2"]
		fimg[0].header["CRPIX1"] = nax1 / 2
		fimg[0].header["CRPIX2"] = nax2 / 2
	
		ra = fimg[0].header["RA"]; dec = fimg[0].header["DEC"]
		coo = astropy.coordinates.ICRS("%sh%sm%ss %sd%sm%ss" % (ra.split(":")[0],
	    	                           ra.split(":")[1], ra.split(":")[2],
	        	                       dec.split(":")[0], dec.split(":")[1],
	            	                   dec.split(":")[2]))
		fimg[0].header["CRVAL1"] = coo.ra.deg
		fimg[0].header["CRVAL2"] = coo.dec.deg
	
		fimg[0].header["CD1_1"] = -LMIPIXSCALE / 3600.0
		fimg[0].header["CD1_2"] = 0.0
		fimg[0].header["CD2_1"] = 0.0
		fimg[0].header["CD2_2"] = LMIPIXSCALE / 3600.0
		del fimg[0].header["CTYPE1U"]
		del fimg[0].header["CRPIX1U"]
		del fimg[0].header["CRVAL1U"]
		del fimg[0].header["CD1_1U"]
		del fimg[0].header["CFINT1"]
		del fimg[0].header["CTYPE2U"]
		del fimg[0].header["CRPIX2U"]
		del fimg[0].header["CRVAL2U"]
		del fimg[0].header["CD2_2U"]
		del fimg[0].header["CD1_2U"]
		del fimg[0].header["CD2_1U"]
		del fimg[0].header["CFINT2"]
		
		if os.path.exists("%s%s" % (wpre, image)):
			os.remove("%s%s" % (wpre, image))
		
		fimg.writeto("%s%s" % (wpre, image))
	
		os.system("python %s %s%s" % (ASTROMCMD, wpre, image))  
	
		# If successful, rename image
		if os.path.exists("a%s%s" % (wpre, image)):
			shutil.move("a%s%s" % (wpre, image), "%s%s" % (wpre, image))
	
	return
   
	
###########################################################################

def lmi_flat(image, fpre="f", fkey="FILTER", flatpre="Flat"):

	'''Flat-field an LMI image.'''
	
	# First need to identify correct flat
	hdr = pyfits.getheader(image)
	if not os.path.exists("%s-%s.fits" % (flatpre, hdr[fkey])):
		print "Error: No Flat for specified filter: %s" % hdr[fkey]
		return
		
	ccdproc = iraf.ccdred.ccdproc
	ccdproc.ccdtype = ""
	ccdproc.noproc = no
	ccdproc.fixpix = no
	ccdproc.overscan = no
	ccdproc.trim = no
	ccdproc.zerocor = no
	ccdproc.darkcor = no
	ccdproc.flatcor = yes
	ccdproc.illumcor = no
	ccdproc.fringecor = no
	ccdproc.readcor = no
	ccdproc.scancor = no
	ccdproc.readaxis = 'line'
	ccdproc.fixfile = ""
	ccdproc.flat = "%s-%s.fits" % (flatpre, hdr[fkey])
	ccdproc.interactive = no
	
	if os.path.exists("%s%s" % (fpre, image)):
		os.remove("%s%s" % (fpre, image))
	ccdproc(images=image, output="%s%s" % (fpre, image))
	
	return
			
###########################################################################

def lmi_debias(image, bpre="b", bfile="Bias.fits", clobber=globclob):

	'''Subtract bias from LMI image.'''
	
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
	ccdproc.zero = bfile
	ccdproc.interactive = no
	ccdproc.function = "legendre"
	ccdproc.order = 1
	ccdproc.sample = "*"
	ccdproc.naverage = 1
	ccdproc.niterate = 1
	ccdproc.low_reject = 3.0
	ccdproc.high_reject = 3.0
	ccdproc.grow = 0
	
	if os.path.exists("%s%s" % (bpre, image)):
		os.remove("%s%s" % (bpre, image))
	ccdproc(images=image, output="%s%s" % (bpre, image))
	
	return
	
###########################################################################

def lmi_coadd(object, filt):

    '''Coadd all images of a given object in a specified filter.'''
    
    # Find the appropriate images
    ims = []
    allims = glob.glob("wfbplmi.????.fits")
    for im in allims:
    	h = pyfits.open(im)
    	if (h[0].header["OBJECT"] == object) and (h[0].header["FILTER"] == filt):
    		ims.append(im)    
    
    if len(ims) == 0:
    	print "No images to coadd; Exiting!"
    	return
    	
    # Setup configuration files
    o1 = open("daofind.param", "w")
    o1.write("NUMBER\nXWIN_IMAGE\nYWIN_IMAGE\nMAG_AUTO\nFLAGS\nA_IMAGE\nB_IMAGE\n")
    o1.write("ELONGATION\nFWHM_IMAGE\nCLASS_STAR\nXWIN_WORLD\nYWIN_WORLD\n")
    o1.write("ERRAWIN_IMAGE\nERRBWIN_IMAGE\nERRTHETAWIN_IMAGE\nERRAWIN_WORLD\n")
    o1.write("ERRBWIN_WORLD\nERRTHETAWIN_WORLD\nFLUX_AUTO\nFLUX_RADIUS\nFLUXERR_AUTO")
    o1.close()
    
    o2 = open("default.conv", "w")
    o2.write("CONV NORM\n# 5x5 convolution mask of a gaussian PSF with FWHM = 3.0 pixels.\n0.092163 0.221178 0.296069 0.221178 0.092163\n0.221178 0.530797 0.710525 0.530797 0.221178\n0.296069 0.710525 0.951108 0.710525 0.296069\n0.221178 0.530797 0.710525 0.530797 0.221178\n0.092163 0.221178 0.296069 0.221178 0.092163")
    o2.close()
	    
    # First run sextractor
    for im in ims:
    	os.system("sex -CATALOG_NAME %s.cat -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME daofind.param -DETECT_THRESH 3.0 -ANALYSIS_THRESH 3.0 -GAIN_KEY GAIN -PIXEL_SCALE 0 %s" % (im[:-5], im))
    	
    # Run scamp for alignment
    imlist = " ".join(ims)
    os.system("scamp -CHECKPLOT_DEV NULL %s" % imlist.replace(".fits", ".cat"))
    
    # Run swarp for coaddition
    os.system("swarp -IMAGEOUT_NAME %s-%s.fits %s" % (object, filt, imlist))
    
    # Remove junk files
    os.system("rm *.head")
    
    return
			
			          
