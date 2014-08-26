"""
	Purpose:	this is a collection of preprocessing functions for use with data from RATIR.

	Usage:
		1)	enter python or ipython environment
		2)	can load all functions using:
			- "from rat_preproc import *" if you want to call functions using just the function's name
			- "import rat_preproc as rp" if you want to call functions using rp.function_name( args )

	Notes:
		* small memory leak in ratdisp and ratdisp_calib

	Future Improvements:
		- 
"""

import os
from fnmatch import fnmatch
import astropy.io.fits as pf
import numpy as np
import matplotlib.pylab as pl
from scipy.ndimage.interpolation import zoom
import shutil
from glob import glob
import gc

import astro_functs as af # contains basic functions and RATIR constants
from astro_functs import show_list # allow user to call show_list without "af." prefix

# Preprocessing constants
FITS_IN_KEY = lambda n: 'IMCMB{:03}'.format(int(n)) # function to make FITS keywords to store file names of combined frames

"""
	Translated from ratlist.sh by John Capone (jicapone@astro.umd.edu).

	Purpose:			To create lists of detected fits files in each subdir of the specified working dir. Each camera has it's own list file.

	Input:
		workdir:		directory with subdirectories containing data (ex. visit numbers)
		cams:			lists created for specified camera numbers

	Usage:
		1)	enter python or ipython environment
		2)	load function -> 'from rat_preproc import ratlist'
		3)	run function -> 'ratlist( workdir = 'path/to/data/', cams = [#,#,...] )'
			- since default workdir is '.', this argument can be ignored if you are in the same directory as the data
		4)	creates lists specifying frames for each camera

	Notes:
		- added option to specify working directory
		- fixed matching so only raw files are listed

	Future Improvements:
		- 
"""
def ratlist( workdir='.', cams=[0,1,2,3] ):

	# move to working directory
	start_dir = os.getcwd()
	os.chdir( workdir ) # move into working dir
	
	d = os.getcwd().split('/')[-1] # name of current directory
	print "\n* * * creating lists of frames in {} * * *".format(d)
	for i in range(len(cams)):
		flist = open( 'C{}_{}.list'.format( cams[i], d ), 'w' )
		matches = glob( '*C{}?.fits'.format(cams[i]) )
		for f in matches: flist.write( f + '\n' ) # add detected file to correct list file
		flist.close() # close list file

	os.chdir( start_dir ) # move back to starting dir

"""
	Translated from plotratir_flat.pro by John Capone (jicapone@astro.umd.edu).
	Modified by Vicki Toy 8/15/14

	Purpose:		display calibration images for verification by user

	Input:
		ftype:		type of frames. defaults for af.FLAT_NAME
		workdir:	directory where function is to be executed
		cams:		camera numbers.  all by default
		auto:		automated selection of frames.  if ftype is af.BIAS_NAME, select all.  if ftype is af.FLAT_NAME, select non-saturated frames with sufficient counts.
		amin:		minimum median value for automated selection as fraction of saturation value
		amax:		maximum median value for automated selection as fraction of saturation value

	Usage:
		1)	enter python or ipython environment
		2)	load function -> 'from rat_preproc import ratdisp_calib'
		3)	run function -> 'ratdisp_calib( ftype = bias or flat name, workdir = 'path/to/data/', cams = [#,#,...] )'
			- since default workdir is '.', this argument can be ignored if you are in the same directory as the data
		4)	select which frames you would like to use to create master calibration frames

	Notes:
		- added option to specify working directory
		- added error handling for if a camera list is missing
		- camera argument can now be int
				- ccd lists are now by filter rather than camera number
		- automated bias and flat frame selection can be done relative to the detectors' staturation points

	Future Improvements:
		- better way to do cameras with multiple filters (current way uses camera # -2, won't work for other instruments)
"""
def ratdisp_calib( ftype=af.FLAT_NAME, workdir='.', cams=[0,1,2,3], overwrite=True, auto=False, amin=0.1, amax=0.8 ):

	# check for non-list camera argument
	if type(cams) is not list:
		cams = [cams] # convert to list
	
	pl.ion() # pylab in interactive mode

	# move to working directory
	start_dir = os.getcwd()
	os.chdir( workdir )
	
	d = os.getcwd().split('/')[-1] # name of current directory
	if not auto:
		print "\n* * * displaying {} frames in {} for selection * * *".format(ftype, d)
	else:
		print "\n* * * automatically selecting {} frames in {} * * *".format(ftype, d)

	# delete existing calibration lists
	fntemp = glob( ftype+'*.list' )
	for fn in fntemp: os.remove( fn )
	

	
	# work on FITs files for specified cameras
	for cam_i in cams:

		# print current camera number
		print "\n* * * CAMERA {} * * *".format( cam_i )
		if af.CAM_BIAS[cam_i] == False and ftype is af.BIAS_NAME:
			print "\t* Warning: Camera C{} does not have {} frames.  Skipping...".format( cam_i, af.BIAS_NAME )
			continue
			
		fn_list = '{}_{}.list'.format( af.CAM_NAMES[cam_i], d )
		try:
			fin = open( fn_list, 'r' ) # open list of FITs files
		except IOError:
			print "Warning: {} not found.  Skipping camera {}.".format( fn_list, cam_i )
		else:

			# look at FITs data sequentially
			count = -1
			for line in fin:
				count = count + 1
					
				fits_fn = line.rstrip() # current fits file name with return removed
				fits_id = fits_fn.split('.')[0] # fits file name with extention removed
				print '{}'.format( fits_fn )

				# open data
				hdulist = pf.open( fits_fn )
				im = hdulist[0].data
				h  = hdulist[0].header
					
				if af.CAM_SPLIT[cam_i]:
					#Overwrites list files that will contain new information
					if count == 0:
						fout = [open( '{}_{}.list'.format( ftype, af.SPLIT_FILTERS[cam_i-2] ), 'w' ), open( '{}_{}.list'.format( ftype, af.SPLIT_FILTERS[cam_i] ), 'w' )] # create list files for new img FITs files (e, w)	
						for f in fout: f.close()
					im1 = im[af.SLICES[af.SPLIT_FILTERS[cam_i-2]]]
					m1  = np.median( im1 )
					s1  = af.robust_sigma( im1 )
					im2 = im[af.SLICES[af.SPLIT_FILTERS[cam_i]]]
					m2  = np.median( im2 )
					s2  = af.robust_sigma( im2 )
					print '\t* Median of left side is {} counts'.format( m1 )
					print '\t* Median of right side is {} counts'.format( m2 )	
									
				else:
					im1 = im[af.SLICES['C'+str(cam_i)]]
					m   = np.median( im1 )
					s   = af.robust_sigma( im1 )
					print '\t* Median is {} counts.'.format( m )								
					print '\t* Filter used: {}'.format( h['FILTER'] )	
					
					#Overwrites list files that will contain new information
					if ftype is af.FLAT_NAME and count == 0:
						fout = open( '{}_{}.list'.format( ftype, h['FILTER'] ), 'w' ) # append if this filter's list exists	
						fout.close()
					if ftype is af.BIAS_NAME and count == 0:
						fout = open( '{}_{}.list'.format( ftype, cam_i), 'w' ) # append if this filter's list exists	
						fout.close()											
				if auto:
					if ftype is af.BIAS_NAME:	
						# all bias frames are selected
						imfits = '{}_{}.fits'.format( fits_id, ftype)
						pf.writeto( imfits, im1, header=h, clobber=True ) # save object frame
						fout = open( '{}_{}.list'.format( ftype, cam_i ), 'a' ) # append if this filter's list exists
						fout.write( '{}_{}\n'.format( fits_id, ftype ) ) # write new file name to list
						fout.close()
													
					if ftype is af.FLAT_NAME:
							
						sat_pt = af.CAM_SATUR[cam_i](h['SOFTGAIN'])
						vmin = amin * sat_pt; vmax = amax * sat_pt							
									
						if af.CAM_SPLIT[cam_i]:
							fout = [open( '{}_{}.list'.format( ftype, af.SPLIT_FILTERS[cam_i-2] ), 'a' ), open( '{}_{}.list'.format( ftype, af.SPLIT_FILTERS[cam_i] ), 'a' )] # create list files for new img FITs files (e, w)	
							# check whether median values are in specified range
							# bottom side
							if m1 > vmin and m1 < vmax:
								print "\t* Bottom side selected."
								imfits_1 = '{}_{}_{}.fits'.format( fits_id, ftype, af.SPLIT_FILTERS[cam_i-2] )
								h['FILTER'] = af.SPLIT_FILTERS[cam_i-2]
								pf.writeto( imfits_1, im1, header=h, clobber=True ) # save left frame
								fout[0].write( '{}_{}_{}\n'.format( fits_id, ftype, af.SPLIT_FILTERS[cam_i-2] ) ) # write new file name to list
								fout[0].close()
							else:
								if m1 < vmin:
									print "\t* Bottom side rejected:\tUNDEREXPOSED."
								else:
									print "\t* Bottom side rejected:\tSATURATED."

							# top side
							if m2 > vmin and m2 < vmax:
								print "\t* Top side selected."
								imfits_2 = '{}_{}_{}.fits'.format( fits_id, ftype, af.SPLIT_FILTERS[cam_i] )
								h['FILTER'] = af.SPLIT_FILTERS[cam_i]
								pf.writeto( imfits_2, im2, header=h, clobber=True ) # save object frame
								fout[1].write( '{}_{}_{}\n'.format( fits_id, ftype, af.SPLIT_FILTERS[cam_i] ) ) # write new file name to list
								fout[1].close()
							else:
								if m2 < vmin:
									print "\t* Top side rejected:\tUNDEREXPOSED."
								else:
									print "\t* Top side rejected:\tSATURATED."
						#Not split frame			
						else:
							
							# check whether median value is in specified range
							if m > vmin and m < vmax:
								print "\t* Frame selected."
								imfits = '{}_{}.fits'.format( fits_id, ftype)
								pf.writeto( imfits, im1, header=h, clobber=True ) # save object frame
								fout = open( '{}_{}.list'.format( ftype, h['FILTER'] ), 'a' ) # append if this filter's list exists
								fout.write( '{}_{}\n'.format( fits_id, ftype ) ) # write new file name to list
								fout.close()
							else:
								if m < vmin:
									print "\t* Frame rejected:\tUNDEREXPOSED."
								else:
									print "\t* Frame rejected:\tSATURATED."

				# display image and prompt user
				else:
					if af.CAM_SPLIT[cam_i]:
						pl.subplot(211)
						axim = pl.imshow( im1, vmin=m1-5*s1, vmax=m1+5*s1, origin='lower', cmap=pl.cm.gray )
						axim.axes.set_xticks([])
						axim.axes.set_yticks([])
						pl.title( r"Median = {}, $\sigma$ = {:.1f}".format( int(m1), s1 ) )
						pl.subplot(212)
						axim = pl.imshow( im2, vmin=m2-5*s2, vmax=m2+5*s2, origin='lower', cmap=pl.cm.gray )
						axim.axes.set_xticks([])
						axim.axes.set_yticks([])
						pl.title( r"Median = {}, $\sigma$ = {:.1f}".format( int(m1), s1 ) )						
					else:
						axim = pl.imshow( im1, vmin=m-5*s, vmax=m+5*s, origin='lower', cmap=pl.cm.gray )
						axim.axes.set_xticks([])
						axim.axes.set_yticks([])						
						pl.title( r"Median = {}, $\sigma$ = {:.1f}".format( int(m), s ) )
							
					# query user until valid response is provided
					valid_entry = False
					while not valid_entry:
						user = raw_input("\nType Y for YES, N for NO, Q for QUIT: ")
								
						if user.lower() == 'y':
									
							if af.CAM_SPLIT[cam_i]:
								imfits_left = '{}_{}_{}.fits'.format( fits_id, ftype, af.SPLIT_FILTERS[cam_i-2] )
								h['FILTER'] = af.SPLIT_FILTERS[cam_i-2]
								pf.writeto( imfits_left, imleft, header=h, clobber=True ) # save object frame
								fout[0].write( '{}_{}_{}\n'.format( fits_id, ftype, af.SPLIT_FILTERS[cam_i-2] ) ) # write new file name to list
								fout[0].close()
									
								imfits_right = '{}_{}_{}.fits'.format( fits_id, ftype, af.SPLIT_FILTERS[cam_i] )
								h['FILTER'] = af.SPLIT_FILTERS[cam_i]
								pf.writeto( imfits_right, imright, header=h, clobber=True ) # save object frame
								fout[1].write( '{}_{}_{}\n'.format( fits_id, ftype, af.SPLIT_FILTERS[cam_i] ) ) # write new file name to list
								fout[1].close()
							else:
								imfits = '{}_{}.fits'.format( fits_id, ftype )
								fout = open( '{}_{}.list'.format( ftype, h['FILTER'] ), 'a' ) # append if this filter's list exists
								fout.write( '{}_{}\n'.format( fits_id, ftype ) ) # write new file name to list
								fout.close()
									
							valid_entry = True
								
						elif user.lower() == 'q': # exit function
							print "Exiting..."
							os.chdir( start_dir ) # move back to starting directory
							pl.close('all') # close image to free memory
							return
								
						elif user.lower() != 'n': # invalid case
							print "'{}' is not a valid entry.".format( user )
								
						else: # 'N' selected, skip
							valid_entry = True

				hdulist.close() # close FITs file
				pl.close('all') # close image to free memory

			# close files
			fin.close()
	# move back to starting directory
	os.chdir( start_dir )

"""
	Translated from plotratir.pro by John Capone (jicapone@astro.umd.edu).
	Modified by Vicki Toy 8/14/14

	Purpose:	display RATIR images for verification by user

	Input:
		workdir:	directory where function is to be executed
		targetdir:	directory where selected frames and lists are output
		cams:		camera numbers.  all by default
		auto:		select all science frames

	Usage:
		1)	enter python or ipython environment
		2)	load function -> 'from rat_preproc import ratdisp'
		3)	run function -> 'ratdisp( workdir = 'path/to/data/', targetdir = 'path/to/new data/', cams = [#,#,...] )'
			- since default workdir is '.', this argument can be ignored if you are in the same directory as the data
		4)	select which frames you would like to use

	Notes:
		- added option to specify working directory
		- added error handling for if a camera list is missing
		- camera argument can now be int
		- added target directory option.  sky and object frames are written to the same dir.
		- prompts user for overwrite
				- ccd lists are now by filter rather than camera number
		- added GAIN and SATURATE keywords to headers
		- added automated option
		- removed frame rotation and added WCS keywords

	Future Improvements:
		- automation of frame selection
			- view 20ish automatically selected frames at a time
		- better way to do cameras with multiple filters (current way uses camera # -2, won't work for other instruments)
"""
def ratdisp( workdir='.', targetdir='.', cams=[0,1,2,3], auto=False ):

	# check for non-list camera argument
	if type(cams) is not list:
		cams = [cams] # convert to list
	
	pl.ion() # pylab in interactive mode

	# move to working directory
	start_dir = os.getcwd()
	os.chdir( workdir )
	
	d = os.getcwd().split('/')[-1] # name of current directory
	if not auto:
		print "\n* * * displaying science frames in {} for selection * * *".format(d)
	else:
		print "\n* * * selecting all science frames in {} * * *".format(d)

	# make target directory if it does not exist
	if not os.path.exists( targetdir ):
		print "Creating target directory: ", targetdir
		os.makedirs( targetdir )
	# warn user if previous files may be overwritten
	else:
		print "Warning: Target directory exists. Existing files will be overwritten."
		resp = raw_input( "Proceed? (y/n): " )
		if resp.lower() != 'y':
			print "Exiting..."
			os.chdir( start_dir ) # move back to starting directory
			return
	
	# remove tailing / from target directory name if present
	if targetdir[-1] == '/':
		targetdir = targetdir[:-1]

	# delete existing object and sky lists
	fntemp = glob( '{}/{}*.list'.format( targetdir, af.OBJ_NAME ) )
	for fn in fntemp: os.remove( fn )

	# work on FITs files for specified cameras
	for cam_i in cams:
		
		# print current camera number
		print "\n* * * CAMERA {} * * *".format( cam_i )
		
		fn_list = '{}_{}.list'.format( af.CAM_NAMES[cam_i], d )
		try:
			fin = open( fn_list, 'r' ) # open list of FITs files
		except IOError:
			print "Warning: {} not found.  Skipping camera {}.".format( fn_list, cam_i )
		else:
				
			# look at FITs data sequentially
			for line in fin:

				fits_fn = line.rstrip() # current fits file name with return removed
				fits_id = fits_fn.split('.')[0] # fits file name with extention removed
				print '\n{}'.format( fits_fn )
					
				# open data
				hdulist = pf.open( fits_fn )
				im = hdulist[0].data
				h = hdulist[0].header

				# check for required header keywords
				if 'PRPSLID' in h:
					prpslid = h['PRPSLID']
				else:
					"ERROR: ratdisp - PRPSLID not found in fits header."
					os.chdir( start_dir ) # move back to starting directory
					pl.close('all') # close image to free memory
					return -1
					
				if 'VSTID' in h:
					vstid = h['VSTID']
				else:
					"ERROR: ratdisp - VSTID keyword not found in fits header."
					os.chdir( start_dir ) # move back to starting directory
					pl.close('all') # close image to free memory
					return -1
					
				if af.CENTER_KEY in h:
					center = h[af.CENTER_KEY].split('center')[0]
				else:
					"ERROR: ratdisp - {} keyword not found in fits header.".format( af.CENTER_KEY )
					os.chdir( start_dir ) # move back to starting directory
					pl.close('all')
					return	
								
				targname = '{}-vis{}'.format( prpslid, vstid )				
									
				# get image statistics
				if af.CAM_SPLIT[cam_i]:
					im1 = im[af.SLICES[af.SPLIT_FILTERS[cam_i-2]]]
					m1  = np.median( im1 )
					s1  = af.robust_sigma( im1 )
					im2 = im[af.SLICES[af.SPLIT_FILTERS[cam_i]]]
					m2  = np.median( im2 )
					s2  = af.robust_sigma( im2 )
				else:
					im1 = im[af.SLICES['C'+str(cam_i)]]				
					m = np.median( im1 )
					s = af.robust_sigma( im1 )
						
				# display image and prompt user
				if not auto:
				
					if af.CAM_SPLIT[cam_i]:
					
						pl.figure(figsize=(10,10))
						pl.subplot(211)
						axim = pl.imshow( im1, vmin=m1-5*s1, vmax=m1+5*s1, origin='lower', cmap=pl.cm.gray )
						axim.axes.set_xticks([])
						axim.axes.set_yticks([])
						pl.title( r"{} band".format( af.SPLIT_FILTERS[cam_i-2] ) + ': ' + "Median = {}, $\sigma$ = {:.1f}".format( int(m1), s1 ) )
						pl.subplot(212)
						axim = pl.imshow( im2, vmin=m2-5*s2, vmax=m2+5*s2, origin='lower', cmap=pl.cm.gray )
						axim.axes.set_xticks([])
						axim.axes.set_yticks([])
						pl.title( r"{} band".format( af.SPLIT_FILTERS[cam_i] ) + ': ' + "Median = {}, $\sigma$ = {:.1f}".format( int(m2), s2 ) )
					else:					
						axim = pl.imshow( im1, vmin=m-5*s, vmax=m+5*s, origin='lower', cmap=pl.cm.gray )
						axim.axes.set_xticks([])
						axim.axes.set_yticks([])
						pl.title( r"{} band: Median = {}, $\sigma$ = {:.1f}".format( h['FILTER'], int(m), s ) )
				
				if af.CAM_SPLIT[cam_i]:
					if center.count( af.SPLIT_FILTERS[cam_i] ) != 0:
						print "\t* The target is focused on the {} filter.".format( af.SPLIT_FILTERS[cam_i] )
					elif center.count( af.SPLIT_FILTERS[cam_i-2] ) != 0:
						print "\t* The target is focused on the {} filter.".format( af.SPLIT_FILTERS[cam_i-2] )
					else:
						print "\t* Warning: The target is NOT focused on an H2RG filter. The target is focused on the {} filter.".format( center )
					### ###
				else:
					# print filter name
					print '\t* Filter used: {}'.format( h['FILTER'] )

				# query user until valid response is provided
				valid_entry = False
				while not valid_entry:
						
					# either select all if auto, or have user select
					if auto:
						user = 'y'
					else:
						user = raw_input("\nType Y for YES, N for NO, Q for QUIT: ")
					
					if user.lower() == 'y' and af.CAM_SPLIT[cam_i]:	
						if center.count( af.SPLIT_FILTERS[cam_i] ) != 0:
							direction = 't'
						elif center.count( af.SPLIT_FILTERS[cam_i-2] ) != 0:
							direction = 'b'
						else:
							print  "\t* Warning: Skipping frame not centered on H2RG filter."
							user = 'n'
							direction = ''				
						
					if user.lower() == 'y':
							
						# set keyword values
						h['CAMERA']   = cam_i
						h['TARGNAME'] = targname
						h['PIXSCALE'] = af.CAM_PXSCALE[cam_i]
						h['WAVELENG'] = af.CAM_WAVE[cam_i]
						h['GAIN']     = (af.CAM_GAIN[cam_i]( h['SOFTGAIN'] ), 'in electrons/DN')
						h['SATURATE'] = (af.CAM_SATUR[cam_i]( h['SOFTGAIN'] ), 'in electrons/DN')
						h['CRPIX1']   = af.CAM_X0[cam_i]
						h['CRPIX2']   = af.CAM_Y0[cam_i]
						h['CTYPE1']   = 'RA---TAN'
						h['CTYPE2']   = 'DEC--TAN'
						h['CD1_1']    =  -af.CAM_SECPIX1[cam_i]*np.cos(af.CAM_THETA[cam_i]*np.pi/180.0)/3600.
						h['CD2_1']    =   af.CAM_SECPIX1[cam_i]*np.sin(af.CAM_THETA[cam_i]*np.pi/180.0)/3600.
						h['CD1_2']    =   af.CAM_SECPIX2[cam_i]*np.sin(af.CAM_THETA[cam_i]*np.pi/180.0)/3600.
						h['CD2_2']    =   af.CAM_SECPIX2[cam_i]*np.cos(af.CAM_THETA[cam_i]*np.pi/180.0)/3600.  
						h['CRVAL1']   =  h[af.RA_KEY]  - af.APOFFS[h[af.CENTER_KEY]][0]/60.0 + h[af.OFFRA_KEY] #includes aperture offsets and target offsets (ie. dithering)
						h['CRVAL2']   =  h[af.DEC_KEY] - af.APOFFS[h[af.CENTER_KEY]][1]/60.0 + h[af.OFFDEC_KEY]

						if af.CAM_SPLIT[cam_i]:
							for key in af.H2RG_ASTR:
								h[key] = af.H2RG_ASTR[key]
							
							if direction.lower() == 'b':
								f_img = cam_i-2
								f_sky = cam_i
							else:
								f_img = cam_i
								f_sky = cam_i-2
						
						if af.CAM_SPLIT[cam_i]:
							imfits = '{}/{}_{}_{}.fits'.format( targetdir, fits_id, af.OBJ_NAME, af.SPLIT_FILTERS[f_img] )
							h['FILTER'] = af.SPLIT_FILTERS[f_img]
							im_img = im[af.SLICES[h['FILTER']]]
							h['NAXIS1'] = af.SLICES[h['FILTER']][1].stop - af.SLICES[h['FILTER']][1].start
							h['NAXIS2'] = af.SLICES[h['FILTER']][0].stop - af.SLICES[h['FILTER']][0].start
							h['CRPIX1'] = af.CAM_X0[cam_i] - af.SLICES[h['FILTER']][1].start
							h['CRPIX2'] = af.CAM_Y0[cam_i] - af.SLICES[h['FILTER']][0].start
							pf.writeto( imfits, im_img, header=h, clobber=True ) # save object frame
							
							# filter side with sky, now saved as object, but different list to keep track
							skyfits = '{}/{}_{}_{}.fits'.format( targetdir, fits_id, af.OBJ_NAME, af.SPLIT_FILTERS[f_sky] )
							h['FILTER'] = af.SPLIT_FILTERS[f_sky]
							im_sky = im[af.SLICES[h['FILTER']]]
							h['NAXIS1'] = af.SLICES[h['FILTER']][1].stop - af.SLICES[h['FILTER']][1].start #Repeat incase filter sizes different
							h['NAXIS2'] = af.SLICES[h['FILTER']][0].stop - af.SLICES[h['FILTER']][0].start
							h['CRPIX1'] = af.CAM_X0[cam_i] - af.SLICES[h['FILTER']][1].start
							h['CRPIX2'] = af.CAM_Y0[cam_i] - af.SLICES[h['FILTER']][0].start
							pf.writeto( skyfits, im_sky, header=h, clobber=True ) # save sky frame
							valid_entry = True
						else:
							imfits = '{}/{}_{}_{}.fits'.format( targetdir, fits_id, af.OBJ_NAME, cam_i )
							im_img = im[af.SLICES['C'+str(cam_i)]]
							h['NAXIS1'] = af.SLICES['C'+str(cam_i)][1].stop - af.SLICES['C'+str(cam_i)][1].start
							h['NAXIS2'] = af.SLICES['C'+str(cam_i)][0].stop - af.SLICES['C'+str(cam_i)][0].start
							h['CRPIX1'] = af.CAM_X0[cam_i] - af.SLICES['C'+str(cam_i)][1].start
							h['CRPIX2'] = af.CAM_Y0[cam_i] - af.SLICES['C'+str(cam_i)][0].start
							pf.writeto( imfits, im_img, header=h, clobber=True )
							valid_entry = True						
						
					elif user.lower() == 'q': # exit function
						print "Exiting..."
						os.chdir( start_dir ) # move back to starting directory
						pl.close('all') # close image to free memory
						return
						
					elif user.lower() != 'n': # invalid case
						print "'{}' is not a valid entry.".format( user )
						
					else: # 'N' selected, skip
						valid_entry = True

				hdulist.close() # close FITs file
				pl.close('all') # close image to free memory

			# close files
			fin.close()
	# move back to starting directory
	os.chdir( start_dir )

"""
	Written by John Capone (jicapone@astro.umd.edu).

	Purpose:		bias list comes directly from ratlist, so the list name must be changed to work with mkmaster. this function is called from mkmaster and should not be called by the user.

	Input:
		cams:		list of camera numbers
		workdir:	directory where function is to be executed

	Usage:
		1)	not intended for useage by user.  should be called by mkmaster()

	Notes:
		- 

	Future Improvements:
		- 
"""
def _prep_bias_list( cams, workdir='.' ):

	# move to working directory
	start_dir = os.getcwd()
	os.chdir( workdir )

	# rename bias list
	for cam in cams:
		start_fn = 'C{}_{}.list'.format( cam, os.getcwd().split('/')[-1] )
		end_fn = '{}_{}.list'.format( af.BIAS_NAME, cam )

		if os.path.exists( end_fn ): # do nothing if list already exists
			os.chdir( start_dir ) # move back to starting directory
			return 
		elif os.path.exists( start_fn ):
			shutil.copy( start_fn, end_fn )
		else:
			d = os.getcwd().split('/')[-1] # name of current directory
			print "Warning: list of bias frames for camera {} was not found in {}.".format( cam, d )

	# move back to starting directory
	os.chdir( start_dir )

"""
	Written by John Capone (jicapone@astro.umd.edu).

	Purpose:		make master bias and master flat frames
					* currently no outlier rejection

	Input:
		mtype:		type of master frame. should be either af.FLAT_NAME or af.BIAS_NAME
		bands:		list of photometric bands or 'ALL'.  if 'ALL', will search for lists and create master frames for all lists.
		workdir:	directory where function is to be executed
		fmin:		minimum number of files needed to make a master frame

	Usage:
		1)	enter python or ipython environment
		2)	load function -> 'from rat_preproc import mkmaster'
		3)	run function -> 'mkmaster( mtype = bias or flat name, workdir = 'path/to/data/', bands = [0,1, and/or H2RG filer names] )'
			- since default workdir is '.', this argument can be ignored if you are in the same directory as the data
			- bands uses camera numbers for the CCDs since the filters change.  H2RG filters are static, so the names specified in the configuration file are used.
		4)	create a master configuration file using the frames previously selected by ratdisp_calib()

	Notes:
		- checks that data being combined used same filter
		- added filter keyword master frame's header
		- added fmin parameter.  allows user to abort if fewer files are found to make a master frame.
		- can now set bands to 'ALL' rather than inputing a list of desired bands.  in this case, script will look for available lists and make a master for each.

	Future Improvements:
		- may want to be generalized for use with dark current frames
		- 
"""
def mkmaster( mtype, bands='ALL', workdir='.', fmin=5 ):

	# move to working directory
	start_dir = os.getcwd()
	os.chdir( workdir )

	# check for valid mtype
	if mtype not in [af.FLAT_NAME, af.BIAS_NAME]:
		print "Error: valid arguments for mtype are", af.FLAT_NAME, "and", af.BIAS_NAME, ". Exiting..."
		os.chdir( start_dir ) # move back to starting directory
		return

	# if making masters for all bands, search for lists
	if isinstance(bands, str):
		if bands.upper() == 'ALL':
			ltemp = glob( '{}_*.list'.format(mtype) )
			bands = []
			for f in ltemp: bands.append( f.split('_')[1].split('.')[0] )
	elif type(bands) is not list:
		print "Error: bands must either be \"ALL\" or a list of desired filters."
		os.chdir( start_dir ) # move back to starting directory
		return

	sorttype = 'BAND'

	if mtype is af.BIAS_NAME:
		_prep_bias_list( cams=bands, workdir='.' )
		sorttype = 'CAMERA' 
		

	d = os.getcwd().split('/')[-1] # name of current directory
	print "\n* * * making master {} frame in {} * * *".format( mtype, d )

	# work on FITs files for specified photometric bands
	for band in bands:

		# print current band
		print "\n* * * {} {} * * *".format( band, sorttype )

		flist = open( '{}_{}.list'.format( mtype, band ), 'r' )
		hdu = pf.PrimaryHDU()
		data_arr = []
		filter_arr = [] # to check that all frames used the same filter
		i = 0
		fs = flist.readlines()
		cont = True
		if len(fs) < fmin:
			if len(fs) == 0:
				print 'No frames available to make master {} for {} {}.'.format( mtype, band, sorttype.lower() )
				cont = False
			else:
				temp = raw_input( "Only {} frames available to make master {} for {} {}.  Continue? (y/n): ".format( len(fs), mtype, band, sorttype.lower() ) )
			
				if temp.lower() == 'y' or temp.lower() == 'yes':
					cont = True
				else:
					cont = False
		if cont:
			for f in fs:
				fits_fn = f.rstrip() # current fits file name with return removed
				fits_id = fits_fn.split('.')[0] # fits file name with extention removed
				fn = fits_id + '.fits'
				hdu.header[FITS_IN_KEY(i)] = fn # add flat fn to master flat header
				print fn
				hdulist = pf.open( fn )
				if mtype is af.FLAT_NAME: # normalize flat frames
					data_arr.append(hdulist[0].data/np.median(hdulist[0].data))
				else:
					data_arr.append(hdulist[0].data)
				filter_arr.append(hdulist[0].header['FILTER'])
				i += 1
			filt0 = None
			for filt in filter_arr:
				if filt0 is None:
					filt0 = filt
				elif filt != filt0:
					print "Error: cannot combine frames with different filters.  Exiting..."
					os.chdir( start_dir ) # move back to starting directory
					return
			hdu.header['FILTER'] = filt0 # add filter keyword to master frame
			hdu.header['CAMERA'] = band  # add camera keyword to master frame

			data_arr = np.array( data_arr )
			
			master = af.imcombine( data_arr, type='median' )
			dtemp = master.astype(np.float)
			if mtype is af.FLAT_NAME:
				hdu.data = dtemp/np.median(dtemp) # normalize master flat
			else:
				hdu.data = dtemp
			hdulist = pf.HDUList( [hdu] )

			hdulist.writeto( '{}_{}.fits'.format( mtype, band ), clobber=True )
		else:
			print "Skipping {}...".format( band )

	# move back to starting directory
	os.chdir( start_dir )