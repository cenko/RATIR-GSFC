"""
Translated from plotratir_flat.pro by John Capone (jicapone@astro.umd.edu).

Notes:
	- frame rotation needs to be fixed
	- added option to specify working directory
	- added error handling for if a camera list is missing
	- camera argument can now be int

Usage:
	1)	enter python or ipython environment
	2)	load function -> 'from ratdisp_calib import ratdisp_calib'
	3)	run function -> 'ratdisp_calib( workdir = 'path/to/data/', cams=[#,#,...] )'
		- since default workdir is '.', this argument can be ignored if you are in the same directory as the data
	4)	select which flats you would like to use for 
"""

import os
import astropy.io.fits as pf
import numpy as np
import matplotlib.pylab as pl
from scipy.ndimage.interpolation import zoom
import astro_functs as af # contains basic functions and RATIR constants

# CONSTANTS
ZOOM_LVL = 0.5

# 270 w/o transpose for ZY and 0 w/ JH

"""
Purpose:		display calibration images for verification by user
Input:
	ftype:		type of frames. defaults for af.FLAT_NAME
	workdir:	directory where function is to be executed
	cams:		camera numbers.  all by default
"""
def ratdisp_calib( ftype=af.FLAT_NAME, workdir='.', cams=[0,1,2,3] ):

	print "* * * displaying {} frames for selection * * *".format(ftype)

	# check for non-list camera argument
	if type(cams) is not list:
		cams = [cams] # convert to list
	
	pl.ion() # pylab in interactive mode

	# move to working directory
	start_dir = os.getcwd()
	os.chdir( workdir )
	
	d = os.getcwd().split('/')[-1] # name of current directory

	# work on FITs files for specified cameras
	for cam_i in cams:
		
		# CCDs
		if cam_i in [0,1]:
			
			fn_list = '{}_{}.list'.format( af.CAM_NAMES[cam_i], d )
			try:
				fin = open( fn_list, 'r' ) # open list of FITs files
			except IOError:
				print "Warning: {} not found.  Skipping camera {}.".format( fn_list, cam_i )
			else:
				fout = open( '{}_{}.list'.format( ftype, cam_i ), 'w' ) # create list file for new FITs files
				
				# look at FITs data sequentially
				for line in fin:
					
					fits_fn = line.rstrip() # current fits file name with return removed
					fits_id = fits_fn.split('.')[0] # fits file name with extention removed
					print fits_fn
					hdulist = pf.open( fits_fn )
					im = hdulist[0].data
					h = hdulist[0].header
					
					# no rotation needed (why is C0 rotated by 270 deg in IDL code?)
					im = np.rot90( im, af.CAM_ROTAT[cam_i] )
					
					# Let user determine if the image is good or not
					implot = zoom( im, ZOOM_LVL ) # change image size for display
					m = np.median( im )
					s = af.robust_sigma( im )
					print 'Median is {} counts.'.format( m )
					
					# display image and prompt user
					axim = pl.imshow( implot, vmin=m-5*s, vmax=m+5*s, origin='lower', cmap=pl.cm.gray )
					pl.title( r"Median = {}, $\sigma$ = {:.1f}".format( int(m), s ) )
					
					# query user until valid response is provided
					valid_entry = False
					while not valid_entry:
						
						direction = raw_input("Type Y for YES, N for NO, Q for QUIT: ")
						
						if direction.lower() == 'y':
							fout.write( fits_id + '\n' ) # write new file name to list
							valid_entry = True
						
						elif direction.lower() == 'q': # exit function
							print "Exiting..."
							os.chdir( start_dir ) # move back to starting directory
							return
						
						elif direction.lower() != 'n': # invalid case
							print "'{}' is not a valid entry.".format( direction )
						
						else: # 'N' selected, skip
							valid_entry = True
				
				# close files
				fin.close()
				fout.close()
			
		# H2RGs
		if cam_i in [2,3]:
			
			fn_list = '{}_{}.list'.format( af.CAM_NAMES[cam_i], d )
			try:
				fin = open( fn_list, 'r' ) # open list of FITs files
			except IOError:
				print "Warning: {} not found.  Skipping camera {}.".format( fn_list, cam_i )
			else:
				fout = [open( '{}_{}.list'.format( ftype, af.H2RG_FILTERS[cam_i-2] ), 'w' ), open( '{}_{}.list'.format( ftype, af.H2RG_FILTERS[cam_i] ), 'w' )] # create list files for new img FITs files (e, w)
				
				# look at FITs data sequentially
				for line in fin:
					
					fits_fn = line.rstrip() # current fits file name with return removed
					fits_id = fits_fn.split('.')[0] # fits file name with extention removed
					print fits_fn
					hdulist = pf.open( fits_fn )
					im = hdulist[0].data
					h = hdulist[0].header
					
					# rotate FITs data to align with filters (Z/J in E-left, Y/H in W-right)
					im = np.rot90( im, af.CAM_ROTAT[cam_i] )
					if cam_i == 3:
						im = np.flipud( im ) # flip C3 about y axis
					
					# Let user determine if the image is good or not
					implot = zoom( im, ZOOM_LVL ) # change image size for display
					imleft = im[af.H2RG_SLICES[cam_i-2]]
					mleft = np.median( imleft )
					sleft = af.robust_sigma( imleft )
					imright = im[af.H2RG_SLICES[cam_i]]
					mright = np.median( imright )
					sright = af.robust_sigma( imright )
					print 'Median of left side is {} counts'.format( mleft )
					print 'Median of right side is {} counts'.format( mright )
					
					# display image and prompt user			***** use subplots to display filter windows *****
					pl.subplot(121)
					pl.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.0, hspace=0.2)
					axim = pl.imshow( imleft, vmin=mleft-5*sleft, vmax=mleft+5*sleft, origin='lower', cmap=pl.cm.gray )
					axim.axes.set_xticks([axim.axes.get_xlim()[0]])
					axim.axes.set_xticklabels(['E'])
					axim.axes.set_yticks(axim.axes.get_ylim())
					axim.axes.set_yticklabels(['S','N'])
					pl.title( r"Median = {}, $\sigma$ = {:.1f}".format( int(mleft), sleft ) )
					pl.subplot(122)
					axim = pl.imshow( imright, vmin=mright-5*sright, vmax=mright+5*sright, origin='lower', cmap=pl.cm.gray )
					axim.axes.set_xticks([axim.axes.get_xlim()[1]])
					axim.axes.set_xticklabels(['W'])
					axim.axes.set_yticks([])
					pl.title( r"Median = {}, $\sigma$ = {:.1f}".format( int(mright), sright ) )
					
					# query user until valid response is provided, Z & J face the EASTERN part of the field
					valid_entry = False
					while not valid_entry:
						
						direction = raw_input("Type Y for YES, N for NO, Q for QUIT: ")
						
						if direction.lower() == 'y':
							
							imfits_left = '{}_{}_{}.fits'.format( fits_id, af.FLAT_NAME, af.H2RG_FILTERS[cam_i-2] )
							im_left = im[af.H2RG_SLICES[cam_i-2]]
							h['FILTER'] = af.H2RG_FILTERS[cam_i-2]
							pf.writeto( imfits_left, im_left, header=h, clobber=True ) # save object frame
							fout[0].write( '{}_{}_{}\n'.format( fits_id, af.FLAT_NAME, af.H2RG_FILTERS[cam_i-2] ) ) # write new file name to list
							
							imfits_right = '{}_{}_{}.fits'.format( fits_id, af.FLAT_NAME, af.H2RG_FILTERS[cam_i] )
							im_right = im[af.H2RG_SLICES[cam_i]]
							h['FILTER'] = af.H2RG_FILTERS[cam_i]
							pf.writeto( imfits_right, im_right, header=h, clobber=True ) # save object frame
							fout[1].write( '{}_{}_{}\n'.format( fits_id, af.FLAT_NAME, af.H2RG_FILTERS[cam_i] ) ) # write new file name to list
							
							valid_entry = True
						
						elif direction.lower() == 'q': # exit function
							print "Exiting..."
							os.chdir( start_dir ) # move back to starting directory
							return
						
						elif direction.lower() != 'n': # invalid case
							print "'{}' is not a valid entry.".format( direction )
						
						else: # 'N' selected, skip
							valid_entry = True
				
				# close files
				fin.close()
				for f in fout: f.close()

	# move back to starting directory
	os.chdir( start_dir )
