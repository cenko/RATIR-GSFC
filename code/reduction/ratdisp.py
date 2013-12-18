"""
Translated from plotratir.pro by John Capone (jicapone@astro.umd.edu).

Notes:
	- frame rotation needs to be fixed
	- added option to specify working directory
	- added error handling for if a camera list is missing
	- camera argument can now be int
	- added target directory option.  sky and object frames are written to the same dir.
	- prompts user for overwrite
"""

import os
import astropy.io.fits as pf
import numpy as np
import matplotlib.pylab as pl
from scipy.ndimage.interpolation import zoom
import astro_functs as af # contains basic functions and RATIR constants

# CONSTANTS
ZOOM_LVL = 0.5

# display RATIR images for verification by user
def ratdisp( workdir='.', targetdir='.', cams=[0,1,2,3] ):
	
	print "* * * displaying science frames for selection * * *"

	# check for non-list camera argument
	if type(cams) is not list:
		cams = [cams] # convert to list
	
	pl.ion() # pylab in interactive mode

	# move to working directory
	start_dir = os.getcwd()
	os.chdir( workdir )
	
	d = os.getcwd().split('/')[-1] # name of current directory

	# make target directory if it does not exist
	if not os.path.exists( targetdir ):
		print "Creating target directory: ", targetdir
		os.makedirs( targetdir )
	else:
		print "Warning: Target directory exists. Existing files will be overwritten."
		resp = raw_input( "Proceed? (y/n): " )
		if resp.lower() != 'y':
			print "Exiting..."
			return
	
	# remove tailing / from target directory name if present
	if targetdir[-1] == '/':
		targetdir = targetdir[:-1]

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
				fout_img = open( '{}/{}_C{}.list'.format( targetdir, af.OBJ_NAME, cam_i ), 'w' ) # create list file for new FITs files
				fout_sky = open( '{}/{}_C{}.list'.format( targetdir, af.SKY_NAME, cam_i ), 'w' ) # create list file for new FITs files
				
				# look at FITs data sequentially
				for line in fin:
					fits_fn = line.rstrip() # current fits file name with return removed
					fits_id = fits_fn.split('.')[0] # fits file name with extention removed
					print fits_fn
					hdulist = pf.open( fits_fn )
					im = hdulist[0].data
					h = hdulist[0].header
					if 'PRPSLID' in h:
						prpslid = h['PRPSLID']
					else:
						"ERROR: ratdisp - PRPSLID not found in fits header."
						os.chdir( start_dir ) # move back to starting directory
						return
					if 'VSTID' in h:
						vstid = h['VSTID']
					else:
						"ERROR: ratdisp - VSTID keywork not found in fits header."
						os.chdir( start_dir ) # move back to starting directory
						return
					targname = '{}-vis{}'.format( prpslid, vstid )
					targname_sky = '{}-sky'.format( targname )
					
					# rotate image or not
					im = np.rot90( im, af.CAM_ROTAT[cam_i] )
					
					# Let user determine if the image is good or not
					implot = zoom( im, ZOOM_LVL ) # change image size for display
					m = np.median( im )
					s = af.robust_sigma( im )
					
					# display image and prompt user
					axim = pl.imshow( implot, vmin=m-5*s, vmax=m+5*s, origin='lower', cmap=pl.cm.gray )
					axim.axes.set_xticks(axim.axes.get_xlim())
					axim.axes.set_xticklabels(['E','W'])
					axim.axes.set_yticks(axim.axes.get_ylim())
					axim.axes.set_yticklabels(['S','N'])
					pl.title( r"Median = {}, $\sigma$ = {:.1f}".format( int(m), s ) )

					# query user until valid response is provided
					valid_entry = False
					while not valid_entry:
						direction = raw_input("Type Y for YES, N for NO, Q for QUIT: ")
						if direction.lower() == 'y':
							h['WAVELEN'] = 'OPT'
							
							# object frame
							imfits = '{}/{}_{}_{}.fits'.format( targetdir, fits_id, af.OBJ_NAME, cam_i )
							h['TARGNAME'] = targname
							hdulist.writeto( imfits, clobber=True ) # save object frame
							fout_img.write( fits_id + '\n' ) # write new file name to list
							
							# sky frame
							skyfits = '{}/{}_{}_{}.fits'.format( targetdir, fits_id, af.SKY_NAME, cam_i )
							h['TARGNAME'] = targname_sky
							hdulist.writeto( skyfits, clobber=True ) # save sky frame
							fout_sky.write( fits_id + '\n' ) # write new file name to list
							
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
				fout_img.close()
				fout_sky.close()
		
		# H2RGs
		if cam_i in [2,3]:
			fn_list = '{}_{}.list'.format( af.CAM_NAMES[cam_i], d )
			try:
				fin = open( fn_list, 'r' ) # open list of FITs files
			except IOError:
				print "Warning: {} not found.  Skipping camera {}.".format( fn_list, cam_i )
			else:
				fout_img = [open( '{}/{}_{}.list'.format( targetdir, af.OBJ_NAME, af.H2RG_FILTERS[cam_i-2] ), 'w' ), open( '{}/{}_{}.list'.format( targetdir, af.OBJ_NAME, af.H2RG_FILTERS[cam_i] ), 'w' )] # create list files for new img FITs files (e, w)
				fout_sky = [open( '{}/{}_{}.list'.format( targetdir, af.SKY_NAME, af.H2RG_FILTERS[cam_i-2] ), 'w' ), open( '{}/{}_{}.list'.format( targetdir, af.SKY_NAME, af.H2RG_FILTERS[cam_i] ), 'w' )] # create list files for new sky FITs files (e, w)
				
				# look at FITs data sequentially
				for line in fin:
					fits_fn = line.rstrip() # current fits file name with return removed
					fits_id = fits_fn.split('.')[0] # fits file name with extention removed
					print fits_fn
					hdulist = pf.open( fits_fn )
					im = hdulist[0].data
					h = hdulist[0].header
					if 'PRPSLID' in h:
						prpslid = h['PRPSLID']
					else:
						"ERROR: ratdisp - PRPSLID not found in fits header."
						os.chdir( start_dir ) # move back to starting directory
						return
					if 'VSTID' in h:
						vstid = h['VSTID']
					else:
						"ERROR: ratdisp - VSTID keywork not found in fits header."
						os.chdir( start_dir ) # move back to starting directory
						return
					targname = '{}-vis{}'.format( prpslid, vstid )
					targname_sky = '{}-sky'.format( targname )
					
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
					
					# display image and prompt user
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

					# query user until valid response is provided
					valid_entry = False
					while not valid_entry:
						direction = raw_input("Type E for EAST, W for WEST, N for NEXT, Q for QUIT: ")
						if direction.lower() == 'e' or direction.lower() == 'w': # selected
							h['WAVELEN'] = 'IR'
							
							if direction.lower() == 'e':
								f_img = cam_i-2
								f_sky = cam_i
							else:
								f_img = cam_i
								f_sky = cam_i-2
							
							# filter side with object
							imfits = '{}/{}_{}_{}.fits'.format( targetdir, fits_id, af.OBJ_NAME, af.H2RG_FILTERS[f_img] )
							im_img = im[af.H2RG_SLICES[f_img]]
							h['NAXIS1'] = af.H2RG_SLICES[f_img][0].stop - af.H2RG_SLICES[f_img][0].start
							h['NAXIS2'] = af.H2RG_SLICES[f_img][1].stop - af.H2RG_SLICES[f_img][1].start
							h['FILTER'] = af.H2RG_FILTERS[f_img]
							h['TARGNAME'] = targname
							pf.writeto( imfits, im_img, header=h, clobber=True ) # save object frame
							fout_img[0].write( '{}_{}_{}\n'.format( fits_id, af.OBJ_NAME, af.H2RG_FILTERS[f_img] ) ) # write new file name to list
							
							# filter side with sky
							skyfits = '{}/{}_{}_{}.fits'.format( targetdir, fits_id, af.SKY_NAME, af.H2RG_FILTERS[f_sky] )
							im_sky = im[af.H2RG_SLICES[f_sky]]
							h['NAXIS1'] = af.H2RG_SLICES[f_sky][0].stop - af.H2RG_SLICES[f_sky][0].start # repeat incase filter sizes differ
							h['NAXIS2'] = af.H2RG_SLICES[f_sky][1].stop - af.H2RG_SLICES[f_sky][1].start
							h['FILTER'] = af.H2RG_FILTERS[f_sky]
							h['TARGNAME'] = targname_sky
							pf.writeto( skyfits, im_sky, header=h, clobber=True ) # save sky frame
							fout_sky[0].write( '{}_{}_{}\n'.format( fits_id, af.SKY_NAME, af.H2RG_FILTERS[f_sky] ) ) # write new file name to list
							
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
				for f in fout_img: f.close()
				for f in fout_sky: f.close()
	
	# move back to starting directory
	os.chdir( start_dir )
