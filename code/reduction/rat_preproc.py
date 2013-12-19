"""
	Purpose:	this is a collection of preprocessing functions for use with data from RATIR.

	Usage:
		1)	enter python or ipython environment
		2)	can load all functions using:
			- "from rat_preproc import *" if you want to call functions using just the function's name
			- "import rat_preproc as rp" if you want to call functions using rp.function_name( args )

	Notes:
		-

	Future Improvements:
		-
"""

import os
from fnmatch import fnmatch
import astropy.io.fits as pf
import numpy as np
import matplotlib.pylab as pl
from scipy.ndimage.interpolation import zoom
import astro_functs as af # contains basic functions and RATIR constants
import shutil
from glob import glob

# Preprocessing constants
ZOOM_LVL = 0.5 # image zoom for display to user in ratdisp() and ratdisp_calib()
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

	Future Improvements:
		- 
"""
def ratlist( workdir='.', cams=[0,1,2,3] ):

	# move to working directory
	start_dir = os.getcwd()
	os.chdir( workdir ) # move into working dir
	
	d = os.getcwd().split('/')[-1] # name of current directory
	print "* * * creating lists of frames in {} * * *".format(d)
	fns = []
	for i in cams:
		fns.append( open( 'C{}_{}.list'.format( i, d ), 'w' ) ) # create list file for each camera
	for f in os.listdir( '.' ):
		for i in range(len(cams)):
                        if fnmatch( f, '*C{}*.fits'.format(cams[i]) ):
				fns[i].write( f + '\n' ) # add detected file to correct list file
	os.chdir( start_dir ) # move back to starting dir

"""
	Translated from plotratir_flat.pro by John Capone (jicapone@astro.umd.edu).

	Purpose:		display calibration images for verification by user

	Input:
		ftype:		type of frames. defaults for af.FLAT_NAME
		workdir:	directory where function is to be executed
		cams:		camera numbers.  all by default

	Usage:
		1)	enter python or ipython environment
		2)	load function -> 'from rat_preproc import ratdisp_calib'
		3)	run function -> 'ratdisp_calib( ftype = bias or flat name, workdir = 'path/to/data/', cams = [#,#,...] )'
			- since default workdir is '.', this argument can be ignored if you are in the same directory as the data
		4)	select which frames you would like to use to create master calibration frames

	Notes:
		- frame rotation needs to be fixed
			- 270 w/o transpose for ZY and 0 w/ JH?
		- added option to specify working directory
		- added error handling for if a camera list is missing
		- camera argument can now be int
                - ccd lists are now by filter rather than camera number

	Future Improvements:
		- option of automated frame selection
			- save order for future automation
		- put preprocessing functions in single script
"""
def ratdisp_calib( ftype=af.FLAT_NAME, workdir='.', cams=[0,1,2,3] ):

	# check for non-list camera argument
	if type(cams) is not list:
		cams = [cams] # convert to list
	
	pl.ion() # pylab in interactive mode

	# move to working directory
	start_dir = os.getcwd()
	os.chdir( workdir )
	
	d = os.getcwd().split('/')[-1] # name of current directory
	print "* * * displaying {} frames in {} for selection * * *".format(ftype, d)

        # delete existing calibration lists
        fntemp = glob( ftype+'*.list' )
        for fn in fntemp: os.remove( fn )

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

				# look at FITs data sequentially
				for line in fin:
					
					fits_fn = line.rstrip() # current fits file name with return removed
					fits_id = fits_fn.split('.')[0] # fits file name with extention removed
                                        print '\n{}'.format( fits_fn )
					hdulist = pf.open( fits_fn )
					im = hdulist[0].data
					h = hdulist[0].header
					
					# no rotation needed (why is C0 rotated by 270 deg in IDL code?)
					im = np.rot90( im, af.CAM_ROTAT[cam_i] )
					
					# Let user determine if the image is good or not
					implot = zoom( im, ZOOM_LVL ) # change image size for display
					m = np.median( im )
					s = af.robust_sigma( im )
                                        print '\t* Median is {} counts.'.format( m )
					
                                        # print filter name
                                        print '\t* Filter used: {}'.format( h['FILTER'] )

					# display image and prompt user
					axim = pl.imshow( implot, vmin=m-5*s, vmax=m+5*s, origin='lower', cmap=pl.cm.gray )
					pl.title( r"Median = {}, $\sigma$ = {:.1f}".format( int(m), s ) )
					
					# query user until valid response is provided
					valid_entry = False
					while not valid_entry:
						
                                                direction = raw_input("\nType Y for YES, N for NO, Q for QUIT: ")
						
						if direction.lower() == 'y':
                                                        fout = open( '{}_{}.list'.format( ftype, h['FILTER'] ), 'a' ) # append if this filter's list exists
							fout.write( fits_id + '\n' ) # write new file name to list
                                                        fout.close()
							valid_entry = True
						
						elif direction.lower() == 'q': # exit function
							print "Exiting..."
							os.chdir( start_dir ) # move back to starting directory
							return
						
						elif direction.lower() != 'n': # invalid case
							print "'{}' is not a valid entry.".format( direction )
						
						else: # 'N' selected, skip
							valid_entry = True

                                        hdulist.close() # close FITs file

				# close files
				fin.close()
			
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
                                        print '\n{}'.format( fits_fn )
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
                                        print '\t* Median of left side is {} counts'.format( mleft )
                                        print '\t* Median of right side is {} counts'.format( mright )
					
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
						
                                                direction = raw_input("\nType Y for YES, N for NO, Q for QUIT: ")
						
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

                                        hdulist.close() # close FITs file

				# close files
				fin.close()
				for f in fout: f.close()

	# move back to starting directory
	os.chdir( start_dir )

"""
	Translated from plotratir.pro by John Capone (jicapone@astro.umd.edu).

	Purpose:	display RATIR images for verification by user

	Input:
		workdir:	directory where function is to be executed
		targetdir:	directory where selected frames and lists are output
		cams:		camera numbers.  all by default

	Usage:
		1)	enter python or ipython environment
		2)	load function -> 'from rat_preproc import ratdisp'
		3)	run function -> 'ratdisp( workdir = 'path/to/data/', targetdir = 'path/to/new data/', cams = [#,#,...] )'
			- since default workdir is '.', this argument can be ignored if you are in the same directory as the data
		4)	select which frames you would like to use

	Notes:
		- frame rotation needs to be fixed
		- added option to specify working directory
		- added error handling for if a camera list is missing
		- camera argument can now be int
		- added target directory option.  sky and object frames are written to the same dir.
		- prompts user for overwrite
                - ccd lists are now by filter rather than camera number

	Future Improvements:
		- automation of frame selection
			- use previous E/W pattern provided by user?
			- could provide this in a calibration file
"""

def ratdisp( workdir='.', targetdir='.', cams=[0,1,2,3] ):

	# check for non-list camera argument
	if type(cams) is not list:
		cams = [cams] # convert to list
	
	pl.ion() # pylab in interactive mode

	# move to working directory
	start_dir = os.getcwd()
	os.chdir( workdir )
	
	d = os.getcwd().split('/')[-1] # name of current directory
	print "* * * displaying science frames in {} for selection * * *".format(d)

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
        fntemp = glob( '{}/{}*.list'.format( targetdir, af.SKY_NAME ) )
        for fn in fntemp: os.remove( fn )

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
				
				# look at FITs data sequentially
				for line in fin:
					fits_fn = line.rstrip() # current fits file name with return removed
					fits_id = fits_fn.split('.')[0] # fits file name with extention removed
                                        print '\n{}'.format( fits_fn )
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

                                        # print filter name
                                        print '\t* Filter used: {}'.format( h['FILTER'] )

					# query user until valid response is provided
					valid_entry = False
					while not valid_entry:
                                                direction = raw_input("\nType Y for YES, N for NO, Q for QUIT: ")
						if direction.lower() == 'y':
                                                        h['WAVELENG'] = 'OPT'
							
							# object frame
							imfits = '{}/{}_{}_{}.fits'.format( targetdir, fits_id, af.OBJ_NAME, cam_i )
							h['TARGNAME'] = targname
							hdulist.writeto( imfits, clobber=True ) # save object frame
                                                        fout = open( '{}/{}_{}.list'.format( targetdir, af.OBJ_NAME, h['FILTER'] ), 'a' ) # append if this filter's list exists
                                                        fout.write( fits_id + '\n' ) # write new file name to list
                                                        fout.close()

							# sky frame
							skyfits = '{}/{}_{}_{}.fits'.format( targetdir, fits_id, af.SKY_NAME, cam_i )
							h['TARGNAME'] = targname_sky
							hdulist.writeto( skyfits, clobber=True ) # save sky frame
                                                        fout = open( '{}/{}_{}.list'.format( targetdir, af.SKY_NAME, h['FILTER'] ), 'a' ) # append if this filter's list exists
                                                        fout.write( fits_id + '\n' ) # write new file name to list
                                                        fout.close()
							
							valid_entry = True
						
						elif direction.lower() == 'q': # exit function
							print "Exiting..."
							os.chdir( start_dir ) # move back to starting directory
							return
						
						elif direction.lower() != 'n': # invalid case
							print "'{}' is not a valid entry.".format( direction )
						
						else: # 'N' selected, skip
							valid_entry = True

                                        hdulist.close() # close FITs file

				# close files
				fin.close()
		
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
                                        print '\n{}'.format( fits_fn )
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
                                                direction = raw_input("\nType E for EAST, W for WEST, N for NEXT, Q for QUIT: ")
						if direction.lower() == 'e' or direction.lower() == 'w': # selected
                                                        h['WAVELENG'] = 'IR'
							
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

                                        hdulist.close() # close FITs file

				# close files
				fin.close()
				for f in fout_img: f.close()
				for f in fout_sky: f.close()
	
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
                bands:		photometric bands
		workdir:	directory where function is to be executed

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

	Future Improvements:
		- may want to be generalized for use with dark current frames
		- need to add filter name to master's header
"""
def mkmaster( mtype, bands, workdir='.' ):

	# move to working directory
	start_dir = os.getcwd()
	os.chdir( workdir )

	# check for valid mtype
	if mtype not in [af.FLAT_NAME, af.BIAS_NAME]:
		print "Error: valid arguments for mtype are", af.FLAT_NAME, "and", af.BIAS_NAME, ". Exiting..."
		os.chdir( start_dir ) # move back to starting directory
		return

	if mtype is af.BIAS_NAME:
		_prep_bias_list( cams=bands, workdir=workdir )

	d = os.getcwd().split('/')[-1] # name of current directory
	print "* * * making master {} frame in {} * * *".format( mtype, d )

        # work on FITs files for specified photometric bands
	for band in bands:
		flist = open( '{}_{}.list'.format( mtype, band ), 'r' )
		hdu = pf.PrimaryHDU()
		data_arr = []
		filter_arr = [] # to check that all frames used the same filter
		i = 0
		for f in flist.readlines():
			fits_fn = f.rstrip() # current fits file name with return removed
			fits_id = fits_fn.split('.')[0] # fits file name with extention removed
			fn = fits_id + '.fits'
			hdu.header[FITS_IN_KEY(i)] = fn # add flat fn to master flat header
			print fn
			hdulist = pf.open( fn )
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
		data_arr = np.array( data_arr )
		master = af.imcombine( data_arr, type='median' )
		dtemp = master.astype(np.float)
		if mtype is af.BIAS_NAME:
			hdu.data = dtemp
		else:
			hdu.data = dtemp/np.median(dtemp) # scale is median
		hdulist = pf.HDUList( [hdu] )

		hdulist.writeto( '{}_{}.fits'.format( mtype, band ), clobber=True )

	# move back to starting directory
	os.chdir( start_dir )