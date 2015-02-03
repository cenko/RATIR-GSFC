"""
    Purpose:    this is a collection of preprocessing functions for use with data from RATIR.

    Usage:
        1)  enter python or ipython environment
        2)  can load all functions using:
            - "from rat_preproc import *" if you want to call functions using just the function's name
            - "import rat_preproc as rp" if you want to call functions using rp.function_name(args)

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
import matplotlib.gridspec as gridspec
from scipy.ndimage.interpolation import zoom
import shutil
from glob import glob
import gc
import datetime
import pickle

import astro_functs as af # contains basic functions and RATIR constants
from astro_functs import show_list # allow user to call show_list without "af." prefix

# Preprocessing constants
FITS_IN_KEY = lambda n: 'IMCMB{:03}'.format(int(n)) # function to make FITS keywords to store file names of combined frames

"""
    Translated from plotratir_flat.pro by John Capone (jicapone@astro.umd.edu).
    Modified by Vicki Toy 8/15/14

    Purpose:        display calibration images for verification by user

    Input:
        ftype:      type of frames. defaults for af.FLAT_NAME
        workdir:    directory where function is to be executed
        cams:       camera numbers.  all by default
        auto:       automated selection of frames.  if ftype is af.BIAS_NAME, select all.  if ftype is af.FLAT_NAME, select non-saturated frames with sufficient counts.
        amin:       minimum median value for automated selection as fraction of saturation value
        amax:       maximum median value for automated selection as fraction of saturation value
        save_select:save dictionary of selected frames to python pickle file

    Usage:
        1)  enter python or ipython environment
        2)  load function -> 'from rat_preproc import choose_calib'
        3)  run function -> 'file_dict = choose_calib(ftype = bias or flat name, workdir = 'path/to/data/', cams = [#,#,...])'
            - since default workdir is '.', this argument can be ignored if you are in the same directory as the data
        4)  select which frames you would like to use to create master calibration frames
        5)  call mkmaster using the dictionary returned by this function

    Notes:
        - added option to specify working directory
        - added error handling for if a camera list is missing
        - camera argument can now be int
            - ccd lists are now by filter rather than camera number
        - automated bias and flat frame selection can be done relative to the detectors' staturation points
            * this is dangerous when selecting flats!!!
        - removed the need for list files.  uses Python dictionaries instead.
            - the returned dictionary MUST be stored to create master frames using mkmaster()

    Future Improvements:
        - better way to do cameras with multiple filters (current way uses camera # -2, won't work for other instruments)
"""
def choose_calib(ftype, workdir='.', cams=[0,1,2,3], auto=False, amin=0.3, amax=0.7, save_select=True):

    if auto and (ftype is af.FLAT_NAME):
        temp = raw_input(af.bcolors.WARNING+"Warning: automated selection of flats is not recommended! Continue? (y/n): "+af.bcolors.ENDC)
        if (temp.lower() != 'y') or (temp.lower() != 'yes'):
            af.print_bold("Exiting...")
            return

    # check for non-list camera argument
    if type(cams) is not list:
        cams = [cams] # convert to list

    # check for non-integer camera designators
    if type(cams[0]) is not int:
        af.print_err("Error: cameras must be specified by an integer. Exiting...")
        return
    
    pl.ion() # pylab in interactive mode

    # move to working directory
    start_dir = os.getcwd()
    os.chdir(workdir)
    
    d = os.getcwd().split('/')[-1] # name of current directory
    if not auto:
        af.print_head("\nDisplaying {} frames in {} for selection:".format(ftype, d))
    else:
        af.print_head("\nAutomatically selecting {} frames in {}:".format(ftype, d))

    # dictionary to store selected fits files by camera or filter
    fits_list_dict = {}

    # open figure for images if not auto
    if not auto:
        fig = pl.figure(figsize=(8,5))

    # work on FITs files for specified cameras
    for cam_i in cams:

        # print current camera number
        af.print_under("\n{:^50}".format('CAMERA {}'.format(cam_i)))
        if af.CAM_BIAS[cam_i] == False and ftype is af.BIAS_NAME:
            af.print_warn("Warning: Camera C{} does not have {} frames.  Skipping...".format(cam_i, af.BIAS_NAME))
            continue
        
        # find raw files of selected type for this camera
        fits_list = glob('????????T??????C{}{}.fits'.format(cam_i, 'b' if ftype is af.BIAS_NAME else 'f'))
        if len(fits_list) == 0:
            af.print_warn("Warning: no fits files found.  Skipping camera {}.".format(cam_i))
            continue

        # look at FITs data sequentially
        for fits_fn in fits_list:
            
            fits_id = fits_fn.split('.')[0] # fits file name with extention removed
            print '{}'.format(fits_fn)

            # open data
            hdulist = pf.open(fits_fn)
            im = hdulist[0].data
            h  = hdulist[0].header

            # get detector's saturation level
            sat_pt = af.CAM_SATUR[cam_i](h['SOFTGAIN'])
            
            if af.CAM_SPLIT[cam_i]:
                im1 = im[af.SLICES[af.SPLIT_FILTERS[cam_i-2]]]
                m1  = np.median(im1)
                s1  = af.robust_sigma(im1)
                im2 = im[af.SLICES[af.SPLIT_FILTERS[cam_i]]]
                m2  = np.median(im2)
                s2  = af.robust_sigma(im2)
                print '\t* Median of left side is {} counts ({:.0%} of saturation level).'.format(m1, float(m1)/sat_pt)
                print '\t* Median of right side is {} counts ({:.0%} of saturation level).'.format(m2, float(m2)/sat_pt)
            else:
                if (h['FILTER'] not in af.RAT_FILTERS) and (ftype is af.FLAT_NAME):
                    af.print_warn("Warning: invalid filter detected.  Skipping {} band.".format(h['FILTER']))
                    continue
                im1 = im[af.SLICES['C'+str(cam_i)]]
                m  = np.median(im1)
                s  = af.robust_sigma(im1)
                print '\t* Median is {} counts ({:.0%} of saturation level).'.format(m, float(m)/sat_pt)
                if ftype is af.FLAT_NAME:
                    print '\t* Filter used: {}'.format(h['FILTER'])

            if auto:
                if ftype is af.BIAS_NAME:   
                    # all bias frames are selected
                    if fits_list_dict.has_key('C{}'.format(cam_i)):
                        fits_list_dict['C{}'.format(cam_i)].append(fits_fn)
                    else:
                        fits_list_dict['C{}'.format(cam_i)] = [fits_fn]
                                                
                if ftype is af.FLAT_NAME:

                    vmin = amin * sat_pt; vmax = amax * sat_pt
                                
                    if af.CAM_SPLIT[cam_i]:
                        # check whether median values are in specified range
                        # bottom side
                        if m1 > vmin and m1 < vmax:
                            af.print_blue("\t* Bottom side selected.")
                            imfits_1 = '{}_{}.fits'.format(fits_id, af.SPLIT_FILTERS[cam_i-2])
                            h['FILTER'] = af.SPLIT_FILTERS[cam_i-2]
                            if os.path.exists(imfits_1): os.remove(imfits_1) # delete old copy
                            pf.writeto(imfits_1, im1, header=h, clobber=True) # save left frame
                            if fits_list_dict.has_key(af.SPLIT_FILTERS[cam_i-2]):
                                fits_list_dict[af.SPLIT_FILTERS[cam_i-2]].append(imfits_1)
                            else:
                                fits_list_dict[af.SPLIT_FILTERS[cam_i-2]] = [imfits_1]
                        else:
                            if m1 < vmin:
                                af.print_warn("\t* Bottom side rejected:\tUNDEREXPOSED.")
                            else:
                                af.print_warn("\t* Bottom side rejected:\tSATURATED.")

                        # top side
                        if m2 > vmin and m2 < vmax:
                            af.print_blue("\t* Top side selected.")
                            imfits_2 = '{}_{}.fits'.format(fits_id, af.SPLIT_FILTERS[cam_i])
                            h['FILTER'] = af.SPLIT_FILTERS[cam_i]
                            if os.path.exists(imfits_2): os.remove(imfits_2) # delete old copy
                            pf.writeto(imfits_2, im2, header=h, clobber=True) # save object frame
                            if fits_list_dict.has_key(af.SPLIT_FILTERS[cam_i]):
                                fits_list_dict[af.SPLIT_FILTERS[cam_i]].append(imfits_2)
                            else:
                                fits_list_dict[af.SPLIT_FILTERS[cam_i]] = [imfits_2]
                        else:
                            if m2 < vmin:
                                af.print_warn("\t* Top side rejected:\tUNDEREXPOSED.")
                            else:
                                af.print_warn("\t* Top side rejected:\tSATURATED.")
                    
                    #Not split frame            
                    else:
                        # check whether median value is in specified range
                        if m > vmin and m < vmax:
                            af.print_blue("\t* Frame selected.")
                            if fits_list_dict.has_key(h['FILTER']):
                                fits_list_dict[h['FILTER']].append(fits_fn)
                            else:
                                fits_list_dict[h['FILTER']] = [fits_fn]
                        else:
                            if m < vmin:
                                af.print_warn("\t* Frame rejected:\tUNDEREXPOSED.")
                            else:
                                af.print_warn("\t* Frame rejected:\tSATURATED.")

            # display image and prompt user
            else:
                if af.CAM_SPLIT[cam_i]:
                    # show top frame
                    ax1 = fig.add_subplot(221)
                    z1, z2 = af.zscale(im1)
                    ax1.imshow(im1, vmin=z1, vmax=z2, origin='lower', cmap=pl.cm.gray)
                    ax1.set_xticks([])
                    ax1.set_yticks([])
                    ax1.set_title(r"Median = {}, $\sigma$ = {:.1f}".format(int(m1), s1))
                    ax1.set_xlabel(r"Median is {:.0%} of saturation level.".format(float(m1)/sat_pt))
                    # show pixel distribution
                    axhist = fig.add_subplot(222)
                    axhist.hist(im1.flatten(), bins=50, normed=True, log=True, range=(0, sat_pt))
                    axhist.set_xlim((0, sat_pt))
                    axhist.set_xticks([0, 0.5*sat_pt, sat_pt])
                    axhist.set_xticklabels(['0%', '50%', '100%'])
                    axhist.set_yticks([])
                    axhist.grid()
                    axhist.set_title("Pixel distribution")
                    axhist.set_ylabel("Log Scale")
                    # show bottom frame
                    ax2 = fig.add_subplot(223)
                    z1, z2 = af.zscale(im2)
                    ax2.imshow(im2, vmin=z1, vmax=z2, origin='lower', cmap=pl.cm.gray)
                    ax2.set_xticks([])
                    ax2.set_yticks([])
                    ax2.set_title(r"Median = {}, $\sigma$ = {:.1f}".format(int(m1), s1))     
                    ax2.set_xlabel(r"Median is {:.0%} of saturation level.".format(float(m2)/sat_pt))
                    # show pixel distribution
                    axhist = fig.add_subplot(224)
                    axhist.hist(im2.flatten(), bins=50, normed=True, log=True, range=(0, sat_pt))
                    axhist.set_xlim((0, sat_pt))
                    axhist.set_xticks([0, 0.5*sat_pt, sat_pt])
                    axhist.set_xticklabels(['0%', '50%', '100%'])
                    axhist.set_yticks([])
                    axhist.grid()
                    axhist.set_title("Pixel distribution")
                    axhist.set_ylabel("Log Scale")
                    fig.subplots_adjust(wspace=0.1, hspace=0.45)
                else:
                    # show frame
                    z1, z2 = af.zscale(im1)
                    ax = fig.add_subplot(121)
                    ax.imshow(im1, vmin=z1, vmax=z2, origin='lower', cmap=pl.cm.gray)
                    ax.set_xticks([])
                    ax.set_yticks([])                        
                    ax.set_title(r"Median = {}, $\sigma$ = {:.1f}".format(int(m), s))
                    ax.set_xlabel(r"Median is {:.0%} of saturation level.".format(float(m)/sat_pt))
                    # show pixel distribution
                    axhist = fig.add_subplot(122)
                    axhist.hist(im1.flatten(), bins=50, normed=True, log=True, range=(0, sat_pt))
                    axhist.set_xlim((0, sat_pt))
                    axhist.set_xticks([0, 0.5*sat_pt, sat_pt])
                    axhist.set_xticklabels(['0%', '50%', '100%'])
                    axhist.set_yticks([])
                    axhist.grid()
                    axhist.set_title("Pixel distribution")
                    axhist.set_ylabel("Log Scale")
                fig.canvas.draw()
                        
                # query user until valid response is provided
                valid_entry = False
                while not valid_entry:
                    user = raw_input("\nType Y for YES, N for NO, Q for QUIT: ")
                            
                    if user.lower() == 'y':
                                
                        if af.CAM_SPLIT[cam_i]:
                            imfits_1 = '{}_{}.fits'.format(fits_id, af.SPLIT_FILTERS[cam_i-2])
                            h['FILTER'] = af.SPLIT_FILTERS[cam_i-2]
                            if os.path.exists(imfits_1): os.remove(imfits_1) # delete old copy
                            pf.writeto(imfits_1, im1, header=h, clobber=True) # save object frame
                            if fits_list_dict.has_key(af.SPLIT_FILTERS[cam_i-2]):
                                fits_list_dict[af.SPLIT_FILTERS[cam_i-2]].append(imfits_1)
                            else:
                                fits_list_dict[af.SPLIT_FILTERS[cam_i-2]] = [imfits_1]
                                
                            imfits_2 = '{}_{}.fits'.format(fits_id, af.SPLIT_FILTERS[cam_i])
                            h['FILTER'] = af.SPLIT_FILTERS[cam_i]
                            if os.path.exists(imfits_2): os.remove(imfits_2) # delete old copy
                            pf.writeto(imfits_2, im2, header=h, clobber=True) # save object frame
                            if fits_list_dict.has_key(af.SPLIT_FILTERS[cam_i]):
                                fits_list_dict[af.SPLIT_FILTERS[cam_i]].append(imfits_2)
                            else:
                                fits_list_dict[af.SPLIT_FILTERS[cam_i]] = [imfits_2]
                        else:
                            if ftype is af.FLAT_NAME:
                                fl_key = h['FILTER']
                            else:
                                fl_key = 'C{}'.format(cam_i)
                            if fits_list_dict.has_key(fl_key):
                                fits_list_dict[fl_key].append(fits_fn)
                            else:
                                fits_list_dict[fl_key] = [fits_fn]
                                
                        valid_entry = True
                            
                    elif user.lower() == 'q': # exit function
                        af.print_bold("Exiting...")
                        os.chdir(start_dir) # move back to starting directory
                        pl.close('all') # close image to free memory
                        return
                            
                    elif user.lower() == 'n': # 'N' selected, skip
                        valid_entry = True
                            
                    else: # invalid case
                        af.print_warn("'{}' is not a valid entry.".format(user))

            if not auto:
                fig.clear() # clear image
            hdulist.close() # close FITs file
        
    if not auto:
        pl.close('all') # close image to free memory

    if auto:
        af.print_head("\nDisplaying automatically selected {} frames:".format(ftype))
        af.show_list(fits_list_dict)

    if save_select:
        dt = datetime.datetime.now()
        fnout = '{}_'.format(ftype)+dt.isoformat().split('.')[0].replace('-','').replace(':','')+'.p' # python pickle extension
        af.print_head("\nSaving selection dictionary to {}".format(fnout))
        pickle.dump( fits_list_dict, open( fnout, 'wb' ) ) # save dictionary to pickle

    os.chdir(start_dir) # move back to starting directory

    return fits_list_dict

"""
    Translated from plotratir.pro by John Capone (jicapone@astro.umd.edu).
    Modified by Vicki Toy 8/14/14

    Purpose:    display RATIR images for verification by user

    Input:
        workdir:    directory where function is to be executed
        targetdir:  directory where selected frames and lists are output
        cams:       camera numbers.  all by default
        auto:       select all science frames
        save_select:save dictionary of selected frames to python pickle file

    Usage:
        1)  enter python or ipython environment
        2)  load function -> 'from rat_preproc import choose_science'
        3)  run function -> 'choose_science(workdir = 'path/to/data/', targetdir = 'path/to/new data/', cams = [#,#,...])'
            - since default workdir is '.', this argument can be ignored if you are in the same directory as the data
        4)  select which frames you would like to use

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
def choose_science(workdir='.', targetdir='.', cams=[0,1,2,3], auto=False, save_select=True):

    # check for non-list camera argument
    if type(cams) is not list:
        cams = [cams] # convert to list

    # check for non-integer camera designators
    if type(cams[0]) is not int:
        af.print_err("Error: cameras must be specified by an integer. Exiting...")
        return
    
    pl.ion() # pylab in interactive mode

    # move to working directory
    start_dir = os.getcwd()
    os.chdir(workdir)
    
    d = os.getcwd().split('/')[-1] # name of current directory
    if not auto:
        af.print_head("\nDisplaying science frames in {} for selection:".format(d))
    else:
        af.print_head("\nSelecting all science frames in {}:".format(d))

    # dictionary to store selected fits files by camera or filter
    fits_list_dict = {}
    
    # remove tailing / from target directory name if present
    if targetdir[-1] == '/':
        targetdir = targetdir[:-1]

    # make target directory if it does not exist
    if not os.path.exists(targetdir):
        af.print_blue("Creating target directory: {}".format(targetdir))
        os.makedirs(targetdir)
    # warn user if previous files may be overwritten
    else:
        af.print_warn("Warning: Target directory exists. Existing files will be overwritten.")
        resp = raw_input("Proceed? (y/n): ")
        if resp.lower() != 'y':
            af.print_bold("Exiting...")
            os.chdir(start_dir) # move back to starting directory
            return
        else:
            shutil.rmtree(targetdir)
            os.makedirs(targetdir)

    # open figure for images if not auto
    if not auto:
        fig = pl.figure(figsize=(8,5))

    # work on FITs files for specified cameras
    for cam_i in cams:
        
        # print current camera number
        af.print_under("\n{:^50}".format('CAMERA {}'.format(cam_i)))
        
        # find raw files of selected type for this camera
        fits_list = glob('????????T??????C{}o.fits'.format(cam_i))
        if len(fits_list) == 0:
            af.print_warn("Warning: no fits files found.  Skipping camera {}.".format(cam_i))
            continue

        # look at FITs data sequentially
        for fits_fn in fits_list:

            fits_id = fits_fn.split('.')[0] # fits file name with extention removed
            print '{}'.format(fits_fn)
                
            # open data
            hdulist = pf.open(fits_fn)
            im = hdulist[0].data
            h = hdulist[0].header

            # check for required header keywords
            if 'PRPSLID' in h:
                prpslid = h['PRPSLID']
            else:
                af.print_err("ERROR: ratdisp - PRPSLID not found in fits header.")
                os.chdir(start_dir) # move back to starting directory
                pl.close('all') # close image to free memory
                return -1
                
            if 'VSTID' in h:
                vstid = h['VSTID']
            else:
                af.print_err("ERROR: ratdisp - VSTID keyword not found in fits header.")
                os.chdir(start_dir) # move back to starting directory
                pl.close('all') # close image to free memory
                return -1
                
            if af.CENTER_KEY in h:
                center = h[af.CENTER_KEY].split('center')[0]
            else:
                af.print_err("ERROR: ratdisp - {} keyword not found in fits header.".format(af.CENTER_KEY))
                os.chdir(start_dir) # move back to starting directory
                pl.close('all')
                return  
                            
            targname = '{}-vis{}'.format(prpslid, vstid)              
                                
            # get image statistics
            if af.CAM_SPLIT[cam_i]:
                im1 = im[af.SLICES[af.SPLIT_FILTERS[cam_i-2]]]
                m1  = np.median(im1)
                s1  = af.robust_sigma(im1)
                im2 = im[af.SLICES[af.SPLIT_FILTERS[cam_i]]]
                m2  = np.median(im2)
                s2  = af.robust_sigma(im2)
            else:
                im1 = im[af.SLICES['C'+str(cam_i)]]             
                m = np.median(im1)
                s = af.robust_sigma(im1)
                    
            # display image and prompt user
            if not auto:
            
                if af.CAM_SPLIT[cam_i]:

                    ax1 = fig.add_subplot(211)
                    z1, z2 = af.zscale(im1)
                    ax1.imshow(im1, vmin=z1, vmax=z2, origin='lower', cmap=pl.cm.gray)
                    ax1.set_xticks([])
                    ax1.set_yticks([])
                    ax1.set_title(r"{} band".format(af.SPLIT_FILTERS[cam_i-2]) + ': ' + "Median = {}, $\sigma$ = {:.1f}".format(int(m1), s1))

                    ax2 = fig.add_subplot(212)
                    z1, z2 = af.zscale(im2)
                    ax2.imshow(im2, vmin=z1, vmax=z2, origin='lower', cmap=pl.cm.gray)
                    ax2.set_xticks([])
                    ax2.set_yticks([])
                    ax2.set_title(r"{} band".format(af.SPLIT_FILTERS[cam_i]) + ': ' + "Median = {}, $\sigma$ = {:.1f}".format(int(m2), s2))
                
                else:

                    ax = fig.add_subplot(111)
                    z1, z2 = af.zscale(im1)
                    ax.imshow(im1, vmin=z1, vmax=z2, origin='lower', cmap=pl.cm.gray)
                    ax.set_xticks([])
                    ax.set_yticks([])
                    ax.set_title(r"{} band: Median = {}, $\sigma$ = {:.1f}".format(h['FILTER'], int(m), s))
                
                fig.canvas.draw()
            
            if af.CAM_SPLIT[cam_i]:
                if center.count(af.SPLIT_FILTERS[cam_i]) != 0:
                    print "\t* The target is focused on the {} filter.".format(af.SPLIT_FILTERS[cam_i])
                elif center.count(af.SPLIT_FILTERS[cam_i-2]) != 0:
                    print "\t* The target is focused on the {} filter.".format(af.SPLIT_FILTERS[cam_i-2])
                else:
                    af.print_warn("\t* Warning: The target is NOT focused on an H2RG filter. The target is focused on the {} filter.".format(center))
            else:
                # print filter name
                print '\t* Filter used: {}'.format(h['FILTER'])

            # query user until valid response is provided
            valid_entry = False
            while not valid_entry:
                    
                # either select all if auto, or have user select
                if auto:
                    user = 'y'
                else:
                    user = raw_input("\nType Y for YES, N for NO, Q for QUIT: ")
                
                if user.lower() == 'y' and af.CAM_SPLIT[cam_i]: 
                    if center.count(af.SPLIT_FILTERS[cam_i]) != 0:
                        direction = 't'
                    elif center.count(af.SPLIT_FILTERS[cam_i-2]) != 0:
                        direction = 'b'
                    else:
                        af.print_warn("\t* Warning: Skipping frame not centered on H2RG filter.")
                        user = 'n'
                        direction = ''

                    
                if user.lower() == 'y':
                        
                    # set keyword values
                    h['CAMERA']   = cam_i
                    h['TARGNAME'] = targname
                    h['PIXSCALE'] = af.CAM_PXSCALE[cam_i]
                    h['WAVELENG'] = af.CAM_WAVE[cam_i]
                    h['GAIN']     = (af.CAM_GAIN[cam_i](h['SOFTGAIN']), 'in electrons/DN')
                    h['SATURATE'] = (af.CAM_SATUR[cam_i](h['SOFTGAIN']), 'in electrons/DN')
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
                        
                        imfits = '{}/{}_{}_{}.fits'.format(targetdir, fits_id, af.OBJ_NAME, af.SPLIT_FILTERS[f_img])
                        h['FILTER'] = af.SPLIT_FILTERS[f_img]
                        im_img = im[af.SLICES[h['FILTER']]]
                        h['NAXIS1'] = af.SLICES[h['FILTER']][1].stop - af.SLICES[h['FILTER']][1].start
                        h['NAXIS2'] = af.SLICES[h['FILTER']][0].stop - af.SLICES[h['FILTER']][0].start
                        h['CRPIX1'] = af.CAM_X0[cam_i] - af.SLICES[h['FILTER']][1].start
                        h['CRPIX2'] = af.CAM_Y0[cam_i] - af.SLICES[h['FILTER']][0].start
                        pf.writeto(imfits, im_img, header=h, clobber=True) # save object frame
                        if fits_list_dict.has_key(af.SPLIT_FILTERS[f_img]):
                            fits_list_dict[af.SPLIT_FILTERS[f_img]].append(imfits)
                        else:
                            fits_list_dict[af.SPLIT_FILTERS[f_img]] = [imfits]
                        
                        # filter side with sky, now saved as object, but different list to keep track
                        skyfits = '{}/{}_{}_{}.fits'.format(targetdir, fits_id, af.OBJ_NAME, af.SPLIT_FILTERS[f_sky])
                        h['FILTER'] = af.SPLIT_FILTERS[f_sky]
                        im_sky = im[af.SLICES[h['FILTER']]]
                        h['NAXIS1'] = af.SLICES[h['FILTER']][1].stop - af.SLICES[h['FILTER']][1].start #Repeat incase filter sizes different
                        h['NAXIS2'] = af.SLICES[h['FILTER']][0].stop - af.SLICES[h['FILTER']][0].start
                        h['CRPIX1'] = af.CAM_X0[cam_i] - af.SLICES[h['FILTER']][1].start
                        h['CRPIX2'] = af.CAM_Y0[cam_i] - af.SLICES[h['FILTER']][0].start
                        pf.writeto(skyfits, im_sky, header=h, clobber=True) # save sky frame
                        if fits_list_dict.has_key(af.SPLIT_FILTERS[f_sky]):
                            fits_list_dict[af.SPLIT_FILTERS[f_sky]].append(skyfits)
                        else:
                            fits_list_dict[af.SPLIT_FILTERS[f_sky]] = [skyfits]

                        valid_entry = True

                    else:
                        
                        imfits = '{}/{}_{}_{}.fits'.format(targetdir, fits_id, af.OBJ_NAME, cam_i)
                        im_img = im[af.SLICES['C'+str(cam_i)]]
                        h['NAXIS1'] = af.SLICES['C'+str(cam_i)][1].stop - af.SLICES['C'+str(cam_i)][1].start
                        h['NAXIS2'] = af.SLICES['C'+str(cam_i)][0].stop - af.SLICES['C'+str(cam_i)][0].start
                        h['CRPIX1'] = af.CAM_X0[cam_i] - af.SLICES['C'+str(cam_i)][1].start
                        h['CRPIX2'] = af.CAM_Y0[cam_i] - af.SLICES['C'+str(cam_i)][0].start
                        pf.writeto(imfits, im_img, header=h, clobber=True)
                        if fits_list_dict.has_key(h['FILTER']):
                            fits_list_dict[h['FILTER']].append(imfits)
                        else:
                            fits_list_dict[h['FILTER']] = [imfits]

                        valid_entry = True                      
                    
                elif user.lower() == 'q': # exit function
                    af.print_bold("Exiting...")
                    os.chdir(start_dir) # move back to starting directory
                    pl.close('all') # close image to free memory
                    return
                
                elif user.lower() == 'n': # 'N' selected, skip
                    valid_entry = True
                        
                else: # invalid case
                    af.print_warn("'{}' is not a valid entry.".format(user))

            if not auto:
                fig.clear() # clear image
            hdulist.close() # close FITs file

    if not auto:
        pl.close('all') # close image to free memory

    if auto:
        af.print_head("\nDisplaying automatically selected science frames:")
        af.show_list(fits_list_dict)

    if save_select:
        dt = datetime.datetime.now()
        fnout = 'object_'+dt.isoformat().split('.')[0].replace('-','').replace(':','')+'.p' # python pickle extension
        af.print_head("\nSaving selection dictionary to {}".format(fnout))
        pickle.dump( fits_list_dict, open( fnout, 'wb' ) ) # save dictionary to pickle

    os.chdir(start_dir) # move back to starting directory

    return fits_list_dict

"""
    Written by John Capone (jicapone@astro.umd.edu).

    Purpose:        make master bias and master flat frames
                    * currently no outlier rejection other than median combine

    Input:
        fn_dict:    dictionary output by choose_calib() containing organized fits file names
        mtype:      type of master frame. should be either af.FLAT_NAME or af.BIAS_NAME
        fmin:       minimum number of files needed to make a master frame

    Usage:
        1)  follow directions for choose_calib()
        2)  make sure to assign the output from choose_calib() to a variable!
        3)  run function -> 'mkmaster(fn_dict=output from choose_calib(), mtype = bias or flat name)'
        4)  create a master configuration file using the frames previously selected by ratdisp_calib()

    Notes:
        - checks that data being combined used same filter
        - added filter keyword master frame's header
        - added fmin parameter.  allows user to abort if fewer files are found to make a master frame.

    Future Improvements:
        - may want to be generalized for use with dark current frames
        - 
"""
def mkmaster(fn_dict, mtype, fmin=5):

    # check for valid mtype
    if mtype not in [af.FLAT_NAME, af.BIAS_NAME]:
        af.print_err("Error: valid arguments for mtype are", af.FLAT_NAME, "and", af.BIAS_NAME, ". Exiting...")
        return

    bands = fn_dict.keys()

    sorttype = 'BAND'

    if mtype is af.BIAS_NAME:
        sorttype = 'CAMERA' 
        
    d = os.getcwd().split('/')[-1] # name of current directory
    af.print_head("\nMaking master {} frame in {}:".format(mtype, d))

    # work on FITs files for specified photometric bands
    for band in bands:

        # print current band
        af.print_under("\n{:^50}".format('{} {}'.format(band, sorttype)))

        hdu = pf.PrimaryHDU()
        data_arr = []
        filter_arr = [] # to check that all frames used the same filter
        i = 0
        fs = fn_dict[band]
        cont = True
        if len(fs) < fmin:
            if len(fs) == 0:
                af.print_err('No frames available to make master {} for {} {}.'.format(mtype, band, sorttype.lower()))
                cont = False
            else:
                temp = raw_input(af.bcolors.WARNING+"Only {} frames available to make master {} for {} {}.  Continue? (y/n): ".format(len(fs), mtype, band, sorttype.lower())+af.bcolors.ENDC)
            
                if temp.lower() == 'y' or temp.lower() == 'yes':
                    cont = True
                else:
                    cont = False
        if cont:
            for fn in fs:
                hdu.header[FITS_IN_KEY(i)] = fn # add flat fn to master flat header
                print fn
                hdulist = pf.open(fn)
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
                    af.print_err("Error: cannot combine frames with different filters.  Exiting...")
                    return
            hdu.header['FILTER'] = filt0 # add filter keyword to master frame
            if band in af.RAT_FILTERS[:3]:
                hdu.header['CAMERA'] = 0  # add camera keyword to master frame
            elif band is af.RAT_FILTERS[3]:
                hdu.header['CAMERA'] = 1  # add camera keyword to master frame
            elif band in af.RAT_FILTERS[4:6]:
                hdu.header['CAMERA'] = 2  # add camera keyword to master frame
            elif band in af.RAT_FILTERS[6:]:
                hdu.header['CAMERA'] = 3  # add camera keyword to master frame

            data_arr = np.array(data_arr)
            
            master = af.imcombine(data_arr, type='median')
            dtemp = master.astype(np.float)
            if mtype is af.FLAT_NAME:
                hdu.data = dtemp/np.median(dtemp) # normalize master flat
            else:
                hdu.data = dtemp
            hdulist = pf.HDUList([hdu])

            hdulist.writeto('{}_{}.fits'.format(mtype, band), clobber=True)
        else:
            af.print_warn("Skipping {}...".format(band))