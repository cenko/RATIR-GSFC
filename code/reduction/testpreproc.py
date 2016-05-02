"""
    Purpose:    this is a collection of preprocessing functions for use with data from RATIR.

    Usage:
        1)  enter python or ipython environment
        2)  can load all functions using:
            - "from preproc import *" if you want to call functions using just the function's name
"""

import os
import itertools
from fnmatch import fnmatch
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import shutil
from glob import glob
import gc
import datetime
import pickle
import sys


# installed modules
import astropy.io.fits as pf
import matplotlib.pylab as pl
import numpy as np
from scipy.ndimage.interpolation import zoom

# custom modules/functions
from zscale import zscale

#import ratir_functs as ratir # contains basic functions and RATIR constants
#import lmi_functs as lmi     # contains basic functions and LMI constants
#from astro_functs import show_list # allow user to call show_list without "af." prefix
import astro_functs as af

from instrument_class import *

# Preprocessing constants
FITS_IN_KEY = lambda n: 'IMCMB{:03}'.format(int(n)) # function to make FITS keywords to store file names of combined frames


def choose_calib(instrument, ftype, workdir='.', cams=[0,1,2,3], auto=False, reject_sat=True, amin=0.2, amax=0.8, save_select=True, figsize=(8,5)):
    """
    NAME:
        choose_calib
    PURPOSE:
        Either auto-select or display calibration images for user verification
    INPUT:
        ftype       - type of calibration frames (ex. 'flat', 'bias', 'dark')
        workdir     - directory where function executed
        cams        - camera numbers (default is all)
        auto        - automated frame selection. If 'bias', will select all, if 'flat' 
                      will select non-saturated frames with sufficient counts
        reject_sat  - reject frames with saturated pixels
        amin        - minimum fraction of saturation value for median (automated)
        amax        - maximum fraction of saturation value for median (automated)
        save_select - save dictionary of selected frames to python pickle file
    EXAMPLE:
        file_dict = choose_calib(ftype = bias, dark or flat name, workdir = 'path/to/data/', cams = [#,#,...])
        *** call mkmaster using dictionary or pickle file ***
    FUTURE IMPROVEMENTS:
        Better way to do cameras with split filter (currently uses camera # -2 this won't
        work for any other instrument)
    """

    if instrument == 'ratir': 
        instrum = ratir

    if auto and (ftype is instrum.flatname):
        temp = raw_input(af.bcolors.WARNING+"Warning: automated selection of flats is not recommended! Continue? (y/n): "+af.bcolors.ENDC)
        if (temp.lower() != 'y') and (temp.lower() != 'yes'):
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
        fig = pl.figure(figsize=figsize)

    # work on FITs files for specified cameras
    for cam_i in cams:

        # print current camera number
        af.print_under("\n{:^50}".format('CAMERA {}'.format(cam_i)))

        # check for valid calibration request
        if instrum.cam_bias[cam_i] == False and ftype is instrum.biasname:
            af.print_warn("Warning: Camera C{} does not have {} frames.  "+
                "Skipping...".format(cam_i, instrum.biasname))
            continue
        if instrum.cam_dark[cam_i] == False and ftype is instrum.darkname:
            af.print_warn("Warning: Camera C{} does not have {} frames.  "+
                "Skipping...".format(cam_i, instrum.darkname))
            continue
          
        # find raw files of selected type for this camera
        fits_list = glob('????????T??????C{}{}.fits'.format(cam_i,
            instrum.ftype_post[ftype]))
        
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
            if instrument == 'ratir':
                sat_pt = instrum.cam_satur[cam_i](h['SOFTGAIN'])
            else:
                sat_pt = h['SATUR']
                
            if reject_sat:
                if np.any(im == sat_pt):
                    af.print_warn("Warning: saturated pixels in frame.  Skipping frame {}.".format(fits_fn))
                    continue

            # check if instrument camera is split, if it is make sure correct specs being used
            if instrum.split == True:
                if instrum.CAM_SPLIT[cam_i]:
                    print '\t* Split filter used'
                else:
                    if (h['FILTER'] not in instrum.filters) and (ftype is instrum.flatname):
                        af.print_warn("Warning: invalid filter detected.  Skipping {} band.".format(h['FILTER']))
                        continue
            
                    if ftype is instrum.flatname:
                        print '\t* Filter used: {}'.format(h['FILTER'])
            
            # return quick image summary    
            [im1, m, s, sfrac] = image_summary(im, sat_pt, cam_i, instrum, split=instrum.CAM_SPLIT[cam_i])
            
            isplit = False
            if instrum.split == True:
                if instrum.CAM_SPLIT[cam_i] == True: 
                    isplit = True 
                    [m1,m2] = m; [s1,s2] = s; [im1,im2] = im1; [sfrac1, sfrac2] = sfrac  
            
            # if auto select then find values with correct ranges
            if auto:

                # all bias and dark frames are selected
                if ftype in [instrum.biasname, instrum.darkname]:
                    addtodict(dict=fits_list_dict, key='C{}'.format(cam_i), value=fits_fn)

                # flats are selected based on median value
                elif ftype is instrum.flatname:

                    vmin = amin * sat_pt; vmax = amax * sat_pt
                    
                    if isplit == True:
                            
                        # check whether median values are in specified range
                        # bottom side
                        if m1 > vmin and m1 < vmax:
                            print '\t* Filter used: {}'.format(instrum.SLICE_FILTERS['C{}a'.format(cam_i)])
                            af.print_blue("\t* Bottom side selected.")
                            imfits_1 = savefile(fits_id, im1, instrum.SLICE_FILTERS['C{}a'.format(cam_i)], h)
                            addtodict(dict=fits_list_dict, 
                                key=instrum.SLICE_FILTERS['C{}a'.format(cam_i)], value=imfits_1)

                        else:
                            if m1 < vmin:
                                af.print_warn("\t* Bottom side rejected:\tUNDEREXPOSED.")
                            else:
                                af.print_warn("\t* Bottom side rejected:\tSATURATED.")

                        # top side
                        if m2 > vmin and m2 < vmax:
                            print '\t* Filter used: {}'.format(instrum.SLICE_FILTERS['C{}b'.format(cam_i)])
                            af.print_blue("\t* Top side selected.")
                            imfits_2 = savefile(fits_id, im2, instrum.SLICE_FILTERS['C{}b'.format(cam_i)], h)
                            addtodict(dict=fits_list_dict, 
                                key=instrum.SLICE_FILTERS['C{}b'.format(cam_i)], value=imfits_2)
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
                            addtodict(dict=fits_list_dict, 
                                key=h['FILTER'], value=fits_fn)

                        else:
                            if m < vmin:
                                af.print_warn("\t* Frame rejected:\tUNDEREXPOSED.")
                            else:
                                af.print_warn("\t* Frame rejected:\tSATURATED.")

            # display image and prompt user
            else:
                
                if isplit == True:
                
                    if (sfrac1 < amin) or (sfrac1 > amax) or (sfrac2 < amin) or (sfrac2 > amax):
                        af.print_warn("Warning: median value outside specified range of {:.0%} - {:.0%} of saturation value in frame.  Skipping frame {}.".format(amin, amax, fits_fn))
                        continue
                    # show top frame
                    ax1 = fig.add_subplot(221)
                    plot_params_calib(ax1, im1, m1, s1, sat_pt, hist=False)
                    
                    # show pixel distribution
                    axhist = fig.add_subplot(222)
                    plot_params_calib(axhist, im1, m1, s1, sat_pt, hist=True)
                    
                    # show bottom frame
                    ax2 = fig.add_subplot(223)
                    plot_params_calib(ax2, im2, m2, s2, sat_pt, hist=False)
                    
                    # show pixel distribution
                    axhist = fig.add_subplot(224)
                    plot_params_calib(axhist, im2, m2, s2, sat_pt, hist=True)
                    fig.subplots_adjust(wspace=0.1, hspace=0.45)
                    
                else:
                    if (sfrac < amin) or (sfrac > amax):
                        af.print_warn("Warning: median value outside specified range of {:.0%} - {:.0%} of saturation value in frame.  Skipping frame {}.".format(amin, amax, fits_fn))
                        continue

                    # show frame
                    ax = fig.add_subplot(121)
                    plot_params_calib(ax, im1, m, s, sat_pt, hist=False)
                    
                    # show pixel distribution
                    axhist = fig.add_subplot(122)
                    plot_params_calib(axhist, im1, m, s, sat_pt, hist=True)
                    
                fig.canvas.draw()
                        
                # query user until valid response is provided
                valid_entry = False
                while not valid_entry:

                    user = raw_input("\nType Y for YES, N for NO, Q for QUIT: ")
                            
                    if user.lower() == 'y':
                                
                        if isplit == True:
                            imfits_1 = savefile(fits_id, im1,
                                instrum.SLICE_FILTERS['C{}a'.format(cam_i)], h)
                            addtodict(dict=fits_list_dict, 
                                key=instrum.SLICE_FILTERS['C{}a'.format(cam_i)],
                                value=imfits_1)
                                
                            imfits_2 = savefile(fits_id, im2,
                                instrum.SLICE_FILTERS['C{}b'.format(cam_i)], h)
                            addtodict(dict=fits_list_dict, 
                                key=instrum.SLICE_FILTERS['C{}b'.format(cam_i)],
                                value=imfits_2)

                        else:
                            if ftype is instrum.flatname:
                                fl_key = h['FILTER']
                            else:
                                fl_key = 'C{}'.format(cam_i)
                            addtodict(dict=fits_list_dict, key=fl_key, value=fits_fn)                            
                                
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

            if not auto: fig.clear() # clear image
            hdulist.close() # close FITs file

    if auto:
        af.print_head("\nDisplaying automatically selected {} frames:".format(ftype))
        af.show_list(fits_list_dict)
    else:
        pl.close('all') # close image to free memory

    if save_select:
        dt = datetime.datetime.now()
        fnout = '{}_'.format(ftype)+dt.isoformat().split('.')[0].replace('-','').replace(':','')+'.p' # python pickle extension
        af.print_head("\nSaving selection dictionary to {}".format(fnout))
        pickle.dump( fits_list_dict, open( fnout, 'wb' ) ) # save dictionary to pickle

    os.chdir(start_dir) # move back to starting directory

    return fits_list_dict
    
    
def image_summary(im, sat_pt, cam_i, instrum, split=False):
    """
    NAME:
        image_summary
    PURPOSE:
        Calculate median, robust scatter, and fraction of saturation point for slice.  If
        split array then will output information for both sides of array
    INPUTS:
        im     - data 
        sat_pt - saturation point
        cam_i  - camera that is being used
        af     - module that contains instrument specific information about cameras and split
        split  - is camera split with filters?
    """
    
    if split == True:
        im1 = im[instrum.slice['C{}a'.format(cam_i)]]
        m1  = np.median(im1)
        s1  = af.robust_sigma(im1)
        sfrac1 = float(m1)/sat_pt
        im2 = im[instrum.slice['C{}b'.format(cam_i)]]
        m2  = np.median(im2)
        s2  = af.robust_sigma(im2)
        sfrac2 = float(m2)/sat_pt
        print '\t* Median of left side is {} counts ({:.0%} of saturation level).'.format(m1, sfrac1)
        print '\t* Median of right side is {} counts ({:.0%} of saturation level).'.format(m2, sfrac2)
        
        return [[im1,im2], [m1,m2], [s1,s2],[sfrac1,sfrac2]]
        
    else:
        im1 = im[instrum.slice['C{}'.format(cam_i)]]
        m  = np.median(im1)
        s  = af.robust_sigma(im1)
        sfrac = float(m)/sat_pt
        print '\t* Median is {} counts ({:.0%} of saturation level).'.format(m, sfrac)

        return [im1,m,s,sfrac]
        
def addtodict(dict=None, key=None, value=None):
    """
    NAME:
        addtodict
    PURPOSE:
        Adds (key, value) to dictionary by initializing or appending
    INPUTS:
        dict  - dictionary to add to
        key   - key to add to dictionary
        value - value of key to add to dictionary
    EXAMPLE:
        addtodict(dict=dictionary, key='test', value='testing')
    """
    
    try:
        dict[key].append(value)
    except:
        dict[key] = [value]
        
def savefile(file, im, filter, h):
    """
    NAME:
        savefile
    PURPOSE:
        Adds filter keyword and saves file (usually for split filter cameras)
    INPUT:
        file   - filename without extension
        im     - data
        filter - filter name to add
        h      - header
    EXAMPLE:
        savefile('test', data, 'Z', header)
    """
    newfile = '{}_{}.fits'.format(file, filter)
    h['FILTER'] = filter
    if os.path.exists(newfile): os.remove(newfile) # delete old copy
    pf.writeto(newfile, im, header=h, clobber=True) # save object frame
    return newfile
    
def plot_params_calib(ax, im, m, s, sat_pt, hist=False):
    """
    NAME:
        plot_params_calib
    PURPOSE:
        Plots calibration files (image and histogram) for user selection
    INPUTS:
        ax     - plot reference that we will be using to plot images
        im     - data to plot
        m      - median
        s      - robust scatter
        sat_pt - saturation point
    EXAMPLE:
        plot_params_calib(ax, im, m, s, sat_pt, hist=False)
    NOTE:
        Will not plot unless you have a show() or something to display    
    """
    
    z1, z2 = af.zscale(im)
    if z2 <= z1:
        z1 = m1 - s1; z2 = m1 + s1
          
    if hist:
        ax.hist(im.flatten(), bins=50, normed=True, log=True, range=(0, sat_pt))
        ax.set_xlim((0, sat_pt))
        ax.set_xticks([0, 0.5*sat_pt, sat_pt])
        ax.set_xticklabels(['0%', '50%', '100%'])
        ax.set_yticks([])
        ax.grid()
        ax.set_title("Pixel distribution")
        ax.set_ylabel("Log Scale")
    else:
        ax.imshow(im, vmin=z1, vmax=z2, origin='lower', cmap=pl.cm.gray, interpolation='none')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(r"Median = {}, $\sigma$ = {:.1f}".format(int(m), s))
        ax.set_xlabel(r"Median is {:.0%} of saturation level.".format(float(m)/sat_pt))    
        

def plot_params_science(ax, im, filter, h, central=False, split=False, 
                        window_zoom=4, calibrate=False):

    disp_im = np.copy(im)
    if calibrate:
        if split:
            disp_im = np.copy(im)/mflat_data
        else:
            disp_im = (np.copy(im)-mbias_data-mdark_data*h['EXPTIME'])/mflat_data
    z1, z2 = af.zscale(disp_im)
    ax.set_xticks([])
    ax.set_yticks([])
    xm, ym = np.array(disp_im.shape, dtype=float)/2
    xr = xm/float(window_zoom); yr = ym/float(window_zoom)
    
    if central:    
        ax.imshow(disp_im[xm-xr:xm+xr,ym-yr:ym+yr], vmin=z1, vmax=z2, origin='lower', 
            cmap=pl.cm.gray, interpolation='none')
        ax.contour(disp_im[xm-xr:xm+xr,ym-yr:ym+yr], levels=[z2], 
            origin='lower', colors='r')
        ax.set_title("Zoomed region")
        
    else:
        ax.imshow(disp_im, vmin=z1, vmax=z2, origin='lower', 
            cmap=pl.cm.gray, interpolation='none')
        ax.contour(disp_im, levels=[z2], origin='lower', colors='r')
        ax.add_patch(Rectangle((ym-yr, xm-xr), 2*yr, 2*xr, ec='b', fc='none', lw=2))
        ax.set_title(r"{} band".format(filter))