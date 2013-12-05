"""
Translated from plotratir.pro by John Capone (jicapone@astro.umd.edu).
"""

import numpy as np
import os
import fnmatch
import astropy.io.fits as pf
from scipy.ndimage.interpolation import zoom
from scipy.misc import bytescale
import pylab as pl
import scipy as sp
from astropy import wcs

import time

def plotratir():
    ''' 
    Save a Matplotlib figure as an image without borders or frames.
       Args:
            fileName (str): String that ends in .png etc.

            fig (Matplotlib figure instance): figure you want to save as the image
        Keyword Args:
            orig_size (tuple): width, height of the original image used to maintain 
            aspect ratio.
    '''
    def SaveFigureAsImage(fileName,fig=None,**kwargs):
        fig_size = fig.get_size_inches()
        w,h = fig_size[0], fig_size[1]
        fig.patch.set_alpha(0)
        if kwargs.has_key('orig_size'): # Aspect ratio scaling if required
            w,h = kwargs['orig_size']
            w2,h2 = fig_size[0],fig_size[1]
            fig.set_size_inches([(w2/w)*w,(w2/w)*h])
            fig.set_dpi((w2/w)*fig.get_dpi())
        a=fig.gca()
        a.set_frame_on(False)
        a.set_xticks([]); a.set_yticks([])
        pl.axis('off')
        pl.xlim(0,h); pl.ylim(w,0)
        fig.savefig(fileName, transparent=True, bbox_inches='tight', \
                            pad_inches=0)

    # general error to raise within plotratir.py
    class plotratirError( Exception ):
        def __init__(self, value):
            self.value = value
        def __str__(self):
            return repr(self.value)

    # find index of xarr and yarr with minimum RSS distance from x and y
    #   Note: needs updated implementation optimized for numpy
    def nearest( x, y, xarr, yarr, mindist ):
        index = -1
        imindist = mindist
        if np.shape(xarr) != np.shape(yarr):
            raise plotratirError( "xarr and yarr must have equal dimensions" )    
        for i in range( np.size(xarr) ):
            dist = np.sqrt( (x - xarr[i])**2. + (y - yarr[i])**2. )
            if dist < imindist:
                imindist = dist
                index = i
        return index

    # make a circle for identifying sources in images
    def circle( xcenter, ycenter, radius ):
        points = np.linspace( 0., 2.*np.pi, 100 )
        x = xcenter + radius * np.cos(points)
        y = ycenter + radius * np.sin(points)
        return np.transpose([x,y])

    # arrays for plotting?
    #   Note: clean up later
    filters = ['r','i','z','y','J','H']
    arr_size = 10000
    plotra = np.zeros(arr_size)
    plotdec = np.zeros(arr_size)
    plotrmag = np.zeros(arr_size)
    plotrmagerr = np.zeros(arr_size)
    plotimag = np.zeros(arr_size)
    plotimagerr = np.zeros(arr_size)
    plotzmag = np.zeros(arr_size)
    plotzmagerr = np.zeros(arr_size)
    plotymag = np.zeros(arr_size)
    plotymagerr = np.zeros(arr_size)
    plotJmag = np.zeros(arr_size)
    plotJmagerr = np.zeros(arr_size)
    plotHmag = np.zeros(arr_size)
    plotHmagerr = np.zeros(arr_size)
    plotcat = np.zeros(arr_size)
    plotfilt = np.zeros(arr_size)

    # returns files in directory "loc" which start with prefix and end with postfix
    def get_files( selection, loc='.' ):
        matches = []
        for files in os.listdir(loc):
            if fnmatch.fnmatch( files, selection ):
                matches.append(files)
        return matches

    # retrieve detection files
    prefchar = 'coadd'
    wildcharimg = '?????-????-?????_?'
    zffiles = get_files( prefchar + wildcharimg + '.crop.multi.fits' )

    time1 = time.clock()

    # find overlapping stars   ***   SLOW   ***
    for i in range(np.size(zffiles)):
        cfilter = zffiles[i].split('_')[1].split('.')[0] # extract filter label from file name
        if cfilter == 'Z' or cfilter == 'Y':
            cfilter = cfilter.lower()
        ifile = 'finalphot' + cfilter + '.am'
        s_id, x, y, ra, dec, mag, magerr = np.loadtxt(ifile, unpack=True)
        if i == 0:
            plotra[0:np.size(ra)] = ra
            plotdec[0:np.size(dec)] = dec
            exec 'plot' + cfilter + 'mag[0:np.size(mag)] = mag'
            exec 'plot' + cfilter + 'magerr[0:np.size(magerr)] = magerr'
            # IDL ;plotcat(0:n_elements(cat)-1)=cat
            nstars = np.size(ra)
        else:
            compra = ra
            compdec = dec
            compmag = mag
            compmagerr = magerr
            # IDL ;compcat=cat
            
            count = 0
            for j in range(np.size(compra)):
                smatch = nearest( compra[j]*np.cos(compdec[j]*np.pi/180.), compdec[j], plotra*np.cos(plotdec*np.pi/180.), plotdec, mindist=1./3600. )
                if smatch >= 0:
                    exec 'plot' + cfilter + 'mag[smatch] = compmag[j]'
                    exec 'plot' + cfilter + 'magerr[smatch] = compmagerr[j]'
                    # IDL ;plotcat(smatch(0))=1
                    count += 1
                else:
                    plotra[nstars] = compra[j]
                    plotdec[nstars] = compdec[j]
                    # IDL ;IF compcat(k) EQ 1 THEN plotcat(nstars) = 1
                    exec 'plot' + cfilter + 'mag[nstars] = compmag[j]'
                    exec 'plot' + cfilter + 'magerr[nstars] = compmagerr[j]'
                    nstars += 1
                    count += 1

    time2 = time.clock()
    print (time2 - time1)

    # output detections
    np.savetxt( 'finalmags.txt', np.array([plotra, plotdec, plotrmag, plotrmagerr, plotimag, plotimagerr, plotzmag, plotzmagerr, plotymag, plotymagerr, plotJmag, plotJmagerr, plotHmag, plotHmagerr]).T ) # IDL ;plotcat

    # plot each image with circles on star identification
    realdetections = np.where( (plotrmagerr > 0) & (plotimagerr > 0) )[0]
    imgarr = []
    cfilterarr = []
    for i in range(np.size(zffiles)):
        ifile = zffiles[i]
        hdulist = pf.open(ifile)
        h = hdulist[0].header
        img = hdulist[0].data
        if i == 0:
            refh = h
        im_size = np.shape(img)
        scalefactor = 1.
        img = zoom( img, 1./scalefactor, order=0 )
        imgarr.append(img)
        cfilter = zffiles[i].split('_')[1].split('.')[0] # extract filter label from file name
        if cfilter == 'Z' or cfilter == 'Y':
            cfilter = cfilter.lower()
        cfilterarr.append(cfilter)

    imgarr = np.array(imgarr)
    cfilterarr = np.array(cfilterarr)
    def find_where( arr, val ):
        arrc = np.copy(arr)
        pos = []
        while np.any( arrc == val ):
            pos.append(np.argmax( arrc == val ))
            arrc[pos[-1]] = '\0'
        if np.size(pos) == 0:
            return -1
        elif np.size(pos) == 1:
            return pos[0]
        else:
            return np.array(pos)

    cr = find_where(cfilterarr, 'r')
    ci = find_where(cfilterarr, 'i')
    cz = find_where(cfilterarr, 'z')
    cy = find_where(cfilterarr, 'y')
    cH = find_where(cfilterarr, 'H')
    cJ = find_where(cfilterarr, 'J')

    if cJ >= 0 and cH >= 0:
        r = imgarr[cJ,:,:] * 0.5 + imgarr[cH,:,:] * 0.5
    if cJ >= 0 and cH < 0:
        r = imgarr[cJ,:,:]
    if cH >= 0 and cJ < 0:
        r = imgarr[cH,:,:]
    if cH < 0 and cJ < 0:
        r = 0

    if cz >= 0 and cy >= 0:
        g = imgarr[cz,:,:] * 0.5 + imgarr[cy,:,:] * 0.5
    if cz >= 0 and cy < 0:
        g = imgarr[cz,:,:]
    if cy >= 0 and cz < 0:
        g = imgarr[cy,:,:]
    if cy < 0 and cz < 0:
        g = 0

    if cr >= 0 and ci >= 0:
        b = imgarr[cr,:,:] * 0.5 + imgarr[ci,:,:] * 0.5
    if cr >= 0 and ci < 0:
        b = imgarr[cr,:,:]
    if ci >= 0 and cr < 0:
        b = imgarr[ci,:,:]
    if ci < 0 and cr < 0:
        b = 0

    red = r
    green = g
    blue = b

    mi = np.min(blue)
    ma = np.max(blue) * 0.01 + mi
    if np.size(blue) > 1:
        blue = bytescale( blue, 0, 8, 250 )

    mi = np.min(green)
    ma = np.max(green) * 0.005 + mi
    if np.size(green) > 1:
        green = bytescale( green, 0, 8, 250 )

    mi = np.min(red)
    ma = np.max(red) * 0.004 + mi
    if np.size(red) > 1:
        red = bytescale( red, 0, 8, 250 )

    if np.size(red) > 1:
        im_size = np.shape(red)
    elif np.size(green) > 1:
        im_size = np.shape(green)
    elif np.size(blue) > 1:
        im_size = np.shape(blue)

    #height = 8 # height of figure in inches
    #pl.figure( 0, figsize=((float(im_size[1])/float(im_size[0]))*height, height) )

    def bytearr( x, y, z ):
        return np.zeros((x,y,z)).astype(np.uint8)
    color = bytearr( im_size[0], im_size[1], 3 )
    if np.size(red) > 1:
        color[:,:,0] = red * 0.5
    if np.size(green) > 1:
        color[:,:,1] = green * 0.5
    if np.size(blue) > 1:
        color[:,:,2] = blue * 0.5
    color = color[::-1,:]
    fig = pl.figure('color image')
    pl.axis('off')
    pl.imshow( color, interpolation='None', origin='lower' )
    sp.misc.imsave( 'color.png', color )

    for i in range(np.size(zffiles)):
        ifile = zffiles[i]
        hdulist = pf.open(ifile)
        h = hdulist[0].header
        img = hdulist[0].data
        cfilter = cfilterarr[i]
        scale = bytescale(img, 0, 10, 255)
        dpi = 72. # px per inch
        figsize = (np.array(img.shape)/dpi)[::-1]
        fig = pl.figure(i)
        pl.imshow( scale, interpolation='None', cmap=pl.cm.gray, origin='lower' )
        xlims = pl.xlim()
        ylims = pl.ylim()
        # Parse the WCS keywords in the primary HDU
        w = wcs.WCS(h)
        world = np.transpose([plotra, plotdec])
        pixcrd = w.wcs_world2pix(world, 1)
        fs = 12
        fw = 'normal'
        lw = 1
        for j in realdetections:
            ctemp = circle( pixcrd[j][0]/scalefactor, pixcrd[j][1]/scalefactor, 10 ).T
            pl.plot( ctemp[0], ctemp[1], c='#00ff00', lw=lw )
            if pixcrd[j][0]/scalefactor+40 < xlims[1] and pixcrd[j][1]/scalefactor+20 < ylims[1]:
                pl.text( pixcrd[j][0]/scalefactor+15, pixcrd[j][1]/scalefactor, `j`, color='#00ff00', fontsize=fs, fontweight=fw )
            else:
                pl.text( pixcrd[j][0]/scalefactor-30, pixcrd[j][1]/scalefactor-10, `j`, color='#00ff00', fontsize=fs, fontweight=fw )
        pl.text( 0.2*xlims[1], 0.9*ylims[1], cfilter+'-Band', color='r', fontsize=fs, fontweight=fw )
        a = pl.gca()
        a.set_frame_on(False)
        a.set_xticks([]); a.set_yticks([])
        pl.axis('off')
        pl.xlim(xlims)
        pl.ylim(ylims)
        fig.set_size_inches(figsize[0], figsize[1])
        ofile = ifile.split('.')[0] + '.png'
        pl.savefig( ofile, bbox_inches='tight', pad_inches=0, transparent=True, dpi=dpi )