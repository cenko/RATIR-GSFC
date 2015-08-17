import os
import glob
import pyfits as pf
import numpy as np
import datetime
import astrometrystats as astst
import cosmics

def pipeprepare(filename, outname=None, biasfile=None, darkfile=None, verbose=1):

    """
    NAME:
        pipeprepare
    PURPOSE:
        Adds additional header keywords needed later on in the pipeline and removes 
        unnecessary header keywords by looking through a list of mandatory keywords.  
        Also runs bias and dark subtraction for filters with an existing master bias/dark 
        (CCDs).  The prepared images are written to disk with outname
    INPUTS:
        filename - name of FITS file, or array of filenames, or file w/list of filenames 
    OPTIONAL KEYWORDS:
        outname  - specify output file to write to disk
        biasfile - full name (including path) of master bias
        darkfile - full name (including path) of master dark
        verbose  - print out comments        
    EXAMPLE:
        pipeprepare(filename, outname=outname, biasfile=biasfile, darkfile=darkfile, verbose=1)
    DEPENDENCIES:
        autoproc_depend.pipeprepare()
    FUTURE IMPROVEMENTS:
        Need to check what additional keywords need to propagate, and check if values 
        that are set with magic numbers can be set from existing keywords.  
        Change airmass keywords for RIMAS
    """
    
    # ------ Process input filenames(s) ------
    
    # Check for empty filename
    if len(filename) == 0:
        print 'No filename specified'
        return

    # If string, check if a file of items or wildcards
    # otherwise store all files
    if isinstance(filename,str):
        fileext = os.path.splitext(filename)[1][1:]
        
        files = [filename]
        
        if fileext in ['cat', 'lis', 'list', 'txt']:
            f = open(filename,'r')
            files = f.read().splitlines() 
            f.close()            
                      
        if '?' in filename or '*' in filename:
            files = glob.glob(filename)
            if len(files) == 0: 
                print 'Cannot find any files matching ', filename
                return
    else:
        files = filename


    # ------ Read data and process header information ------        
    for file in files: 
        f = pf.open(file)
        head = f[0].header
        data = f[0].data
        f.close()
        
        # Grabs starting airmass, can alternatively use ETROBAM (ending time observed airmass)
        # and saves as AIRMASS.  CHANGE FOR RIMAS
        head['AIRMASS'] = head['STROBAM']
        head['DATE-OBS'] = head['SDATE']
        
        mandatorykey = ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2',
 					    'HISTORY','DATE-OBS','EXPOSURE','EXPTIME','INSTRUME',
 					    'ORIGIN','LATITUDE','LONGITUD',
 					    'CCD_TYPE','CCD_SER','SATURATE','BINNING','BINY','BINX',
 					    'WAVELENG','TARGNAME','CAMERA','UTC','UT','OBJECT','PIXSCALE',
 					    'SUN_ALT','SMNSP','CD1_1','CD1_2','CD2_1','CD2_2',
 					    'CRPIX1','CRPIX2','CRVAL1','CRVAL2','CTYPE1','CTYPE2', 
 					    'SOFTGAIN','FILTER','AVERAGE','STDEV','GAIN','AIRMASS','CCD_NAME', 
 					    'PV1_1','PV2_1','PV1_17','PV2_17','PV1_19','PV2_19','PV1_21','PV2_21',
 					    'PV1_31','PV2_31','PV1_33','PV2_33','PV1_35','PV2_35','PV1_37','PV2_37']
 					    
 	    # Finds list of unnecessary keywords, then deletes extraneous
        newhead = head
        for oldkey in head.keys():
            if oldkey not in mandatorykey:
                del newhead[oldkey]
        
        # If biasfile keyword set subtract master bias from current file with given master bias file
        # If they are not the same size, quit program without saving with preparation prefix (will not move
        # on in following processing steps)
        if biasfile != None:
            bf = pf.open(biasfile)
            bias = bf[0].data
            bf.close()
            
            if np.shape(data) != np.shape(bias):
                
                print file + ' could not be bias subtracted because it is not the same' +\
                             ' size as the master bias, remove file to avoid confusion'
                return
            
            if verbose > 0: print '    bias subtracting'
            
            newdata = data - bias
            
            # If darkfile keyword set subtract master dark from current file with given master dark file
            # If they are not the same size, quit program without saving with preparation prefix (will not move
            # on in following processing steps)
            if darkfile != None:
                ff = pf.open(darkfile)
                dark = ff[0].data * newhead['EXPTIME']
                ff.close()
                
                if np.shape(data) != np.shape(dark):
                    print ' '
                    print file + ' could not be dark subtracted because it is not the same' +\
                                 ' size as the master dark, remove file to avoid confusion'
                    return  
                          
                if pipevar['verbose'] > 0: print '    dark subtracting'
                
                newdata = newdata - dark
            else:
                print file, 'could not be dark subtracted because the master dark file was not provided'
        else:
            newdata = data
        
        # Write changes to disk
        pf.writeto(outname, newdata, newhead)
        
        if verbose > 0: print file, '-> ', outname
        
def flatpipeproc(filename, flatname, flatminval=0, flatmaxval=0):

    """
    NAME:
        flatpipeproc
    PURPOSE:
        Checks if flat is same size as data, then divides for correct filter
    INPUTS:
        filename - name of FITS file, or array of filenames, or file w/list of filenames 
        flatname - name of FITS master flat file
    OPTIONAL KEYWORDS:
        flatminval - if not set to 0 below this value will set to NaNs
        flatmaxval - if not set to 0 above this value will set to NaNs
    EXAMPLE:
        flatpipeproc(filename, flatname, flatminval=0.3)
    """

    # ------ Process input filenames(s) ------
    
    # Check for empty filename
    if len(filename) == 0:
        print 'No filename specified'
        return

    # If string, check if a file of items or wildcards
    # otherwise store all files
    if isinstance(filename,str):
        fileext = os.path.splitext(filename)[1][1:]
        
        files = [filename]
        
        if fileext in ['cat', 'lis', 'list', 'txt']:
            f = open(filename,'r')
            files = f.read().splitlines() 
            f.close()            
                      
        if '?' in filename or '*' in filename:
            files = glob.glob(filename)
            if len(files) == 0: 
                print 'Cannot find any files matching ', filename
                return
    else:
        files = filename
        
    f = pf.open(flatname)
    flat = f[0].data
    f.close()
    
    med = np.median(flat)
    if (med < 0.5) or (med > 2.0): print 'Warning: flat is not normalized to one'
    
    for file in files:
        f = pf.open(file)
        data = f[0].data
        head = f[0].header
        f.close()
        
        if np.shape(data) != np.shape(flat):
            print file + ' could not be dark subtracted because it is not the same' +\
                         ' size as the master dark, remove file to avoid confusion'
            return  
        
        # Set values too low/high to NaNs
        if flatminval > 0:
            flat[flat < flatminval] = float('NaN')
        goodsignal = np.where(flat-1.0 < 0.1)
   
        if flatmaxval > 0:
            flat[flat > flatminval] = float('NaN')
        
        # Divides out flattened field and adds keywords to header to show change
        fdata = data / flat         
        
        head['FLATFLD'] = flatname
        skycts = np.median(fdata[goodsignal])        
        head['SKYCTS']  = (skycts, 'Sky counts')
        
        try:
            head['CTRATE'] = (skycts/head['EXPTIME'], 'Sky counts per second')
        except:
            print 'No EXPTIME keyword'
            
        date = datetime.datetime.now().isoformat()
        head.add_history('Processed by flatproc ' + date)
        
        fileroot = os.path.basename(file)
        filedir  = os.path.dirname(file)
        outnameim = filedir + '/f' + fileroot
        
        pf.writeto(outnameim, fdata, head)
        
def skypipecombine(filelist, outfile, filt, pipevar, removeobjects=None, 
    objthresh=6, algorithm='median', trimlo=None, trimhi=None, mincounts=1, 
    maxcounts=55000, satlevel=30000, type=None):
    
    """
    NAME:
        skypipecombine
    PURPOSE:
        Create sigma clipped median sky flat.  Scales each file based on the overall 
        sigma clipped median, then removes objects selected with sextractor (uses flux 
        fraction radius) in each file. Removes saturated pixels.  Calculates sigma clipped 
        median of each pixel and saves anything with non-finite values (saturated or 
        source) to the median of the entire frame.  Save with outfile name.
    INPUTS:
        filelist - files to be processed
        outfile  - name for output fits file
        filt	 - filter of files
        pipevar  - pipeline parameters in dictionary   
    OPTIONAL KEYWORDS:
        removeobjects 	- specifies if you want objects removed
        objthresh		- sets sigma in removeobjects (default is 6)
        algorithm       - algorithm to solve (mean or median, default is median)
        trimlo			- trim off bottom of data in mean algorithm mode (default is 25%)
        trimhi			- rim off top of data in mean algorithm mode (default is 25%)
        mincounts		- sets minimum counts allowed (default is 1)
        maxcounts		- sets maximum counts allowed (default is 55000)
        satlevel		- sets saturation level (default is 30000)
        type			- sets 'SKYTYPE' keyword in header of outfile to this string
    EXAMPLE:
        skypipecombine(filelist, 'sky-filt.fits', removeobjects=True, type='sky')
    DEPENDENCIES:
        medclip, Sextractor
    FUTURE IMPROVEMENTS:
        medclip slow find faster solution
        Need to take saturation level from header?
        Saved header is from middle file, maybe use blank?
    """
    
    # Sets defaults for trimming (25% of list)
    if trimlo != None: (len(filelist)+1)/4
    if trimhi != None: trimlo
    
    # If given list, then grab all filenames, saved to files
    if len(filelist) == 1:
        f = open(filelist,'r')
        files = f.read().splitlines() 
        f.close()
    else:
        files = filelist
    
    nfiles = len(files)
    nmid = len(files)/2
    
    # Read in middle file and initialize arrays
    f = pf.open(files[nmid])
    data_m = f[0].data
    head_m = f[0].header
    f.close()
    
    nx = head_m['NAXIS1']
    ny = head_m['NAXIS2']
    
    data     = np.zeros((nfiles, ny, nx)) + float('NaN')
    skymeds  = []
    usefiles = []
    
    z = 0
    # For each file and make sure size matches middle file, calculate sigma clipped 
    # median (3sig, 6 iter), then if within counts limit save data into 3d data cube 
    # and save clipped median into skymed, and mark file as usable
    # Increment z by one when this is true
    for file in files:
        f = pf.open(file)
        data_i = f[0].data
        head_i = f[0].header
        f.close()        
        
        inx = head_i['NAXIS1']
        iny = head_i['NAXIS2']        

        if (inx != nx) or (iny != ny):
            print 'File ' + file + ' has wrong dimensions ('+str(inx)+ \
                  ' x '+ str(iny)+'; should have '+str(nx)+' x '+str(ny)+')'
              
        # Perform 3 sigma clipped median and save to inmeds
        inmed, instd = medclip(data_i, clipsig=3, maxiter=3)
        
        
        # If median is within limits save data, otherwise exclude files
        if inmed >= mincounts and inmed <= maxcounts:
            if pipevar['verbose'] > 0:
                print file + ' ('+str(inmed) + ' counts/pix)'
            
            skymeds  += [inmed]
            usefiles += [file]
            data[z, :,:] = data_i
            z += 1
        else:
            if inmed < mincounts:
                print file + ' (' + str(inmed) + ' counts/pix) - too few counts; excluding'
            if inmed > maxcounts:
                print file + ' (' + str(inmed) + ' counts/pix) - too many counts; excluding'            
        
    if z < 2:
        print 'ERROR - Not enough counts to make a flat with these data!'
        return
    
    # Median of sigma clipped medians
    medsky = np.median(skymeds)
    
    # Scale each file by median of sigma clipped medians divided by median of data
    # Corrects for each flat's changing sky background
    for f in np.arange(z-1):
        factor = medsky / skymeds[f]
        data[f,:,:] = data[f,:,:]*factor
    
    # Removes extraneous indexes in data for skipped files
    if z != nfiles: data = data[0:z,:,:]
    
    # Removes objects from field by calculating iterative median sigma clipping 
    # (5 sigma, 5 iter) and using the calculated stddev to remove 6sigma (or non-default 
    # object threshold) data from the median along with values above the saturation limit.
    # Find 3sigma clipped median of each pixel from remaining values or mean of middle 50%
    
    if removeobjects != None:
        if pipevar['verbose'] > 0: print '  Identifying objects...'
                
        for f in np.arange(z-1):
            
            indata = data[f,:,:]
            
            # Set sources above objthresh  limit to NaN
            datamed, datastd = medclip(indata, clipsig=5, maxiter=5)
            sourcepixels = np.where( abs(indata-datamed) >= objthresh*datastd)
            
            if len(sourcepixels[0]) > 0:
                indata[sourcepixels] = float('NaN')
            
            satpixels = np.where( indata >= satlevel )
            
            if len(satpixels[0]) > 0:
                indata[satpixels] = float('NaN')
            
            data[f,:,:] = indata
        
        reflat = np.zeros((ny, nx)) + float('NaN')
            
        # If algorithm set to median, find 3 sigma clipped median of each pixel 
        # excluding NaN values (which are eventually set to median)
        if algorithm == 'median':
            if pipevar['verbose'] > 0: print '  Median-combining...'
            
            for y in np.arange(ny):
                for x in np.arange(nx):
                    vector = data[:,y,x]
                    temp = np.isfinite(vector)
                    
                    if len(vector[temp]) == 2:
                        reflat[y,x] = np.median(vector[temp])
                        continue                    
                        
                    if len(vector[temp]) < 2: continue
                    
                    me, st = medclip(vector[temp], clipsig=3, maxiter=5)
                    reflat[y,x] = me

            # Replace bad pixels with median of entire sky
            good = np.isfinite(reflat)
            allmed = np.median(reflat[good])
            bad = ~good # Opposite of boolean array good
            reflat[bad] = allmed
        
        # If algorithm set to mean, takes mean of trimmed sorted values. Default is to 
        # trim 25% off top and bottom, if not enough good data, set trimming to 0
        if algorithm == 'mean':
        
            if pipevar['verbose'] > 0: print '  Combining via trimmed mean...'
            
            for y in np.arange(ny):
                for x in np.arange(nx):
                    slice = data[:,y,x]
                    good = np.isfinite(slice)
                    
                    cslice = slice[good]
                    ctgood = len(cslice)
                    
                    if ctgood == 0:
                        reflat[y,x] = 1
                    
                    itrimlo = trimlo
                    itrimhi = trimhi
                    
                    while ctgood-itrimlo-itrimhi < 1:
                        itrimlo = max(itrimlo - 1, 0)
                        itrimhi = max(itrimhi - 1, 0)
                        
                    cslice = np.sort(cslice)
                    cslice = cslice[itrimlo:ctgood-itrimhi]
                    reflat[y,x] = np.mean(cslice)
        
        # Adds header information to signify what files we used 
        for f in np.arange(z-1):
            head_m['SKY'+str(f)] = usefiles[f]
            
        
        if type != None:
            head_m['SKYTYPE'] = type
        
        date = datetime.datetime.now().isoformat()
        head_m.add_history('Processed by skypipecombine ' + date) 
        
        if pipevar['verbose'] > 0: print '  Median-combining...'
        
        pf.writeto(outfile, reflat, head_m)       
                    
def skypipeproc(filename, flatname, outfile, flatminval=None, flatmaxval=None):

    """
    NAME:
        skypipeproc
    PURPOSE:
        Subtracts sky flat from data and then subtracts median of that from remaining data. 
    INPUTS:
    	filename - file or list of files to be sky subtracted
    	flatname - sky flat fits file 
    	outfile  - name of output file
    OPTIONAL KEYWORDS:
        flatminval - minimum required value in flat (default for skycts calculation is 0.1)
        flatmaxval - maximum required value in flat
    EXAMPLE:
        skypipeproc(filename, flatname, outfile)        
    """
    
    # ------ Process input filenames(s) ------
    
    # Check for empty filename
    if len(filename) == 0:
        print 'No filename specified'
        return

    # If string, check if a file of items or wildcards
    # otherwise store all files
    if isinstance(filename,str):
        fileext = os.path.splitext(filename)[1][1:]
        
        files = [filename]
        
        if fileext in ['cat', 'lis', 'list', 'txt']:
            f = open(filename,'r')
            files = f.read().splitlines() 
            f.close()            
                      
        if '?' in filename or '*' in filename:
            files = glob.glob(filename)
            if len(files) == 0: 
                print 'Cannot find any files matching ', filename
                return
    else:
        files = filename  
    
    # Open flat    
    f = pf.open(flatname)
    flat = f[0].data  
    f.close()
    
    med = np.median(flat)
    
    # For each input file check if same size as flats (required). If there is a minimum 
    # or maximum flat value set, forces values outside of that range to NaN. Use finite 
    # values above 0.1 to determine skycounts, and subtract flat along with median of 
    # flattened data. Saves to new fits file
    for file in files:
        f = pf.open(file)
        data = f[0].data
        head = f[0].header
        f.close()
        
        if np.shape(data) != np.shape(flat):
            print file + ' could not be flat subtracted because it is not the same' +\
                         ' size as the master flat, remove file to avoid confusion'
            return
                
        if flatminval != None:
            w = np.where(flat < flatminval)  
            if len(w[0]) != 0:
                flat[w] = float('NaN')
            goodsignal = np.where(np.logical_and(flat >= flatminval, np.isfinite(flat)))
        else:
            goodsignal = np.where(np.logical_and(flat >= 0.1, np.isfinite(flat)))

        if flatmaxval != None:
            w = np.where(flat > flatminval) 
            if len(w[0]) != 0:
                flat[w] = float('NaN')            
                
        # Scale skyflat, subtract scaled skyflat, and subtract median of subsequent flat 
        # subtracted data. Calculate skycounts from data (above minimum, or 
        # by default above 0.1)
        flattmp = np.median(flat[np.isfinite(flat)])
        imgtmp  = np.median(data)
        
        scalefr = imgtmp/flattmp
        fdata   = data - scalefr * flat
        
        tmp     = np.median(fdata)
        fdata   = fdata - tmp

        skycts  = np.median(fdata[goodsignal])
        
        # Adds header keywords to denote new median counts and file we used to flatfield
        head['SFLATFLD'] = flatname
        head['SKYCTS']   = (skycts, 'Sky counts')
        
        try:
            head['CTRATE'] = (skycts/head['EXPTIME'], 'Sky counts per second')
        except:
            print 'No EXPTIME keyword'        

        date = datetime.datetime.now().isoformat()
        head.add_history('Processed by skypipeproc ' + date)
        
        pf.writeto(outfile, fdata, head) 


def cosmiczap(filename, outname, sigclip=6.0, maxiter=3, verbose=True):

    """
    NAME:
        cosmiczap
    PURPOSE:
        Removes cosmic rays using Laplacian cosmic ray identification written in cosmics.py 
    INPUTS:
    	filename - file or list of files to be cosmic ray zapped
    	outfile  - name of output file
    OPTIONAL KEYWORDS:
        sigclip  - sigma to clip
        maxiter  - maximum number of times to iterate loop
        verbose  - quiet?
    EXAMPLE:
        cosmiczap(filename, outname)
    DEPENDENCIES:
        cosmic.py (described in http://arxiv.org/pdf/1506.07791v3.pdf)  
    FUTURE IMPROVEMENTS:
        Read readnoise from header?    
    """
    
    data, head = cosmics.fromfits(filename, verbose=False)
    
    gain = head['GAIN']
    c = cosmics.cosmicsimage(data, gain=gain, readnoise=18, sigclip=sigclip,
        sigfrac = 0.4, objlim = 5.0, verbose=False)
    
    tot = c.run(maxiter=maxiter, verbose=False)
    
    head['NPZAP'] = (tot, "Num. of pixels zapped by cosmiczap")
    date = datetime.datetime.now().isoformat()
    head.add_history('Processed by cosmiczap ' + date)    
    
    if verbose: print '  Zapped %d total affected pixels (%.3f%% of total)' \
                      %(tot,tot*100.0/np.size(data))
    
    cosmics.tofits(outname, c.cleanarray, head, verbose=False)


def medclip(indata, clipsig=3.0, maxiter=5, verbose=0):

    """
    NAME:
        medclip
    PURPOSE:
        Median iterative sigma-clipping
    INPUTS:
        indata - array to be clipped
    OPTIONAL KEYWORDS:
        clipsig - sigma to clip around
        maxiter - maximum number of times to clip
        verbose - allow to print messages
    EXAMPLE:
        med, sigma = medclip(indata, sigma=5.0)
    """
    
    # Flatten array
    skpix = indata.reshape( indata.size, )
 
    ct = indata.size
    iter = 0
    numrej = len(skpix)
    ndata  = len(skpix)
    
    while (iter < maxiter) and (numrej > min(ndata*0.01, 50)):
        lastct = ct
        medval = np.median(skpix)
        sig = np.std(skpix)
        wsm = np.where( abs(skpix-medval) < clipsig*sig )
        ct = len(wsm[0])
        if ct > 0:
            skpix = skpix[wsm]
 
        numrej = abs(ct - lastct)
        if ct <=2: return 'Too few remaining'
        iter += 1
 
    med   = np.median( skpix )
    sigma = np.std( skpix )
 
    if verbose:
        print '%.1f-sigma clipped median' % (clipsig)
        print 'Mean computed in %i iterations' % (iter)
        print 'Mean = %.6f, sigma = %.6f' % (med, sigma)
 
    return med, sigma
