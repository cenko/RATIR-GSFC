import glob
import os
import pyfits as pf
import numpy as np
import autoproc_depend as apd

inpipevar = {'autoastrocommand':'autoastrometry', 'getsedcommand':'get_SEDs', 
			'sexcommand':'sex' , 'swarpcommand':'swarp' , 'rmifiles':0,  
			'prefix':'', 'datadir':'' , 'imworkingdir':'' , 'overwrite':0 , 'verbose':1, 
			'flatfail':'' , 'catastrofail':'' , 'relastrofail':'' , 'fullastrofail':'' ,
			'pipeautopath':'' , 'refdatapath':'', 'defaultspath':'' }

def autopipedefaults(pipevar=inpipevar):

    """
    NAME:
        autopipedefaults
    PURPOSE:
        Sets commonly used variables for pipeautoproc to use throughout each step
        Uses pipeautoproc.par to set variables, otherwise set to default values
        Saves in a dictionary
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro, 
                    but can be set to default) 
    EXAMPLE:
        autopipedefaults(pipevar=inpipevar)
    """

    print 'Setting pipeline parameters (DEFAULTS)'
    
    path = os.path.dirname(__file__)
        
    pipevar['pipeautopath'] = path
    sfile = path+'/pipeautoproc.par'
    
    if os.path.isfile(sfile):
        f = open(sfile,'r')
        for line in f.readlines():
            line = ''.join(line.split())
            colpos = line.find(':')
            
            keyword = line[0:colpos]
            value = line[colpos+1:]
            pipevar[keyword] = value
        f.close()

    if pipevar['refdatapath'] == '': 
        pipevar['refdatapath'] = pipevar['pipeautopath']+'/refdata'

    if pipevar['defaultspath'] == '': 
        pipevar['defaultspath'] = pipevar['pipeautopath']+'/defaults'

    if pipevar['imworkingdir'] != '' and not(os.path.exists(pipevar['imworkingdir'])): 
        print 'Creating imaging working directory: ',  pipevar['imworkingdir']
        os.makedirs(pipevar['imworkingdir'])
        
def autopipeprepare(pipevar=inpipevar):

    """
    NAME:
        autopipeprepare
    PURPOSE:
        Runs pipeprepare on every valid file and saves files with prefix 'p'.  Changes 
        header with more manageable keywords and does bias/dark subtraction if bias/dark 
        master exists (compares header keywords in files and bias/dark master)
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro, 
                   but can be set to default)
    EXAMPLE:
        autopipeprepare(pipevar=inpipevar)
    DEPENDENCIES:
        autoproc_depend.pipeprepare()
    FUTURE IMPROVEMENTS:
        MOST REFER TO pipeprepare.pro: Need to check what additional keywords need to 
        propagate, and check if values that are set with magic numbers can be set from 
        existing keywords.
        Check pipeprepare for RIMAS, RATIR, or VLT/VT to see changes that 
        need to be made for RIMAS pipeline
    """
    
    print 'PREPARE'
    
    # Looks for existing files in given data directory using prefix
    files = glob.glob(pipevar['datadir'] + pipevar['prefix'] + '*.fits')
    pfiles = glob.glob(pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
    
    if len(files) == 0:
        print 'Did not find any files! Check your data directory path!'
        return
    
    if pipevar['verbose']: print 'Found', len(files), 'files'
    
    # Finds any master bias files and filter name from header keyword
    # Assumes camera name is in header under CAMERA
    biasfiles = glob.glob(pipevar['imworkingdir'] + 'bias*')
    biascamera = []
    if len(biasfiles) > 0:
        for bfile in biasfiles:
            head = pf.getheader(bfile)  
            camera = head['CAMERA'] 
            biascamera += [camera]
            
    # Finds any master dark files and filter name from header keyword
    # Assumes camera name is in header under CAMERA
    darkfiles = glob.glob(pipevar['imworkingdir'] + 'dark*')
    darkcamera = []
    if len(darkfiles) > 0:
        for dfile in darkfiles:
            head = pf.getheader(dfile)  
            camera = head['CAMERA'] 
            darkcamera += [camera]     
        
    # For each file (that doesn't have an existing p file or can be overwritten), 
    # run pipeprepare on it with output file being saved into the imworkingdir, 
    # will run bias subtraction if bias master available (checks based on how bias 
    # file and data file are named
    for file in files:
        print file
        fileroot = os.path.basename(file)
        outnameim = pipevar['imworkingdir'] + 'p' + fileroot
        
        head = pf.getheader(file)
        camera = head['CAMERA']
        
        try:
            bcamloc  = biascamera.index(str(camera))    
            biasfile = biasfiles[bcamloc]      
        except:
            biasfile = None
        
        try:
            dcamloc  = darkcamera.index(str(camera))    
            darkfile = darkfiles[bcamloc]      
        except:
            darkfile = None

        if (outnameim not in pfiles) or (pipevar['overwrite'] != 0):
            apd.pipeprepare(file, outname=outnameim, biasfile=biasfile, darkfile=darkfile,
                             verbose = pipevar['verbose'])
        else:
            print 'Skipping prepare. File already exists'
            
def autopipeimflatten(pipevar=inpipevar):
    
    """
    NAME:
        autopipeflatten
    PURPOSE:
        Flatten data using flat with matching filter name
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro, 
                   but can be set to default)
    EXAMPLE:
        autopipeflatten(pipevar=inpipevar)
    DEPENDENCIES:
        autoproc_depend.flatpipeproc()
    """
    
    print 'FLATTEN'
    
    # Finds prepared files and checks to see if there are any existing flattened files
    # Find flats in imworkingdir with name flat somewhere in a fits file name
    files  = glob.glob(pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
    ffiles = glob.glob(pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
    flats  = glob.glob(pipevar['imworkingdir'] + '*flat*.fits')

    if len(files) == 0:
        print 'Did not find any files! Check your data directory path!'
        return
    
    # If there are flats, then grab the filter from each of them, 
    # otherwise end program
    flatfilts = []
    if len(flats) > 0:
        for flat in flats:
            head = pf.getheader(flat)
            filter = head['FILTER']
            flatfilts += [filter]
    else:
        print 'No flats found for any filter!'
        return
    
    # Create outfile name and check to see if outfile already exists.  If it doesn't or
    # overwrite enabled then take filter from file and find where the flat filter matches
    # If no flats match filter, store in pipevar.flatfail, otherwise run flatproc on file
    for file in files:
        print file
        fileroot = os.path.basename(file)
        outnameim = pipevar['imworkingdir'] + 'f' + fileroot
        
        if (outnameim not in ffiles) or (pipevar['overwrite'] != 0):
            head = pf.getheader(file)
            filter = head['FILTER']

            try:
                flatfileno = flatfilts.index(filter)
            except:
                print 'Flat field not found for '+ file +' (filter='+filter+')'
                pipevar['flatfail'] += ' ' + file
                continue
            
            flatfile = flats[flatfileno]
        
            if pipevar['verbose']: print 'Flattening', file, 'using', flatfile
            
            apd.flatpipeproc(file, flatfile, flatminval=0.3)
        
        else:
            print 'Skipping flatten. File already exists'
            
    # If remove intermediate files keyword set, delete p(PREFIX)*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        
def autopipemakesky(pipevar=inpipevar):
    """
    NAME:
        autopipemakesky
    PURPOSE:
        Combine sky flats based on filter type (sigma clipping for sources)
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro, 
                   but can be set to default)
    EXAMPLE:
        autopipemakesky(pipevar=inpipevar)
    DEPENDENCIES:
        astroproc_depend.skypipecombine, astroproc_depend.medclip, Sextractor
    FUTURE IMPROVEMENTS:
        skypipecombine slow, find better algorithm
    """   
    
    print 'MAKE SKY'
    
    # Copies necessary parameter file for sextractor if not in current working directory
    if not os.path.isfile('source.param'): 
        os.system('cp ' + pipevar['defaultspath'] + '/source.param .')
    if not os.path.isfile('sex_source.config'): 
        os.system('cp ' + pipevar['defaultspath'] + '/sex_source.config .')
    if not os.path.isfile('sex.conv'): 
        os.system('cp ' + pipevar['defaultspath'] + '/sex.conv .')
    if not os.path.isfile('defaulf.nnw'): 
        os.system('cp ' + pipevar['defaultspath'] + '/default.nnw .')
        
    # Finds files with given prefix
    files  = glob.glob(pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')

    if len(files) == 0:
        print 'Did not find any files! Check your data directory path!'
        return
    
    filters = []
    for file in files:
        head = pf.getheader(file)
        filter = head['FILTER']
        filters += filter
    filters = np.array(filters)
    
    # Unique list of filters
    filterlist = set(filters)
    
    # For each unique filter, combine sky files using skycombine if more than 2 files
    # Otherwise return list of unprocessed files
    for filt in filterlist:
        skyflats = np.where(filters == filt)
        outflatname = pipevar['imworkingdir'] + 'sky-' + filt + '.fits'
        
        if len(skyflats[0]) >= 2:
        
            if os.path.isfile(outflatname) and pipevar['overwrite'] == 0:
                print 'Skipping makesky for '+filt+'. File already exists'
                continue

            files = np.array(files)
            if pipevar['verbose']:
                print filt, '-band sky flats.'
                print files[skyflats]
            
                apd.skypipecombine(files[skyflats], outflatname, file,
                    pipevar, removeobjects=True, type='sky')
        else:
            print 'Unable to produce a flat field for this setting: ' + filt
            print 'Will not be able to further process ' + str(len(skyflats)) + \
                  ' image(s) without a flat from another source:'

            for i in np.arange(len(skyflats[0])):
                print '    ' + files[skyflats[i]]           

    # If remove intermediate files keyword set, delete p(PREFIX)*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        
def autopipeskysub(pipevar=inpipevar):
    """
    NAME:
        autopipeskysub
    PURPOSE:
        Subtracts master sky flat from data and subtracts median.
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro, 
                   but can be set to default)
    EXAMPLE:
        autopipeskysub(pipevar=inpipevar)
    DEPENDENCIES:
        autoproc_depend.skypipeproc
    """
    
    print 'SKY-SUBTRACT'
    
    # Find data that needs to be sky subtracted
    files  = glob.glob(pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
    sfiles = glob.glob(pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')
    
    if len(files) == 0:
        print 'Did not find any files! Check your data directory path!'
        return
    
    skys = glob.glob(pipevar['imworkingdir'] + '*sky-*.fits')
    
    if len(skys) == 0:
        print 'No master sky files found, cannot sky subtract'
        return
        
    # Find the associated filter of each master skyflat
    skyfilts = []
    for sky in skys:
        head = pf.getheader(sky)
        filter = head['FILTER']
        skyfilts += [filter]
    
    # For each file if output files don't exist or override set check if we have master 
    # skyflat for filter, sky subtract if it exists using skypipeproc
    for file in files:

        fileroot = os.path.basename(file)
        outfile = pipevar['imworkingdir'] + 's' + fileroot    
    
        if os.path.isfile(outfile) and pipevar['overwrite'] == 0:
            print 'Skipping sky subtraction for '+file+'. File already exists'
            continue
            
        head = pf.getheader(file)
        filter = head['FILTER'] 

        # Find corresponding master skyflat
        try:
            skyloc  = skyfilts.index(filter)    
            skyfile = skys[skyloc]      
        except:
            print 'Sky field not found for ', file
            pipevar['flatfail'] += ' ' + file
            continue
        
        if pipevar['verbose'] > 0:
            print 'Sky subtracting', file, 'using', skyfile
        
        apd.skypipeproc(file, skyfile, outfile)


    # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
    # and sky-*.fits files
    if pipevar['rmifiles'] != 0:
        
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')


def autopipecrcleanim(pipevar=inpipevar):
    
    """
    NAME:
        autopipecrcleanim
    PURPOSE:
        Removes cosmic rays
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro, 
                   but can be set to default)
    EXAMPLE:
        autopipecrcleanim(pipevar=inpipevar)
    DEPENDENCIES:
        autoproc_depend.cosmiczap
    FUTURE IMPROVEMENTS:
        Slow, alter cosmics.py?
        Get readnoise from header
    """
    
    print 'CRCLEAN'

    # Find data that needs to be cosmic ray zapped
    files  = glob.glob(pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')
    zfiles = glob.glob(pipevar['imworkingdir'] + 'zsfp' + pipevar['prefix'] + '*.fits')
    
    if len(files) == 0:
        print 'Did not find any files! Check your data directory path!'
        return
 
    # For each file check that objects meet count limits and exposure time
   	# (i.e. short exposure time with lot of counts will be ignored), also targets that are
   	# calibration files will be ignored.
   	# Run cosmiczap on the files and have output files be 'z'+file plus weight files
    
    for file in files:
    
        head = pf.getheader(file)
                
        try:
            target  = head['TARGNAME']
        except:
            print 'Requires header keywords: TARGNAME. Check file.'
            continue
        
        if 'flat' in target.lower(): continue
        if 'twilight' in target.lower(): continue

        fileroot = os.path.basename(file)
        outfile = pipevar['imworkingdir'] + 'z' + fileroot 
          
        if os.path.isfile(outfile) and pipevar['overwrite'] == 0:
            print 'Skipping crzap for '+file+'. File already exists'
            continue 
                
        if pipevar['verbose'] > 0:
            print 'Cleaning cosmic rays from', file      
        
        # Runs cosmics.py
        apd.cosmiczap(file, outfile, sigclip=6.0, maxiter=3, verbose=pipevar['verbose'])
        
        
    # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
    # sky-*.fits, sfp(PREFIX)*.fits files
    if pipevar['rmifiles'] != 0:
        
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')