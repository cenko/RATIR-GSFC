RATIR-GSFC
==========

GSFC/UMd RATIR Pipeline Repository

Directories:
code/ 
images/
sandbox/

Outline
-------

1. [SETUP](#setup)

2. [REDUCTION](#redux)

3. [PHOTOMETRY](#photo)

1. Setup {#setup}
--------
1.1 Requires (most installable using Macports):

	-IDL
	-Python
	-SExtractor
	-SWarp
	-cdsclient package

1.2 Run startup.sh in Unix shell to source code

1.3 Alter pipeautoproc.par to point to full path of autoastrometry.py on individual computer (under code/reduction/ratauto)



2. Reduction {#redux}
------------
2.1 Run preprocessing scripts (JOHN PUT INSTRUCTIONS IN)

2.2 Run ratautoproc.pro in data directory (in directory above processed data and reduction folder)

	IDL
	ratautoproc, datadir='raw/', redo=1
	
	datadir specifies where processed data is stored (should be in lower directory), will run all data in this directory  
	Will save to specified directory (imworkingdir) in pipeautoproc.par
	redo keyword overwrites previously reduction processed files
	Additional keywords:
		start, stop, step, only (allows you to run particular steps if you don't want to run full reduction)
		nocrclean (if set skips cosmic ray cleaning)	
	
	Runs these steps in this order unless specified:	
		steps = ['prepare', 'flatten', 'makesky', 'skysub', 'crclean', 'astrometry', 'stack']
		
	1) Prepare (autopipeprepare.pro, pipeprepare.pro)
		Header information manipulation and bias subtraction for images with master bias files
	
	2) Flatten (autopipeimflatten.pro, flatpipeproc.pro)
		Divides master flat for each filter
		
	3) Makesky (autopipemakesky.pro, skypipecombine.pro)
		Creates master sky by removing sources (using outlier rejection) and then median iterative sigma clipping for each pixel
		
	4) Skysub (autopipeskysub.pro, skypipeproc.pro)
		Subtracts sky, then subtracts median, then adds 1000
	
	5) Crclean (autopipecrcleanim.pro, pzap_perley.pro)
		Cosmic ray cleaning, details fuzzy
		
	6) Astrometry (autopipeastrometry.pro, vlt_autoastrometry.py)
		Fixes WCS coordinates by using pair-distance matching and asterism matching
	
	7) Stack (autopipestack.pro)
		Uses SWarp to stack images with same filter

3. Photometry {#photo}
-------------
3.1 Run autoredux.py in photometry folder

Can automatically run full photometry reduction using autoredux.py
Needs to run inside directory with coadd*.fits files:

	python
	import autoredux
	autoredux.autoredux()

autoredux.py runs the following programs in this order:

	1) Identify point sources in individual frames after cropping (automatic crop
    	using weight files, can run manual crop using manualcrop keyword).
		icoords.py
		Output: initially cropped files (.crop.fits)
		  	coords files for w/ RA and DEC identified by sextractor

	2) Run photometry and calculate zero points
		calcoff.py

	3) Run association program to make a master star list
		assoc.py

	4) run final photometry
		finalphot.py

	5) create webpage
		plotratir.py



The sandbox was meant as a place to test new code developement.  Note that any changes made to this directory will not be reflected in the repository.

