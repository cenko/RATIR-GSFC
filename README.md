RATIR-GSFC
==========

GSFC/UMd RATIR Pipeline Repository

Directories:

* code/

* images/

* sandbox/

Outline
-------

1. [Setup](https://github.com/cenko/RATIR-GSFC#1-setup)

2. [Reduction](https://github.com/cenko/RATIR-GSFC#2-reduction)

3. [Photometry](https://github.com/cenko/RATIR-GSFC#3-photometry)

1. Setup
--------

1. Dependencies:

	- IDL

	- Python

		+ **ADD PACKAGE DEPENDENCIES?**

	- SExtractor

	- SWarp

	- cdsclient package  

	Most can be installed using Macports.  

2. Run *startup.sh* in Unix shell to source code:

	```bash
	source startup.sh
	```

3. Alter *pipeautoproc.par* to point to full path of *autoastrometry.py* on individual computer (under code/reduction/ratauto)  

2. Reduction
------------

### 2.1 Run preprocessing scripts

1. Enter python environment.

2. Load the preprocessing commands (this will also load astro_functs.py as af):

	```python
	In [1]: from rat_preproc import *
	```

3. Create lists of fits files for specified cameras (must be done for each directory containing FITs files you'll be using):

	```python
	In [2]: ratlist( workdir = 'path/to/FITS/files/', cams = [0,1,2,3] )
	```

4. Select calibration frames you want to use:

	```python
	In [3]: ratdisp_calib( ftype=af.FLAT_NAME or af.BIAS_NAME, workdir='path/to/FITS/flats/', cams=[0,1,2,3], auto=True, amin=0.1, amax=0.8 )
	```
	
5. Select science frames you want to use:

	```python
	In [4]: ratdisp( workdir='path/to/FITS/files/', targetdir='path/to/new/FITS/files/', cams=[0,1,2,3], auto=True )
	```
	
6. Make master bias or flat frame:

	```python
	In [5]: mkmaster( af.BIAS_NAME or af.FLAT_NAME, bands='ALL', workdir='.', fmin=5 )
	```
	
More detailed instructions can be found in *reduction_instructions.rtf* or code comments in *rat_preproc.py* and *astro_functs.py*.

### 2.2 Run *ratautoproc.pro*

This IDL script should be run from the data directory (the directory above processed data and reduction folder)

```IDL
ratautoproc, datadir='raw/', redo=1
```
	
* *datadir* specifies where processed data is stored (should be in lower directory), will run all data in this directory  

* Will save to specified directory (imworkingdir) in pipeautoproc.par

* *redo* keyword overwrites previously reduction processed files

* Additional keywords:

	- start

	- stop

	- step

	- only (allows you to run particular steps if you don't want to run full reduction)

	- nocrclean (if set skips cosmic ray cleaning)	
	
* Runs these steps in this order unless specified:	

	- steps = ['prepare', 'flatten', 'makesky', 'skysub', 'crclean', 'astrometry', 'stack']
	
Description of each step:

1. Prepare (*autopipeprepare.pro*, *pipeprepare.pro*)

	- Header information manipulation and bias subtraction for images with master bias files
	
2. Flatten (*autopipeimflatten.pro*, *flatpipeproc.pro*)

	- Divides master flat for each filter
		
3. Makesky (*autopipemakesky.pro*, *skypipecombine.pro*)

	- Creates master sky by removing sources (using outlier rejection) and then median iterative sigma clipping for each pixel
		
4. Skysub (*autopipeskysub.pro*, *skypipeproc.pro*)

	- Subtracts sky, then subtracts median, then adds 1000
	
5. Crclean (*autopipecrcleanim.pro*, *pzap_perley.pro*)

	- Cosmic ray cleaning, details fuzzy
		
6. Astrometry (*autopipeastrometry.pro*, *vlt_autoastrometry.py*)

	- Fixes WCS coordinates by using pair-distance matching and asterism matching
	
7. Stack (*autopipestack.pro*)

	- Uses SWarp to stack images with same filter

3. Photometry
-------------

### 3.1 Run autoredux.py

Run this python scrip from the photometry folder.

Can automatically run full photometry reduction using *autoredux.py*  
Needs to run inside directory with coadd*.fits files:

```python
import autoredux
autoredux.autoredux()
```

*autoredux.py* runs the following programs in this order:

1. Identify point sources in individual frames after cropping (automatic crop

	- using weight files, can run manual crop using manualcrop keyword).

	- icoords.py

	- Output: initially cropped files (.crop.fits)

		- coords files for w/ RA and DEC identified by sextractor

2. Run photometry and calculate zero points

	- calcoff.py

3. Run association program to make a master star list

	- assoc.py

4. run final photometry

	- finalphot.py

5. create webpage

	- plotratir.py

The sandbox was meant as a place to test new code developement.  Note that any changes made to this directory will not be reflected in the repository.
