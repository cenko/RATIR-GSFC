RATIR-GSFC
==========

GSFC/UMd RATIR Pipeline Repository

**This is a development pipeline that is no longer being supported. 
  The most up-to-date version is now being maintained on 
  https://github.com/vickitoy/photometry_pipeline**

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

    - SExtractor (v2.19.5)

    - SWarp (v2.38.0)
    
    - Scamp (v2.0.1) (For installing with macports dependencies see http://geha.commons.yale.edu/resources/installing-photometry-tools-on-a-recent-mac/, may need to replace deprecated PLplot functions with find and replace)
    
    - Missfits (v2.8.0)

    - cdsclient package  

    Most can be installed using Macports.  

2. Run *startup.sh* in Unix shell to source code:

    ```bash
    source startup.sh
    ```

3. Alter *pipeautoproc.par* to point to full path of *autoastrometry.py* on individual computer (under code/reduction/ratauto)  

4. **For configuring for other instruments see README.txt under code/**

2. Reduction
------------

### 2.1 Run preprocessing scripts

1. Enter python environment.

2. Load the preprocessing commands (this will also load astro_functs.py as af):

    ```python
    In [1]: from preproc import *
    ```

3. Select calibration frames you want to use:

    ```python
    In [2]: calibration_dict = choose_calib( instrument, 
                                             ftype,
                                             workdir='path/to/FITS/flats/',
                                             cams=[0,1,2,3],
                                             auto=False,
                                             reject_sat=True,
                                             amin=0.2, amax=0.8,
                                             save_select=True )
    ```
    
    #### Input
    - *instrument*:
        - string name of instrument to process
        - instruments defined in specific_instruments.py (ex. 'ratir')  
    - *ftype*:
        - the type of calibration (af.FLAT_NAME or af.BIAS_NAME)
    - *workdir*:
        - defaults to current directory
    - *cams*:
        - selects cameras by number
    - *auto*:
        - function will automatically select frames based on median values
    - *reject_sat*:
        - frames are rejected if they contain **any** saturated pixels
    - *amin/amax*:
        - values are fractions of a detector's staturation value
        - only frames with median values in this range can be selected
    - *save_select*:
        - save python dictionary of selected frames to a python pickle if True
    
    #### Return
    - returns a python dictionary of selected frames.  this dictionary is used by *mkmaster*.

    **WARNING:** RATIR flats are often bad even when the median value is in the acceptable range.  Auto mode is only recommended for bias frame selection.
    
4. Select science frames you want to use:

    ```python
    In [3]: science_dict = choose_science( instrument,
                                           workdir='path/to/FITS/files/',
                                           targetdir='path/to/new/FITS/files/',
                                           cams=[0,1,2,3],
                                           auto=True,
                                           save_select=True.
                                           calibrate=False )
    ```

    #### Input
    - *instrument*:
        - string name of instrument to process
        - instruments defined in specific_instruments.py (ex. 'ratir')  
    - *workdir*:
        - defaults to current directory
    - *cams*:
        - selects cameras by number
    - *auto*:
        - function will select all science frames if True
    - *save_select*:
        - save python dictionary of selected frames to a python pickle if True
    - *calibrate*
        - uses calibration files (needs bias, dark, and flat) to display images with
          less noise
    
    #### Return
    - returns a python dictionary of selected frames.

    **WARNING:** When auto is True, all science frames are selected.  Since the telescope occasionally has tracking issues, it is recommended to check all frames.
    
5. Make master bias or flat frame:

    ```python
    In [4]: mkmaster( instrument, calibration_dict, ftype, fmin=5 )
    ```
    
    #### Input
    - *instrument*:
        - string name of instrument to process
        - instruments defined in specific_instruments.py (ex. 'ratir')  
    - *calibration_dict*:
        - python dictionary created by *choose_calib*
        - can also provide the file name of a pickled dictionary to re-load it from a previous session
    - *ftype*:
        - type of master calibration frame being created (af.FLAT_NAME or af.BIAS_NAME)
    - *fmin*:
        - the minimum number of calibration frames allowed for a given camera or band
    
More detailed instructions can be found in code comments in *preproc.py* and *specific_instruments.py*.

### 2.2 Run *ratautoproc.pro* OR run *autoproc.py*

Either script should be run from the data directory (the directory above processed data and reduction folder).  
(Now written in Python with same keywords. IDL version will not be maintained after 9/1/2015)


```IDL>
ratautoproc, datadir='raw/', imdir='reduced/', redo=1
```

or

```python>
autoproc(datadir='raw/', imdir='reduced/', redo=1)
```
    
* *datadir* specifies where processed data is stored (should be in lower directory), will run all data in this directory  

* Will save to specified directory (imworkingdir) in pipeautoproc.par (can specify with imdir)

* Need to move master flats and master biases into imworkingdir for bias subtraction and flat fielding

* *redo* keyword overwrites previously reduction processed files

* Additional keywords:

    - start

    - stop

    - step

    - only (allows you to run particular steps if you don't want to run full reduction)
    
    - nomastersky (if set does not create master sky and subtracts median of entire frame)

    - nocrclean (if set skips cosmic ray cleaning)  
    
    - quiet (mainly silent operation unless errors)
    
    - rmifiles (removes intermediate files)
    
* Runs these steps in this order unless specified:  

    - steps = ['prepare', 'flatten', 'makesky', 'skysub', 'crclean', 'astrometry', 'stack']
    
Description of each step:

1. Prepare (Python: *autoproc_steps.py*, *autoproc_depend.py*, IDL: *autopipeprepare.pro*, *pipeprepare.pro*)

    - Header information manipulation and bias subtraction for images with master bias files
    
2. Flatten (Python: *autoproc_steps.py*, *autoproc_depend.py*, IDL: *autopipeimflatten.pro*, *flatpipeproc.pro*)

    - Divides master flat for each filter
        
3. Makesky (Python: *autoproc_steps.py*, *autoproc_depend.py*, IDL: *autopipemakesky.pro*, *skypipecombine.pro*)

    - Creates master sky by removing sources (using outlier rejection) and then median iterative sigma clipping for each pixel
        
4. Skysub (Python: *autoproc_steps.py*, *autoproc_depend.py*, IDL: *autopipeskysub.pro*, *skypipeproc.pro*)

    - Subtracts sky, then subtracts median, then adds 1000
    
5. Crclean (Python: *autoproc_steps.py*, *autoproc_depend.py*, *cosmics.py*, IDL: *autopipecrcleanim.pro*, *pzap_perley.pro*)

    - Cosmic ray cleaning
        
6. Astrometry (Python: *autoproc_steps.py*, *autoproc_depend.py*, *vlt_autoastrometry.py*, IDL: *autopipeastrometry.pro*, *vlt_autoastrometry.py*)

    - Fixes WCS coordinates by using pair-distance matching and asterism matching (also runs Scamp for extra astrometry correction)
    
7. Stack (Python: *autoproc_steps.py*, *autoproc_depend.py*, IDL: *autopipestack.pro*)

    - Uses SWarp to stack images with same filter and calculates both relative and absolute zeropoint info (absolute only for coadded)

3. Photometry
-------------

### 3.1 Run autoredux.py

Move desired coadd files ('coadd????.fits', 'coadd???.weight.fits') into a separate photometry folder.
Run this python script from the photometry folder.

Can automatically run full photometry reduction using *autoredux.py*  
Needs to run inside directory with coadd*.fits files:

```python
from autoredux import *
autoredux()
```

*autoredux.py* runs the following programs in this order:

1.  Creates same sampling and crop for all files, creates multicolor image and using multicolor image to find all sources, then finds photometry
    of these sources for each (resampled) file using sextractor (with corrected zeropoint).  

    - photom.py

    - Output: aperture photometry

        - *.am (absolute magnitude) files for w/ RA and DEC identified by sextractor

2. create webpage

    - plotratir.py

The sandbox was meant as a place to test new code development.  Note that any changes made to this directory will not be reflected in the repository.
