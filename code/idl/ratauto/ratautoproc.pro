; NAME:
;   ratautoproc.pro                               (release version 2012.07.24)
;
; PURPOSE:
;   Fully-automated reduction of LRIS imaging and (2D) spectroscopy.
;
; CALLING SEQUENCE:
;   ratautoproc, [settings]
;
; INPUTS (all optional):
;   datadir      - Location of raw data (current directory if unspecified)
;   modestr      - Mode to process (imaging or spectroscopy)
;   camerastr    - Camera to process (red or blue)
;   chipstr      - Chip to process (left or right)
;   start        - Start with this step, skipping previous ones
;   stop         - End with this step, skipping subsequent ones
;   only         - Do only this step
;   step         -  (completely identical to only, takes precedence)
;   chip         - Process only one chip (left or right)
;   redo         - Repeat step(s), overwriting any existing files
;   nocrclean    - Do not zap cosmic rays
;
; ADDITIONAL OPTIONS:
;   If any of the following files are found in the directory where lrisautoproc
;   is run, it will change the default behavior.
;
;   lrisautoproc.par - Contains various defaults (mainly directory paths)
;   catalog.txt - Source positional catalog.  Use to auto-correct the object
;        name in the header based on position (independent of TARGNAME.)
;   brokenlamps.txt - List of arc lamps that are not functioning.
;   tracepos.txt - List of y-coordinates for user-specified trace extraction
;        of certain images.
;   seeing.txt - A file with only two numbers, the minimum and maximum seeing
;        experienced in arcseconds.
;   badflattargets.txt - List of fields (TARGNAMEs) not to use when 
;        constructing supersky flats.
;   (backgrounds.txt) - Background size to subtract for each imaging field, or 
;        a list of fields for which scalar background subtraction is needed.
;
;
;  (Options in parenthesis above are NOT YET ENABLED or have not been tested.)
;
;
; COMMENTS:
;   This code is meant to be fully automated, producing science-quality
;   output with a single command line with zero intervention.  Reduces
;   both LRIS-B and LRIS-R 2/3 data (LRIS-R1 will be enabled soon) in
;   either imaging or spectroscopy mode.
; 
;   The production of images is quite high-level and includes photometric 
;   calibration of the field (although the accuracy of this has not been
;   robustly tested).  Full spectroscopic extraction and reduction has recently
;   been added but has not been tested to the same extent as the imaging
;   pipeline.
;
;   The program works in a series of steps following standard CCD reduction
;   techniques, automatically recognizing calibration files, matching
;   calibrations together and science frames with appropriate calibrations,
;   etc.  Users looking for more control over the reductions can run each
;   step individually with the "step" command.
;
;   If the reduction is interrupted, the program will resume without
;   incident (although any failed steps may be repeated); each task is 
;   run independently of all others.
;   
;   The program tries very hard to correct for observer mistakes, such
;   as inconsistent windowing of the data.  But it's not perfect.
;   If problems occur, generally the easiest fix is to delete any 
;   offending files (images that are difficult to distinguish whether they 
;   are twilight flats can cause various problems, for example) and rerun.
;   More significant problems generally require direct modification of the
;   code, which is still in production.  See 'Troubleshooting' below for
;   more information.
;
;   Filenames for the input raw images are expected to be in 
;   [b|r]YYMMDD_NNNN.fits format for LRIS-B and LRIS-R, respectively.
;   They can be either in the working directory or in a different directory
;   specified by datadir.
; 
;   The code runs in a series of steps in the following order:
;
;   1. "Prepare" the data by converting the multi-header extension FITS files
;      to a standard single frame, correcting/flagging known pixel defects,
;      cropping unused area of the chip, bias-subtracting using the overscan,
;      and adding extra information to the header.  The conversion/bias 
;      correction employs a modified version of readmhdu.pro.  The X and Y axes
;      for spectroscopic frames are transposed for easier viewing of traces on 
;      standard monitors and to match 1D spectrum plots (wavelength=x-axis).
;      Output: p[b|r]*.fits (written to ./imredux/ and ./spredux/ by default.)
;
;   2. Create flat-fields.  The program searches through all science exposures
;      and determines what flat fields are needed, then attempts to make the
;      best possible flat for each case (generally, the order of preference is
;      supersky > twilight > dome > internal halogen, but which type of flat
;      is actually made depends on the number of frames available and their
;      quality.)  In the case of imaging flats, stars are identified and
;      masked out.
;      Output:  fp[b|r]*.fits
;
;   3. Flat-correct data.  Divide each image by the flatfield.  A more
;      refined cropping is also done at this stage, depending on the placement
;      of the image area during the run (which is variable.)
;   
;      (Fringe correction routines are not needed for LRIS since 2009, but
;      were available for earlier versions of this code and will be re-enabled
;      in the near future.)
;
;   4. Split the data into separate right and left chips.  (Because of the
;      physical gap and rotation between the LRIS chips, they must be treated
;      independently after this point.)
;      Output:  fp[b|r]*[l|r].fits
;
;   5. Removes cosmic rays, using the independent routines pzap.pro and
;      pzapspec.pro.  See those programs for more information.  
;      This can be a time-consuming process.
;      Output: zfp[b|r]*[l|r].fits
;
;   At this point different procedures are applied for imaging and spectroscopy.
;   For imaging:
;  
;   6. Solve astrometry of the field against the best available online catalog
;      (SDSS or USNO-B1.0), using the independent autoastrometry.py code.
;      (Also requires sextractor to be installed.)
;      One image is solved relative to a reference catalog, then the
;      remaining images are solved against the first.  NO distortion corrections
;      are currently applied.
;      Output: azfp[b|r]*[l|r].fits
;   
;   7. Produce photometric solution of the field.  The GSFC IDLastro aper.pro
;      program is used to measure photometry of point sources in each field
;      for comparison of different exposures (on the same target) to establish 
;      the relative zeropoints.  
;      Output: [object].cal.list and other text files
;  
;   8. Stack exposures.  A weighted, masked median is performed for each field.
;      Requires swarp to be installed and properly linked.
;      Output: coadd[object].[filter].[l|r|]fits
;
;
; 
; EXAMPLES:
;   1.  In a directory containing a full night worth of LRIS data, enter
;
;       IDL> lrisautoproc
;
;       This will automatically execute all of the steps above.
;       Depending on the quantity of data from your run and the speed of your
;       computer, the reduction of data from an entire night takes 20 minutes
;       up to about two hours.
;
;
;   2.  Fully reduce only spectroscopy, and only the right half of the red CCD:
;
;       IDL> lrisautoproc, mode='s', camera='red', chip='r'
;
;
;   3.  Run only the "prepare" step (bias subtraction/formatting), on data 
;       stored in a separate directory:
;
;       IDL> lrisautoproc, step='prepare', datadir='/scr3/user/lris/20120501/'
;
;
;   4.  For live-pipeline usage, I recommend placing all the data in a
;       subdirectory 'raw/', then running lrisautoproc one step down as 
;       follows:
;
;       IDL> lrisautoproc, datadir='raw/', chip='r'
;
;       This will save some processing time that would otherwise be spent on
;       e.g. cleaning cosmic rays from the left chip.  The left chip is
;       marginally useful for imaging (more stars for calibration) but of
;       no value for spectroscopy (except for rare cases of highly extended or
;       widely separated objects.)
;
;
;
;
;
;   INSTALLATION:
;
;   Create the subdirectory 'lrisauto' somewhere on your hard drive (probably
;   in your IDL directory), and unpack the contents of this installation file
;   there (e.g.,  tar -xvf lrisauto.tar.gz).  You will need to tell IDL about
;   the existence of this new directory by editing your idl_startup file
;   (which in turn is specified by .idlenv or in your .bashrc or .cshrc file)
;   to add the string ":+/path/to/lrisauto:+/path/to/lrisauto/dependencies/"
;   to the IDL_PATH argument.  The GSFC IDLastro routines must also be
;   installed (and visible within IDL_PATH); see http://idlastro.gsfc.nasa.gov/
;   for programs and instructions.
;
;   In order to process images, you will also need to have autoastrometry,
;   swarp and sextractor installed (this requirement will eventually be removed
;   via a simplification of the astrometric solver method).  If the latter two
;   cannot simply be called via "swarp" and "sex", you will need to edit the 
;   file lrisautoproc.par AND also edit the dependencies/autoastrometry.py file
;   to indicate the actual commands to call these routines in the global
;   variables at the top of the code.  The standard UNIX routine wget is also 
;   used to download star catalogs.
;
;   For a Caltech astronomy cluster matchine, it is recommended to simply point
;   your IDL path variable to the latest release version in Daniel Perley's 
;   home directory (see citsetup.txt for more info.)  This will ensure you 
;   always have the most up-to-date version.
;
;
;
;
;   TROUBLESHOOTING:
;
;   While this pipeline is designed to deal with all possible observing
;   circumstances (including many common mistakes), much more testing and 
;   development will be required before this ideal is fully reached.
;   Despite best efforts, the program may crash if it encounters an 
;   unanticipated situation or has problems accomplishing its goals and is
;   unable to proceed.  If you encounter problems, try e-mailing Daniel Perley
;   (dperley@astro.caltech.edu) for assistance, after checking the below.
;
;   If the pipeline crashes during operation:
;
;   * Did it report a 'file not found' error?  
;       If so, a subroutine, configuration file, or catalog it needed was not 
;       found at the expected location.  Generally this is an installation 
;       problem, but you can inspect the paths in lrisautoproc.par to make sure
;       that the directories and files specified there actually exist.  
;  
;   * Did an external call (autoastrometry, swarp, sex) fail?
;       Check that you have specified the paths correctly and actually have the
;       required software installed for astrometry/coadding (see "installation", 
;       above.)  Alternatively if you don't care about reducing images, specify
;       mode='s' at runtime to process only the spectroscopy, which uses only
;       IDL routines and requires no external packages (except IDLastro).
;   
;   * Did it encounter some other error while processing a file?  
;       Check which file the program was working on when it crashed.  If the
;       file is not essential, try deleting it and re-running the pipeline
;       starting with the current step (or delete ALL of the file's precursors
;       with the same file number and rerun the pipeline.)  Please report the
;       error (and perhaps include a copy of the deleted file) so the problem
;       can be identified and fixed in future releases.
;   
;
;   If processing completed, but the results are problematic:
;
;   * Did it return without actually processing any data? 
;       If this is the first run of the pipeline on a night, make sure that
;       it is in the current working directory or that you have correctly 
;       pointed to the directory containing raw data with the "datadir" 
;       keyword.  If you are re-doing a step, it will not overwrite existing
;       files by default.  Set the /redo keyword to overwrite files from a
;       previously-attempted step (be sure to set "start" or "step" unless
;       you want to restart the pipeline from the beginning in this case.)
;       Finally, the pipeline will not construct flats it doesn't need, so
;       if no science exposures exist yet it will not undertake steps beyond
;       'prepare'.
;
;   * Do the reduced images look bad?
;       The most common source of this is bad flat-fielding.  Especially if
;       twilight flats were not taken, the routine may have attempted to
;       construct a supersky flat, which is difficult to do with Keck due
;       to bright stellar halos and (especially) galaxies.  There is a way 
;       to specify fields not to use in supersky flats (badflattargets.txt)
;       Fields with bright galaxies are also susceptible to sky over-
;       subtraction.  Sky subtraction can be turned off for these fields
;       (currently manually by modifying the swarp file, but a semi-automated
;       method will be coming soon.)  Other times, a bad image can result from
;       unavoidable instrument issues such as scattered light from a bright
;       star outside the field.
;
;   * Is the cosmic ray removal malfunctioning?
;       In some circumstances bright nebular emission lines can be zapped, 
;       especially in conditions of good seeing and binned data.  The cosmic
;       ray cleaner settings can be modified internally if this is a problem.
;       Check both the 1D reduction and the 2D subtracted frames to be sure
;       cosmic ray cleaning is not causing problems if nebular emission lines
;       are important to your science.
;       The imaging cosmic ray clean threshold is set relatively high to avoid
;       zapping stellar cores (cosmic rays are removed by median combining
;       later anyway so it's not critical), but if you see this occurring
;       please report it immediately.
;
;   * Is the zeropoint calibration of the image unreasonable?
;       The imaging zeropoint-calibration has not been well-tested yet, so 
;       please report this problem so this part of the pipeline can be 
;       improved.
;
;   * Is the wavelength calibration of a spectrum off?
;       Wavelength calibration is probably the most fragile part of the
;       pipeline as there is not yet a careful check on the accuracy of the
;       wavelength solution, and gratings the pipeline hasn't been tested on
;       are not guaranteed to produce good results.  If the wavelength
;       calibration looks suspect, load the relevant *.wavsol.ps plot (in 
;       addition, the f*.wavsol.ps plots contain the master arc solutions)
;       and inspect.  If the wavelength solution is bad, check the input
;       arc file for problems (perhaps it wasn't warmed up?).  If one of the
;       arc lamps was broken (typically the Ne lamp), you can inform the
;       pipeline of this by placing a file containing the abbreviation of 
;       the broken arc element in brokenarcs.txt (in your main reduction 
;       directory).  If all else fails, you can try to solve the arc yourself
;       and place the polynomial terms in the format "arcname.sol".
;       Note that it is sometimes possible that the grating got 'stuck' and
;       misreported its position, causing the arc and a file to look like
;       they are the same setup (in the header) but actually have a very
;       different wavelength solution.  If the actual central wavelength is
;       more than 200 Angstroms from the keyword value, this will cause a
;       failure in the pipeline unless the WAVELEN keyword in the affected
;       files is edited to the correct value.
;
;   * Were some files skipped?
;       If the pipeline encountered a non-fatal problem processing an
;       individual image or spectrum (such as a failed wavelength calibration
;       or inability to flatfield) then it will simply not process that file 
;       any further.  For most such cases a summary of the problems will be
;       printed out at the end of processing.  If a file is not being
;       processed and you do not see it in the final summary, you can simply
;       rerun the pipeline (without deleting any files and without setting the
;       redo flag) and it will try to repeat any failed steps of this nature. 
;
;   * Did it extract the wrong object, or no object at all?
;       It can't always tell which source you're interested in - it guesses the
;       brightest one in the appropriate part of the chip, but sometimes that
;       is not a valid assumption (and this can cause particular problems if
;       different objects are "brightest" on the red and blue sides!)  If the 
;       slit was moved away from the slit-B position, the source could also have
;       moved out of the pipeline's default search radius.  To tell the
;       program which object you care about, create the file "tracepos.txt"
;       in the main working directory and, for each file you want to modify the
;       tracing settings for, include a line of the following format:
;         szfpr120101_0100r.fits   30-40
;       That is, the filename of the sky-subtracted 2D spectrum and the zone to
;       restrict the trace search.  This refers to the location in the center
;       of the CCD.  Note that you have to do this for every exposure of an
;       object.  The units are in image pixels from the bottom of the chip.
;
;   * Does the flux go to zero (or even negative) for much of the trace?
;       If the object is faint or has a neighbor that is marginally resolved at
;       some wavelengths, the tracing may have failed.  You can try to restrict
;       the trace zone with tracepos.txt (above). If this doesn't work,
;       automated extraction will have to wait for further releases of the
;       pipeline.
;
;   * Does the output spectrum have many nonexistent absorption lines?
;       A common problem for observations of distant transients arises when an
;       offset host galaxy is contained within the background aperture and its 
;       emission lines are subtracted from the transient, creating false 
;       absorption lines or absorption+emission features.  Eventually the 
;       background will be chosen automatically to exclude neighboring objects,
;       but in the meantime you can manually specify the background region 
;       using tracepos.txt (see "Did it extract the wrong object?", above) by 
;       adding the minimum and maximum *distance* after the tracing region:
;         szfpr120101_0100r.fits   30-40   15-23
;       This will use a background region centered 19 pixels from the trace 
;       center with a diameter of 8 pixels.
;       (Another possiblility is that cosmic rays in the background aperture
;       missed by the cosmic-ray cleaner.  Please report such 2D images so the
;       cosmic ray settings can be fine-tuned to avoid this situation.)
;
;   The risk of encountering pipeline issues can be significantly reduced
;       by following standard observing procedures, below.  
;
;
;
;   OBSERVING TIPS:
;
;   Spectroscopy -
;
;   Binning:  Do not bin the red side by more than 2x1 for low-resolution 
;      spectroscopy.  Further binning causes severe saturation issues on bright
;      arc lines needed for automated wavelength calibration, and makes 
;      cosmic-ray identification more difficult.  If seeing conditions are
;      very good, spatial binning is not recommended.  In either case, be sure
;      to get calibrations in the same binning setting(s) as your science data
;      for each setup.
;
;   Arcs:  Please turn the mercury lamp on during red-side arc exposures (in
;      addition to the standard neon and argon lamps), and during setup tune 
;      the central wavelength to be sure it appears on the red side.  Having 
;      this line helps the wavelength solution across the D560 dichroic.
;
;   Flats:  Halogen flats are very poor: there is signficant spatial banding
;      and a large overall gradient across the y-axis, and significant
;      scattered red light reaches the blue camera.  The intensity of the lamp
;      is also not constant between exposures.  Only in desparate straits
;      should the halogen lamps be relied upon.  Spectroscopic dome flats are
;      much superior:  turn on the spectral dome lamp and get flats as you 
;      would for imaging.  Note that the UV cannot be pixel-calibrated with 
;      either the halogen OR with dome flats - the lamp is not hot enough. 
;      (In principle twilight spectroscopic flats can accomplish this, and 
;      while not supported by the pipeline yet, this is recommended if precise
;      flat-fielding is essential to your UV science.)
;
;   Standards:  The list of recognized spectrophotometric flux-calibration 
;      standards can be viewed in the 'refdata' subdirectory as the file
;      specstandards.dat, or online at:
;      http://astro.caltech.edu/~dperley/programs/lrisauto/specstandards.txt
;      Standards not on this list will not be recognized as standards (if you
;      have a favorite standards not on this list, please suggest it.)  Ideally
;      the same standard should be observed twice at significantly different 
;      airmasses in each setup.
;      It is recommended to take arcs alongside observations of bright 
;      standards where a short integration time is employed.  This will ensure
;      the standard has a good wavelength solution (no sky-line refinement is
;      possible for the shortest exposures: those less than about 10 seconds in
;      dark conditions, more in bright conditions or twilight).  This will help 
;      ensure accurate flux-calibration and (once implemented) telluric-line
;      removal.
;
;
;   Imaging -
;
;   Flats:  In general, twilight flats produce the most reliable results, so
;      try to get them if possible.  Observing in the U-band REQUIRES twilight
;      flats since there are not enough dark-sky counts and the dome lamp is
;      not hot enough to emit in the UV.  (Twilight flats do seem to create
;      artifacts with the glue defects in some filters, however.)
;
;   Standards:  If you are imaging any uncalibrated fields and conditions are 
;      photometric, try to get at least two measurements of 
;      Landolt photometric standard fields (preferably two observations of 
;      the same standard field) per filter at significantly different 
;      airmasses.  This can be done during twilight, even as part of getting 
;      sky flats.  (Landolt standards are better than using SDSS to determine
;      the telescope zeropoints because there is minimal contamination from 
;      nearby objects in short exposures of bright standards.)
;
;   Science fields:  Always specify that you want the source placed at the 
;      "LRIS-B imaging position" to avoid the chip gap.  A test exposure is
;      usually not necessary if the operator is trustworthy, but a short
;      exposure (< ~30s) will help significantly with calibration later if a 
;      moderately deep (e.g. SDSS) calibration of the field is not available.
;      When dithering, try to make sure that every image in the
;      dither sequence is at a different vertical position (as seen on the
;      live DS9 display) by e.g. dithering 1" vertically to accompany each
;      otherwise horizontal dither. LRIS-R3 has some bad columns (horizontal
;      lines in the viewer orientation) which one should avoid keeping over 
;      the same horizontal location (helpful, since these columns actually 
;      fluctuate in intensity and do not always flat-field out.)
;
;   Organization -
; 
;   Use save_state and restore_state (NOT the GUI) to store your LRIS 
;      configurations.  This includes the CCD parameters and will ensure that
;      you do not forget to chanage the binning and windowing which are NOT
;      included in the XLRIS save/load options.
;
;   To check what calibrations (etc.) you have or need, a useful tool is the
;      lriscat.py code, available at: 
;         http://www.astro.caltech.edu/~dperley/programs/lriscat.py
;      Run this at the command-line (python lriscat.py *.fits) to produce an
;      immediate catalog table containing all relevant LRIS settings.
;
;
;
; -----------------------------------------------------------------------------





; still to do:
; archival response curves and airmass adjustment
; telluric correction
; gain correction
; simplify matching by creating readlrisheaders.pro that creates a structure array
; check that USNO catalogs are loaded OK
; check for looking for "nickel" catalogs
; need an option for skipping spectroscopic flats entirely
; use system variables to make life easier
; faster CR cleaner for blue (and old red): done?
; old red side re-enable
; better text output; error/warning log written to disk
; background subtraction - if there is a neighbor it will be messed up.
; consistent aperture size for repeated observations?  Or, stack in 2D space?
; photometry table - need to be able to update on reruns, and deal with un-photometry-able fields
; automated halo removal for bright stars?
; option to make flats that aren't needed if files exist (beginning of night, compare dome vs. sky, etc.)
; during flux calibration:  re-normalize consistently for proper airmass interpolation/extrapolation.
; split r and l into multi-extension files in prepare and enable ccd/ files to recognize this.
; not sure if bpm is working, still a bad line on red side (when binned, anyway)
;  this appears to be due to movement of the weakened columns.  Make sure this isn't a crop problem and 
;  if not try to do something to detect this...
; interpolate standards in log space
; need better sky subtraction for twilight blue spectra (try a poly model with lines removed)
; future options
;  gratings=
;  filters=
;  files=
; sky subtraction of galaxies still causes issues, see 120427_0114
; still some tracing issues, see 120427_0111
; matching of apertures for multiple reobservations and red/blue...?
; simple median 'counts' from lrisprepare not appropriate for some gratings/dichroics
; need a new method for the slit-averaged flats - maybe linear extrapolation fit
; specify files=100-102, or filter='B', or grating='400/8000', etc. - partially done
; standard quality - red vs. blue, actual line recognition

; probably want a global arc solution, rather than each exposure (more robust) - done
; need sky-line refined wavelength solution - done
; deal with flats taken in wrong binning mode - done
; recognize mira files - done
; need to recognize saturated standards - done
; cleanup of autoastrometry files - done


; -------------

pro autolrismakesky, chip=chip, camera=camera, outpipevar=outpipevar, inpipevar=inpipevar

	;Setup pipeline variables that carry throughout the pipeline
	if keyword_set(inpipevar) then begin
		pipevar = inpipevar
		print, 'Using provided pipevar'
	endif else begin
		pipevar = {autoastrocommand:'autoastrometry' , sexcommand:'sex' , swarpcommand:'swarp' , $
					datadir:'' , imworkingdir:'' , overwrite:0 , modestr:'',$
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:'' }
	endelse
 
    prefchar = '2'

    wildcharsky = '?????????????????_sky_?'
    files = findfile(pipevar.imworkingdir+'fp'+prefchar+wildcharsky+'.fits')
    isdomeflat = intarr(n_elements(files))
    isskyflat  = intarr(n_elements(files))
    isscience  = intarr(n_elements(files))
    filters =  strarr(n_elements(files))  
    dichs =  strarr(n_elements(files))    
    wins =  strarr(n_elements(files))     
    binnings = strarr(n_elements(files))  
    counts = fltarr(n_elements(files))    
    targets = strarr(n_elements(files))   
    exps = fltarr(n_elements(files))      
    for f = 0, n_elements(files)-1 do begin
       if files[f] eq '' then continue
       h = headfits(files[f], /silent)

       filters[f] = clip(sxpar(h, 'FILTER'))
       targets[f] = clip(sxpar(h,'TARGNAME'))
       isscience[f]=1           ;all images are sky frames here
       isskyflat[f]=1           ;all images are sky frames here
    endfor
    sciindex = where(isscience, ct)
    if ct eq 0 then return
    scifilterlist  = unique(filters[sciindex])

    for i = 0, n_elements(scifilterlist)-1 do begin
       filt = scifilterlist[i]
       flatsuccess = 0
       ; Regular sky flat.
       if flatsuccess eq 0 then begin
          skyflats = where(isskyflat and filters eq filt, ctsky)
          outflatname = pipevar.imworkingdir+'sky-'+filt+'.fits'
          if ctsky ge 2 then begin
             print, filt, '-band sky flats.'
             if file_test(outflatname) and pipevar.overwrite eq 0 then continue ; flat exists already
             skycombine, files[skyflats], outflatname, /removeobjects, type='sky'
             flatsuccess = 1
          endif
       endif

       if flatsuccess eq 0 then begin
          print, 'Unable to produce a flat field for this setting: ', filt, '-band / ', 'D'+dich,' '+binstr
          noproc = where(isscience and filters eq filt and dichs eq dich and binnings eq bin, ctnoproc)
          print, 'Will not be able to further process ', clip(ctnoproc), ' images without a flat from another source:'
          for j = 0, ctnoproc-1 do print, '   ', files[noproc[j]], exps[noproc[j]], '  ', targets[noproc[j]]
          
       endif

    endfor
    
	outpipevar = pipevar
end

; -------------------------
pro autolrisimflatten, chip=chip, camera=camera, outpipevar=outpipevar, inpipevar=inpipevar

	;Setup pipeline variables that carry throughout the pipeline
	if keyword_set(inpipevar) then begin
		pipevar = inpipevar
		print, 'Using provided pipevar'
	endif else begin
		pipevar = {autoastrocommand:'autoastrometry' , sexcommand:'sex' , swarpcommand:'swarp' , $
					datadir:'' , imworkingdir:'' , overwrite:0 , modestr:'',$
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:'' }
	endelse

   prefchar = '2'
   fchip = strmid(chip,0,1)

   files = findfile(pipevar.imworkingdir+'p'+prefchar+pipevar.wildchar+fchip+'.fits')
   ffiles = findfile(pipevar.imworkingdir+'fp'+prefchar+pipevar.wildchar+fchip+'.fits')
   if n_elements(files) eq 1 and files[0] eq '' then return
   flats = findfile(pipevar.imworkingdir+'*flat*.fits')
   flatfilts = strarr(n_elements(flats))
   flatchips = strarr(n_elements(flats))
   flatdichs = strarr(n_elements(flats))
   flatwins = strarr(n_elements(flats))
   flatbins = strarr(n_elements(flats))
   if flats[0] ne '' then begin
     for f = 0, n_elements(flats)-1 do begin
       h = headfits(flats[f], /silent)
       filter = clip(sxpar(h, 'FILTER'))
       flatfilts[f] = filter
     endfor
   endif else begin
     nopflats = 1
     ; Don't actually raise the error until we know flats were necessary.
  endelse
   
    for f = 0, n_elements(files)-1 do begin
      if files[f] eq '' then continue
      outfile = fileappend(files[f], 'f')
      match = where(outfile eq ffiles, ct) ; check if output file exists
      if ct eq 0 or pipevar.overwrite then begin               
         h = headfits(files[f], /silent)
         counts = sxpar(h, 'SKYCTS') > sxpar(h,'COUNTS')
         exptime = sxpar(h, 'ELAPTIME')
         azimuth = sxpar(h,'AZ')
         elevation = sxpar(h,'EL')
         domeazimuth = sxpar(h,'DOMEPOSN')
         target = sxpar(h,'TARGNAME')
         filter = clip(sxpar(h, 'FILTER'))
         binning='1'
         flatfileno = where(flatfilts eq filter, ct) ;all we care about is the filter for RATIR

         if ct eq 0 then begin
            flatfilenoothbin = where(flatfilts eq filter,ctwobin) ;all we care about is the filter for RATIR
            ; Currently doesn't support flat-fielding with a different binning until I can be sure the case where both windowing and binning
            ; are wrong are dealt with
            if ctwobin gt 0 then binstr = ', '+repstr(binning,',','x') else binstr = ''
            print, 'Flat field not found for ', removepath(files[f]), ' (', filter+', '+'D'+binstr,')' ; , ' / ', win
            pipevar.flatfail = pipevar.flatfail +' '+ files[f]
            continue
         endif
         flatfile = flats[flatfileno[0]]
         print, 'Flattening ', removepath(files[f]), ' using ', removepath(flatfile)
         flatproc, files[f], flatfile, flatminval=0.3, crop='auto'  ; add autocropping of low-signal zones
      endif
    endfor
    
	outpipevar = pipevar
 end

; -------------------------
pro autolrisskysub, chip=chip, camera=camera, outpipevar=outpipevar, inpipevar=inpipevar

	;Setup pipeline variables that carry throughout the pipeline
	if keyword_set(inpipevar) then begin
		pipevar = inpipevar
		print, 'Using provided pipevar'
	endif else begin
		pipevar = {autoastrocommand:'autoastrometry' , sexcommand:'sex' , swarpcommand:'swarp' , $
					datadir:'' , imworkingdir:'' , overwrite:0 , modestr:'',$
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:'' }
	endelse
	
   prefchar = '2'
   fchip = strmid(chip,0,1)
   wildcharimg = '?????????????????_img_?'

   files = findfile(pipevar.imworkingdir+'fp'+prefchar+wildcharimg+'.fits')
   sfiles = findfile(pipevar.imworkingdir+'sfp'+prefchar+wildcharimg+'.fits')
   if n_elements(files) eq 1 and files[0] eq '' then return
   skys = findfile(pipevar.imworkingdir+'*sky-*.fits')
   skyfilts = strarr(n_elements(skys))
   skycams = strarr(n_elements(skys))
   skychips = strarr(n_elements(skys))
   skydichs = strarr(n_elements(skys))
   skywins = strarr(n_elements(skys))
   skybins = strarr(n_elements(skys))

   if skys[0] ne '' then begin
     for f = 0, n_elements(skys)-1 do begin
       h = headfits(skys[f], /silent)
       filter = clip(sxpar(h, 'FILTER'))
       skyfilts[f] = filter
     endfor
   endif else begin
     nopskys = 1
     ; Don't actually raise the error until we know skys were necessary.
   endelse
   for f = 0, n_elements(files)-1 do begin
      if files[f] eq '' then continue
      outfile = fileappend(files[f], 's')
      match = where(outfile eq sfiles, ct) ; check if output file exists
      if ct eq 0 or pipevar.overwrite then begin               
         h = headfits(files[f], /silent)
         camera = sxpar(h,'WAVELENG')
         camera = strcompress(camera,/remove_all)
         counts = sxpar(h, 'SKYCTS') > sxpar(h,'COUNTS')
         exptime = sxpar(h, 'ELAPTIME')
         azimuth = sxpar(h,'AZ')
         elevation = sxpar(h,'EL')
         domeazimuth = sxpar(h,'DOMEPOSN')
         target = sxpar(h,'TARGNAME')

         filter = clip(sxpar(h, 'FILTER'))
         binning='1'
         skyfileno = where(skyfilts eq filter, ct) ;all we care about is the filter for RATIR

         if ct eq 0 then begin
            skyfilenoothbin = where(skyfilts eq filter,ctwobin) ;all we care about is the filter for RATIR
            ; Currently doesn't support flat-fielding with a different binning until I can be sure the case where both windowing and binning
            ; are wrong are dealt with
            if ctwobin gt 0 then binstr = ', '+repstr(binning,',','x') else binstr = ''
            print, 'Sky field not found for ';, removepath(files[f]), ' (', filter+', '+'D'+dich+ binstr,')' ; , ' / ', win
stop
            pipevar.flatfail = pipevar.flatfail +' '+ files[f]
            continue
         endif
         skyfile = skys[skyfileno[0]]
         print, 'Sky Subtracting ', removepath(files[f]), ' using ', removepath(skyfile)
         if strcmp(camera,'OPT') then skyprocopt, files[f], skyfile else skyproc, files[f], skyfile

      endif
      skipskysub:
    endfor
    
	outpipevar = pipevar

end

; -------------------------
pro autolriscrcleanim, chip=chip, camera=camera, outpipevar=outpipevar, inpipevar=inpipevar

	;Setup pipeline variables that carry throughout the pipeline
	if keyword_set(inpipevar) then begin
		pipevar = inpipevar
		print, 'Using provided pipevar'
	endif else begin
		pipevar = {autoastrocommand:'autoastrometry' , sexcommand:'sex' , swarpcommand:'swarp' , $
					datadir:'' , imworkingdir:'' , overwrite:0 , modestr:'',$
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:'' }
	endelse

   prefchar = '2'
   chipchar = strmid(chip,0,1)

   wildcharsky = '?????????????????_img_?'
   files = choosefiles(prefchar+pipevar.wildchar+'.fits',pipevar.imworkingdir+'isfp',pipevar.imworkingdir+'sfp')
   if pipevar.overwrite eq 0 then files = unmatched(files,'z')
   for f = 0, n_elements(files)-1 do begin
      if files[f] eq '' then continue
      h = headfits(files[f])
      counts = sxpar(h,'COUNTS')
      exptime = sxpar(h,'ELAPTIME')
      target = sxpar(h,'TARGNAME')
      if counts gt 40000. then continue
      if counts gt 30000. and exptime lt 10. then continue
      if counts gt 20000. and exptime gt 5. and exptime lt 10. then continue
      if target eq '' then continue   ; blank target is certainly not a science object.
      if strpos(target,'flat') ge 0 then continue
      if strpos(target,'sky') ge 0 then continue
      if strpos(target,'twilight') ge 0 then continue

      print, 'Cleaning cosmic rays from ', removepath(files[f])
      zeal = 0.85
      if camera eq 'blue' then zeal = 0.75
      if (camera eq 'red' and pipevar.lrisversion eq 1) then usamp=1 else usamp=0
      if camera eq 'red' and pipevar.lrisversion gt 1 then zeal = 0.6 ;evidently should be the new default
      if n_elements(setzeal) eq 1 then zeal = setzeal
      slashpos = strpos(files[f],'/')
      dir = strmid(files[f],0,slashpos+1)
      outname = dir + 'z' + strmid(files[f],slashpos+1)
      if file_test(outname) and pipevar.overwrite eq 0 then continue
      pzap_perley, files[f], /weight, zeal=zeal, usamp=usamp, /quiet
   endfor

	outpipevar = pipevar

end

; -------------------------
pro autolrisastrometry, camera=camera, chip=chip, outpipevar=outpipevar, inpipevar=inpipevar

	;Setup pipeline variables that carry throughout the pipeline
	if keyword_set(inpipevar) then begin
		pipevar = inpipevar
		print, 'Using provided pipevar'
	endif else begin
		pipevar = {autoastrocommand:'autoastrometry' , sexcommand:'sex' , swarpcommand:'swarp' , $
					datadir:'' , imworkingdir:'' , overwrite:0 , modestr:'',$
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:'' }
	endelse
	
    achip = strmid(chip,0,1)
    prefchar = '2'
    wildcharimg = '?????????????????_img_?'
    zffiles = choosefiles(prefchar+wildcharimg+'.fits',pipevar.imworkingdir+'zisfp',pipevar.imworkingdir+'zsfp',$
                                                          pipevar.imworkingdir+'isfp',pipevar.imworkingdir+'sfp')
    
    filetargets = strarr(n_elements(zffiles))
    fileexposures = strarr(n_elements(zffiles))
    filecounts = fltarr(n_elements(zffiles))
    filefilt = strarr(n_elements(zffiles))

    for f = 0, n_elements(zffiles)-1 do begin
       if zffiles[f] eq '' then continue
       h = headfits(zffiles[f], /silent)
       filetargets[f] = repstr(repstr(strtrim(sxpar(h,'TARGNAME'),2),' ', '_'),'/','_')  ; spaces cause barfing in filenames
       filefilt[f] = sxpar(h,'FILTER')
       fileexposures[f] = string(sxpar(h,'EXPTIME'))
       filecounts[f] = sxpar(h,'COUNTS')
    endfor

    ; Make a reference catalog using a representative image out of an image block (several images of the same field)

    targets = unique(filetargets)
    filters = unique(filefilt)

    for t = 0, n_elements(targets)-1 do begin
       for f = 0, n_elements(filters)-1 do begin

          targetfiles = where(zffiles eq targets[t])
          refcatfile = strcompress(pipevar.imworkingdir+targets[t]+'.'+filters[f]+'.cat',/remove_all)

          if file_test(refcatfile) and pipevar.overwrite eq 0 then continue
          thistarget = where(filetargets eq targets[t] and filefilt eq filters[f])
          maxexp = max(fileexposures(thistarget))
          minctrate = min(filecounts[thistarget]/fileexposures[thistarget])
          
          if maxexp lt 5 or minctrate gt 5000 then begin
             print, targets[t], ' is a shallow field (standard or sky): not making a reference catalog.'
             continue
          endif
          
          imagesthistarg = zffiles[thistarget]
          refimagename = imagesthistarg[n_elements(imagesthistarg)/2]

          print
          print, 'Making reference catalog for ', targets[t], ' using ', refimagename
          print, pipevar.autoastrocommand+' '+refimagename;+' -upa 2 -q'
          spawn, pipevar.autoastrocommand+' '+refimagename;+' -upa 2 -q'

          outfile = fileappend(refimagename,'a')
          if file_test(outfile) eq 0 then begin
             catastrofail = catastrofail +' '+ refimagename
             print, 'WARNING - astrometry on the reference image was unsuccessful!'
             print
          endif else begin
             print, pipevar.autoastrocommand+' '+outfile+' -n '+refcatfile + ' -x 55000 -q'
             spawn, pipevar.autoastrocommand+' '+outfile+' -n '+refcatfile + ' -x 55000 -q'
             print
          endelse
      ;   (respecifying px and pa should no longer be necessary with the new splitlrisred and splitlris blue which add astrometry.)
       endfor
    endfor

    if pipevar.overwrite eq 0 then zffiles = unmatched(zffiles,'a')  ; catalogs should always be the same; only do this now

    ; Use the reference catalog to do a more precise relative astrometric solution
    for f = 0, n_elements(zffiles)-1 do begin
       if zffiles[f] eq '' then continue
       outfile = fileappend(zffiles[f],'a')
       if file_test(outfile) and pipevar.overwrite eq 0 then continue
       h = headfits(zffiles[f], /silent)
       pa = strtrim(string((float(sxpar(h, 'ROTPOSN'))+0) MOD 360),2)
       exptime = sxpar(h,'ELAPTIME')
       counts = sxpar(h,'COUNTS')
       filt = sxpar(h,'FILTER')
       if counts/exptime gt 5000. or exptime lt 5 then begin 
         ; only do catalog astrometry for short exposures (no catalog).
          if counts gt 30000 and exptime lt 10 then begin
             print, zffiles[f], ' is in bright twilight; not solving astrometry.'
             continue
          endif
          print, zffiles[f], ' is a twilight/standard frame, solving astrometry directly against a catalog.'
          print, pipevar.autoastrocommand+' '+zffiles[f];+' -upa 2 -q'
          spawn, pipevar.autoastrocommand+' '+zffiles[f];+' -upa 2 -q'
          print
          if file_test(outfile) eq 0 then pipevar.fullastrofail = pipevar.fullastrofail +' '+ zffiles[f] 
       endif else begin
         ; use the short exposure first, but fall back on direct catalog.
          targname = repstr(strtrim(sxpar(h,'TARGNAME'),2),' ', '_') ; spaces cause barfing in filenames
          if targname eq '' then continue
          if strpos(targname,'flat') ge 0 then continue
          refcatfile = strcompress(pipevar.imworkingdir+targname+'.'+filt+'.cat',/remove_all)
          
          if file_test(refcatfile) then begin
             print, targname
             print, pipevar.autoastrocommand+' '+zffiles[f]+' -c '+refcatfile;+' -upa 2' + ' -x 55000 -q' ;-upa 0.1
             spawn, pipevar.autoastrocommand+' '+zffiles[f]+' -c '+refcatfile;+' -upa 2' + ' -x 55000 -q' ;-upa 0.1
          endif else begin
             print, 'No reference catalog '+refcatfile+' exists for this field.'
          endelse
          if file_test(outfile) eq 0 then begin
             print, 'Refined astrometry of ', zffiles[f], ' was not successful.  Trying direct astrometry:'
             print, pipevar.autoastrocommand+' '+zffiles[f];+' -upa 2 -q'
             spawn, pipevar.autoastrocommand+' '+zffiles[f];+' -upa 2 -q'

             if file_test(outfile) then pipevar.relastrofail  = pipevar.relastrofail +' '+ zffiles[f] $
             else pipevar.fullastrofail = pipevar.fullastrofail +' '+ zffiles[f]
          endif
       endelse

    endfor

	outpipevar = pipevar
	
end



; ------------------------

pro autolrisphotometry, camera=camera, chip=chip, outpipevar=outpipevar, inpipevar=inpipevar

	;Setup pipeline variables that carry throughout the pipeline
	if keyword_set(inpipevar) then begin
		pipevar = inpipevar
		print, 'Using provided pipevar'
	endif else begin
		pipevar = {autoastrocommand:'autoastrometry' , sexcommand:'sex' , swarpcommand:'swarp' , $
					datadir:'' , imworkingdir:'' , overwrite:0 , modestr:'',$
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:'' }
	endelse

  ; need to allow doing this on only one chip.
   lrisautofieldphot,blue=bl,red=re,chip=ch

	outpipevar = pipevar
	
end

; -------------------------
pro autolrisstack, camera=camera, chip=chip, outpipevar=outpipevar, inpipevar=inpipevar

	;Setup pipeline variables that carry throughout the pipeline
	if keyword_set(inpipevar) then begin
		pipevar = inpipevar
		print, 'Using provided pipevar'
	endif else begin
		pipevar = {autoastrocommand:'autoastrometry' , sexcommand:'sex' , swarpcommand:'swarp' , $
					datadir:'' , imworkingdir:'' , overwrite:0 , modestr:'',$
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:'' }
	endelse

  if file_test('default.swarp') eq 0 then $
     spawn, pipevar.swarpcommand+' -d > default.swarp'

  if n_elements(chip) eq '' then chip = '*'
  chipchar = strmid(chip,0,1)
  prefchar = '2'
  wildcharimg = '?????????????????_img_?'
  azffiles = findfile(pipevar.imworkingdir+'a*'+prefchar+wildcharimg+'.fits')
  realfiles = where(azffiles ne '', ct)
   if ct eq 0 then return ; can't stack if there are no astrometry files
   if ct ge 1 then azffiles = azffiles[realfiles]
   ;chipchar = 'lr'  ; this is a hack to combine both sides

   camver = camera
   if camver eq 'red' then camver = camver + clip(pipevar.lrisversion)
   if file_test(pipevar.imworkingdir+'autophotsummaryflux.txt') then begin
      readcol,pipevar.imworkingdir+'autophotsummaryflux.txt',pfile,pexp,pfilt,pair,dum,pdmag,pfluxratio,pseeing,format='a,i,a,f,a,f,f,f,f',/silent
      photodata = replicate({filename:'',dmag:0.,fluxratio:0.,seeing:0.},n_elements(pfile))
      photodata.filename = pfile
      photodata.dmag = pdmag
      photodata.fluxratio = pfluxratio
      photodata.seeing = pseeing
   endif

   filetargets = strarr(n_elements(azffiles))
   fileexposures = fltarr(n_elements(azffiles))
   filefilters = strarr(n_elements(azffiles))
   filesatval = fltarr(n_elements(azffiles))
   fileskyval = fltarr(n_elements(azffiles))
   fileairval = fltarr(n_elements(azffiles))
   fileseeingpix = fltarr(n_elements(azffiles))
   filefluxratio = fltarr(n_elements(azffiles))
   datestr = ''

   for f = 0, n_elements(azffiles)-1 do begin
     h = headfits(azffiles[f], /silent)
     if f eq 0 then datestr = sxpar(h,'DATE')
     filetargets[f] = repstr(strtrim(sxpar(h,'TARGNAME'),2),' ', '_')  ; spaces cause barfing in filenames
     fileexposures[f] = 60.
     filefilters[f] = string(sxpar(h,'FILTER'))
     filesatval[f] = sxpar(h,'SATURATE')
     fileskyval[f] = sxpar(h,'COUNTS')
     fileairval[f] = sxpar(h,'AIRMASS')
     if n_elements(photodata) gt 0 then begin
       pmatch = where(azffiles[f] eq photodata.filename, ct)
       if ct gt 0 then begin
         fileseeingpix[f] = photodata[pmatch].seeing
         filefluxratio[f] = photodata[pmatch].fluxratio
       endif
     endif
   endfor
   targets = unique(filetargets)

   for t = 0, n_elements(targets)-1 do begin
     target = targets[t]
     thistarget = where(filetargets eq target, cttarg)
     if cttarg eq 0 then continue
     thistargetfilts = unique(filefilters[thistarget])
     for l = 0, n_elements(thistargetfilts)-1 do begin
        filter = thistargetfilts[l]
        thistargfilt = where(filetargets eq target and filefilters eq filter, cttargfilt)
        if cttargfilt eq 0 then continue
        expthresh = (median(fileexposures[thistargfilt]) * 0.5) > (max(fileexposures[thistargfilt]) * 0.19) ; formerly median*0.8
        stacki = where(filetargets eq target and filefilters eq filter, ctstack)
        if ctstack eq 0 then continue
        stacklist = azffiles[stacki]
        stackexps = fileexposures[stacki]
        medianexp = median(stackexps)
        medair = median(fileairval[stacki])
        minair = min(fileairval[stacki])
        maxair = max(fileairval[stacki])
        minseeingpix = min(fileseeingpix[stacki])
        maxseeingpix = max(fileseeingpix[stacki])
        medseeingpix = median(fileseeingpix[stacki])
        minfluxratio = min(filefluxratio[stacki])
        maxfluxratio = max(filefluxratio[stacki])
        medfluxratio = median(filefluxratio[stacki])
        totalexp = total(stackexps)
        nstack = n_elements(stacklist)
        if nstack gt 1 then stdfluxratio = stdev(filefluxratio[stacki]) else stdfluxratio = 0

        outfile       = pipevar.imworkingdir + 'coadd' + strtrim(target,2) +'_'+ strtrim(filter,2) + '.fits'
        outweightfile = pipevar.imworkingdir + 'coadd' + strtrim(target,2) +'_'+ strtrim(filter,2) + '.weight.fits'
        stackcmd = pipevar.swarpcommand+' '    ;'swarp '
        for s = 0, n_elements(stacklist)-1 do begin
           if s eq 0 then stackcmd = stackcmd + stacklist[s] 
           if s gt 0 then stackcmd = stackcmd + ',' + stacklist[s] 
        endfor
        for s = 0, n_elements(stacklist)-1 do begin
           weightfilename = strmid(stacklist[s],0,strlen(stacklist[s])-5-0)  + '.weight.fits'
           weightexists = file_test(weightfilename)
           if weightexists eq 0 then begin
              slashpos = strpos(weightfilename,'/',/reverse_search)
              weightfilenameinit = pipevar.imworkingdir + strmid(weightfilename,slashpos+2)  
              weightexists = file_test(weightfilenameinit)
              if weightexists then begin
                  print, 'mv '+weightfilenameinit+' '+weightfilename
                  spawn, 'mv '+weightfilenameinit+' '+weightfilename
              endif
           endif
        endfor

        stackcmd = stackcmd + ' -IMAGEOUT_NAME ' + outfile + ' -WEIGHTOUT_NAME ' + outweightfile
        if nstack gt 1 then begin ; used to be 3...
          stackcmd = stackcmd + ' -WEIGHT_TYPE MAP_WEIGHT'
        endif else begin
          stackcmd = stackcmd + ' -WEIGHT_TYPE BACKGROUND'
        endelse
        stackcmd = stackcmd + ' -FSCALE_DEFAULT '
        for s = 0, n_elements(stacklist)-1 do stackcmd = stackcmd + strtrim(string(medianexp/stackexps[s]),2) + ','   ; multiplicative factor...
        stackcmd = stackcmd + ' -COPY_KEYWORDS OBJECT,TARGNAME,TELESCOP,FILTER,DICHNAME,SLITNAME,'+$
                               'GRISNAME,GRANAME,GRANGLE,INSTRUME,UT,UTC,MJD_OBS,MJD-OBS,ST,'+$
                               'DATE,DATE-OBS,AIRMASS,AZ,RA,DEC,EL,HA,TELFOCUS,ROTPOSN,DOMEPOSN,CURRINST,OBSERVER,'+$
                               'FLATFLD,FLATTYPE'

        if (file_test(outfile) eq 0) or pipevar.overwrite then begin
           if nstack eq 1 then $ ; used to be 3 for some bizarre reason?
              print, 'Warning - only ', clip(nstack), ' exposures; not flagging bad pixels.'
           print, 'Stacking ', target, ':'
           for s = 0, n_elements(stacklist)-1 do print, stacklist[s], '  ', strtrim(string(stackexps[s]),2) + 's'
           print, stackcmd
           spawn, stackcmd  ; do the stack
           data = mrdfits(outfile, 0, h, /silent)
           sxaddpar, h, 'DATE', datestr
           sxaddpar, h, 'NSTACK', nstack
           sxaddpar, h, 'AIRMASS', medair, 'Median exposure airmass'
           sxaddpar, h, 'AIRMIN', minair, 'Minimum exposure airmass'
           sxaddpar, h, 'AIRMAX', maxair, 'Maximum exposure airmass'
           sxaddpar, h, 'EXPOSURE', medianexp, 'Effective rescaled exposure time'
           sxaddpar, h, 'TOTALEXP', totalexp, 'Total summed integration time'
           sxaddpar, h, 'MAXEXP', max(stackexps), 'Length of longest exposure'
           sxaddpar, h, 'MINEXP', min(stackexps), 'Length of shortest exposure'
           sxaddpar, h, 'SATURATE', min(filesatval[stacki]-fileskyval[stacki])
           sxaddpar, h, 'MEDSKY', median(fileskyval[stacki], /even)
           if minseeingpix gt 0 then begin
              sxaddpar, h, 'SEEING', medseeingpix, 'Median seeing in pixels'
              sxaddpar, h, 'SEEMIN', minseeingpix
              sxaddpar, h, 'SEEMAX', maxseeingpix
              if camera eq 'red' and pipevar.lrisversion eq 1 then medseeingarcsec = medseeingpix*0.210 else medseeingarcsec = medseeingpix*0.135
              sxaddpar, h, 'SEEARCSC', medseeingarcsec, 'Median seeing in arcsec'
           endif
           if minfluxratio gt 0 then begin
              sxaddpar, h, 'TRANSRAT', maxfluxratio/minfluxratio, 'Transmission ratio of max/min images'
              sxaddpar, h, 'TRANSSTD', stdfluxratio, 'Std. dev. of relative transmission'
           endif
           for s = 0, n_elements(stacklist)-1 do  sxaddpar, h, 'IMAGE'+strtrim(string(s),2), stacklist[s]
           for s = 0, n_elements(stacklist)-1 do  sxaddpar, h, 'IMEXP'+strtrim(string(s),2), stackexps[s]
           get_date, now
           sxdelpar, h, ['SOFTNAME','SOFTVERS','SOFTDATE','SOFTAUTH','SOFTINST','AUTHOR']
           hist = 'Processed by lrisautoproc '+now
           sxaddhist, hist, h, /COMMENT
           mwrfits, data, outfile, h, /create ; add this header info

           ; In the future, there will be a final alignment step.

        endif
     endfor
  endfor

  qq = 0
  if qq eq 0 then begin; chip eq 'b' and (camera ne 'red' and lrisversion gt 1) then begin 
     ; Combine left and right sides (or multiple exposure blocks, if necessary) for COADDS

     coaddfiles = findfile(pipevar.imworkingdir+'coadd*.fits')

     if coaddfiles ne [''] then coaddfiles = coaddfiles[where(coaddfiles ne '')] ;else continue

     ; match every r with an l and coadd
     for f = 0, n_elements(coaddfilesr)-1 do begin
        filename = coaddfiles[f]
        if ct eq 1 then begin
           outfile = coaddfiles[f]
           if file_test(outfile) eq 0 or pipevar.overwrite then begin
              weightname = repstr(filename, '.fits','.weight.fits') 
              print
              stackcmd = stackcmd + ' -IMAGEOUT_NAME ' + outfile
              stackcmd = stackcmd + ' -WEIGHT_TYPE MAP_WEIGHT'
              stackcmd = stackcmd + ' -COPY_KEYWORDS OBJECT,TARGNAME,TELESCOP,FILTER,DICHNAME,'+$
                                  'SLITNAME,GRISNAME,GRANAME,GRANGLE,INSTRUME,CAMERA,UTC,MJD_OBS,ST,'+$
                                  'DATE,DATE-OBS,AIRMASS,AZ,RA,DEC,EL,HA,TELFOCUS,ROTPOSN,DOMEPOSN,CURRINST,OBSERVER,'+$
                                  'DOMEPOSN,CURRINST,ADCWAVE0,ADCWAVE1,AMPMODE,CCDGAIN,CCDSPEED,NUMAMPS,'+$
                                  'FLATFLD,FLATTYPE,EXPOSURE,TOTALEXP,MINEXP,MAXEXP'+$
                                  'NSTACK,MEDSKY,AIRMASS,AIRMIN,AIRMAX,SEEING,SEEMIN,SEEMAX,SEEARCSEC,TRANSRAT,TRANSSTD'
              print, stackcmd
              spawn, stackcmd
              data = mrdfits(outfile, 0, h)
              hr = headfits(rfilename)
              hl = headfits(lfilename)
              for i = 0, 100 do begin
                rimagekeyn = sxpar(hr,'IMAGE'+clip(i))
                if clip(rimagekeyn) ne '0' then sxaddpar, h, 'RIMAGE'+clip(i), clip(rimagekeyn)
              endfor
              for i = 0, 100 do begin
                limagekeyn = sxpar(hl,'IMAGE'+clip(i))
                if clip(limagekeyn) ne '0' then sxaddpar, h, 'LIMAGE'+clip(i), clip(limagekeyn)
              endfor
              sxdelpar, h, ['SOFTNAME','SOFTVERS','SOFTDATE','SOFTAUTH','SOFTINST','AUTHOR']
              get_date, now
              hist = 'Processed by lrisautoproc '+now
              sxaddhist, hist, h, /COMMENT
              mwrfits, data, outfile, h, /create
           endif
        endif
     endfor


     ; For individual images (standards, found if exptime < 30 sec)
     
     azffiles = findfile(pipevar.imworkingdir+'a*f*'+prefchar+pipevar.wildchar+'r.fits')
     realfiles = where(azffiles ne '', ct)

     ; Combine l and r images.
     doindlrcombine = 0
     if doindlrcombine and (ct gt 0) then begin
      for f = 0, n_elements(azffiles)-1 do begin
        infiler = azffiles[f]
        h = headfits(infiler)
        exptime = sxpar(h, 'ELAPTIME')
        filen = strmid(infiler,strpos(infiler,'r.fits')-4,4)
        target = repstr(clip(sxpar(h,'TARGNAME')),' ', '_')
        if exptime le 30 then begin
           infilel = repstr(infiler,'r.','l.')
           outfile = target + '_' + filen + 'o.fits'
           if file_test(infilel) eq 0 then continue                    ; mate doesn't exist
           if file_test(outfile) eq 1 and pipevar.overwrite eq 0 then continue  ; already stacked these two

           print
           print, 'Preparing to combine ' + infiler + ' and ' + infilel
           print, pipevar.autoastrocommand+' '+infiler+' -q'
           spawn, pipevar.autoastrocommand+' '+infiler+' -q'
           print, pipevar.autoastrocommand+' '+infilel+' -q'
           spawn, pipevar.autoastrocommand+' '+infilel+' -q'

           stackcmd = pipevar.swarpcommand+' ' + extractpath(infiler)+'a' + removepath(infiler) + ',' + extractpath(infilel) + 'a' + removepath(infilel)
           print
           print, 'Combining ' + extractpath(infiler)+'a'+removepath(infiler) + ' and a' + extractpath(infilel)+'a'+removepath(infilel)
           stackcmd = stackcmd + ' -IMAGEOUT_NAME ' + outfile ; no weighting for this one
           stackcmd = stackcmd + ' -VERBOSE_TYPE QUIET' ; fast, so hide this
           stackcmd = stackcmd + ' -COPY_KEYWORDS OBJECT,TARGNAME,TELESCOP,FILTER,'+filtkey+',DICHNAME,'+$
                               'SLITNAME,GRISNAME,GRANAME,GRANGLE,INSTRUME,CAMERA,UTC,MJD_OBS,ST,'+$
                               'DATE-OBS,AIRMASS,AZ,RA,DEC,EL,HA,TELFOCUS,ELAPTIME,'+$
                               'DOMEPOSN,CURRINST,ADCWAVE0,ADCWAVE1,AMPMODE,CCDGAIN,CCDSPEED,NUMAMPS,'
           print, stackcmd
           spawn, stackcmd
           print, ' -> ', outfile
           print
        endif
      endfor
     endif
   endif

	outpipevar = pipevar
	
end



; -------------------------

; need to allow the user to completely ignore the left chip at all stages by specifying the chip.
; need to restore the gain correction.
; need to do something about when crashes, leaves you in imredux (check if you are already in the imredux directory)

pro ratautoproc, datadirectory=datadirectory, modestr=modestr, camerastr=camerastr, chipstr=chipstr, start=start, stop=stop, only=only, step=step, nocrclean=nocrclean, redo=redo
;   modestr      - Mode (imaging or spectroscopy)
;   camerastr    - Camera to process (red or blue)
;   chipstr      - Chip to process (left or right)
;   start        - Start with this step, skipping previous ones
;   stop         - End with this step, skipping subsequent ones
;   only         - Do only this step
;   step         -  (completely identical to only, takes precedence)
;   redo         - Overwrite any existing products
;   nocrclean    - Do not zap cosmic rays

!quiet = 1

; Load default parameters and interpret user arguments.

close, /all

	;VLT REMOVE WILDCHAR and LRISVERSION WHEN FULLY COMPLETE
	
	pipevar = {autoastrocommand:'autoastrometry' , sexcommand:'sex' , swarpcommand:'swarp' , $
					datadir:'' , imworkingdir:'' , overwrite:0 , modestr:'',$
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:'', wildchar: '?????????????????_???_?', lrisversion:3 }
	pipevar.modestr='im'
	
	if keyword_set(redo) then pipevar.overwrite=1
	if n_elements(datadirectory) gt 0 then pipevar.datadir = datadirectory		

	autopipedefaults, outpipevar=pipevar, inpipevar=pipevar
	modes = pipevar.modestr
	
nocrclean = keyword_set(nocrclean)

cd, current=pwd
dirtree = strsplit(pwd,'/',/extract,count=nd)
lastdir = dirtree[nd-1]
if lastdir ne '' then lastdir += '/'  ; whether or not a slash is on the end is very
                                      ; confusing, need to rethink this.
if lastdir eq 'imredux/' then begin
     print, 'Currently in a reduction subdirectory.'
     print, 'Type cd.. and rerun.'
  ; could reinterpret this as run with the mode set to whatever this directory is...
     return
endif

cameras=['blue']                ;placeholder

; --- Process chip options
if n_elements(chipstr) eq 0 then chipstr = 'r'
chipstr = strlowcase(chipstr)
if chipstr eq 'lr' or chipstr eq 'l,r' then chips = ['l', 'r']
if chipstr eq 'rl' or chipstr eq 'r,l' or chipstr eq 'both' or chipstr eq 'b' then chips = ['r', 'l']
if n_elements(chips) eq 0 then begin  
   if strmid(chipstr,0,1) eq 'l' then chips = ['l']
   if strmid(chipstr,0,1) eq 'r' then chips = ['r']
endif
if n_elements(chips) eq 0 then begin
  print, 'Cannot recognize chip request: ', chipstr
  return
endif

; --- Process step options
steps = ['prepare', 'makeflat', 'flatten', 'makesky', 'skysub', 'crclean', 'skysubtract', 'astrometry', 'photometry', 'stack']
if n_elements(start) gt 0 then begin
   w = (where(steps eq start, ct)) [0]
   if ct eq 0 then begin
      print, "Invalid starting step '", start, "
      print, "Must be one of: ", steps
      return
   endif
   steps = steps[w:*]
endif
if n_elements(stop) gt 0 then begin
   w = (where(steps eq stop, ct)) [0]
   if ct eq 0 then begin
      print, "Invalid stopping step '", stop, "'
      print, "Must be one of: ", steps
      if n_elements(start) gt 0 then print, 'Note that start is also set.'
      return
   endif
   steps = steps[0:w]
endif
if n_elements(step) gt 0 then only = step
if n_elements(only) gt 0 then begin
   w = (where(steps eq only, ct)) [0]
   if ct eq 0 then begin
      print, "Invalid step '", only, "'
      print, "Must be one of: ", steps
      if n_elements(start) gt 0 then print, 'Note that start is also set.'
      if n_elements(stop) gt 0 then print, 'Note that stop is also set.'
      return
   endif
   steps = steps[w]
endif

; For imaging, check if autoastrometry, sextractor, swarp are installed and functioning

if total(modes eq 'im') ge 1 then begin
   if total(steps eq 'astrometry') gt 0 or total(steps eq 'photometry') gt 0 or total(steps eq 'stack') gt 0 then begin
      if file_test('temp.txt') then spawn, 'rm -f temp.txt'
      cmd = getsexpath()+'sex -d '+' > temp.txt'
      spawn, cmd
      c = countlines('temp.txt')
      if c eq 0 then begin
         print, 'Error: Sextractor is not installed or not configured.'
         print, '   Cannot run image alignment steps.'
         print, "   Configure or set mode='s' or stop='crclean'"
         return
      endif
   endif

   if total(steps eq 'stack') gt 0 then begin
      if file_test('temp.txt') then spawn, 'rm -f temp.txt'
      cmd = pipevar.swarpcommand+' -d > temp.txt'
      spawn, cmd
      c = countlines('temp.txt')
      if c eq 0 then begin
         print, 'Error: Swarp is not installed or not configured.'
         print, '   Cannot run image coadds.'
         print, "   Configure or set mode='s' or stop='photometry'"
         return
      endif
   endif

   if total(steps eq 'astrometry') or total(steps eq 'stack') gt 0 then begin
      if file_test('temp.txt') then spawn, 'rm -f temp.txt'
      cmd = pipevar.autoastrocommand+' > temp.txt'
      spawn, cmd
      c = countlines('temp.txt')
      if c eq 0 then begin
         print, 'Error: Autoastrometry is not installed or not configured.'
         print, '   Cannot run image alignment steps.'
         print, "   Configure or set mode='s' or stop='crclean' "
         return
      endif
   endif
endif

lrisversion = 3

if keyword_set(nocrclean) eq 0 then flag = 1

nsteps = n_elements(steps)
nmodes = n_elements(modes)
ncameras = n_elements(cameras)
nchips = n_elements(chips)

for istep = 0, nsteps-1 do begin
   instep = steps[istep]

   for icam = 0, ncameras-1 do begin
      camera = cameras[icam]
      ca = strmid(cameras[icam],0,2)

      if instep eq 'prepare' then begin
         ; the mode setting is passed on to the routine itself.
         autopipeprepare, outpipevar=pipevar, inpipevar=pipevar
      endif

      for imode = 0, nmodes-1 do begin
         mo = strmid(modes[imode],0,2)
         if mo eq 'im' then begin   
            if instep eq 'flatten'  then autolrisimflatten,  cam=camera, chip='', outpipevar=pipevar, inpipevar=pipevar
            if instep eq 'makesky' then autolrismakesky,cam=camera, chip='', outpipevar=pipevar, inpipevar=pipevar
            if instep eq 'skysub'  then autolrisskysub,  cam=camera, chip='', outpipevar=pipevar, inpipevar=pipevar
         endif

         for ichip = 0, nchips-1 do begin
            ch = strmid(chips[ichip],0,1)

            if nocrclean eq 0 and instep eq 'crclean' then begin
               if mo eq 'im' then autolriscrcleanim,    cam=camera, chip=ch, outpipevar=pipevar, inpipevar=pipevar
            endif 

            if mo eq 'im' then begin   
               if instep eq 'astrometry' then autolrisastrometry, chip=ch, cam=camera, outpipevar=pipevar, inpipevar=pipevar
               if camera eq 'blue' then bl = 1 else bl = 0
               if camera eq 'red' then  re = 1 else re = 0
               if n_elements(chips) eq 1 then begin
                  if instep eq 'photometry' then autolrisphotometry, chip=ch, camera=camera, outpipevar=pipevar, inpipevar=pipevar
               endif else begin
                  if instep eq 'photometry' and ichip eq 1 then autolrisphotometry, camera=camera, outpipevar=pipevar, inpipevar=pipevar
               endelse
               if instep eq 'stack' then autolrisstack, chip=ch, cam=camera, outpipevar=pipevar, inpipevar=pipevar
            endif
         endfor
      endfor ; mode
   endfor ; camera
endfor ; step

print

if strlen(pipevar.flatfail) gt 0 then begin
	print
  	print, 'Unable to flat-field the following images:'
  	ffailfile = strsplit(pipevar.flatfail, /extract)
  	for f = 0, n_elements(ffailfile)-1 do begin
  		print, ffailfile[f]
  	endfor
endif

nafail = strlen(pipevar.relastrofail) + strlen(pipevar.fullastrofail) + strlen(pipevar.catastrofail)
if nafail gt 0 then begin
	print
	
	if strlen(pipevar.catastrofail) gt 0 then begin
		print, 'Unable to produce astrometric catalogs for the following reference images:'
		cafailfile = strsplit(pipevar.catastrofail, /extract)
		for f=0, n_elements(cafailfile)-1 do begin
			print, cafailfile[f]
		endfor
	endif
			
	if strlen(pipevar.relastrofail) gt 0 then begin	
		print, 'Relative astrometry failed for the following images, but absolute was successful:'
		rafailfile = strsplit(pipevar.relastrofail, /extract)
		for f=0, n_elements(rafailfile)-1 do begin
			print, rafailfile[f]
		endfor
	endif
			
	if strlen(pipevar.fullastrofail) gt 0 then begin
		print, 'All astrometry failed for the following images (not stacked):'
		fafailfile = strsplit(pipevar.fullastrofail, /extract)
		for f=0, n_elements(fafailfile)-1 do begin
			print, fafailfile[f]
		endfor
	endif	
endif	

print, 'Processing complete.'
!quiet = 0

; cleanup- should ultimately put these things in imredux

if file_test('temp*.*') gt 0 then spawn, 'rm -f temp*.*'
if file_test('det.*') gt 0 then spawn, 'rm -f det.*'
if file_test('cat.*') gt 0 then spawn, 'rm -f cat.*'

end

 
