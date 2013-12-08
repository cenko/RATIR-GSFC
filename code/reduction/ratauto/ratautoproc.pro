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
;   pipeautoproc.par - Contains various defaults (mainly directory paths)
;   catalog.txt - Source positional catalog.  Use to auto-correct the object
;        name in the header based on position (independent of TARGNAME.)
;   seeing.txt - A file with only two numbers, the minimum and maximum seeing
;        experienced in arcseconds.
;
;
;
; COMMENTS:
;   This code is meant to be fully automated, producing science-quality
;   output with a single command line with zero intervention.  Reduces
;   RATIR data.
; 
;   The production of images is quite high-level and includes photometric 
;   calibration of the field (although the accuracy of this has not been
;   robustly tested).
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
;   2*.fits format. They can be either in the working directory or 
;	in a different directory specified by datadir.
; 
;   The code runs in a series of steps in the following order:
;
;   1. "Prepare" the data by converting the multi-header extension FITS files
;      to a standard single frame, correcting/flagging known pixel defects,
;      cropping unused area of the chip, bias-subtracting using the overscan,
;      and adding extra information to the header.  The conversion/bias 
;      correction employs a modified version of readmhdu.pro.  
;      Output: p[b|r]*.fits (written to ./imredux/ by default.)
;
;   2. (TAKEN OUT, assume that have master flat in format flat_(filter).fits in imredux/ folder)
;	   Create flat-fields.  The program searches through all science exposures
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
;      of the image area during the run (which is variable.) Assumes master flat
;	   (ex. flat_H.fits in imredux/ folder)
;
;   4. Removes cosmic rays, using the independent routines pzap.pro and
;      pzapspec.pro.  See those programs for more information.  
;      This can be a time-consuming process.
;      Output: zfp[b|r]*[l|r].fits
;  
;   5. Solve astrometry of the field against the best available online catalog
;      (SDSS or USNO-B1.0), using the independent vlt_autoastrometry.py code.
;      (Also requires sextractor to be installed.)
;      One image is solved relative to a reference catalog, then the
;      remaining images are solved against the first.  NO distortion corrections
;      are currently applied.
;      Output: azfp[b|r]*[l|r].fits
;   
;   6. (REMOVED) Produce photometric solution of the field.  The GSFC IDLastro aper.pro
;      program is used to measure photometry of point sources in each field
;      for comparison of different exposures (on the same target) to establish 
;      the relative zeropoints.  
;      Output: [object].cal.list and other text files
;  
;   7. Stack exposures.  A weighted, masked median is performed for each field.
;      Requires swarp to be installed and properly linked.
;      Output: coadd[object].[filter].[l|r|]fits
;
;
; 
; EXAMPLES:
;   1.  In a directory containing a full night worth of RATIR data, enter
;
;       IDL> ratautoproc
;
;       This will automatically execute all of the steps above.
;       Depending on the quantity of data from your run and the speed of your
;       computer, the reduction of data from an entire night takes 20 minutes
;       up to about two hours.
;
;
;   2.  Run only the "prepare" step (bias subtraction/formatting), on data 
;       stored in a separate directory:
;
;       IDL> ratautoproc, mode='im', step='prepare', datadir='raw/'
;
;
;   3.  For live-pipeline usage, I recommend placing all the data in a
;       subdirectory 'raw/', then running ratautoproc one step down as 
;       follows:
;
;       IDL> ratautoproc, datadir='raw/', mode='im'
;
;       This will save some processing time.
;
;
;
;
;
;   INSTALLATION:
;
;   Create the subdirectory 'ratauto' somewhere on your hard drive (probably
;   in your IDL directory), and unpack the contents of this installation file
;   there (e.g.,  tar -xvf ratauto.tar.gz).  You will need to tell IDL about
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
;   file ratautoproc.par AND also edit the dependencies/vlt_autoastrometry.py file
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
; check that USNO catalogs are loaded OK
; check for looking for "nickel" catalogs
; faster CR cleaner for blue (and old red): done?
; better text output; error/warning log written to disk
; background subtraction - if there is a neighbor it will be messed up.
; consistent aperture size for repeated observations?  Or, stack in 2D space?
; photometry table - need to be able to update on reruns, and deal with un-photometry-able fields
; automated halo removal for bright stars?
; not sure if bpm is working, still a bad line on red side (when binned, anyway)
;  this appears to be due to movement of the weakened columns.  Make sure this isn't a crop problem and 
;  if not try to do something to detect this...
; interpolate standards in log space
; sky subtraction of galaxies still causes issues, see 120427_0114
; still some tracing issues, see 120427_0111
; matching of apertures for multiple reobservations and red/blue...?
; ------------------------

; -------------------------
; need to restore the gain correction.

;   start        - Start with this step, skipping previous ones
;   stop         - End with this step, skipping subsequent ones
;   only         - Do only this step
;   step         -  (completely identical to only, takes precedence)
;   redo         - Overwrite any existing products
;   nocrclean    - Do not zap cosmic rays
pro ratautoproc, datadirectory=datadirectory, start=start, stop=stop, only=only, step=step, nocrclean=nocrclean, redo=redo

	!quiet = 1
	close, /all
	
	; Load default parameters and interpret user arguments.
	pipevar = {autoastrocommand:'autoastrometry' , sexcommand:'sex' , swarpcommand:'swarp' , $
					datadir:'' , imworkingdir:'' , overwrite:0 ,$
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:''}
	
	if keyword_set(redo) then pipevar.overwrite=1
	if n_elements(datadirectory) gt 0 then pipevar.datadir = datadirectory		

	autopipedefaults, outpipevar=pipevar, inpipevar=pipevar

	;Step options
	steps = ['prepare', 'flatten', 'makesky', 'skysub', 'crclean', 'astrometry', 'stack']

	;If start is specified, truncate steps to start at specified step.  If invalid step end program with error
	if n_elements(start) gt 0 then begin
   		w = (where(steps eq start, ct)) [0] 	
   		if ct eq 0 then begin
      		print, "Invalid starting step '", start, "
      		print, "Must be one of: ", steps
      		return
   		endif 		
   		steps = steps[w:*]
	endif
	
	;If stop is specified, truncate steps to end at specified step.  If invalid step end program with error
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
	
	;If step specified set only to specified step
	if n_elements(step) gt 0 then only = step

	;If only specified (including step), set steps to run as specified step.  If invalid step end program with error
	if n_elements(only) gt 0 then begin
   		w = (where(steps eq only, ct)) [0]
   		
   		if ct eq 0 then begin
      		print, "Invalid step '", only, "'
      		print, "Must be one of: ", steps
      		if n_elements(start) gt 0 then print, 'Note that start is also set.'
      		if n_elements(stop)  gt 0 then print, 'Note that stop is also set.'
      		return
   		endif
   		
   		steps = steps[w]
	endif
	
	nocrclean = keyword_set(nocrclean)
	
	;Check if autoastrometry, sextractor, swarp are installed and functioning before running steps
   	if total(steps eq 'astrometry') gt 0 or total(steps eq 'photometry') gt 0 or total(steps eq 'stack') gt 0 then begin
      	if file_test('temp.txt') then spawn, 'rm -f temp.txt'
      	cmd = getsexpath()+'sex -d '+' > temp.txt'
      	spawn, cmd
		c = countlines('temp.txt')
      		
      	if c eq 0 then begin
     		print, 'Error: Sextractor is not installed or not configured.'
     		print, '   Cannot run image alignment steps.'
         	print, "   Configure or stop='crclean'"
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
         	print, "   Configure or stop='photometry'"
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
         	print, "   Configure or stop='crclean' "
         	return
      	endif
   	endif

	;Runs each processing step specified in the correct order (crclean is optional)
	for istep = 0, n_elements(steps)-1 do begin
   		instep = steps[istep]

      	if instep eq 'prepare' then autopipeprepare, outpipevar=pipevar, inpipevar=pipevar      		 
        if instep eq 'flatten' then autopipeimflatten, outpipevar=pipevar, inpipevar=pipevar
        if instep eq 'makesky' then autopipemakesky,   outpipevar=pipevar, inpipevar=pipevar
        if instep eq 'skysub'  then autopipeskysub,    outpipevar=pipevar, inpipevar=pipevar
    	
    	if nocrclean eq 0 and instep eq 'crclean' then begin
    		autopipecrcleanim, outpipevar=pipevar, inpipevar=pipevar
        endif 
            		 
        if instep eq 'astrometry' then autopipeastrometry, outpipevar=pipevar, inpipevar=pipevar
        if instep eq 'stack' then autopipestack, outpipevar=pipevar, inpipevar=pipevar
	endfor

	;Prints the files that were not flat fielded due to problems with file
	if strlen(pipevar.flatfail) gt 0 then begin
		print
  		print, 'Unable to flat-field the following images:'
  		ffailfile = strsplit(pipevar.flatfail, /extract)
  		for f = 0, n_elements(ffailfile)-1 do begin
  			print, ffailfile[f]
  		endfor
	endif

	;Prints the files that were not astrometry corrected due to problems with the file 
	;(specifies catalog, relative, and absolute failures)
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

	;Remove any files that were created during the reduction process
	if file_test('temp*.*') gt 0 then spawn, 'rm -f temp*.*'
	if file_test('det.*')   gt 0 then spawn, 'rm -f det.*'
	if file_test('cat.*')   gt 0 then spawn, 'rm -f cat.*'

end