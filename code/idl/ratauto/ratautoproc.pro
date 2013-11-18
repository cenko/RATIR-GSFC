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
;   continuous   - Keep running, assimilating new data as it appears
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
;  (nofringe)    - Skip fringe production and correction for LRIS-R1
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

; -------------------------
pro autolrisdefaults
   common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite, fileseq, filterstr, gratingstr, targetstr
   common lrisfail, flatfail, catastrofail, relastrofail, fullastrofail, extractfail, wavsolfail, wavsolwarn, fluxcalfail
   common lrisconfig, lrisautopath, refdatapath, defaultspath

   paths = strsplit(!path,':',/extract)
   for p = -1, n_elements(paths)-1 do begin   ; find the directory this file is in
      if p ge 0 then path = paths[p] $
                     else path = ''
      sfile = path+'/ratautoproc.pro'
      n = countfiles(sfile)
      if n gt 0 then lrisautopath = path
      if n gt 0 then break
   endfor
   if n eq 0 then begin
      print, 'Problem finding lrisautoproc directory.  This error should not occur.'
   endif
   for p = -1, 1 do begin
      if p eq -1 then path = '.'
      if p eq  0 then path = lrisautopath
      if p eq  1 then path = lrisautopath + '/defaults'
      sfile = path+'/lrisautoproc.par'       ; get the configuration file
      n = countfiles(sfile)
      if n gt 0 then begin
         openr, 5, sfile
         iline = ''
         while not eof(5) do begin
            readf, 5, iline
            colonpos = strpos(iline,':')
            if colonpos lt 1 or colonpos ge strlen(iline) then continue
            varname = clip(strmid(iline,0,colonpos))
            varset = clip(strmid(iline,colonpos+1))
            if varname eq 'refdatadir' then refdatapath = varset
            if varname eq 'defaultsdir' then defaultspath = varset
            if varname eq 'autoastrocommand' then autoastrocommand = varset
            if varname eq 'swarpcommand' then swarpcommand = varset
            if varname eq 'sexcommand' then sexcommand = varset
            if varname eq 'datadir' and n_elements(datadir) eq 0 then datadir = varset ; command line overrides file value
            if varname eq 'imworkingdir' then imworkingdir = varset
            if varname eq 'spworkingdir' then spworkingdir = varset
            if varname eq 'imoutputdir'  then imoutputdir = varset
            if varname eq 'spoutputdir'  then spoutputdir = varset
         endwhile
         break ; found the file: stop looping
      endif
   endfor

   ; use internal defaults to set anything still unspecified
   if n_elements(refdatapath) eq 0 then refdatapath = lrisautopath+'/refdata'
   if n_elements(defaultspath) eq 0 then defaultspath = lrisautopath+'/defaults'
   if n_elements(autoastrocommand) eq 0 then autoastrocommand = 'autoastrometry'
   if n_elements(swarpcommand) eq 0 then swarpcommand = 'swarp'
   if n_elements(sexcommand) eq 0 then sexcommand = 'sex'
   if n_elements(datadir) eq 0 then datadir = '' 
   if n_elements(imworkingdir) eq 0 then imworkingdir = ''
   if n_elements(spworkingdir) eq 0 then spworkingdir = ''
   if n_elements(imfinaldir) eq 0 then imfinaldir = ''
   if n_elements(spfinaldir) eq 0 then spfinaldir = ''
   if n_elements(wildchar) eq 0 then wildchar = '?????????????????_???_?'
   if n_elements(overwrite) eq 0 then overwrite = 0

   setswarppath, extractpath(swarpcommand)
   setsexpath, extractpath(sexcommand)

   ; these must always start blank
   flatfail = ['']
   catastrofail = ['']
   relastrofail = ['']
   fullastrofail = ['']
   extractfail = ['']
   wavsolfail = ['']
   wavsolwarn = ['']
   fluxcalfail = ['']

  if strlen(imworkingdir) gt 0 and file_test(imworkingdir) eq 0 then begin
     print, 'Creating imaging working directory ', imworkingdir
     spawn, 'mkdir '+imworkingdir
     check = file_test(imworkingdir)
     if check eq 0 then begin 
       print, 'WARNING - Failed to create image working directory!'
       print, '          Using current directory.  Check lrisautoproc.par inputs.'
     endif
  endif
  if strlen(spworkingdir) gt 0 and file_test(spworkingdir) eq 0 then begin
     print, 'Creating spectroscopy working directory ', spworkingdir
     spawn, 'mkdir '+spworkingdir
     check = file_test(spworkingdir)
     if check eq 0 then begin 
       print, 'WARNING - Failed to create spectrum working directory!'
       print, '          Using current directory.  Check lrisautoproc.par inputs.'
     endif
  endif

end

; -------------------------
pro autolrisprepare, modestr=modestr, camstr=camstr
   common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite, fileseq, filtestr, gratingstr, targetstr
   common lrisconfig, lrisautopath, refdatapath, defaultspath

  ; Prepare images for processing (debias, crop, fix)
  ; Current version always merges the two sides together, leaving splitting to later tasks.

  if n_elements(modestr) gt 0 then begin
     if modestr eq 'is' or modestr eq 'i,s' or modestr eq 'im,sp'  then mode = ''
     if modestr eq 'si' or modestr eq 's,i' or modestr eq 'sp,im'  then mode = ''
     if n_elements(mode) eq 0 then mode = strmid(modestr,0,1)
  endif else begin
     mode = ''
  endelse

  prefchar = '2'

  files = findfile(datadir+prefchar+wildchar+'.fits')
  pfiles = [findfile(imworkingdir+'p'+prefchar+wildchar+'.fits'), $
            findfile(spworkingdir+'p'+prefchar+wildchar+'.fits')]

  if datadir ne '' then  begin
     print, 'Looking for raw data at: ', datadir+prefchar+wildchar+'.fits'

     if n_elements(files) gt 0 then begin
        print, 'Found ', clip(n_elements(files)), ' files'
     endif else begin
        print, 'Did not find any files!  Check your data directory path!'
     endelse
  endif

  namefixfiles = ['']
  catalogfile = 'catalog.txt' ; current directory
  if file_test(catalogfile) then namefixfiles = [namefixfiles, catalogfile]
  landoltfieldposfile = refdatapath+'/landoltpos.txt'
  if file_test(landoltfieldposfile) then namefixfiles = [namefixfiles, landoltfieldposfile]
  specstandardposfile = refdatapath+'/specstandards.dat'
  if file_test(specstandardposfile) then namefixfiles = [namefixfiles, specstandardposfile]
  if n_elements(namefixfiles) gt 1 then namefixfiles = namefixfiles[1:*] else delvarx, namefixfiles

  for f = 0, n_elements(files)-1 do begin
     if files[f] eq '' then continue
     if n_elements(fileseq) gt 0 then if total(filenum(files[f]) eq fileseq) eq 0 then continue
    
     slashpos = strpos(files[f],'/',/reverse_search)    
     fileroot = strmid(files[f],slashpos+1)
     outnameim = imworkingdir + 'p' + fileroot
     outnamesp = spworkingdir + 'p' + fileroot
     matchi = where(outnameim eq pfiles, ct1)
     matchs = where(outnamesp eq pfiles, ct2)
     ct = ct1+ct2

     if ct eq 0 or overwrite gt 0 then begin
        filemode = 'i'          ;placeholder
        outname=outnameim

print,outname
        lrisprepare, files[f], gain=[1,1,1,1], flag=flag, /merge, crop='auto', $
          outname=outname, namefixfiles=namefixfiles
     endif
  endfor

end

; -------------------------
pro autolrismakeimflat, chip=chip, camera=camera

   common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite, fileseq, filterstr, gratingstr, targetstr

    if strlen(chip) gt 1 then fchip = strmid(chip,0,1) else fchip = chip

    prefchar = strmid(camera,0,1)

    files = findfile(imworkingdir+'p'+prefchar+wildchar+fchip+'.fits')
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
       if n_elements(fileseq) gt 0 then if total(filenum(files[f]) eq fileseq) eq 0 then continue
       h = headfits(files[f], /silent)

       slitname = strtrim(sxpar(h, 'SLITNAME'),2)
       if slitname ne 'direct' then continue
       filters[f] = clip(sxpar(h, 'FILTER'))
       counts[f] =  sxpar(h,'COUNTS')
       dichs[f] = clip(sxpar(h,'DICHNAME'))
       wins[f] = clip(sxpar(h,'WINDOW'))
       if clip(sxpar(h,'PANE')) ne '0' then wins[f] = clip(sxpar(h,'PANE')) else wins[f] = strmid(wins[f],2)
       binnings[f] = repstr(clip(sxpar(h,'BINNING')),',','x')
       targets[f] = clip(sxpar(h,'TARGNAME'))
       trapdoor = clip(sxpar(h,'TRAPDOOR'))

       exptime = sxpar(h, 'ELAPTIME')
       exps[f] = exptime
       azimuth = sxpar(h,'AZ')
       elevation = sxpar(h,'EL')
       domeazimuth = sxpar(h,'DOMEPOSN')

       domerelaz = domeazimuth - azimuth
       while domerelaz lt -180.0 do domerelaz += 360.
       while domerelaz ge  180.0 do domerelaz -= 360.
       if trapdoor eq 'open' then begin  ;abs(domerelaz-90) lt 1 or abs(domerelaz+90) lt 1))
         if (abs(elevation-45) lt 0.5 and abs(domerelaz) gt 3) or strpos(targets[f],'cass') ge 0 then begin
             isdomeflat[f] = 1  ;if counts[f] gt 7000 and counts[f] lt 40000 then 
         endif else begin
             if ((counts[f] gt 9000 and camera eq 'blue') or (counts[f] gt 12000 and camera eq 'red')) and counts[f] lt 43000 and exps[f] lt 60 then isskyflat[f] = 1
             if abs(domerelaz) lt 5.0 and ((counts[f] lt 40000) and (counts[f] le 15000 or (filters[f] eq 'I' and counts[f] le 30000) or (filters[f] eq 'RG850') or (exps[f] ge 60))) then isscience[f] = 1
         endelse
       endif

    endfor

    sciindex = where(isscience, ct)
    if ct eq 0 then return
    scifilterlist  = unique(filters[sciindex])
    for i = 0, n_elements(scifilterlist)-1 do begin
      filt = scifilterlist[i]

      dichlist = unique(dichs[sciindex[where(filt eq filters[sciindex])]])
      for d = 0, n_elements(dichlist)-1 do begin
        dich = dichlist[d]
        if n_elements(dichlist) gt 1 then dichstr = dichlist[d] else dichstr = ''

        binlist = unique(binnings[sciindex[where(filt eq filters[sciindex] and dich eq dichs[sciindex])]])
        for b = 0, n_elements(binlist)-1 do begin
          bin = binlist[b]
          if n_elements(binlist) gt 1 then binstr = 'b'+binlist[b] else binstr = ''

          outflatname = imworkingdir+'lris'+strmid(camera,0,1)+'flat'+filt+dichstr+binstr+fchip+'.fits'
          if file_test(outflatname) and overwrite eq 0 then continue ; flat exists already

          flatsuccess = 0

          ; Good supersky flat.
          ctsupersky = 0
          if (lrisversion eq 1 and camera eq 'red' and (filt eq 'I' or filt eq 'RG850' or filt eq 'GG570')) ne 1 then begin ; don't supersky fringed frames
           ; First, see if we have enough frames to do a supersky flat.
            poss = where(filters eq filt and dichs eq dich and binnings eq bin and isdomeflat eq 0 and counts gt 600 and exps gt 10., ctposs)
            if ctposs gt 0 then begin ;otherwise usually means some weird situation like the dichroic info is missing, etc.
               medcounts = median(counts[poss])
               if camera eq 'red' then superskymincount = 1000 > medcounts/3.
               if camera eq 'blue' then superskymincount = 500 > medcounts/3.
               superskymaxcount = (medcounts * 5) < (30000. > medcounts*1.2) ; 
               superskyminexp = 10
               superskyflats = where(filters eq filt and dichs eq dich and binnings eq bin and isdomeflat eq 0 and $
                                     counts gt superskymincount and counts lt superskymaxcount and exps ge superskyminexp, ctsupersky)
 
               if file_test('badflattargets.txt') eq 0 then begin
                  badflattargets = [''] 
               endif else begin
                  nbf = countlines ('badflattargets.txt')
                  openr, 6, 'badflattargets.txt'
                  badflattargets = strarr(nbf)
                  readf, 6, badflattargets
                  close, 6
                  notblank = where(badflattargets ne '', ct)
                  if ct gt 0 then badflattargets = badflattargets[where(badflattargets[notblank])]
               endelse


               for bad = 0, n_elements(badflattargets)-1 do begin
                  w = where(targets[superskyflats] ne badflattargets[bad], ctsupersky)
                  superskyflats = superskyflats[w]
               endfor
  
               if ctsupersky gt 1 then superskyfields = unique(targets[superskyflats]) else superskyfields = ['']
               nsuperskyfields = n_elements(superskyfields)
               ; this criterion is a bit loose and depends on the number of very bright stars in the image and how big the dither steps were.
               ; eventually, should make flatcombine able to find bright (ultra-saturated) objects and add halo masks, and make
               ; lrisautoproc treat big dithers like new fields.
               if (ctsupersky ge 21 and nsuperskyfields ge 5) or (ctsupersky ge 25 and nsuperskyfields ge 4) or $
                 (ctsupersky gt 29 and nsuperskyfields gt 3) or ctsupersky ge 33 then begin
                  print, filt, '-band / ', 'D'+dich,' '+binstr+': ', strtrim(string(ctsupersky),2), ' supersky flats (', clip(nsuperskyfields), ' unique fields)'
                  flatcombine, files[superskyflats], outflatname, /removeobjects, type='supersky', xnormkey='AMPR1'
                  flatsuccess = 1
               endif
            endif
         endif

         ; Regular sky flat.
         if flatsuccess eq 0 then begin
            skyflats = where(isskyflat and filters eq filt and dichs eq dich and binnings eq bin, ctsky)
            if ctsky ge 5 then begin
               print, filt, '-band / ', 'D'+dich,' '+binstr+': ', strtrim(string(ctsky),2), ' sky flats.'
               flatcombine, files[skyflats], outflatname, /removeobjects, type='sky', normzone=normzone, xnormkey='AMPR1'
               flatsuccess = 1
            endif
         endif

        ; Dome flat.
        if flatsuccess eq 0 then begin
           if strmid(camera,0,1) eq 'b' then mincounts=5000. else mincounts=7000.
           domeflats = where(isdomeflat and filters eq filt and dichs eq dich and binnings eq bin and $
                             counts gt mincounts and counts lt 39000., ctdome)
           if ctdome ge 3 then begin
              print, filt, '-band / ', 'D'+dich,' '+binstr+': ', strtrim(string(ctdome),2), ' dome flats'
              flatcombine, files[domeflats], outflatname, type='dome', normzone=normzone, xnormkey='AMPR1'
              flatsuccess = 1
           endif
        endif

        ; marginal sky flat
        if flatsuccess eq 0 then begin
            print, 'WARNING - insufficient images for a good flat.  ('+clip(ctdome)+' dome, '+clip(ctsky)+' sky, '+clip(ctsupersky)+' science).'
            if ctsky ge 3 then begin
               print, filt, '-band / ', 'D'+dich,' '+binstr+': ', strtrim(string(ctsky),2), ' sky flats (WARNING - marginal!)'
               flatcombine, files[skyflats], outflatname, /removeobjects, type='sky', xnormkey='AMPR1'
               flatsuccess = 1
            endif
        endif

        ; marginal supersky flat
        if flatsuccess eq 0 then begin
            if ctsupersky ge 5 then begin
               print, filt, '-band / ', 'D'+dich,' '+binstr+': ', strtrim(string(ctsupersky),2), ' supersky flats  (WARNING - EXTREMELY marginal)'
               print, 'Flat field for this filter/dichroic combination is likely to be extremely poor.'
               flatcombine, files[superskyflats], outflatname, /removeobjects, type='supersky', xnormkey='AMPR1'
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
     endfor  
  endfor


end


; -------------

pro autolrismakesky, chip=chip, camera=camera

    common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite, fileseq, filterstr, gratingstr, targetstr
  
    prefchar = '2'

    wildcharsky = '?????????????????_sky_?'
    files = findfile(imworkingdir+'fp'+prefchar+wildcharsky+'.fits')
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
       if n_elements(fileseq) gt 0 then if total(filenum(files[f]) eq fileseq) eq 0 then continue
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
          outflatname = imworkingdir+'sky-'+filt+'.fits'
          if ctsky ge 2 then begin
             print, filt, '-band sky flats.'
             if file_test(outflatname) and overwrite eq 0 then continue ; flat exists already
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
       
                                ;endfor
                                ;endfor  
    endfor

end

; -------------------------
pro autolrisimflatten, chip=chip, camera=camera

   common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite, fileseq, filterstr, gratingstr, targetstr
   common lrisfail, flatfail, catastrofail, relastrofail, fullastrofail, extractfail, wavsolfail, wavsolwarn, fluxcalfail

   prefchar = '2'
   fchip = strmid(chip,0,1)

   files = findfile(imworkingdir+'p'+prefchar+wildchar+fchip+'.fits')
   ffiles = findfile(imworkingdir+'fp'+prefchar+wildchar+fchip+'.fits')
   if n_elements(files) eq 1 and files[0] eq '' then return
   flats = findfile(imworkingdir+'*flat*.fits')
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
      if n_elements(fileseq) gt 0 then if total(filenum(files[f]) eq fileseq) eq 0 then continue
      outfile = fileappend(files[f], 'f')
      match = where(outfile eq ffiles, ct) ; check if output file exists
      if ct eq 0 or overwrite then begin               
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
            flatfail = [flatfail, files[f]]
            continue
         endif
         flatfile = flats[flatfileno[0]]
         print, 'Flattening ', removepath(files[f]), ' using ', removepath(flatfile)
         flatproc, files[f], flatfile, flatminval=0.3, crop='auto'  ; add autocropping of low-signal zones
      endif
    endfor

 end

; -------------------------
pro autolrisskysub, chip=chip, camera=camera

   common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite, fileseq, filterstr, gratingstr, targetstr
   common lrisfail, flatfail, catastrofail, relastrofail, fullastrofail, extractfail, wavsolfail, wavsolwarn, fluxcalfail

   prefchar = '2'
   fchip = strmid(chip,0,1)
   wildcharimg = '?????????????????_img_?'

   files = findfile(imworkingdir+'fp'+prefchar+wildcharimg+'.fits')
   sfiles = findfile(imworkingdir+'sfp'+prefchar+wildcharimg+'.fits')
   if n_elements(files) eq 1 and files[0] eq '' then return
   skys = findfile(imworkingdir+'*sky-*.fits')
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
      if n_elements(fileseq) gt 0 then if total(filenum(files[f]) eq fileseq) eq 0 then continue
      outfile = fileappend(files[f], 's')
      match = where(outfile eq sfiles, ct) ; check if output file exists
      if ct eq 0 or overwrite then begin               
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
            flatfail = [flatfail, files[f]]
            continue
         endif
         skyfile = skys[skyfileno[0]]
         print, 'Sky Subtracting ', removepath(files[f]), ' using ', removepath(skyfile)
         if strcmp(camera,'OPT') then skyprocopt, files[f], skyfile else skyproc, files[f], skyfile

      endif
      skipskysub:
    endfor

end

; -------------------------
pro autolriscrcleanim, chip=chip, camera=camera

   common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite, fileseq, filterstr, gratingstr, targetstr

   prefchar = '2'
   chipchar = strmid(chip,0,1)

   wildcharsky = '?????????????????_img_?'
   files = choosefiles(prefchar+wildchar+'.fits',imworkingdir+'isfp',imworkingdir+'sfp')
   if overwrite eq 0 then files = unmatched(files,'z')
   for f = 0, n_elements(files)-1 do begin
      if files[f] eq '' then continue
      if n_elements(fileseq) gt 0 then if total(filenum(files[f]) eq fileseq) eq 0 then continue
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
      if (camera eq 'red' and lrisversion eq 1) then usamp=1 else usamp=0
      if camera eq 'red' and lrisversion gt 1 then zeal = 0.6 ;evidently should be the new default
      if n_elements(setzeal) eq 1 then zeal = setzeal
      slashpos = strpos(files[f],'/')
      dir = strmid(files[f],0,slashpos+1)
      outname = dir + 'z' + strmid(files[f],slashpos+1)
      if file_test(outname) and overwrite eq 0 then continue
      pzap_perley, files[f], /weight, zeal=zeal, usamp=usamp, /quiet
   endfor

end

pro autolrisextract, camera=camera, chip=chip

   common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite, fileseq, filterstr, gratingstr, targetstr
   common lrisconfig, lrisautopath, refdatapath, defaultspath
   common lrisfail, flatfail, catastrofail, relastrofail, fullastrofail, extractfail, wavsolfail, wavsolwarn, fluxcalfail

   camchar = strmid(camera,0,1)
   chipchar = strmid(chip,0,1)

   if n_elements(badarclist) eq 0 then begin
     badarclist = ['']
     if file_test('badarc.txt') then begin
       nbarc = countlines ('badflattargets.txt')
       openr, 6, 'badarc.txt'
       badarcs = strarr(nbf)
       badarcstr = ''
       while not eof(6) do begin
          readf, 6, inline
          badarcstr += (inline + ' ')
       endwhile
       close, 6 
       badarclist = strlowcase(strsplit(badarcstr,/extract))
       notblank = where(badarcs ne '', ct)
     endif
   endif

   ; search for arcs
   allfiles = choosefiles(camchar+wildchar+chipchar+'.fits',spworkingdir+'fp')
   nf = n_elements(allfiles)
   isarc = intarr(nf)
   arcairmass = fltarr(nf)
   arcgr = strarr(nf)
   arcbinning = strarr(nf)
   arcwindow = strarr(nf)
   arccwav = fltarr(nf)
   arcdich = strarr(nf)
   arclamps = strarr(nf)
   arcnlamps = intarr(nf)
   arcra = fltarr(nf)
   arcdec = fltarr(nf)
   for f = 0, nf-1 do begin
      if allfiles[f] eq '' then return
      filename = allfiles[f]
      h = headfits(filename)
      lampstr = clip(sxpar(h,'LAMPS'))


      lamps = fix(strsplit(lampstr,',',/extract))
      fixlampstr, lampstr, badarclist
      arclamps[f] = lampstr
      Hg = lamps[0]
      Neon = lamps[1]
      Ar = lamps[2]
      Cd = lamps[3] ; is this actually the order?
      Zn = lamps[4]
      hal = lamps[5]


      exptime = sxpar(h,'EXPTIME')
      if camchar eq 'b' and Hg + Cd + Zn ge 2 and exptime lt 10. then isarc[f] = 1
      if camchar eq 'r' and Neon + Ar ge 1 and exptime lt 10. then isarc[f] = 1
      if camchar eq 'b' then arcnlamps[f] = Hg + Cd + Zn
      if camchar eq 'r' then arcnlamps[f] = Hg + Neon + Ar
      arcairmass[f] = sxpar(h,'AIRMASS')
      arcgr[f] = clip(sxpar(h,'GRNAME'))
      arcbinning[f] = clip(sxpar(h,'BINNING'))
      arcwindow[f] = clip(sxpar(h,'PANE'))
      if camchar eq 'r' then arccwav[f] = sxpar(h,'WAVELEN') 
      if camchar eq 'b' then begin
        if arcgr[f] eq '400/3400' then arccwav[f] = 3600. ; 3505??
        if arcgr[f] eq '600/4000' then arccwav[f] = 4305.
      endif      
      arcdich[f] = clip(sxpar(h,'DICHNAME'))
      arcra[f] = ten(sxpar(h,'RA'))*15.
      arcdec[f] = ten(sxpar(h,'DEC'))
   endfor
   w = where(isarc, ct)
   if ct eq 0 then begin
      print, 'Cannot find any decent arcs.  You will need to extract the spectra manually.'
      return
   endif

   ; anytime two arcs were taken in succession with identical parameters, 
   ; remove the first one - it might not be warmed up
   keeparc = intarr(n_elements(w))+1
   for a = 1, n_elements(w)-1 do begin
      i = w[a]
      if isarc[i-1] eq 0 then continue
      if arcgr[i] eq arcgr[i-1] and arcbinning[i] eq arcbinning[i-1] and arcwindow[i] eq arcwindow[i-1] and $
         arccwav[i] eq arccwav[i-1] and arcnlamps[i] eq arcnlamps[i-1] then $
         keeparc[a-1] = 0 ; should be safe
   endfor
   w = w[where(keeparc)]
   arcfiles = allfiles[w]
   arcairmass = arcairmass[w]
   arcgr = arcgr[w]
   arcbinning = arcbinning[w]
   arcwindow = arcwindow[w]
   arccwav = arccwav[w]
   arcnlamps = arcnlamps[w]
   arcra = arcra[w]    ; az and el would actually be more useful
   arcdec = arcdec[w]

   ; ----- If overwrite mode, destroy old arc solutions -----

   if overwrite then begin
      arcsols = findfile(spworkingdir+'*'+camchar+wildchar+chipchar+'.sol', count=nsol)
      for a = 0, nsol-1 do $
        file_delete, arcsols[a]
   endif

   ; -------- Now do the extraction ------------

   
   if file_test('tracepos.txt') gt 0 then begin
      ntr = countlines('tracepos.txt')
      openr, 9, 'tracepos.txt'
      inline = ''
      tr_files = strarr(ntr)
      tr_ymin = -1 + fltarr(ntr)
      tr_ymax = -1 + fltarr(ntr)
      tr_yfix = -1 + fltarr(ntr)
      tr_bakmin = -1 + fltarr(ntr)
      tr_bakmax = -1 + fltarr(ntr)
      tr_tracefile = strarr(ntr)
      tt = 0
      for t = 0, ntr-1 do begin
          readf, 9, inline
          hashpos = strpos(inline,'#')
          if hashpos ge 0 then inline = strmid(inline,hashpos)
          inline = clip(inline)
          if strlen(inline) eq 0 then continue
          inarr = strsplit(inline,/extract)
          if n_elements(inarr) lt 2 then continue
          tr_files[tt] = inarr[0]
          yinfo = inarr[1]
          if strpos(yinfo,'-') eq 0 then begin
            tr_yfix[tt] = float(yinfo)
          endif else begin
            yarr = strsplit(yinfo,'-',/extract)
            tr_ymin[tt] = float(yarr[0])
            tr_ymax[tt] = float(yarr[1])
          endelse
          if n_elements(inarr) gt 2 then begin
            barr = strsplit(inarr[2],'-',/extract)
            tr_bakmin[tt] = float(barr[0])
            tr_bakmax[tt] = float(barr[1])
          endif
          if n_elements(inarr) gt 3 then begin
             tr_tracefile[tt] = strsplit(inarr[3],/extract)
          endif
          tt += 1
      endfor
      if tt ge 1 then tr_files = tr_files[0:tt-1] ; ignore the rest
      close, 9
   endif

   files = choosefiles(camchar+wildchar+chipchar+'.fits',spworkingdir+'szfp',spworkingdir+'sfp')


   if camchar eq 'r' then grkey = 'GRATNAME' else grkey = 'GRISNAME'
   if camchar eq 'r' then filtkey = 'REDFILT' else filtkey = 'BLUFILT'

   keywords = ['FILTER','EXPTIME','OBJECT','TARGNAME','ROTPOSN','DICHNAME','GRNAME',grkey,filtkey,'SLITNAME','INSTRUME','TELESCOP','ELAPTIME','UTC','MJD-OBS','AZ','EL','RA','DEC','EQUINOX','AIRMASS','SUNELEV','MOONELEV','MOONDIST','SLITPA','DATE_BEG','DATE_END','BINNING','CAMERA','FLATFLE','FLATTYPE']
   if camchar eq 'r' then keywords  = [keywords, 'WAVELEN']

   for f = 0, n_elements(files)-1 do begin
      if n_elements(fileseq) gt 0 then if total(filenum(files[f]) eq fileseq) eq 0 then continue
      filename = files[f]
      h = headfits(filename)
      nx = sxpar(h,'NAXIS1')
      ;specbin = sxpar(h,'XBIN')
      spatbin = sxpar(h,'YBIN')
      gr = clip(sxpar(h,'GRNAME'))
      binning = clip(sxpar(h,'BINNING'))
      window = clip(sxpar(h,'PANE'))
      airmass = sxpar(h,'AIRMASS')
      dispersion = sxpar(h,'DISPERSN')
      target = sxpar(h,'TARGNAME')
      ra = ten(sxpar(h,'RA'))*15.
      dec = ten(sxpar(h,'DEC'))
      if camchar eq 'r' then cwav = sxpar(h,'WAVELEN') ; seems that this correction is necessary?
      if camchar eq 'b' then begin
        if gr eq '400/3400' then cwav = 3600. ; 3505??
        if gr eq '600/4000' then cwav = 4305.
      endif

      outfile = fileappend(repstr(filename,'.fits', '.spec'),'e')

      if overwrite eq 0 and file_test(outfile) then continue

      print, 'Extracting ', removepath(filename), ' ('+clip(target)+')'

      if camchar eq 'r' and spatbin eq 2 then normaltracepos = 76  ; different depending on binning, probably depends on the crop...
      if camchar eq 'r' and spatbin eq 1 then normaltracepos = 92
      if camchar eq 'b' then normaltracepos = 144
      defaultysearch = normaltracepos*[0.6,1.4]/spatbin

      if n_elements(ycent) gt 0 then delvarx, ycent
      if n_elements(bakdist) gt 0 then delvarx, bakdist
      if n_elements(bakwidth) gt 0 then delvarx, bakwidth
      if n_elements(tracefile) gt 0 then delvarx, tracefile
      ysearch = defaultysearch
      if n_elements(tr_files) gt 0 then begin
         t = where(tr_files eq removepath(files[f]), ct)
         print, tr_files , removepath(files[f]), ct
         if ct gt 0 then begin
            t = t[0]
            if tr_yfix[t] ge 0 then ycent = tr_yfix[t]
            if tr_ymin[t] ge 0 and tr_ymax[t] gt 0 then begin
               ysearch = [tr_ymin[t], tr_ymax[t]] 
            endif 
            if tr_bakmin[t] ge 0 and tr_bakmax[t] ge 0 then begin 
               bakdist = (tr_bakmin[t]+tr_bakmax[t])/2.
               bakwidth = (tr_bakmax[t]-tr_bakmin[t])
            endif
            if tr_tracefile[t] ne '' then   tracefile = tr_tracefile[t] 
         endif
      endif

      skylist = refdatapath + '/' + 'skylines.dat'

      if n_elements(arcfile) gt 0 then delvarx, arcfile
      if n_elements(skyfile) gt 0 then delvarx, skyfile
      if n_elements(initsol) gt 0 then delvarx, initsol
      if n_elements(extractarcfile) gt 0 then delvarx, extractarcfile

      if camchar eq 'r' then uncrefwav = 0.3
      if camchar eq 'b' then uncrefwav = 0.3
 

      ; find a matching arc
      matches = where(arcgr eq gr and arcbinning eq binning and arcwindow eq window and abs(arccwav-cwav) lt 10., ct)
      if ct eq 0 then begin
         print, 'WARNING: No arcs were taken in the same configuration as ', removepath(filename), ' ('+clip(sxpar(h,'TARGNAME'))+')'
         print, '     gr='+gr, ' binning='+repstr(binning,',','x'), ' window='+window, ' cwav='+clip(cwav)
         print, 'Cannot wavelength calibrate.  Skipping...'
         continue
      endif else begin
         archere = 0
         w = where(abs(arcra[matches]-ra) lt 0.1 and abs(arcdec[matches]-dec) lt 0.1, ct)
         if ct gt 0 then begin
            w = max(w)
            archere = 1
         endif else begin
            ; get an arc not near horizon, if we can.
            w = where(arcairmass[matches] lt 3.0, ct)
            if ct gt 0 then matches = matches[w]

            ; as many lamps as possible.
            w = where(arcnlamps[matches] eq max(arcnlamps[matches]), ct) 
            matches = matches[w]
 
            ; at this exact central wavelength, if possible
            w = where(abs(arccwav[matches]-cwav) lt 0.1, ct)
            if ct gt 0 then matches = matches[w]

            ; as close to this elevation as possible.
            w = where(abs(arcairmass[matches] - airmass) eq min(abs(arcairmass[matches]-airmass)))
            w = max(w)
         endelse

         arcfile = arcfiles[matches[w]]
         arcair = arcairmass[matches[w]]
         solfile = repstr(arcfile,'.fits','.sol')
         solexists = file_test(solfile)

         if solexists eq 0 or archere then begin
 
            arch = headfits(arcfile)
                                      ;'~/research/redux/arc/lrisall.list'
            arclistfile = refdatapath+'/'+'lrisarcall.list'   
            lampstr = clip(sxpar(arch,'LAMPS'))
            fixlampstr, lampstr, badarclist
            lamps = strsplit(lampstr,',',/extract)
            lampnames = ['']
            if fix(lamps[0]) then lampnames = [lampnames, 'Hg']
            if fix(lamps[1]) then lampnames = [lampnames, 'Ne']
            if fix(lamps[2]) then lampnames = [lampnames, 'Ar']
            if fix(lamps[3]) then lampnames = [lampnames, 'Cd']
            if fix(lamps[4]) then lampnames = [lampnames, 'Zn']
            if n_elements(lampnames) gt 1 then lampnames = lampnames[1:*]
      
            readcol, arclistfile, wav, name, format='f,a'
            keep = intarr(n_elements(name)) + 1
            for l = 0, n_elements(name)-1 do begin
            linenames = strsplit(name[l],'+',/extract) ; signify blends with +
            for b = 0, n_elements(linenames)-1 do begin
               if total(linenames[b] eq lampnames) eq 0 then keep[l] = 0
            endfor
            endfor
            wav = wav[where(keep)]
            arclist = wav

            brightarclistfile = refdatapath+'/'+'lrisarcbright.list'
            readcol, brightarclistfile, wav, name, format='f,a'
            keep = intarr(n_elements(name)) + 1
            for l = 0, n_elements(name)-1 do begin
               linenames = strsplit(name[l],'+',/extract) ; signify blends with +
               for b = 0, n_elements(linenames)-1 do begin
                  if total(linenames[b] eq lampnames) eq 0 then keep[l] = 0
               endfor
            endfor
            wav = wav[where(keep)]
            brightarclist = wav

            lampstr = ''
            for l = 0, n_elements(lampnames)-1 do lampstr += lampnames[l]
            print, '  Calibrating with arc ', arcfile, ' (',lampstr,')'

            if archere eq 0 then begin
               ; Use an arc at some other position.  Extract a generic arc spectrum and save the solution.

               arcdata = mrdfits(arcfile, 0, /silent)
               arcdata = arcdata[*,defaultysearch[0]:defaultysearch[1]]
               farc = total(arcdata,2)
               nx = n_elements(farc)

               checkplotfile = repstr(arcfile,'.fits','.wavsol.ps')
               arcsol = solvearcspec(farc, wav, refwav=arclist, brightrefwav=brightarclist, $
                  cwav=cwav, unccwav=200., dispersion=dispersion, uncdisp=0.075*dispersion, $
                  maxsep=nx*dispersion/2., uncrefwav=uncrefwav, initsol=initsol, nbrightline=20, $
                  nmatch=nmatch, checkplotfile=checkplotfile, maxorder=7)
               checkplotfile = ''

              openw, 9, solfile
              for i = 0, n_elements(arcsol)-1 do printf, 9, arcsol[i]
              close, 9  ; could print out a header (and even search it), but it's not really necessary now

            endif else begin
               ; Use an arc at THIS position.  (Extraction will be done in extractspec, and there is no sky line refinement.)
               print, '  Arc is at target position; skipping sky-line refinement.'
               extractarcfile = arcfile

               skylist = '' ; don't sky-refine if an arc was taken here.
            endelse


         endif else begin
            ; Arc is already solved; just load a saved solution.

            print, '  Wavelength previously solved; loading '+removepath(solfile)
            if solexists then initsol = grabcolumn(solfile,0)

         endelse
      endelse

      if sxpar(h,'EXPTIME') le 11.0 then skylist = ''  ; don't sky-refine short exposures

      if camchar eq 'b' and sxpar(h,'SUNELEV') gt -11 then skylist = '' ; don't sky-refine twilight exposures.  
                                                                        ;This is conservative, better solar-subtraction will enable improvement.
      if camchar eq 'r' and sxpar(h,'SUNELEV') gt -6  then skylist = ''

      if camchar eq 'r' then begin
         maxnudge = 220.
      endif
      if camchar eq 'b' then maxnudge = 3. ; actually 10 more than this, this is for the search

      specextract, filename, ycent, $ 
      ysearch=ysearch, bakdist=bakdist, bakwidth=bakwidth, $
      arcfile=extractarcfile, arclist=arclist, brightarclist=brightarclist, skyfile=1, skylist=skylist, $
      initsol=initsol, cwav=cwav, dispersion=dispersion, unccwav=200., uncdisp=0.075*dispersion, uncrefwav=uncrefwav, $
      outfile=outfile, keywords=keywords, maxorder=7, success=success, maxnudge=maxnudge

      if success lt 0 and success ne -4 then extractfail = [extractfail, filename]
      if success eq -4 then wavsolfail = [wavsolfail, filename]
      if success eq 2 then wavsolwarn = [wavsolwarn, filename]

  endfor

end

; -------------------------

pro autolrisfluxcal, camera=camera, chip=chip

   common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite, fileseq, filterstr, gratingstr, targetstr
   common lrisfail, flatfail, catastrofail, relastrofail, fullastrofail, extractfail, wavsolfail, wavsolwarn, fluxcalfail
   common lrisconfig, lrisautopath, refdatapath, defaultspath

   common standardspos, standardslist, standardsra, standardsdec

   standardpath = refdatapath

   camchar = strmid(camera,0,1)
   chipchar = strmid(chip,0,1)

   files = choosefiles(camchar+wildchar+chipchar+'.spec',spworkingdir+'eszfp',spworkingdir+'esfp')
   nf = n_elements(files)
   if nf eq 1 and files[0] eq '' then return

   names = strarr(nf)
   airmasses = fltarr(nf)
   gratings = strarr(nf)
   dichroics = strarr(nf)
   filters = strarr(nf)
   isfluxcal = intarr(nf)
   slits = strarr(nf)
   nwavs = intarr(nf)

   ; identify flux calibrators
   nstdobs = 0
   for f = 0, n_elements(files)-1 do begin
      if files[f] eq '' then continue
      if n_elements(fileseq) gt 0 then if total(filenum(files[f]) eq fileseq) eq 0 then continue
      filename = files[f]
      h = headspec(filename)

      target = sxpar(h,'TARGNAME')
      object = sxpar(h,'OBJECT')
      rastr = sxpar(h,'RA')
      decstr = sxpar(h,'DEC')
      if strpos(rastr,':') gt 0 then ra = 15*ten(rastr) else ra = float(rastr)
      if strpos(decstr,':') gt 0 then dec = 15*ten(decstr) else ra = float(decstr)

      airmasses[f] = sxpar(h,'AIRMASS')
      gratings[f] = sxpar(h,'GRNAME')
      dichroics[f] = sxpar(h,'DICHNAME')  
      filters[f] = sxpar(h,'FILTER')
      slits[f] = sxpar(h,'SLITNAME')
      sunel = sxpar(h,'SUNELEV')
      nwav = sxpar(h,'NWAV',count=ctnwav)      
      if ctnwav gt 0 and nwav eq 0 then continue
      nsaturate = sxpar(h,'NSATURAT')
      if nsaturate gt 3 then continue
      name = target

      readstandard, target, /check, result=result, path=standardpath+'/'
      if result eq 0 then begin
         readstandard, object, /check, result=result, path=standardpath+'/'
         if result eq 1 then name = object
      endif

      isfluxcal[f] = result
      names[f] = name ; targname generally, but object name if it matched a standard

      ; should also check wavelen, binning?

      xx = 0
      if xx eq 1 and result eq 1 then begin  ; need to only print this when necessary
         if sunel gt -18.0 then sunelstr = fpr(sunel,3.1) else sunelstr = 'night'
         if nstdobs eq 0 then print, '   Standard-star observations are:'
         if nstdobs eq 0 then print, clip('filename', 25), ' ', clip('standard',16), ' ', clip('air',4), ' ', clip('grating',9), ' ', clip('dich',4), ' ', clip('filt',5), ' ', clip('sunel',5)
         print, clip(removepath(files[f]),25), ' ', clip(names[f],16),  ' ', fpr(airmasses[f],1.2), ' ', clip(gratings[f],9), ' ', clip(dichroics[f],4), ' ', clip(filters[f],5), ' ', clip(sunelstr,5)
         nstdobs += 1
      endif
   endfor

   wfluxcal = where(isfluxcal, complement=wscience, ct)
   if ct eq 0 then begin
      print, 'Found no LRIS-'+strupcase(camchar)+' standard-star observations!'
      print, '  If standards were taken, check that ', standardpath, ' is a valid directory containing the standard-star library.'
      return
   endif
   if ct eq nf then begin
      print, 'Found no science observations!'
      return
   endif

   nsci = n_elements(wscience)

   ; read them in
   ; convert to transmission plot

   for s = 0, nsci-1 do begin
      wsci = wscience[s]
      filename = files[wsci]

      outfile = fileappend(filename,'c')

      if file_test(outfile) and overwrite eq 0 then continue
      
      matchstandards = where(filters[wfluxcal] eq filters[wsci] and gratings[wfluxcal] eq gratings[wsci] and $ 
                             dichroics[wfluxcal] eq dichroics[wsci], nmatchstd)
      if nmatchstd eq 0 then begin
         print, 'No matching standard for ', removepath(files[wsci]), ' ('+gratings[wsci]+', '+filters[wsci]+', D'+clip(dichroics[wsci])+')'
         fluxcalfail = [fluxcalfail, files[wsci]]
         continue
      endif

      sci = readspec(files[wsci], h, /timenorm, cutoff=3000.)
      wav = sci.wav
      nwav = sxpar(h,'NWAV',count=ctnwav)
      if ctnwav gt 0 and nwav eq 0 then begin  
         print, filename, ' has no wavelength calibration.'
         continue
      endif

      sciair = airmasses[wsci]
      print, 'Flux-calibrating ', removepath(files[wsci]), ' (',names[wsci],')  airmass=', fpr(sciair,1.4)
      stdair = airmasses[wfluxcal[matchstandards]]

         hydrogen = [6562.80, 4861.34, 4340.48, 4101.75, 3970.09, 3889.07, $
                     3835.40, 3797.92, 3770.65, 3750.17, 3734.39, 3721.96, $
                     3711.99, 3703.88, 3697.17, 3691.58, 3686.85, 3682.83]
         hv = 0.003 + (hydrogen gt 3900)*0.007
         calcium = [3968.5, 3933.7]
         cav = 0.003 + fltarr(n_elements(calcium))
         telluric = [6890., 7270., 7670., 8150., 9060., 9480.]
         tellv = [0.008, 0.02, 0.011, 0.004, 0.013, 0.022]
         lines = [hydrogen, calcium, telluric]
         linev = [hv, cav, tellv]

      if nmatchstd eq 1 then begin
         if abs(sciair-stdair[0]) gt 0.1 then begin
            print, '     Warning: Only one standard-star observation in this configuration.'
            print, '     Airmass:  science ('+fpr(sciair,2.2)+'), standard ('+fpr(stdair[0],2.2)+')'
         endif
         wstd = wfluxcal[matchstandards[0]]

         ; 'std' is the observation of the standard.
         ; 'ref' is the reference of the standard.

         print, '     Calibrator  ', removepath(files[wstd]), ' ', names[wstd], ' @ air=', fpr(airmasses[wstd],1.3)

         std = readspec(files[wstd], hstd, /timenorm, cutoff=3000.)
         stdflux = interpol2(std.flux, std.wav, wav)

         removelines, wav, stdflux, stdflux_nl, lines=lines, vel=linev, /interpol

         readstandard, names[wstd], wav, refflux, unit='flambda/ang', path=standardpath+'/'
         removelines, wav, refflux, refflux_nl, lines=lines, vel=linev, /interpol
 
         response = stdflux_nl / refflux_nl

         sxaddpar, h, 'FLUXCAL1', removepath(files[wstd])
         sxaddpar, h, 'FLUXNAM1', names[wstd]
         sxaddpar, h, 'FLUXAIR1', airmasses[wstd]

      endif else begin
        if nmatchstd eq 2 then begin
            wlo = wfluxcal[matchstandards[(where(stdair eq min(stdair))) [0]]]
            whi = wfluxcal[matchstandards[(where(stdair eq max(stdair))) [0]]]
        endif

        if nmatchstd ge 3 then begin
           below = where(stdair lt sciair, ctbelow)
           above = where(stdair gt sciair, ctabove)
           ; try to bracket the source if possible
           if ctbelow ge 1 and ctabove gt 1 then begin
              wlo = wfluxcal[matchstandards[below[(where(stdair[below] eq max(stdair[below]))) [0]]]]
              whi = wfluxcal[matchstandards[above[(where(stdair[above] eq min(stdair[above]))) [0]]]]
           endif else begin
              wlo = wfluxcal[matchstandards[(where(stdair eq max(stdair))) [0]]]
              whi = wfluxcal[matchstandards[(where(stdair eq min(stdair))) [0]]]
           endelse
        endif
 
        print, '     Calibrators ', removepath(files[wlo]), ' ', names[wlo], ' @ air='+fpr(airmasses[wlo],1.3)
        print, '                 ', removepath(files[whi]), ' ', names[whi], ' @ air='+fpr(airmasses[whi],1.3)

        stdlo = readspec(files[wlo], hstdlo, /timenorm, cutoff=3000.)
        stdhi = readspec(files[whi], hstdhi, /timenorm, cutoff=3000.)    

        stdloflux = interpol2(stdlo.flux, stdlo.wav, wav)
        stdhiflux = interpol2(stdhi.flux, stdhi.wav, wav)

        removelines, wav, stdloflux, stdloflux_nl, lines=lines, vel=linev, /interpolate
        removelines, wav, stdhiflux, stdhiflux_nl, lines=lines, vel=linev, /interpolate

        readstandard, names[wlo], wav, refloflux, unit='flambda/ang', path=standardpath+'/'
        readstandard, names[whi], wav, refhiflux, unit='flambda/ang', path=standardpath+'/'
 
        removelines, wav, refloflux, refloflux_nl, lines=lines, vel=linev, /interpolate
        removelines, wav, refhiflux, refhiflux_nl, lines=lines, vel=linev, /interpolate

        responselo = stdloflux_nl / refloflux_nl
        responsehi = stdhiflux_nl / refhiflux_nl

         ; anything less than 1/10000 of max is read noise or dead, and will mess up log fit - remove these values.
        bad = where(responselo lt max(responselo)*1e-4, ct, complement=good)
        if ct gt 0 then responselo[bad] = min(responselo[good]) 
        bad = where(responsehi lt max(responsehi)*1e-4, ct, complement=good)
        if ct gt 0 then responsehi[bad] = min(responsehi[good])

        normlo = percentile75(responselo)
        normhi = percentile75(responsehi)  
        norm = (normlo + normhi)/2.    ; could add some sort of anticipated air transmission adjustment...?
        responselo *= norm/normlo
        responsehi *= norm/normhi

        r = (airmasses[wsci]-airmasses[wlo])/(airmasses[whi] - airmasses[wlo]) ; airmass ratio: (sci-lo)/(hi-lo)
        off = min(abs(airmasses[wsci]-[airmasses[whi],airmasses[wlo]]))        ; offset from nearest:  abs(sci-nearest)

        if abs(r) lt 3. and (airmasses[whi]-airmasses[wlo] gt 0.08) then begin
          response = responselo + (responsehi-responselo)*r  ; Is this actually a reasonable thing to do?
        endif else begin
          if r lt 0 then response = responselo
          if r gt 1 then response = responsehi
          if r ge 0 and r le 1 then response = (responselo + responsehi)/2.
        endelse

        if abs(r) gt 2. and off gt 0.15 then begin
          print, '    Warning: Object airmass is well outside range of standards:'
          stdairstr = ''
          for ss = 0, n_elements(stdair)-1 do stdairstr += fpr(stdair[ss],2.2) + ' '
          print, '    Airmass:  science ('+fpr(sciair,2.2)+'), standards ('+stdairstr+')'
        endif

        sxaddpar, h, 'FLUXCAL1', files[wlo]
        sxaddpar, h, 'FLUXCAL2', files[whi]
        sxaddpar, h, 'FLUXAIR1', airmasses[wlo]
        sxaddpar, h, 'FLUXAIR2', airmasses[whi]
        sxaddpar, h, 'FLUXNAM1', names[wlo]
        sxaddpar, h, 'FLUXNAM2', names[whi]
      endelse

      bad = where(response lt median(response)*1e-4, ct, complement=good) ; just in case any negatives snuck through
      if ct gt 0 then response[bad] = min(response[good])
      par = [0,0,0]
      n = n_elements(wav)
      c = median(wav)
      r = median(response[5:n_elements(wav)-5])
      for o = 3, 16 do begin
        par = polyfitg((wav[5:n-5]-c)/1e3, alog((response[5:n-5]/r)), o, guess=[par, 0]) ; in theory I shouldn't be using bad regions for this
      endfor
      modelresponse = r*exp(poly((wav-c)/1e3, par))

      wsignal = where(response gt 0.01*max(response))
      outspec = sci[wsignal]

      outspec.flux =  (sci[wsignal].flux) / modelresponse[wsignal]
      outspec.sky  =  (sci[wsignal].sky)  / modelresponse[wsignal]
      outspec.unc  =  (sci[wsignal].unc)  / modelresponse[wsignal]
      outspec.pix  = sci[wsignal].pix
      outspec.resp = modelresponse[wsignal]

      ni = 1
      if nmatchstd ge 2 then ni = 2
      if ni eq 1 then xsize = 600 else xsize = 1000

      window, 0, xsize=xsize, ysize=800
      !p.multi = [2,1,2]

      for i = 1, ni do begin
        if nmatchstd eq 1 then begin
          xl = 0.1
          xr = 0.98
          resp = response
          name = names[wstd]
          file = files[wstd]
        endif
        if nmatchstd ge 2 then begin
          ; for plotting purposes
          if i eq 1 then begin
            xl = 0.1
            xr = 0.49
            std = stdlo
            stdflux_nl = stdloflux_nl
            refflux_nl = refloflux_nl
            wstd = wlo
            resp = responselo
            name = names[wlo]
            file = files[wlo]
          endif else begin
            xl = 0.59
            xr = 0.98
            std = stdhi
            stdflux_nl = stdhiflux_nl
            refflux_nl = refhiflux_nl
            wstd = whi
            resp = responsehi
            name = names[whi]
            file = files[whi]
          endelse
        endif

        xrange = [min(wav), max(wav)]

        ; could iterate this, doing a psopen the second time?

        !p.position = [xl, 0.8, xr, 0.99]
        plot, [0], [0], xrange=xrange, yrange=[0, 1.02*max([refflux_nl])], /xstyle, /ystyle
        if n_elements(origwav) eq 0 then delvarx, origwav
        readstandard, names[wstd], origwav, origflux, unit='flambda/ang', path=standardpath+'/'
        oplot, origwav, origflux, color=511  ;, psym=1, size=3
        oplot, wav, refflux_nl
        xyouts, xr-(xr-xl)*0.35, !p.position[3]-0.02, /norm, name+' (database)'

        !p.position = [xl, 0.5, xr, 0.78]
        plot, [0], [0], yrange=[0, max([std.flux,stdflux_nl])], xrange=xrange, /xstyle
        oplot, std.wav, std.flux, color=511 ; noninterpolated
        oplot, wav, stdflux_nl  ; interpolated
        xyouts, xr-(xr-xl)*0.65, !p.position[3]-0.02, /norm, name+' ('+removepath(file)+')'
        !p.multi = [2,1,2]

        !p.position = [xl, 0.20, xr, 0.48]
        plot, [0], [0], xrange=xrange, yrange=[0, max(response)*1.05], /xsty, /ysty;, /ylog
        if nmatchstd gt 1 then oplot, wav, resp, color=391482905
        oplot, wav, response
        oplot, wav, modelresponse, color=124890
        xyouts, xl+(xr-xl)*0.05, !p.position[3]-0.02, /norm, 'Solved response'

        if i eq ni then begin
        !p.position = [xl, 0.03, xr, 0.18]
        wn = where(outspec.wav gt 3400 and outspec.wav lt 10000 and $
                   outspec.wav lt max(outspec.wav)*0.95 and outspec.wav gt min(outspec.wav) * 1.05)
        plot, outspec.wav, outspec.flux, /xstyle, xrange=xrange, psym=10, yrange=[0,max(outspec[wn].flux)]
        xyouts, xr-(xr-xl)*0.93, !p.position[3]-0.02, /norm, 'Calibrated '+clip(sxpar(h,'TARGNAME'))+' spectrum ( -> '+removepath(outfile)+')'
        endif
        !p.multi = [2,1,2]

      endfor
      !p.multi = 0
      !p.position = 0

      writespec, outfile, outspec, h

   endfor


end

; -----------------------------

pro autolriscombine, camera=camera

   common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite, fileseq, filterstr, gratingstr, targetstr

   common lrisconfig, lrisautopath, refdatapath, defaultspath

   standardpath = refdatapath

   ; automatic combination of 1D spectra of a common source in a common setting.
   ; assumes the central wavelength wasn't wildly adjusted between exposures.
 
   camchar = strmid(camera,0,1)
   chipchar = 'r'

   files = choosefiles(camchar+wildchar+chipchar+'.spec',spworkingdir+'ceszfp',spworkingdir+'cesfp')
   nf = n_elements(files)

   if nf eq 1 and files[0] eq '' then return

   targets = strarr(nf)
   igrs = strarr(nf)
   ; worry about different cwavs?  probably not since it can shift slightly.

   for f = 0, n_elements(files)-1 do begin
     filename = files[f]
     if n_elements(fileseq) gt 0 then if total(filenum(files[f]) eq fileseq) eq 0 then continue
     h = headspec(filename)
     target = clip(sxpar(h,'TARGNAME'))
     if strpos(target,'_S') gt 0 then begin  ; strip the offset star from the target info
        ustriptarget = strmid(target,0,strpos(target,'_', /reverse_search))
        target = ustriptarget
     endif
     targets[f] = target
     inst = clip(sxpar(h,'INSTRUME'))
     gr = clip(sxpar(h,'GRNAME'))
     incamchar = 'u' ; unknown
     if inst eq 'LRIS' then incamchar = 'r'
     if inst eq 'LRISBLUE' then incamchar = 'b'
     igrs[f] = incamchar + gr ; this is sort of obsolete now, but allows this to work with camera='?'
   endfor

   targetlist = unique(targets)
   ntargets = n_elements(targetlist)

   for t = 0, ntargets-1 do begin 
      target = targetlist[t]
      readstandard, target, /check, result=isstandard, path=standardpath+'/'
      if isstandard then continue ; don't stack observations of flux standards

      wtarget = where(targets eq target, nexp)
      
      ; determine intra-color normalization

      igrlist = unique(igrs[wtarget])
      ngr = n_elements(igrlist) 

      minwav = fltarr(ngr)
      maxwav = fltarr(ngr)
      for g = 0, ngr-1 do begin
         igr = igrlist[g]
         cam = strmid(igr,0,1)
         gr  = strmid(igr,1)
               
         outfile =  spworkingdir+'vc'+repstr(target,'_','')+'_'+repstr(igr,'/','-')+'.spec'
         if file_test(outfile) and overwrite eq 0 then continue
            
         wtargigr = where(targets eq target and igrs eq igr, ngexp)

         fluxnormval = fltarr(ngexp)
         totexp = 0.
         print, ' Combining ', format='($,A)'
         for i = 0, ngexp-1 do begin
            f = wtargigr[i]
            print, removepath(files[f])+' ', format='($,A)'
            s = readspec(files[f], h)
            inexp = sxpar(h,'ELAPTIME')

            if i eq 0 then sums = s
            if i gt 0 then begin
               sums = addspec(sums, s, w1=totexp, w2=inexp)
            endif
            totexp += inexp

            ; header info 
            if i eq 0 then begin
               outh = h  
            endif

            if ngexp gt 1 then begin
               sxaddpar, outh, 'APRAD'+clip(i), sxpar(h,'APRAD')
               sxaddpar, outh, 'YEXTRAC'+clip(i), sxpar(h,'YEXTRACT')
            endif

            sxaddpar, outh, 'INFILE'+clip(i), files[f]

         endfor
         sxaddpar, outh, 'DATE_END', sxpar(h, 'DATE_END')
         sxaddpar, outh, 'AIRMASS', mean([sxpar(h, 'AIRMASS') , sxpar(outh,'AIRMASS')])
         sxaddpar, outh, 'AZ', mean([sxpar(h, 'AZ') , sxpar(outh,'AZ')])
         sxaddpar, outh, 'EL', mean([sxpar(h, 'EL') , sxpar(outh,'EL')])
         sxaddpar, outh, 'TOTALEXP', totexp
         print
         print, '           -> ', removepath(outfile)

         writespec, outfile, sums, outh

      endfor

   endfor

end



pro autolrisconnect

   common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite, fileseq, filterstr, gratingstr, targetstr

   files = findfile(spworkingdir+'vc*.spec')
   files = files[sort(files)]
   nf = n_elements(files)

   if nf eq 1 and files[0] eq '' then return

   targets = strarr(nf)

   for f = 0, n_elements(files)-1 do begin
     filename = files[f]
     h = headspec(filename)
     target = clip(sxpar(h,'TARGNAME'))
     if strpos(target,'_') gt 0 then begin
        posstarget = strmid(target,0,strpos(target,'_'))
        if posstarget ne 'HD' and posstarget ne 'BD' and posstarget ne 'PTF' then target = posstarget
     endif
     targets[f] = target 
   endfor

   targetlist = unique(targets)
   ntargets = n_elements(targetlist)

   ncon = 0

   for t = 0, ntargets-1 do begin 
      target = targetlist[t]

      outfile = spfinaldir+'lris'+target+'.spec'
      if file_test(outfile) and overwrite eq 0 then continue

      if ncon eq 0 then begin
          window, 0
          !p.multi = [0,2,4]
      endif
      ncon += 1

      wtarget = where(targets eq target, nexp)

      vcfiles = files[wtarget]

      ; This in theory allows a spectrum in a bunch of pieces to be put together, 
      ;  but in reality it's only relevant for two

      rescales = [1. ] ; not used, could re-renormalize

      print, 'Connecting', format='($,A)'
      saveh = ['']
      cc = 0
      nsg = n_elements(vcfiles)
      if nsg eq 1 then print, ' (only one side present - trivial connection)', format='($,A)'

      for gg = 0, nsg-1 do begin
         s = readspec(vcfiles[gg], h)
         print, ' '+removepath(vcfiles[gg]), format='($,A)'

         if gg eq 0 then begin 
            comb = s
            outh = h  
         endif else begin
            if max(comb.wav) gt max(s.wav) then begin
               store = s
               s = comb   ; red must be first below
               comb = store
            endif
            overlapcenter = mean([max(comb.wav), min(s.wav)])
            overlapwidth =      (max(comb.wav) - min(s.wav))/3.3
            nooverlap = 0
            if overlapwidth lt 10. then begin
               print
               if overlapwidth gt 0 then print, '   Insufficient blue/red overlap!  (Only ', fpr(overlapwidth,2.1), ' Angstrom)', format='($,A)' $
                                    else print, '   No red/blue overlap.', format='($,A)'
               nooverlap = 1
            endif

           if nooverlap eq 0 then begin
            if overlapwidth lt 50. then begin
               print
               print, '   Caution: marginal blue/red overlap.'
            endif
            cdwav = comb.wav - shift(comb.wav, 1)
            cdwav[0] = cdwav[1]            
            cozone = where(comb.wav gt overlapcenter-overlapwidth and comb.wav lt overlapcenter+overlapwidth and finite(comb.flux), ct)
            if ct eq 0 then begin
               conorm = -1
            endif else begin
               conorm = total(comb[cozone].flux * cdwav[cozone]) / total(cdwav[cozone])
             endelse

            dwav = s.wav - shift(s.wav,1)
            dwav[0] = dwav[1]
            ozone = where(s.wav gt overlapcenter-overlapwidth and s.wav lt overlapcenter+overlapwidth and finite(s.flux), ct)
            if ct eq 0 then begin
               onorm = -1
            endif else begin
               onorm = total(s[ozone].flux * dwav[ozone]) / total(dwav[ozone])
            endelse


            if conorm gt 0 and onorm gt 0 then begin
              rescale = conorm/onorm

              print, '(x'+clip(clip(rescale,5))+') ', format='($,A)'
              plot, [0], [0], xrange=overlapcenter+[-1,1]*overlapwidth, yrange=[0,1.05*max([comb[cozone].flux,s[ozone].flux])], /xsty, /ysty, charsize=2
              oplot, comb[cozone].wav, comb[cozone].flux, color=rgb(70,130,255), psym=10
                            ;, yrange=[0, 5*median([comb[ozone].flux,s[ozone].flux])]
              oplot, s[ozone].wav, s[ozone].flux, color=rgb(64,0,0), psym=10
              s.flux *= rescale
              if tag_exist(s, 'sky') then s.sky *= rescale
              if tag_exist(s, 'unc') then s.unc *= rescale
              if tag_exist(s, 'resp') then s.resp *= rescale
              oplot, s[ozone].wav, s[ozone].flux, color=rgb(255,0,0), psym=10

            endif else begin
               print, '  Unable to rescale blue/red spectra; jump may be present.'
               rescale = 1
            endelse
            rescales = [rescales, rescale]

           endif

           comb = connectspec(comb, s, attachwav=overlapcenter)

         endelse


         if sxpar(h,'INSTRUME') eq 'LRISBLUE' then sxaddpar, outh, 'BLUEXPTM', sxpar(h,'TOTALEXP')
         if sxpar(h,'INSTRUME') eq 'LRIS'     then sxaddpar, outh, 'REDEXPTM', sxpar(h,'TOTALEXP')

         sxaddpar, outh, 'INCFILE'+clip(gg), vcfiles[gg]
 
         if gg gt 0 then begin
            for c = 0, n_elements(h)-1 do begin
               nam = clip(strmid(h[c],0,strpos(h[c],'=')))
               val = sxpar(h, nam)
               if strmid(nam,0,6) eq 'INFILE' then begin
                  sxaddpar, outh, 'INFILE'+clip(cc), val, after='INCFILE'+clip(gg)
               endif else begin
                  outval = sxpar(outh,nam,count=ct)
                  if ct eq 0 then begin
                     sxaddpar, outh, nam, val
                  endif 
               endelse
            endfor
         endif

      endfor

      sxdelpar, outh, 'GRNAME' ; no longer relevant for combining both sides

      print
      print, '             -> ', removepath(outfile)
 
      writespec, outfile, comb, outh

   endfor

end


; -------------------------
pro autolrisastrometry, camera=camera, chip=chip

   common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite, fileseq, filterstr, gratingstr, targetstr
   common lrisfail, flatfail, catastrofail, relastrofail, fullastrofail, extractfail, wavsolfail, wavsolwarn, fluxcalfail
   
    achip = strmid(chip,0,1)
    prefchar = '2'
    wildcharimg = '?????????????????_img_?'
    zffiles = choosefiles(prefchar+wildcharimg+'.fits',imworkingdir+'zisfp',imworkingdir+'zsfp',$
                                                          imworkingdir+'isfp',imworkingdir+'sfp')
    
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
          refcatfile = strcompress(imworkingdir+targets[t]+'.'+filters[f]+'.cat',/remove_all)

          if file_test(refcatfile) and overwrite eq 0 then continue
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
          print, autoastrocommand+' '+refimagename;+' -upa 2 -q'
          spawn, autoastrocommand+' '+refimagename;+' -upa 2 -q'

          outfile = fileappend(refimagename,'a')
          if file_test(outfile) eq 0 then begin
             catastrofail = [catastrofail, refimagename]
             print, 'WARNING - astrometry on the reference image was unsuccessful!'
             print
          endif else begin
             print, autoastrocommand+' '+outfile+' -n '+refcatfile + ' -x 55000 -q'
             spawn, autoastrocommand+' '+outfile+' -n '+refcatfile + ' -x 55000 -q'
             print
          endelse
      ;   (respecifying px and pa should no longer be necessary with the new splitlrisred and splitlris blue which add astrometry.)
       endfor
    endfor

    if overwrite eq 0 then zffiles = unmatched(zffiles,'a')  ; catalogs should always be the same; only do this now

    ; Use the reference catalog to do a more precise relative astrometric solution
    for f = 0, n_elements(zffiles)-1 do begin
       if zffiles[f] eq '' then continue
       outfile = fileappend(zffiles[f],'a')
       if file_test(outfile) and overwrite eq 0 then continue
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
          print, autoastrocommand+' '+zffiles[f];+' -upa 2 -q'
          spawn, autoastrocommand+' '+zffiles[f];+' -upa 2 -q'
          print
          if file_test(outfile) eq 0 then fullastrofail = [fullastrofail, zffiles[f]] 
       endif else begin
         ; use the short exposure first, but fall back on direct catalog.
          targname = repstr(strtrim(sxpar(h,'TARGNAME'),2),' ', '_') ; spaces cause barfing in filenames
          if targname eq '' then continue
          if strpos(targname,'flat') ge 0 then continue
          refcatfile = strcompress(imworkingdir+targname+'.'+filt+'.cat',/remove_all)
          
          if file_test(refcatfile) then begin
             print, targname
             print, autoastrocommand+' '+zffiles[f]+' -c '+refcatfile;+' -upa 2' + ' -x 55000 -q' ;-upa 0.1
             spawn, autoastrocommand+' '+zffiles[f]+' -c '+refcatfile;+' -upa 2' + ' -x 55000 -q' ;-upa 0.1
          endif else begin
             print, 'No reference catalog '+refcatfile+' exists for this field.'
          endelse
          if file_test(outfile) eq 0 then begin
             print, 'Refined astrometry of ', zffiles[f], ' was not successful.  Trying direct astrometry:'
             print, autoastrocommand+' '+zffiles[f];+' -upa 2 -q'
             spawn, autoastrocommand+' '+zffiles[f];+' -upa 2 -q'

             if file_test(outfile) then relastrofail  = [relastrofail,  zffiles[f]] $
             else fullastrofail = [fullastrofail, zffiles[f]] 
          endif
       endelse

    endfor

end



; ------------------------

pro autolrisphotometry, camera=camera, chip=chip

  ; need to allow doing this on only one chip.

   common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite, fileseq, filterstr, gratingstr, targetstr

   lrisautofieldphot,blue=bl,red=re,chip=ch

end

; -------------------------
pro autolrisstack, camera=camera, chip=chip

   common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite, fileseq, filterstr, gratingstr, targetstr
   common lrisfail, flatfail, catastrofail, relastrofail, fullastrofail, extractfail, wavsolfail, wavsolwarn, fluxcalfail

  if file_test('default.swarp') eq 0 then $
     spawn, swarpcommand+' -d > default.swarp'

  if n_elements(chip) eq '' then chip = '*'
  chipchar = strmid(chip,0,1)
  prefchar = '2'
  wildcharimg = '?????????????????_img_?'
  azffiles = findfile(imworkingdir+'a*'+prefchar+wildcharimg+'.fits')
  realfiles = where(azffiles ne '', ct)
   if ct eq 0 then return ; can't stack if there are no astrometry files
   if ct ge 1 then azffiles = azffiles[realfiles]
   ;chipchar = 'lr'  ; this is a hack to combine both sides

   camver = camera
   if camver eq 'red' then camver = camver + clip(lrisversion)
   if file_test(imworkingdir+'autophotsummaryflux.txt') then begin
      readcol,imworkingdir+'autophotsummaryflux.txt',pfile,pexp,pfilt,pair,dum,pdmag,pfluxratio,pseeing,format='a,i,a,f,a,f,f,f,f',/silent
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

        outfile       = imworkingdir + 'coadd' + strtrim(target,2) +'_'+ strtrim(filter,2) + '.fits'
        outweightfile = imworkingdir + 'coadd' + strtrim(target,2) +'_'+ strtrim(filter,2) + '.weight.fits'
        stackcmd = swarpcommand+' '    ;'swarp '
        for s = 0, n_elements(stacklist)-1 do begin
           if s eq 0 then stackcmd = stackcmd + stacklist[s] 
           if s gt 0 then stackcmd = stackcmd + ',' + stacklist[s] 
        endfor
        for s = 0, n_elements(stacklist)-1 do begin
           weightfilename = strmid(stacklist[s],0,strlen(stacklist[s])-5-0)  + '.weight.fits'
           weightexists = file_test(weightfilename)
           if weightexists eq 0 then begin
              slashpos = strpos(weightfilename,'/',/reverse_search)
              weightfilenameinit = imworkingdir + strmid(weightfilename,slashpos+2)  
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

        if (file_test(outfile) eq 0) or overwrite then begin
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
              if camera eq 'red' and lrisversion eq 1 then medseeingarcsec = medseeingpix*0.210 else medseeingarcsec = medseeingpix*0.135
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

     coaddfiles = findfile(imworkingdir+'coadd*.fits')

     if coaddfiles ne [''] then coaddfiles = coaddfiles[where(coaddfiles ne '')] ;else continue

     ; match every r with an l and coadd
     for f = 0, n_elements(coaddfilesr)-1 do begin
        filename = coaddfiles[f]
        if ct eq 1 then begin
           outfile = coaddfiles[f]
           if file_test(outfile) eq 0 or overwrite then begin
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
     
     azffiles = findfile(imworkingdir+'a*f*'+prefchar+wildchar+'r.fits')
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
           if file_test(outfile) eq 1 and overwrite eq 0 then continue  ; already stacked these two

           print
           print, 'Preparing to combine ' + infiler + ' and ' + infilel
           print, autoastrocommand+' '+infiler+' -q'
           spawn, autoastrocommand+' '+infiler+' -q'
           print, autoastrocommand+' '+infilel+' -q'
           spawn, autoastrocommand+' '+infilel+' -q'

           stackcmd = swarpcommand+' ' + extractpath(infiler)+'a' + removepath(infiler) + ',' + extractpath(infilel) + 'a' + removepath(infilel)
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


end



; -------------------------

; need to allow the user to completely ignore the left chip at all stages by specifying the chip.
; need to restore the gain correction.
; need to do something about when crashes, leaves you in imredux (check if you are already in the imredux directory)

pro ratautoproc, datadirectory=datadirectory, modestr=modestr, camerastr=camerastr, chipstr=chipstr, files=files, filters=filters, gratings=gratings, targets=targets, start=start, stop=stop, only=only, step=step, nocrclean=nocrclean, redo=redo, nofringe=nofringe, continuous=continuous
;   modestr      - Mode (imaging or spectroscopy)
;   camerastr    - Camera to process (red or blue)
;   chipstr      - Chip to process (left or right)
;   files        - File numbers to process
;   filters      - Filters to process
;   gratings     - Gratings to process
;   targets      - Targets to process
;   start        - Start with this step, skipping previous ones
;   stop         - End with this step, skipping subsequent ones
;   only         - Do only this step
;   step         -  (completely identical to only, takes precedence)
;   redo         - Overwrite any existing products
;   nocrclean    - Do not zap cosmic rays
;   continuous   - Keep running, assimilating new data as it appears, not implemented yet
; XX quick        - Set defaults to: imaging only, right chip only, no cr-cleaning
; XX  nofringe     - Skip fringe production and correction for LRIS-R1

!quiet = 1

; Load default parameters and interpret user arguments.

common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite, fileseq, filterstr, gratingstr, targetstr
common lrisfail, flatfail, catastrofail, relastrofail, fullastrofail, extractfail, wavsolfail, wavsolwarn, fluxcalfail


close, /all

if n_elements(datadirectory) gt 0 then datadir = datadirectory

if keyword_set(quick) then begin
  if n_elements(chipstr) eq 0 then chipstr = 'r'
  if n_elements(modestr) eq 0 then modestr = 'i'
  if n_elements(crclean) eq 0 then nocrclean = 1
endif

if n_elements(files) then fileseq = strseq(files)
if n_elements(filterstr) then filterseq = strsplit(filters,',',/extract)
if n_elements(targetstr) then filterseq = strsplit(targets,',',/extract)
if n_elements(gratingstr) then gratingseq = strsplit(gratings,',',/extract)

redo = keyword_set(redo)
nocrclean = keyword_set(nocrclean)
nofringe = keyword_set(nofringe)

cd, current=pwd
dirtree = strsplit(pwd,'/',/extract,count=nd)
lastdir = dirtree[nd-1]
if lastdir ne '' then lastdir += '/'  ; whether or not a slash is on the end is very
                                      ; confusing, need to rethink this.
if lastdir eq 'imredux/' or lastdir eq 'spredux/' then begin
     print, 'Currently in a reduction subdirectory.'
     print, 'Type cd.. and rerun.'
  ; could reinterpret this as run with the mode set to whatever this directory is...
     return
endif

autolrisdefaults



; --- Process mode options
if n_elements(modestr) eq 0 then modestr = 'is'
modestr = strlowcase(modestr)
if modestr eq 'is' or modestr eq 'i,s' or modestr eq 'im,sp'  then modes = ['im', 'sp']
if modestr eq 'si' or modestr eq 's,i' or modestr eq 'sp,im'  then modes = ['sp', 'im']
if n_elements(modes) eq 0 then begin
   if strmid(modestr,0,1) eq 'i' then modes = ['im']
   if strmid(modestr,0,1) eq 's' then modes = ['sp']
endif
if n_elements(modes) eq 0 then begin
  print, 'Cannot recognize mode request: ', modestr
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
steps = ['prepare', 'makeflat', 'flatten', 'makesky', 'skysub', 'makefringe', 'rmfringe', 'split', 'crclean', 'skysubtract', 'extract', 'fluxcal', 'combine', 'connect', 'astrometry', 'photometry', 'stack']
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
      cmd = swarpcommand+' -d > temp.txt'
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
      cmd = autoastrocommand+' > temp.txt'
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

if redo then overwrite = 1 else overwrite = 0

if keyword_set(nocrclean) eq 0 then flag = 1

nsteps = n_elements(steps)
nmodes = n_elements(modes)
ncameras = n_elements(cameras)
nchips = n_elements(chips)

npipeiter = 1L
if keyword_set(continuous) eq 0 then maxpipeiter = 1L else maxpipeiter = 100000L

while npipeiter le maxpipeiter do begin
for istep = 0, nsteps-1 do begin
   instep = steps[istep]

   for icam = 0, ncameras-1 do begin
      camera = cameras[icam]
      ca = strmid(cameras[icam],0,2)

      if instep eq 'prepare' then begin
         ; the mode setting is passed on to the routine itself.
         autolrisprepare, modestr=modestr, camstr=ca
      endif

      for imode = 0, nmodes-1 do begin
         mo = strmid(modes[imode],0,2)
         if mo eq 'im' then begin   
            if instep eq 'flatten'  then autolrisimflatten,  cam=camera, chip=''
            if instep eq 'makesky' then autolrismakesky,cam=camera, chip=''
            if instep eq 'skysub'  then autolrisskysub,  cam=camera, chip=''
            if lrisversion eq 1 and nofringe eq 0 and cameras[icam] eq 'red' then begin
               if instep eq 'makefringe' then autolrismakefringe, chip=ch, cam=ca
               if instep eq 'rmfringe'   then autolrisrmfringe, chip=ch, cam=ca
            endif
         endif

         for ichip = 0, nchips-1 do begin
            ch = strmid(chips[ichip],0,1)

            if nocrclean eq 0 and instep eq 'crclean' then begin
               if mo eq 'im' then autolriscrcleanim,    cam=camera, chip=ch
               if mo eq 'sp' then autolriscrcleanspec,  cam=camera, chip=ch
            endif 

            if mo eq 'im' then begin   
               if instep eq 'astrometry' then autolrisastrometry, chip=ch, cam=camera
               if camera eq 'blue' then bl = 1 else bl = 0
               if camera eq 'red' then  re = 1 else re = 0
               if n_elements(chips) eq 1 then begin
                  if instep eq 'photometry' then autolrisphotometry, chip=ch, camera=camera
               endif else begin
                  if instep eq 'photometry' and ichip eq 1 then autolrisphotometry, camera=camera
               endelse
               if instep eq 'stack' then autolrisstack, chip=ch, cam=camera
            endif
            if mo eq 'sp' then begin   
               if instep eq 'skysubtract' then autolrisskysubtract, chip=ch, cam=camera
               if instep eq 'extract' and ch eq 'r' then autolrisextract, chip='r', cam=camera
            endif
         endfor
      endfor ; mode
   endfor ; camera
endfor ; step

print

nflatfail = n_elements(flatfail)-1
if nflatfail gt 1 then begin ;element 0 is blank
  print
  print, 'Unable to flat-field the following images:'
  for f = 1, n_elements(flatfail)-1 do begin
    h = headfits(flatfail[f])
    gg = sxpar(h,'GRNAME')
    if clip(gg) eq '0' then gg = sxpar(h,'SLITNAME')
    print, clip(removepath(flatfail[f]),25), string(sxpar(h,'ELAPTIME'),format='(I5)'), ' ', clip(sxpar(h,'FILTER'),5), ' ', clip(sxpar(h,'DICHNAME'),7),' ',clip(gg,12), ' ',clip(sxpar(h,'TARGNAME'),16), ' ', repstr(clip(sxpar(h,'BINNING')),',','x')
  endfor
  flats = [findfile(imworkingdir+'*flat*.fits'),findfile(spworkingdir+'*flat*.fits')]
  nflats = total(flats ne '')
  if nflats gt 0 then begin
    flats = flats[where(flats ne '')]
    print, 'Available flats are:'
    for f = 0, nflats-1 do begin
       h = headfits(flats[f], /silent)
       gg = sxpar(h,'GRNAME')
       if clip(gg) eq '0' then gg = sxpar(h,'SLITNAME')

print, clip(removepath(flats[f]),25), string(sxpar(h,'ELAPTIME'),format='(I5)'), ' ', clip(sxpar(h,'FILTER'),5), ' ', clip(sxpar(h,'DICHNAME'),7),' ',clip(gg,12), ' ',clip(sxpar(h,'TARGNAME'),16), ' ', repstr(clip(sxpar(h,'BINNING')),',','x')

     endfor
   endif else begin
     print, 'No processed flat-fields exist!  Check inputs or copy raw flat-field files to data directory.'
   endelse
endif

nefail = n_elements(extractfail)-1
if nefail ge 1 then begin
  print, 'Unable to model trace and extract 1D spectrum for the following 2D spectra:'
  for f = 1, n_elements(extractfail)-1 do begin
    h = headfits(extractfail[f])
    print, clip(removepath(extractfail[f]),25), rclip(sxpar(h,'ELAPTIME'),5)+'s', ' ', clip(sxpar(h,'GRNAME'),9), ' ', clip(sxpar(h,'TARGNAME'),16)
  endfor
endif
nwfail = n_elements(wavsolfail)-1
if nwfail ge 1 then begin
   print, 'Unable to wavelength-calibrate the following spectra:'
   for f = 1, n_elements(wavsolfail)-1 do begin
      h = headfits(wavsolfail[f])
      print, clip(removepath(wavsolfail[f]),25), rclip(sxpar(h,'ELAPTIME'),5)+'s', ' ', clip(sxpar(h,'GRNAME'),9), ' ', clip(sxpar(h,'TARGNAME'),16)
   endfor
endif
nwwarn = n_elements(wavsolwarn)-1
if nwwarn ge 1 then begin
   print, 'Wavelength solutions for the following spectra are questionable; check sky-line verification plots:'
   for f = 1, n_elements(wavsolwarn)-1 do begin
      h = headfits(wavsolwarn[f])
      print, clip(removepath(wavsolwarn[f]),25), rclip(sxpar(h,'ELAPTIME'),5)+'s', ' ', clip(sxpar(h,'GRNAME'),9), ' ', clip(sxpar(h,'TARGNAME'),16)
   endfor
endif
nfluxfail = n_elements(fluxcalfail)-1
if nfluxfail ge 1 then begin
   print, 'Unable to flux-calibrate the following spectra:'
   for f = 1, n_elements(fluxcalfail)-1 do begin
      h = headspec(fluxcalfail[f])
      print, clip(removepath(fluxcalfail[f]),25), rclip(sxpar(h,'ELAPTIME'),5)+'s', ' ', clip(sxpar(h,'GRNAME'),9), ' ', 'D'+clip(sxpar(h,'DICHNAME'),7), ' ', clip(sxpar(h,'FILTER'),5), ' ', clip(sxpar(h,'TARGNAME'),16)
   endfor
endif

nafail = n_elements(relastrofail)-1 + n_elements(fullastrofail)-1 + n_elements(catastrofail)-1
if nafail ge 1 then begin
  print
  if n_elements(catastrofail) gt 1 then print, 'Unable to produce astrometric catalogs for the following reference images:'
  for f = 1, n_elements(catastrofail)-1 do begin
    h = headfits(catastrofail[f])
    print, clip(removepath(catastrofail[f]),25), string(sxpar(h,'ELAPTIME'),format='(I5)'), ' ', clip(sxpar(h,'FILTER'),4), ' ', clip(sxpar(h,'TARGNAME'),16), ' ', clip(sxpar(h,'COUNTS'))
  endfor
  if n_elements(relastrofail) gt 1 then print, 'Relative astrometry failed for the following images, but absolute was successful:'
  for f = 1, n_elements(relastrofail)-1 do begin
    h = headfits(relastrofail[f])
    print, clip(removepath(relastrofail[f]),25), string(sxpar(h,'ELAPTIME'),format='(I5)'), ' ', clip(sxpar(h,'FILTER'),4), ' ', clip(sxpar(h,'TARGNAME'),16), ' ', clip(sxpar(h,'COUNTS'))
  endfor
  if n_elements(fullastrofail) gt 1 then print, 'All astrometry failed for the following images (not stacked):'
  for f = 1, n_elements(fullastrofail)-1 do begin
    h = headfits(fullastrofail[f])
    print, clip(removepath(fullastrofail[f]),25), string(sxpar(h,'ELAPTIME'),format='(I5)'), ' ', clip(sxpar(h,'FILTER'),4), ' ', clip(sxpar(h,'TARGNAME'),16), ' ', clip(sxpar(h,'COUNTS'))
  endfor
endif

if keyword_set(continuous) then begin
   print, 'Pipeline iteration complete.'
   print, 'Continuous mode is active.'
   print, 'Pausing for 30 seconds, then will loop from the beginning.  Press Ctrl+C to stop.'
   for s = 30, 0, -1 do begin 
      wait, 1.0
   endfor
endif

npipeiter += 1
endwhile


print, 'Processing complete.'
!quiet = 0

; cleanup- should ultimately put these things in imredux

if file_test('temp*.*') gt 0 then spawn, 'rm -f temp*.*'
if file_test('det.*') gt 0 then spawn, 'rm -f det.*'
if file_test('cat.*') gt 0 then spawn, 'rm -f cat.*'

end

 
