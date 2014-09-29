; NAME:
;   ratautoproc.pro                               (release version 2014.09.26)
;
; PURPOSE:
;   Fully-automated reduction of imaging.
;
; CALLING SEQUENCE:
;   ratautoproc, [settings]
;
; INPUTS (all optional):
;   datadir 	- Location of raw data (current directory if unspecified)
;	imdir		- Location of processed data ('imredux' if unspecified)
;   start		- Start with this step, skipping previous ones
;   stop    	- End with this step, skipping subsequent ones
;   only    	- Do only this step
;   step    	-  (completely identical to only, takes precedence)
;   redo    	- Repeat step(s), overwriting any existing files
;   nocrclean	- Do not zap cosmic rays
;	quiet		- (mainly) silent output unless errors
;	rmifiles	- Removes intermediate files
;
; ADDITIONAL OPTIONS:
;   If any of the following files are found in the directory where ratautoproc
;   is run, it will change the default behavior.
;
;   pipeautoproc.par - Contains various defaults (mainly directory paths)
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
;   code, which is still in production. 
;
;   Filenames for the input raw images are expected to be in 
;   2*.fits format. They can be either in the working directory or 
;	in a different directory specified by datadir.
; 
;   The code runs in a series of steps in the following order:
;
;   1. "Prepare" the data by converting the multi-header extension FITS files
;      to a standard single frame, and adding extra information to the header. 
;      Output: p*.fits (written to ./imredux/ by default.)
;
;   2. (TAKEN OUT, assume that have master flat in format flat_(filter).fits in imredux/ folder)
;	   Create flat-fields.  The program searches through all science exposures
;      and determines what flat fields are needed, then attempts to make the
;      best possible flat for each case (generally, the order of preference is
;      supersky > twilight > dome > internal halogen, but which type of flat
;      is actually made depends on the number of frames available and their
;      quality.)  In the case of imaging flats, stars are identified and
;      masked out.
;      Output:  fp*.fits
;
;   3. Flat-correct data.  Divide each image by the flatfield.  A more
;      refined cropping is also done at this stage, depending on the placement
;      of the image area during the run (which is variable.) Assumes master flat
;	   (ex. flat_H.fits in imredux/ folder)
;
;   4. Removes cosmic rays, using the independent routines pzap.pro and
;      pzapspec.pro.  See those programs for more information.  
;      This can be a time-consuming process.
;      Output: zfp*.fits
;  
;   5. Solve astrometry of the field against the best available online catalog
;      (SDSS/2MASS/APASS/USNO-B1.0), using the independent vlt_autoastrometry.py code.
;      (Also requires sextractor to be installed.)
;      Uses two passes of Scamp for a secondary correction.  Scamp accounts for 
;	   distortion.  Uses higher distortion parameters if already supplied distortion keywords.
;      Output: azfp*.fits
;  
;   6. Stack exposures.  A weighted, masked median is performed for each field.
;      Requires swarp to be installed and properly linked.
;      Output: coadd[object].[filter].fits
;
; 
; EXAMPLES:
;   1.  In a directory containing a full night worth of RATIR data, enter
;
;       IDL> ratautoproc
;
;       This will automatically execute all of the steps above.
;
;
;   2.  Run only the "prepare" step, on data stored in a separate directory:
;
;       IDL> ratautoproc, step='prepare', datadir='raw/'
;
;
;   TROUBLESHOOTING:
;
;   If the pipeline crashes during operation:
;  
;   * Did an external call (autoastrometry, swarp, sex) fail?
;       Check that you have specified the paths correctly and actually have the
;       required software installed for astrometry/coadding (see "installation", 
;       above.)
;   
;   * Did it encounter some other error while processing a file?  
;       Check which file the program was working on when it crashed.  If the
;       file is not essential, try deleting it and re-running the pipeline
;       starting with the current step (or delete ALL of the file's precursors
;       with the same file number and rerun the pipeline.) 
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
;
;   * Were some files skipped?
;       If the pipeline encountered a non-fatal problem processing an
;       individual image (such as an inability to flatfield) then it will 
;		simply not process that file any further.  For most such cases a 
;		summary of the problems will be printed out at the end of processing.
;       If a file is not being processed and you do not see it in the final summary,
;       you can simply rerun the pipeline (without deleting any files and without
;       setting the redo flag) and it will try to repeat any failed steps of this nature. 
;- 

pro ratautoproc, datadir=datadir, imdir=imdir, start=start, stop=stop, only=only, step=step,$
	 nocrclean=nocrclean, redo=redo, quiet=quiet,rmifiles=rmifiles

	!quiet = 1
	close, /all
	
	; Load default parameters and interpret user arguments.
	pipevar = {autoastrocommand:'autoastrometry', getsedcommand:'get_SEDs', $
					sexcommand:'sex' , swarpcommand:'swarp' , $
					prefix:'', datadir:'' , imworkingdir:'' , overwrite:0 , verbose:0, rmifiles:0,$
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:'' }
	
	if keyword_set(redo) then pipevar.overwrite=1
	if keyword_set(quiet) then pipevar.verbose = 0 else pipevar.verbose = 1
	if n_elements(datadir) gt 0 then pipevar.datadir = datadir	

	autopipedefaults, outpipevar=pipevar, inpipevar=pipevar	
	if keyword_set(imdir) then pipevar.imworkingdir = imdir
	if keyword_set(rmifiles) then pipevar.rmifiles = 1	

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

      	if instep eq 'prepare' then autopipeprepare,   outpipevar=pipevar, inpipevar=pipevar      		 
        if instep eq 'flatten' then autopipeimflatten, outpipevar=pipevar, inpipevar=pipevar
        if instep eq 'makesky' then autopipemakesky,   outpipevar=pipevar, inpipevar=pipevar
        if instep eq 'skysub'  then autopipeskysub,    outpipevar=pipevar, inpipevar=pipevar
    	
    	if nocrclean eq 0 and instep eq 'crclean' then begin
    		autopipecrcleanim, outpipevar=pipevar, inpipevar=pipevar
        endif 
            		 
        if instep eq 'astrometry' then autopipeastrometry, outpipevar=pipevar, inpipevar=pipevar
        if instep eq 'stack'      then autopipestack,      outpipevar=pipevar, inpipevar=pipevar
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