;+
; NAME:
;	autopipeastrometry
;
; PURPOSE:
;	Calculate astrometry of image files to fix WCS coordinates (shift and rotation) in header.
;	Using fast astrometry solver (vlt_autoastrometry.py) that using pair-distance matching 
;	and asterism matching.  Returns file with corrected WCS coordinates saved as 'a'+fitsfile
;
; OUTPUT:
;	Creates output file with name 'a'+fitsfile for all files that astrometry could correct.
;
; OPTIONAL KEYWORDS:
;	outpipevar - output pipeline parameters
;	inpipevar  - input pipeline parameters (typically set in autopipedefaults.pro or ratautoproc.pro, but can be set to default) 
;
; EXAMPLE:
;	autopipeastrometry, outpipevar=pipevar, inpipevar=pipevar
;
; DEPENDENCIES:
;	vlt_autoastrometry.py (run as an executable, with various calls to other python scripts)
;
; Written by Dan Perley 
; Modified by Vicki Toy 11/19/2013
;
; FUTURE IMPROVEMENTS:
;	prefchar in variable structure? save region and matchline files with different names 
;	to check image after all have run?
;-
pro autopipeastrometry, outpipevar=outpipevar, inpipevar=inpipevar

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
	
	;CHANGE FOR RIMAS
    prefchar = '2'
	;Find files that have been processed giving zapped cosmic rays preference over unzapped
    zffiles = choosefiles(prefchar+'*_img_?.fits',pipevar.imworkingdir+'zsfp',pipevar.imworkingdir+'sfp')
    
    filetargets   = strarr(n_elements(zffiles))
    fileexposures = strarr(n_elements(zffiles))
    filecounts    = fltarr(n_elements(zffiles))
    filefilt      = strarr(n_elements(zffiles))

	;For each file grab header values (target name, filter, exposure time, and counts) and store in arrays
    for f = 0, n_elements(zffiles)-1 do begin
       	if zffiles[f] eq '' then continue
       	h = headfits(zffiles[f], /silent)
       	filetargets[f]   = repstr(repstr(strtrim(sxpar(h,'TARGNAME'),2),' ', '_'),'/','_')  ; spaces cause barfing in filenames
       	filefilt[f]      = sxpar(h,'FILTER')
       	fileexposures[f] = string(sxpar(h,'EXPTIME'))
       	filecounts[f]    = sxpar(h,'COUNTS')
    endfor

    ;Find unique targets and filters
    targets = unique(filetargets)
    filters = unique(filefilt)

    ; Make a reference catalog using a representative image out of an image block (several images of the same field)
	;For each unique target and unique filter
    for t = 0, n_elements(targets)-1 do begin
       	
       	for f = 0, n_elements(filters)-1 do begin
			
			;Skip if cat file exists and overwrite is set to false
          	refcatfile = strcompress(pipevar.imworkingdir+targets[t]+'.'+filters[f]+'.cat',/remove_all)
          	if file_test(refcatfile) and pipevar.overwrite eq 0 then continue
          	
          	;Find which files have the correct target and filter and find maximum
          	;exposure and minimum count rate from saved arrays
          	thistarget = where(filetargets eq targets[t] and filefilt eq filters[f])
          	maxexp 	   = max(fileexposures(thistarget))
          	minctrate  = min(filecounts[thistarget]/fileexposures[thistarget])
          
          	;Skip files that have too shallow of a field
          	if maxexp lt 5 or minctrate gt 5000 then begin
             	print, targets[t], ' is a shallow field (standard or sky): not making a reference catalog.'
             	continue
          	endif
          
          	;Use middle file from list
          	imagesthistarg = zffiles[thistarget]
          	refimagename = imagesthistarg[n_elements(imagesthistarg)/2]
          	
          	;Run astrometry correction on this middle file that is assumed to be representative of the filter
          	print, 'Making reference catalog for ', targets[t], ' using ', refimagename
          	print, pipevar.autoastrocommand+' '+refimagename
          	spawn, pipevar.autoastrocommand+' '+refimagename

			;The new astrometry corrected file should be saved with the same name, but with 'a' prefix
          	outfile = fileappend(refimagename,'a')
          
          	;If astrometry corrected file was not created add to list of failed files
          	;If it was created, run truncated astrometry on corrected file (will just run sextractor
          	;and pull good sources out) using saturation level as input for sextractor
          	if file_test(outfile) eq 0 then begin
             	pipevar.catastrofail = pipevar.catastrofail +' '+ refimagename
             	print, 'WARNING - astrometry on the reference image was unsuccessful!'
          	endif else begin
             	print, pipevar.autoastrocommand+' '+outfile+' -n '+refcatfile + ' -l 55000 -q'
             	spawn, pipevar.autoastrocommand+' '+outfile+' -n '+refcatfile + ' -l 55000 -q'
          	endelse
          	
       	endfor
    endfor

	;If overwrite not set, make sure you don't rewrite astrometry corrected middle file that was processed (wouldn't do much harm)
    if pipevar.overwrite eq 0 then zffiles = unmatched(zffiles,'a')

    ; Use the reference catalog to do a more precise relative astrometric solution
    for f = 0, n_elements(zffiles)-1 do begin
    
       	if zffiles[f] eq '' then continue
       	outfile = fileappend(zffiles[f],'a')
       	if file_test(outfile) and pipevar.overwrite eq 0 then continue
       	h = headfits(zffiles[f], /silent)
       	exptime = sxpar(h,'ELAPTIME')
       	counts = sxpar(h,'COUNTS')
       	filt = sxpar(h,'FILTER')
       	
       	;If count rate is high or exposure time is low then run astrometry as long as there aren't too many counts
       	;combined with a short exposure time (very bright frame).  Check that new astrometry file created (prefix 'a')
       	;and if was not add to failure list.
       	;If count rate too high or exposure time is too low, run astrometry using catalog of sources
       	;created from middle file of filter.  If catalog doesn't work, try direct astrometry
       	if counts/exptime gt 5000. or exptime lt 5 then begin 
         	; only do catalog astrometry for short exposures (no catalog).
          
          	if counts gt 30000 and exptime lt 10 then begin
             	print, zffiles[f], ' is in bright twilight; not solving astrometry.'
             	continue
          	endif
          
          	print, zffiles[f], ' is a twilight/standard frame, solving astrometry directly against a catalog.'
          	print, pipevar.autoastrocommand+' '+zffiles[f]
          	spawn, pipevar.autoastrocommand+' '+zffiles[f]

          	if file_test(outfile) eq 0 then pipevar.fullastrofail = pipevar.fullastrofail +' '+ zffiles[f]
          	
       	endif else begin
         	; use the short exposure first, but fall back on direct catalog.
          	targname = repstr(strtrim(sxpar(h,'TARGNAME'),2),' ', '_') ; spaces cause barfing in filenames
          	if targname eq '' then continue
          	if strpos(targname,'flat') ge 0 then continue
          	refcatfile = strcompress(pipevar.imworkingdir+targname+'.'+filt+'.cat',/remove_all)

          	if file_test(refcatfile) then begin
             	print, targname
             	print, pipevar.autoastrocommand+' '+zffiles[f]+' -c '+refcatfile
             	spawn, pipevar.autoastrocommand+' '+zffiles[f]+' -c '+refcatfile
          	endif else begin
             	print, 'No reference catalog '+refcatfile+' exists for this field.'
          	endelse
          		if file_test(outfile) eq 0 then begin
             		print, 'Refined astrometry of ', zffiles[f], ' was not successful.  Trying direct astrometry:'
             		print, pipevar.autoastrocommand+' '+zffiles[f]
             		spawn, pipevar.autoastrocommand+' '+zffiles[f]
             		if file_test(outfile) then pipevar.relastrofail  = pipevar.relastrofail + ' ' +  zffiles[f] $
             			else pipevar.fullastrofail = pipevar.fullastrofail  + ' ' + zffiles[f] 
          		endif
       	endelse
    endfor
    
    outpipevar = pipevar
end

