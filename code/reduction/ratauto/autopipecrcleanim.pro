;+
; NAME:
;	autopipecrcleanim
;
; PURPOSE:
;	Identify and remove cosmic ray candidates using pzap_perley.pro (Remove cosmic rays 
;	from a 2-D image using model subtraction / percolation.) Saves file as 'z'+file and 
;	also creates weight file
;
; OPTIONAL KEYWORDS:
;	outpipevar - output pipeline parameters
;	inpipevar  - input pipeline parameters (typically set in autopipedefaults.pro or ratautoproc.pro, but can be set to default) 
;
; EXAMPLE:
;	autopipecrcleanim, outpipevar=pipevar, inpipevar=pipevar
;
; DEPENDENCIES:
;	pzap_perley.pro <--- Has not been examined carefully
;
; Written by Dan Perley 
; Modified by Vicki Toy 11/21/2013
;
; FUTURE IMPROVEMENTS:
;	prefchar in variable structure? New zapping function??? 
;-

pro autopipecrcleanim, outpipevar=outpipevar, inpipevar=inpipevar

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

	;CHANGE FOR RIMAS VLT
   	prefchar = '2'

	;Choose files that have been flatfielded and sky subtracted, if overwrite not set
	;then the files are altered to remove any that have already been zapped
   	files = choosefiles(prefchar+'*_img_?.fits',pipevar.imworkingdir+'sfp')
   	if pipevar.overwrite eq 0 then files = unmatched(files,'z')
   
   	;For each file read in the header and check that objects meet count limits and exposure time
   	;(i.e. short exposure time with lot of counts will be ignored), also targets that are on 
   	;calibration files will be ignored.
   	;Then run pzap_perley on the files and have output files be 'z'+file plus weight files
   	for f = 0, n_elements(files)-1 do begin
      	if files[f] eq '' then continue
      	h = headfits(files[f])
      	counts  = sxpar(h,'COUNTS')
      	exptime = sxpar(h,'ELAPTIME')
      	target  = sxpar(h,'TARGNAME')
      	if counts gt 40000. then continue
      	if counts gt 30000. and exptime lt 10. then continue
      	if counts gt 20000. and exptime gt 5. and exptime lt 10. then continue
      	if target eq '' then continue   ; blank target is certainly not a science object.
      	if strpos(target,'flat') ge 0 then continue
      	if strpos(target,'sky') ge 0 then continue
      	if strpos(target,'twilight') ge 0 then continue

		;Cleans cosmic rays from file using hardcoded zeal number in pzap_perley
      	print, 'Cleaning cosmic rays from ', removepath(files[f])
		zeal = 0.75
      	if n_elements(setzeal) eq 1 then zeal = setzeal
      	slashpos = strpos(files[f],'/')
      	dir = strmid(files[f],0,slashpos+1)
      	outname = dir + 'z' + strmid(files[f],slashpos+1)
      	if file_test(outname) and pipevar.overwrite eq 0 then continue
      	pzap_perley, files[f], /weight, zeal=zeal, usamp=usamp, /quiet
   	endfor

	outpipevar = pipevar

end
