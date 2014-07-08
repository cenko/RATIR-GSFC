;+
; NAME:
;	autopipemakesky
;
; PURPOSE:
;	Combine sky flats based on filter type (uses sextractor source identification
;	for 'r' and 'i' - 2x 100% flux fraction radius - and uses sigma clipping for sources for rest of filters)
;
; OPTIONAL KEYWORDS:
;	outpipevar - output pipeline parameters
;	inpipevar  - input pipeline parameters (typically set in autopipedefaults.pro or ratautoproc.pro, but can be set to default) 
;
; EXAMPLE:
;	autopipemakesky, outpipevar=pipevar, inpipevar=pipevar
;
; DEPENDENCIES:
;	skypipecombine_altsex, sigclipstats_vt, djs_iterstat, Sextractor
;
; Written by Dan Perley 
; Modified by Vicki Toy 11/18/2013
;
; FUTURE IMPROVEMENTS:
;	change max iteration for sigma clipping? change default values of skypipecombine?
;	Change source identification to be user modified?
;-

pro autopipemakesky, outpipevar=outpipevar, inpipevar=inpipevar


	;Setup pipeline variables that carry throughout the pipeline
	if keyword_set(inpipevar) then begin
		pipevar = inpipevar
		print, 'Using provided pipevar'
	endif else begin
		pipevar = {autoastrocommand:'autoastrometry', getsedcommand:'get_SEDs', $
					sexcommand:'sex' , swarpcommand:'swarp' , $
					prefix: '', datadir:'' , imworkingdir:'' , overwrite:0 , $
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:'' }
	endelse
	
	
	;Copies necessary parameter file for sextractor if not in current working directory
	if file_test('source.param') eq 0 then spawn, 'cp '+ pipevar.defaultspath +'/source.param .'
	if file_test('sex_source.config') eq 0 then spawn, 'cp '+ pipevar.defaultspath +'/sex_source.config .'
	if file_test('sex.conv') eq 0 then spawn, 'cp '+ pipevar.defaultspath +'/sex.conv .'
	if file_test('default.nnw') eq 0 then spawn, 'cp '+ pipevar.defaultspath +'/default.nnw .'
	
	;Finds files with given prefix
    files = findfile(pipevar.imworkingdir+'fp'+pipevar.prefix+'*sky*.fits')
    
    filters =  strarr(n_elements(files))   
    
    ;Creates filter list from opening each file    
    for f = 0, n_elements(files)-1 do begin
       if files[f] eq '' then continue
       h = headfits(files[f], /silent)
       filters[f] = clip(sxpar(h, 'FILTER'))
    endfor  
    
    ;Unique list of filters
    filterlist  = unique(filters)

	;For each unique filter, combine sky files using skycombine if more than 2 files
	;Else return list of unprocessed files
    for i = 0, n_elements(filterlist)-1 do begin
    	filt = filterlist[i]
       
       	skyflats = where(filters eq filt, ctsky)
       	outflatname = pipevar.imworkingdir+'sky-'+filt+'.fits'
       
       	if ctsky ge 2 then begin
       		print, filt, '-band sky flats.'
            if file_test(outflatname) and pipevar.overwrite eq 0 then continue
            
            print, files[skyflats]
            
            ;COMMENTED OUT 12/30
            ;skypipecombine, files[skyflats], outflatname, /removeobjects, type='sky'
            
            skypipecombine_altsex, files[skyflats], outflatname, filt, pipevar,/removeobjects, type='sky'

        endif else begin
        
          	print, 'Unable to produce a flat field for this setting: ' + filt
          	noproc = where(filters eq filt, ctnoproc)
          	print, 'Will not be able to further process ', clip(ctnoproc), ' images without a flat from another source:'
          	
          	for j = 0, ctnoproc-1 do print, '   '+ files[noproc[j]]
          
       	endelse
    endfor
    
	outpipevar = pipevar

end