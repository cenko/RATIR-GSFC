;+
; NAME:
;	autopipemakesky
;
; PURPOSE:
;	Combine sky flats based on filter type
;
; OPTIONAL KEYWORDS:
;	outpipevar - output pipeline parameters
;	inpipevar  - input pipeline parameters (typically set in autopipedefaults.pro or ratautoproc.pro, but can be set to default) 
;
; EXAMPLE:
;	autopipemakesky, outpipevar=pipevar, inpipevar=pipevar
;
; DEPENDENCIES:
;	skypipecombine, sigclipstats_vt, djs_iterstat
;
; Written by Dan Perley 
; Modified by Vicki Toy 11/18/2013
;
; FUTURE IMPROVEMENTS:
;	prefchar in variable structure?, change max iteration for sigma clipping? change default values of skypipecombine?
;-

pro autopipemakesky, outpipevar=outpipevar, inpipevar=inpipevar

	;Setup pipeline variables that carry throughout the pipeline
	if keyword_set(inpipevar) then begin
		pipevar = inpipevar
		print, 'Using provided pipevar'
	endif else begin
		pipevar = {autoastrocommand:'autoastrometry' , sexcommand:'sex' , swarpcommand:'swarp' , $
					datadir:'' , imworkingdir:'' , overwrite:0 , $
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:'' }
	endelse
	
	;Finds files with given prefix CHANGE FOR RIMAS VLT
	prefchar = '2'

    files = findfile(pipevar.imworkingdir+'fp'+prefchar+'*sky*.fits')
    
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
            
            skypipecombine, files[skyflats], outflatname, /removeobjects, type='sky'

        endif else begin
        
          	print, 'Unable to produce a flat field for this setting: ' + filt
          	noproc = where(filters eq filt, ctnoproc)
          	print, 'Will not be able to further process ', clip(ctnoproc), ' images without a flat from another source:'
          	
          	for j = 0, ctnoproc-1 do print, '   '+ files[noproc[j]]
          
       	endelse
    endfor
    
	outpipevar = pipevar

end