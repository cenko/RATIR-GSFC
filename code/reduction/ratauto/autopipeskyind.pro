pro autopipeskyind, outpipevar=outpipevar, inpipevar=inpipevar

	print, 'MAKE/SUBTRACT SKY INDIVIDUAL'
	
	
	;Setup pipeline variables that carry throughout the pipeline
	if keyword_set(inpipevar) then begin
		pipevar = inpipevar
		if pipevar.verbose gt 0 then print, 'Using provided pipevar'
	endif else begin
		pipevar = {autoastrocommand:'autoastrometry', getsedcommand:'get_SEDs', $
					sexcommand:'sex' , swarpcommand:'swarp' , $
					prefix:'', datadir:'' , imworkingdir:'' , overwrite:0 , verbose:0, $
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:'' }
	
	pipevar = {autoastrocommand:'autoastrometry', getsedcommand:'get_SEDs', $
					sexcommand:'sex' , swarpcommand:'swarp' , $
					prefix:'', datadir:'' , imworkingdir:'' , overwrite:0 , verbose:0, $
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:'' }
					
	autopipedefaults, outpipevar=pipevar, inpipevar=pipevar
	

	endelse
	
	;Copies necessary parameter file for sextractor if not in current working directory
	if file_test('source.param') eq 0 then spawn, 'cp '+ pipevar.defaultspath +'/source.param .'
	if file_test('sex_source.config') eq 0 then spawn, 'cp '+ pipevar.defaultspath +'/sex_source.config .'
	if file_test('sex.conv') eq 0 then spawn, 'cp '+ pipevar.defaultspath +'/sex.conv .'
	if file_test('default.nnw') eq 0 then spawn, 'cp '+ pipevar.defaultspath +'/default.nnw .'
	
	;Finds files with given prefix
    files = findfile(pipevar.imworkingdir+'fp'+pipevar.prefix+'*img_J.fits')
    
    filters =  strarr(n_elements(files))   
    dithers =  strarr(n_elements(files))   
    
    ;Creates filter list from opening each file    
    for f = 0, n_elements(files)-1 do begin
       if files[f] eq '' then continue
       h = headfits(files[f], /silent)
       filters[f] = clip(sxpar(h, 'FILTER'))
       dithers[f] = clip(sxpar(h, 'DITHER'))
    endfor  
    
    ;Unique list of filters
    filterlist  = unique(filters)
    ditherlist  = unique(dithers)

	;For each unique filter, combine sky files using skycombine if more than 2 files
	;Else return list of unprocessed files
    for i = 0, n_elements(filterlist)-1 do begin
        filt = filterlist[i]
        
        for j = 0, n_elements(ditherlist)-1 do begin
    	    dith = ditherlist[j]
    	    
       
       	    skyflats = where((filters eq filt) and (dithers eq dith), ctsky)
       	    if ctsky eq 0 then continue
       	    ;outflatname = pipevar.imworkingdir+'sky-'+filt+'.fits'
            print, '***********'
            print, filt
            print, dith
            print, files[skyflats]
        
            indfiles = files[skyflats]
            numfiles = n_elements(indfiles)
            skyframes = 4
            if dith eq 1 then continue
            for t = 0, n_elements(indfiles)-1 do begin
                outflatname = pipevar.imworkingdir+'sky-'+filt+strcompress(t,/REMOVE_ALL)+'.fits'
                if (t+skyframes) gt numfiles-1 then filelist = indfiles[numfiles-skyframes-1:numfiles-1] else filelist = indfiles[t:t+skyframes]
                skypipecombine_altsex, filelist, outflatname, filt, pipevar,/removeobjects, type='sky'
                stop
                ;skypipeproc, indfiles[t], outflatname
            endfor
        
        
       	;if ctsky ge 2 then begin
       	;	if pipevar.verbose gt 0 then print, filt, '-band sky flats.'
        ;   if file_test(outflatname) and pipevar.overwrite eq 0 then continue       
        ;    if pipevar.verbose gt 0 then print, files[skyflats]
            
        ;    skypipecombine_altsex, files[skyflats], outflatname, filt, pipevar,/removeobjects, type='sky'

        ;endif else begin
        
          	;print, 'Unable to produce a flat field for this setting: ' + filt
          	;noproc = where(filters eq filt, ctnoproc)
          	;print, 'Will not be able to further process ', clip(ctnoproc), ' images without a flat from another source:'
          	
          	;for j = 0, ctnoproc-1 do print, '   '+ files[noproc[j]]
          
       	;endelse
       	endfor
    endfor
    
	outpipevar = pipevar

end