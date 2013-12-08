;+
; NAME:
;	autopipestack
;
; PURPOSE:
;	Stacks images with same target and filter using SWarp.  Save files as coadd*_(FILTER).fits
;
; OPTIONAL KEYWORDS:
;	outpipevar - output pipeline parameters
;	inpipevar  - input pipeline parameters (typically set in autopipedefaults.pro or ratautoproc.pro, but can be set to default) 
;
; EXAMPLE:
;	autopipestack, outpipevar=pipevar, inpipevar=pipevar
;
; DEPENDENCIES:
;	SWarp
;
; Written by Dan Perley 
; Modified by Vicki Toy 12/08/2013
;
; FUTURE IMPROVEMENTS:
;	prefchar in variable structure?, header keywords to keep
;-

pro autopipestack, outpipevar=outpipevar, inpipevar=inpipevar

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
  
  	; if swarp configuration file is not present as 'default.swarp', have swarp output default configuration to this file name
  	if file_test('default.swarp') eq 0 then $
    	spawn, pipevar.swarpcommand+' -d > default.swarp'
 	
 	;Find files that have had astrometry performed on them, stop program if don't exist
 	;VLT CHANGE FOR RIMAS
  	prefchar = '2'
  	azffiles = findfile(pipevar.imworkingdir+'a*'+prefchar+'*_img_?.fits')
  	if azffiles[0] eq '' then return

  	filetargets   = strarr(n_elements(azffiles))
  	fileexposures = fltarr(n_elements(azffiles))
  	filefilters   = strarr(n_elements(azffiles))
  	filesatval    = fltarr(n_elements(azffiles))
  	fileskyval 	  = fltarr(n_elements(azffiles))
  	fileairval    = fltarr(n_elements(azffiles))
  	datestr = ''

	;Grab information in the headers of astrometry corrected file and save to array
  	for f = 0, n_elements(azffiles)-1 do begin
    	h = headfits(azffiles[f], /silent)
    	if f eq 0 then datestr = sxpar(h,'DATE')
    	
    	filetargets[f]   = repstr(strtrim(sxpar(h,'TARGNAME'),2),' ', '_')  ; spaces cause barfing in filenames
    	fileexposures[f] = sxpar(h, 'EXPOSURE')
    	filefilters[f]   = string(sxpar(h,'FILTER'))
    	filesatval[f]    = sxpar(h,'SATURATE')
    	fileskyval[f]    = sxpar(h,'COUNTS')
    	fileairval[f]    = sxpar(h,'AIRMASS')
    
  	endfor
  
  	targets = unique(filetargets)
  	for t = 0, n_elements(targets)-1 do begin
  	
  		;Finds files with same target and the filters associated with this target
    	target     = targets[t]
    	thistarget = where(filetargets eq target, cttarg)
    	if cttarg eq 0 then continue
    	thistargetfilts = unique(filefilters[thistarget])
    
    	;For each filter find files that have the same target and same filter 
    	;and store information on the exposure times and airmass
    	for l = 0, n_elements(thistargetfilts)-1 do begin
      		filter = thistargetfilts[l]
      		stacki = where(filetargets eq target and filefilters eq filter, ctstack)
      		if ctstack eq 0 then continue
      		
      		stacklist = azffiles[stacki]
      		stackexps = fileexposures[stacki]
      		
      		medianexp = median(stackexps)
      		medair    = median(fileairval[stacki])
      		minair    = min(fileairval[stacki])
      		maxair    = max(fileairval[stacki])
      		totalexp  = total(stackexps)
      		nstack    = n_elements(stacklist)

			;Create output variables that will be used by SWarp
      		outfile       = pipevar.imworkingdir + 'coadd' + strtrim(target,2) +'_'+ strtrim(filter,2) + '.fits'
      		outweightfile = pipevar.imworkingdir + 'coadd' + strtrim(target,2) +'_'+ strtrim(filter,2) + '.weight.fits'
      		stackcmd      = pipevar.swarpcommand+' '
      		
      		;Add each file to stack to the stack command separated by a comma
      		;and check if a weight file exists for the file, if it does not, check in imredux
      		;directory.  If the file is imredux/ then move it to current working directory with
      		;same name
      		for s = 0, n_elements(stacklist)-1 do begin
        		if s eq 0 then stackcmd = stackcmd + stacklist[s] else stackcmd = stackcmd + ',' + stacklist[s] 

       			weightfilename = strmid(stacklist[s],0,strlen(stacklist[s])-5-0)  + '.weight.fits'
        		weightexists   = file_test(weightfilename)
        
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

			;Add output names
      		stackcmd = stackcmd + ' -IMAGEOUT_NAME ' + outfile + ' -WEIGHTOUT_NAME ' + outweightfile
      		
      		;If there is more than one image to stack use relative weights, otherwise use background weight
      		if nstack gt 1 then begin
        		stackcmd = stackcmd + ' -WEIGHT_TYPE MAP_WEIGHT'
      		endif else begin
        		stackcmd = stackcmd + ' -WEIGHT_TYPE BACKGROUND'
      		endelse
      		
      		;Sets flux scale to the median exposure divided by the file's exposure (should be 1 since most files have same exposure time)
      		stackcmd = stackcmd + ' -FSCALE_DEFAULT '
      		for s = 0, n_elements(stacklist)-1 do stackcmd = stackcmd + strtrim(string(medianexp/stackexps[s]),2) + ','
      
      		;Keywords to keep CHANGE FOR RIMAS ALTER
      		stackcmd = stackcmd + ' -COPY_KEYWORDS OBJECT,TARGNAME,TELESCOP,FILTER,'+$
                             	'INSTRUME,OBSERVAT,ORIGIN,CCD_TYPE,JD,SOFTGAIN,'+$
                             	'WAVELENG,DATE-OBS,AIRMASS,FLATFLD,FLATTYPE'

			;If output file doesn't exist or overwrite then run SWarp with the stack command
			;and save data from output file with new header keywords
      		if (file_test(outfile) eq 0) or pipevar.overwrite then begin
       			if nstack eq 1 then print, 'Warning - only ', clip(nstack), ' exposures; not flagging bad pixels.'
        		print, 'Stacking ', target, ':'
        		for s = 0, n_elements(stacklist)-1 do print, stacklist[s], '  ', strtrim(string(stackexps[s]),2) + 's'
        		print, stackcmd
        		spawn, stackcmd
        		
        		;Create new header keywords with min, max, and median values (since stacked)
        		data = mrdfits(outfile, 0, h, /silent)
        		sxaddpar, h, 'DATE'    , datestr
        		sxaddpar, h, 'NSTACK'  , nstack
        		sxaddpar, h, 'AIRMASS' , medair, 'Median exposure airmass'
        		sxaddpar, h, 'AIRMIN'  , minair, 'Minimum exposure airmass'
        		sxaddpar, h, 'AIRMAX'  , maxair, 'Maximum exposure airmass'
        		sxaddpar, h, 'EXPOSURE', medianexp, 'Effective rescaled exposure time'
        		sxaddpar, h, 'TOTALEXP', totalexp, 'Total summed integration time'
        		sxaddpar, h, 'MAXEXP'  , max(stackexps), 'Length of longest exposure'
        		sxaddpar, h, 'MINEXP'  , min(stackexps), 'Length of shortest exposure'
        		sxaddpar, h, 'SATURATE', min(filesatval[stacki]-fileskyval[stacki])
        		sxaddpar, h, 'MEDSKY'  , median(fileskyval[stacki], /even)
        		
        		for s = 0, n_elements(stacklist)-1 do  sxaddpar, h, 'IMAGE'+strtrim(string(s),2), stacklist[s]
        		for s = 0, n_elements(stacklist)-1 do  sxaddpar, h, 'IMEXP'+strtrim(string(s),2), stackexps[s]
        
        		get_date, now
        		sxdelpar, h, ['SOFTNAME','SOFTVERS','SOFTDATE','SOFTAUTH','SOFTINST','AUTHOR']
        		hist = 'Processed by ratautoproc '+now
        		sxaddhist, hist, h, /COMMENT
        		mwrfits, data, outfile, h, /create
      		endif
    	endfor
  	endfor

	outpipevar = pipevar
	
end

