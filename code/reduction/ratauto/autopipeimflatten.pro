;+
; NAME:
;	autopipeimflatten
;
; PURPOSE:
;	Flatten data using flat with matching filter name
;
; OPTIONAL KEYWORDS:
;	outpipevar - output pipeline parameters
;	inpipevar  - input pipeline parameters (typically set in autopipedefaults.pro or ratautoproc.pro, but can be set to default) 
;
; EXAMPLE:
;	autopipeimflatten, outpipevar=pipevar, inpipevar=pipevar
;
; DEPENDENCIES:
;	flatpipeproc
;
; Written by Dan Perley 
; Modified by Vicki Toy 11/18/2013
;
; FUTURE IMPROVEMENTS:
;	prefchar in variable structure?
;-

pro autopipeimflatten, outpipevar=outpipevar, inpipevar=inpipevar

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

	;CHANGE FOR RIMAS VLT  	
	prefchar = '2'

	;Finds prepared files and checks to see if there are any existing flattened files
	;Find flats in imworkingdir with name flat somewhere in a fits file name
   	files  = findfile(pipevar.imworkingdir+'p'+prefchar+'*.fits')
   	ffiles = findfile(pipevar.imworkingdir+'fp'+prefchar+'*.fits')
   	flats  = findfile(pipevar.imworkingdir+'*flat*.fits')
   	
   	;End program if there are no prepared files
   	if n_elements(files) eq 1 and files[0] eq '' then return
   	
   	flatfilts = strarr(n_elements(flats))

   	;If there are flats, then grab the filter from each of them, otherwise set keyword that there are no flats
   	if flats[0] ne '' then begin
		for f = 0, n_elements(flats)-1 do begin
       		h = headfits(flats[f], /silent)
       		filter = clip(sxpar(h, 'FILTER'))
       		flatfilts[f] = filter
     	endfor
	endif
   	
   	;Create outfile name and check to see if outfile already exists.  If it doesn't or
   	;overwrite enabled then take filter from file and find where the flat filter matches
   	;If no flats match filter, store it in pipevar.flatfail, otherwise run flatproc.pro on file
	for f = 0, n_elements(files)-1 do begin
	
    	if files[f] eq '' then continue

      	outfile = fileappend(files[f], 'f')
      	match = where(outfile eq ffiles, ct)
      	if ct eq 0 or pipevar.overwrite then begin               
         	h = headfits(files[f], /silent)

         	filter = clip(sxpar(h, 'FILTER'))
         	binning='1'
         	flatfileno = where(flatfilts eq filter, flct)

         	if flct eq 0 then begin

            	print, 'Flat field not found for '+ removepath(files[f])+' (filter='+filter+')'
            	pipevar.flatfail = pipevar.flatfail +' '+ files[f]
            	continue
            	
         	endif
         
         	flatfile = flats[flatfileno[0]]
         	print, 'Flattening ', removepath(files[f]), ' using ', removepath(flatfile)
         	flatpipeproc, files[f], flatfile, flatminval=0.3

     	endif
	endfor
	
	outpipevar = pipevar
 end