;+
; NAME:
;	autopipeprepare
;
; PURPOSE:
;	Runs pipeprepare on every valid file and saves files with prefix 'p'.  Changes header
;	with more manageable keywords and does bias/dark subtraction if bias/dark master exists (compares
;	filter names in FILENAMEs of files and bias/dark master)
;
; OPTIONAL KEYWORDS:
;	outpipevar - output pipeline parameters
;	inpipevar  - input pipeline parameters (typically set in autopipedefaults.pro or ratautoproc.pro, but can be set to default) 
;
; EXAMPLE:
;	autopipeprepare, outpipevar=pipevar, inpipevar=pipevar
;
; DEPENDENCIES:
;	pipeprepare
;
; Written by Dan Perley 
; Modified by Vicki Toy 11/18/2013
; Modified by John Capone 04/22/2015
;
; FUTURE IMPROVEMENTS:
;	MOST REFER TO pipeprepare.pro: Need to check what additional keywords need to propagate, and check if values that are set with 
;	magic numbers can be set from existing keywords.
;	Check pipeprepare for RIMAS, RATIR, or VLT/VT to see changes that need to be made for RIMAS pipeline
;	Bias subtraction based on parameter file?
;-

pro autopipeprepare, outpipevar=outpipevar, inpipevar=inpipevar
	print, 'PREPARE'
	
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
	endelse 	
	
	;Looks for existing files in given data directory using prefix and tells user how many files found	
  	files  = findfile(pipevar.datadir+pipevar.prefix+'*.fits')
  	pfiles = findfile(pipevar.imworkingdir+'p'+pipevar.prefix+'*.fits')
  	
  	if pipevar.datadir ne '' then  begin
  	
     	if pipevar.verbose gt 0 then print, 'Looking for raw data at: ', pipevar.datadir+pipevar.prefix+'*.fits'

     	if n_elements(files) gt 0 then begin
        	if pipevar.verbose gt 0 then print, 'Found ', clip(n_elements(files)), ' files'
     	endif else begin
        	print, 'Did not find any files!  Check your data directory path!'
     	endelse
     	
 	endif

	;Creates list of spectroscopic standard stars
  	namefixfiles = ['']
  	catalogfile = 'catalog.txt'
  	if file_test(catalogfile) then namefixfiles = [namefixfiles, catalogfile]
  	landoltfieldposfile = pipevar.refdatapath+'/landoltpos.txt'
  	if file_test(landoltfieldposfile) then namefixfiles = [namefixfiles, landoltfieldposfile]
  	specstandardposfile = pipevar.refdatapath+'/specstandards.dat'
  	if file_test(specstandardposfile) then namefixfiles = [namefixfiles, specstandardposfile]
  	if n_elements(namefixfiles) gt 1 then namefixfiles = namefixfiles[1:*] else delvarx, namefixfiles
	
	
	;Finds any master bias files and filter name from header keyword
	;MUST CHANGE FOR RIMAS Assumes camera name is in filename before '.'
	biasfiles = findfile(pipevar.imworkingdir+'bias*', count=bct)
	biascamera = []
	if bct gt 0 then begin
	
		for i = 0, n_elements(biasfiles)-1 do begin
			header = headfits(biasfiles[i])
			camera = strmid(biasfiles[i], strpos(biasfiles[i], '.', /reverse_search)-1,1)
			biascamera = [biascamera,camera]
		endfor
		
	endif

	;Finds any master dark files and filter name from header keyword
	;MUST CHANGE FOR RIMAS Assumes camera name is in filename before '.'
	darkfiles = findfile(pipevar.imworkingdir+'dark*', count=dct)
	darkcamera = []
	if dct gt 0 then begin
	
		for i = 0, n_elements(darkfiles)-1 do begin
			header = headfits(darkfiles[i])
			camera = strmid(darkfiles[i], strpos(darkfiles[i], '.', /reverse_search)-1,1)
			darkcamera = [darkcamera,camera]
		endfor
		
	endif

	;For each file (that doesn't have an existing p file or can be overwritten), run pipeprepare on it with
	;output file being saved into the imworkingdir, will run bias subtraction if bias master available (checks
	;based on how bias file and data file are named (ASSUMES THAT FILTER NAME IS 1 CHARACTER BEFORE *.fits or after bias_* / dark_*
  	for f = 0, n_elements(files)-1 do begin
  	
     	if files[f] eq '' then continue
    
     	slashpos = strpos(files[f],'/',/reverse_search)    
     	fileroot = strmid(files[f],slashpos+1)
     	outnameim = pipevar.imworkingdir + 'p' + fileroot
     	matchi = where(outnameim eq pfiles, ct)
     	
		camera = strmid(files[f], strpos(files[f], 'o_', /reverse_search)-1,1)
		bcamloc = where(camera eq biascamera, brealcam)
		
		if brealcam gt 0 then begin
			biasfile = biasfiles[bcamloc]
			dcamloc = where(camera eq darkcamera, drealcam)
			if drealcam gt 0 then begin
				darkfile = darkfiles[dcamloc]
			endif else begin
				print, 'Did not find master dark file for camera ', camera, '!'
				darkfile=''
			endelse
		endif else begin
			biasfile=''
		endelse

     	if ct eq 0 or pipevar.overwrite gt 0 then begin
			pipeprepare, files[f], pipevar, outname=outnameim, namefixfiles=namefixfiles, biasfile=biasfile, darkfile=darkfile
     	endif
     	
  	endfor
  		
	outpipevar = pipevar
end
