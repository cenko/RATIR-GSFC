;+
; NAME:
;	autopipeprepare
;
; PURPOSE:
;	Runs pipeprepare on every valid file and saves files with prefix 'p'
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
;
; FUTURE IMPROVEMENTS:
;	MOST REFER TO pipeprepare.pro: Need to check what additional keywords need to propagate, and check if values that are set with 
;	magic numbers can be set from existing keywords.  Possibly include prefix character into variable structure?
;	Check pipeprepare for RIMAS, RATIR, or VLT/VT to see changes that need to be made for RIMAS pipeline
;-

pro autopipeprepare, outpipevar=outpipevar, inpipevar=inpipevar

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
	
	;Prefix character for data files CHANGE NEEDED FOR RIMAS VT 
	;Looks for existing files in given data directory using prefchar and tells user how many files found
	
  	prefchar = '2'
	
  	files  = findfile(pipevar.datadir+prefchar+'*.fits')
  	pfiles = findfile(pipevar.imworkingdir+'p'+prefchar+'*.fits')

  	if pipevar.datadir ne '' then  begin
  	
     	print, 'Looking for raw data at: ', pipevar.datadir+prefchar+'*.fits'

     	if n_elements(files) gt 0 then begin
        	print, 'Found ', clip(n_elements(files)), ' files'
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
	
	;For each file (that doesn't have an existing p file or can be overwritten), run pipeprepare on it with
	;output file being saved into the imworkingdir
	
  	for f = 0, n_elements(files)-1 do begin
  	
     	if files[f] eq '' then continue
    
     	slashpos = strpos(files[f],'/',/reverse_search)    
     	fileroot = strmid(files[f],slashpos+1)
     	outnameim = pipevar.imworkingdir + 'p' + fileroot
     	matchi = where(outnameim eq pfiles, ct)

     	if ct eq 0 or pipevar.overwrite gt 0 then begin

     		if pipevar.modestr eq 'im' then outname = outnameim
			
			pipeprepare, files[f], outname=outname, namefixfiles=namefixfiles
			
     	endif
     	
  	endfor
  	
	outpipevar = pipevar
end
