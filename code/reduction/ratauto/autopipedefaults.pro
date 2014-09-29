;+
;NAME:
;   autopipedefaults
;PURPOSE:
;   Sets commonly used variables for pipeautoproc to use throughout each step
;   Uses pipeautoproc.par to set variables, otherwise set to default values
;	Saves in a structure
;OPTIONAL KEYWORDS:
;	outpipevar - output pipeline parameters
;	inpipevar  - input pipeline parameters (typically set in autopipedefaults.pro or ratautoproc.pro, but can be set to default) 
;EXAMPLE:
;   autopipedefaults, outpipevar=pipevar, inpipevar=pipevar
;-

pro autopipedefaults, outpipevar=outpipevar, inpipevar=inpipevar

	print, 'Setting pipeline parameters (DEFAULTS)'
	
	;Setup pipeline variables that carry throughout the pipeline
	if keyword_set(inpipevar) then begin
		pipevar = inpipevar
		if pipevar.verbose gt 0 then print, 'Using provided pipevar'
	endif else begin
		pipevar = {autoastrocommand:'autoastrometry', getsedcommand:'get_SEDs', $
					sexcommand:'sex' , swarpcommand:'swarp' , $
					prefix:'', datadir:'' , imworkingdir:'' , overwrite:0 , verbose:0, rmifiles:0,$
					flatfail:'' , catastrofail:'' , relastrofail:'' , fullastrofail:'' , $
					pipeautopath:'' , refdatapath:'', defaultspath:'' }
	endelse

   	paths 	= strsplit(!path,':',/extract)
   	
   	;Finds the directory that 'ratautoproc.pro' is in
   	;saves it to pipevar.pipeautopath  	
   	for p 	= 0, n_elements(paths)-1 do begin
      	n 	= countfiles(paths[p]+'/ratautoproc.pro')
      	if n gt 0 then pipevar.pipeautopath = paths[p]
      	if n gt 0 then break    	
   	endfor
   	
   	;Error message if could not find directory CHANGE TO STOP PROGRAM VLT (change stop?)
   	  	
   	if n_elements(pipevar.pipeautopath) eq 0 then begin
    	print, 'Problem finding pipeautoproc directory in IDL paths. Redefine paths to include appropriate directory'
		stop
   	endif
   	
   	;Finds parameter file and sets common block variables from the parameter file
   	;(checks current directory, pipevar.pipeautopath, pipevar.pipeautopath/defaults in that order)
   	
   	for p 	= -1, 1 do begin
   		if p eq -1 then path = '.'
      	if p eq  0 then path = pipevar.pipeautopath
      	if p eq  1 then path = pipevar.pipeautopath + '/defaults'
      	sfile = path+'/pipeautoproc.par'       ; get the configuration file
      	n = countfiles(sfile)
      	
      	if n gt 0 then begin
    		openr, 5, sfile
         	iline = ''
         	
         	while not eof(5) do begin
            	readf, 5, iline
            	colonpos 	= strpos(iline,':')
            	if colonpos lt 1 or colonpos ge strlen(iline) then continue
            	varname 	= clip(strmid(iline,0,colonpos))
            	varset 		= clip(strmid(iline,colonpos+1))
            	if varname eq 'refdatadir' 		 then pipevar.datapath = varset
            	if varname eq 'defaultsdir' 	 then pipevar.defaultspath = varset
            	if varname eq 'autoastrocommand' then pipevar.autoastrocommand = varset
            	if varname eq 'getsedcommand'	 then pipevar.getsedcommand = varset
            	if varname eq 'swarpcommand' 	 then pipevar.swarpcommand = varset
            	if varname eq 'sexcommand' 		 then pipevar.sexcommand = varset
            	if varname eq 'prefix'			 then pipevar.prefix = varset
            	if varname eq 'datadir' and pipevar.datadir eq '' then pipevar.datadir = varset
            	if varname eq 'imworkingdir' 	 then pipevar.imworkingdir = varset
         	endwhile
         	
         	free_lun,5 
			break
		endif
   	endfor

   	;Use internal defaults to set anything still unspecified
   	
	if pipevar.refdatapath eq '' then pipevar.refdatapath = pipevar.pipeautopath+'/refdata'
   	if pipevar.defaultspath eq '' then pipevar.defaultspath = pipevar.pipeautopath+'/defaults'
	
   	setswarppath, extractpath(pipevar.swarpcommand)
   	setsexpath, extractpath(pipevar.sexcommand)

  	if strlen(pipevar.imworkingdir) gt 0 and file_test(pipevar.imworkingdir) eq 0 then begin
     	print, 'Creating imaging working directory ', pipevar.imworkingdir
     	spawn, 'mkdir '+pipevar.imworkingdir
     	check = file_test(pipevar.imworkingdir)
     	
     	if check eq 0 then begin 
       		print, 'WARNING - Failed to create image working directory!'
       		print, '          Using current directory.  Check pipeautoproc.par inputs.'
     	endif
  	endif

	outpipevar = pipevar

end