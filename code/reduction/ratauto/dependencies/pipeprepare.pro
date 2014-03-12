;+
; NAME: 
;	pipeprepare
;
; PURPOSE: 
;	Adds additional header keywords needed later on in the pipeline
;   and removes unnecessary header keywords by looking through a list
;   of mandatory keywords.  Also runs bias subtraction for filters with
;	an existing master bias (CCDs).  
;
; INPUTS:
;	filename - name of the FITS file to prepare, or array of filenames, or file w/list of filenames
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;	outname  - specify output file to write to disk
;   namefix  - correct names using a catalog
;	biasfile - full name (including path) of master bias
;
; OPTIONAL OUTPUT KEYWORD PARAMETERS:
;   header
; 
; OUTPUTS:
;   No formal outputs. The prepared images are written to disk with outputname
;
; EXAMPLE: 
;	pipeprepare, filename, outname='p'+filename, namefixfiles=namefix.txt, header=readheader
;
; DEPENDENCIES:
;	IDL Astronomy User's Library routines (http://idlastro.gsfc.nasa.gov/) 
;
; Written by Dan Perley 
; Modified by Vicki Toy 11/20/2013
;
; FUTURE IMPROVEMENTS:
;	Need to check what additional keywords need to propagate, and check if values that are set with 
;	magic numbers can be set from existing keywords.  Possibly include prefix character into variable structure?
;	Master bias list in parameter file?
;	Change airmass keywords for RIMAS
;-

;----------------------------------------------------------------------------------------
pro pipeprepare, filename, outname=outname, namefixfiles=namefixfiles, header=header, biasfile=biasfile

; ---------- Process input filename(s), call recursively if necessary----------

	;Check for empty filename	
	if n_elements(filename) eq 0 then begin
   		print, 'No filename specified.'
   		return
	endif

	;Check for array input	
	if n_elements(filename) gt 1 then begin
   		files    = filename
	endif

	;Check for list input, if so read in filenames from file then check for wildcard input	
	if n_elements(filename) eq 1 then begin
   		filearr  = strsplit(filename,'.', /extract)
   		fileroot = filearr[0]
   		fileext  = filearr[1]
   		
   		if fileext eq 'cat' or fileext eq 'lis' or fileext eq 'list' or fileext eq 'txt' then begin
     		readcol, filename, files, format='a'
   		endif
	
   		starpos  = strpos(filename,'*')
   		qpos     = strpos(filename,'?')
   		
   		if starpos ge 0 or qpos ge 0 then begin
      		files = findfile(filename)
      		if n_elements(files) eq 1 then begin
        		if files[0] eq '' then print, 'Cannot find any files matching ', filename
        		return
      		endif
   		endif
	endif

	; Check for array filename, loop recursively if necessary	
	if n_elements(files) gt 1 then begin
  		for f = 0, n_elements(files)-1 do begin
    		pipeprepare, filename, outname=outname, namefixfiles=namefixfiles, header=header, timer=timer
  		endfor
  		return
	endif
	
	print, filename, format='($,A)'
	

; ---------- Read data and process header information ----------

	;Read in the data and header
	array = readfits(filename, header)
	
	;Fixes names if they match positional catalog 	
	if n_elements(namefixfiles) gt 0 then begin
   		namefix, header, namefixfiles
	endif

	;Add additional keywords and remove unnecessary ones
	;These keywords may need to be changed for RIMAS CHANGE VLT                       
	sxaddpar, header, 'OBSERVAT',"SPM", before='TELESCOP'
	
	;Grabs starting airmass, can alternatively use ETROBAM (ending time observed airmass) 
	;and saves as AIRMASS.  CHANGE FOR RIMAS
	am  = sxpar(header, 'STROBAM')
	sxaddpar, header, 'AIRMASS', am
	
   	filter = strtrim(sxpar(header, 'FILTER'))
   	
	;Target requested RA and DEC CHANGE for RIMAS VLT - ending time coordinates
	radeg  = sxpar(header,'ETRRQRA')
	decdeg = sxpar(header,'ETRRQDE')
	
	;These keywords may need to be changed for RIMAS CHANGE VLT
   	IF strcmp(filter,'J') OR strcmp (filter,'Z') then begin
  		sxaddpar, header, 'CRPIX1', 700. 
     	sxaddpar, header, 'CRPIX2', 900.
      	sxaddpar, header, 'CD1_1',  -8.17218901106E-05
      	sxaddpar, header, 'CD1_2',  3.7651099887E-06
      	sxaddpar, header, 'CD2_1',  4.20827225712E-06
      	sxaddpar, header, 'CD2_2',  8.26704009041E-05
   	endif else if strcmp(filter,'H') OR strcmp (filter,'Y') then begin
      	sxaddpar, header, 'CRPIX1', 200.
      	sxaddpar, header, 'CRPIX2', 900.
      	sxaddpar, header, 'CD1_1',  -8.17218901106E-05
      	sxaddpar, header, 'CD1_2',  3.7651099887E-06
      	sxaddpar, header, 'CD2_1',  4.20827225712E-06
      	sxaddpar, header, 'CD2_2',  8.26704009041E-05
   	endif else begin
      	sxaddpar, header, 'CD1_1', -8.80977078079E-05 
      	sxaddpar, header, 'CD1_2',  1.86753101419E-06 
      	sxaddpar, header, 'CD2_1',  1.86716671065E-06
      	sxaddpar, header, 'CD2_2',  8.81208878042E-05 
      	sxaddpar, header, 'CRPIX1', 512.
      	sxaddpar, header, 'CRPIX2', 512.
   	endelse
   		
   	sxaddpar, header, 'CRVAL1', float(radeg), 'RA (deg)'
   	sxaddpar, header, 'CRVAL2', float(decdeg), 'DEC (deg)'
   	sxaddpar, header, 'CTYPE1', 'RA---TAN'
   	sxaddpar, header, 'CTYPE2', 'DEC--TAN'

   	exptime = sxpar(header,'EXPTIME')
   	sxaddpar, header, 'ELAPTIME',exptime

	;List of mandatory header keywords
	mandatorykey = ['SIMPLE','BITPIX','NAXIS','NAXIS1','NAXIS2', $
 					'HISTORY','CTIME','USEC','JD','DATE-OBS','EXPOSURE', $
 					'EXPTIME','INSTRUME','OBSERVAT','TELESCOP','ORIGIN', $
 					'LATITUDE','LONGITUD','CCD_TYPE','CCD_SER','SATURATE',$
 					'RDNOISE','BINNING','BINY','BINX','WAVELENG','TARGNAME',$
 					'CAMERA','UTC','UT','ORIGOBJ','OBJECT','PXSCALE',$
 					'SUN_ALT','SMNSP','CD1_1','CD1_2','CD2_1','CD2_2',$
 					'CRPIX1','CRPIX2','CRVAL1','CRVAL2','CTYPE1','CTYPE2','ELAPTIME', $
 					'SOFTGAIN','FILTER','AVERAGE','STDEV','GAIN','AIRMASS','CCD_NAME']
 
 	;Finds list of unnecessary keywords, then deletes entire list
	unneckey = ''
	newheader = header
		
		
	h = strsplit(header,/EXTRACT)
	unneckey = strarr(n_elements(h))
	
	for hind = 0,n_elements(h)-2 do begin
	
		key = strmid(h[hind,0], 0, 8)
	
		if where(key eq mandatorykey) eq -1 then begin
			unneckey[hind] = key 
		endif
		
	endfor

	sxdelpar, newheader, unneckey

	;If biasfile keyword set subtract master bias from current file with given master bias file
	;If they are not the same size, quick program without saving with preparation prefix (will not move
	;on in following processing steps)
	if keyword_set(biasfile) then begin
		bias = readfits(biasfile)
		if n_elements(array) ne n_elements(bias) then begin
			print, filename + ' could not be bias subtracted because it is not the same size as the master bias, remove file to avoid confusion'
			return
		endif
		print, ' '
		print, '   bias subtracting'
		newdata = array-bias
	endif else begin
		newdata = array
	endelse

	;Write changes to disk
	outfilename = outname
	mwrfits, newdata, outfilename, newheader, /create
	print, ' -> ', outfilename

end