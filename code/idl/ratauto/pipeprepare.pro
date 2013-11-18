;+
; NAME: 
;	pipeprepare
;
; PURPOSE: 
;	Adds additional header keywords needed later on in the pipeline
;   and removes unnecessary header keywords by looking through a list
;   of mandatory keywords.
;
; INPUTS:
;	filename - name of the FITS file to prepare, or array of filenames, or file w/list of filenames
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;	outname  - specify output file to write to disk
;   namefix  - correct names using a catalog
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
; Modified by Vicki Toy 11/18/2013
;
; FUTURE IMPROVEMENTS:
;	Need to check what additional keywords need to propagate, and check if values that are set with 
;	magic numbers can be set from existing keywords.  Possibly include prefix character into variable structure?
;	Check pipeprepare for RIMAS, RATIR, or VLT/VT to see changes that need to be made for RIMAS pipeline
;-

;----------------------------------------------------------------------------------------
pro pipeprepare, filename, outname=outname, namefixfiles=namefixfiles, header=header

; ---------- Process input filename(s), call recursively if necessary----------

	;Check for empty filename	
	if n_elements(filename) eq 0 then begin
   		print, 'No filename specified.'
   		return
	endif

	;Check for array input	
	if n_elements(filename) gt 1 then begin
   		files = filename
	endif

	;Check for list input, if so read in filenames from file then check for wildcard input	
	if n_elements(filename) eq 1 then begin
   		filearr = strsplit(filename,'.', /extract)
   		fileroot = filearr[0]
   		fileext = filearr[1]
   		
   		if fileext eq 'cat' or fileext eq 'lis' or fileext eq 'list' or fileext eq 'txt' then begin
     		readcol, filename, files, format='a'
   		endif
	
   		starpos = strpos(filename,'*')
   		qpos = strpos(filename,'?')
   		
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
	xbin = 1                        
	ybin = 1                        
	pxscale = 0.3
	latitude = 30.7500	

	sxaddpar, header, 'XBIN', xbin, after='BINNING'
	sxaddpar, header, 'YBIN', ybin, after='BINNING'
	sxaddpar, header, 'PXSCALE', pxscale, 'arcsec/pix'
	sxaddpar, header, 'OBSERVAT',"SPM", before='TELESCOP'
	sxaddpar, header, 'RDNOISE', 4, before='BINNING'
	sxaddpar, header, 'AIRMASS', 1.
	sxaddpar, header, 'LTM1_1', 1./xbin
	sxaddpar, header, 'LTM2_2', 1./ybin
	
	;Sun/moon ephemeris
	jd = sxpar(header, 'JD')
	lst = sxpar(header, 'LST')

	;Calculate sun and moon altitude from start of exposure. The middle (or "worst") would be more useful.
	sunpos, jd, sunra, sundec
	hadec2altaz, lst-sunra, sundec, latitude, sunalt, sunaz

	moonpos, jd, moonra, moondec
	hadec2altaz, lst-moonra, moondec, latitude, moonalt, moonaz
	
	;Target requested RA and DEC CHANGE for RIMAS VLT
	radeg = sxpar(header,'ETRRQRA')
	decdeg = sxpar(header,'ETRRQDE')

	;Calculate separation between moon and target
	gcirc, 2, radeg, decdeg, moonra, moondec, moondist
	moonsep  = moondist/3600

	sxaddpar, header, 'SUNELEV', sunalt, (sunalt>0.?'Above horizon':'Below horizon')
	sxaddpar, header, 'MOONELEV', moonalt, (moonalt>0?'Above horizon':'Below horizon')
	sxaddpar, header, 'MOONDIST', moonsep

   	filter = strtrim(sxpar(header, 'FILTER'))

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
 					'FOC_NAME','SCRIPREP','SCRIPT','SCR_COMM','COMM_NUM', $
 					'CCD_TYPE','CCD_SER','SATURATE','RDNOISE','BINNING', $
 					'YBIN','XBIN','BINX','WAVELENG','TARGNAME','CAMERA', $
 					'UTC','UT','ORIGOBJ','OBJECT','PXSCALE','SUNELEV','MOONELEV', $
 					'MOONDIST','CD1_1','CD1_2','CD2_1','CD2_2','CRPIX1', $
 					'CRPIX2','CRVAL1','CRVAL2','CTYPE1','CTYPE2','ELAPTIME', $
 					'SOFTGAIN','LTM1_1','LTM2_2','FILTER']
 
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

	;Write changes to disk
	outfilename = outname
	mwrfits, array, outfilename, newheader, /create
	print, ' -> ', outfilename

end