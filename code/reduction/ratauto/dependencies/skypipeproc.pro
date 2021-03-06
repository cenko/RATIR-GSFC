;+
; NAME:
;	skypipeproc
;
; PURPOSE:
;	Subtracts sky flat from data and then subtracts median of that from remaining data.  
;	Then crops and saves new fits file.
;
; INPUTS:
;	filename - file or list of files to be sky subtracted
;	flatname - sky flat fits file
;
; OPTIONAL KEYWORDS:
;	flatminval - minimum required value in flat (default for skycounts calculation is 0.1)
;	flatmaxval - maximum required value in flat
;
; EXAMPLE:
;	skypipeproc, filename, flatname
;
; Written by Dan Perley 
; Modified by Vicki Toy 8/15/14
;-

pro skypipeproc, filename, flatname, flatminval=flatminval, flatmaxval=flatmaxval

	;If only one object in filename check to see if it's a list. If so grab filenames from list.
	;Otherwise filenames equal the input given.
	if n_elements(filename) eq 1 then begin
  		filearr = strsplit(filename,'.', /extract)
  		fileroot = filearr[0]
  		fileext = filearr[1]

  		if fileext eq 'cat' or fileext eq 'lis' or fileext eq 'list' or fileext eq 'txt' then begin
    		filenames = grabcolumn(filename,0,/str)
  		endif else  begin
    		filenames = filename
  		endelse
  		
	endif else begin
  		filenames = filename
	endelse

	;Open the given flat file and read in as fits.  Find median of flat
	nfile = n_elements(filenames)
	flat = mrdfits(flatname, 0, hflat, /silent)

	med = median(flat)

	;For each input file, read in as fits and check if same size as flats (if it isn't program will
	;stop).  If there is a minimum or maximum flat value set, forces values outside of that range to NaN.  
	;Use finite values above 0.1 to determine skycounts, and subtract out flat along with median of flattened data. 
	;Then saves to new fits file
	for f = 0, nfile-1 do begin
   		data = mrdfits(filenames[f],0,h, /silent)
   		
   		if n_elements(data) ne n_elements(flat) then begin
      		print, 'WARNING - Data and flat field are different dimensions!!'
			stop
   		endif

   		;Set values too low/high to NaNs
   		if n_elements(flatminval) gt 0 then begin
      		w = where(flat lt flatminval, ct, complement=goodsignal)
      		if ct gt 0 then flat[w] = !values.f_nan
   		endif else begin
      		w = where(flat lt 0.1, complement=goodsignal) ; for sky count determination only
   		endelse
   		
   		if n_elements(flatmaxval) gt 0 then begin
      		w = where(pflat gt flatmaxval, ct)
       		if ct gt 0 then flat[w] = !values.f_nan
   		endif
   		
   		;Scale skyflat, then subtract scaled skyflat, and subtract the median of the subsequent flat subtracted data
   		;Then calculate skycounts from data (above a minimum, or by default above 0.1)
   		flattmp = median(flat)
   		imgtmp  = median(data)
   		
   		scalefr = imgtmp/flattmp
   		fdata   = data - scalefr * flat 
   		tmp    = median(fdata)
   		fdata  = fdata - tmp
   		
   		skycts = median(fdata[goodsignal])
   
   		;Adds header keywords to denote new median counts and that we flatfielded with a sky flat
   		sxaddpar, h, 'SFLATFLD', flatname
   		skycts = median(fdata[goodsignal])
   		sxaddpar, h, 'SKYCTS', skycts, 'Sky counts'
   
   		if sxpar(h,'EXPTIME') gt 0 then $
      		sxaddpar, h, 'CTRATE', skycts/sxpar(h,'EXPTIME'), after='SKYCTS', 'Sky counts per second'

   		get_date, now
   		hist = 'Processed by skyproc '+now
   		sxaddhist, hist, h
   
   		lastslashpos = strpos(filenames[f],'/',/reverse_search)
   
   		;Creates outname with an added 's' prefix
   		if lastslashpos gt 0 then $
     		outname = strmid(filenames[f],0,lastslashpos) + '/' + 's' + strmid(filenames[f],lastslashpos+1) $
   		else  $
     		outname = 's'+filenames[f]

   		mwrfits, fdata, outname, h, /create
	endfor

end
