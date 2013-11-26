;+
; NAME:
;	flatpipeproc
;
; PURPOSE:
;	Checks if flat is same size as data, then divides for correct filter
;
; EXAMPLE:
;	flatpipeproc, filename, flatname, flatminval=0.3
;
; Written by Dan Perley 
; Modified by Vicki Toy 11/18/2013
;-

pro flatpipeproc, filename, flatname, flatminval=flatminval, flatmaxval=flatmaxval

	;Looks at extension of file and if it's a list of files, pulls them from file and stores as filenames
	if n_elements(filename) eq 1 then begin
  		filearr = strsplit(filename,'.', /extract)
  		fileext = filearr[1]
	endif
	
  	if fileext eq 'cat' or fileext eq 'lis' or fileext eq 'list' or fileext eq 'txt' then begin
    	filenames = grabcolumn(filename,0,/str)
  	endif else  begin
    	filenames = filename
	endelse
	
	;Reads in flat file and warns you if the flat is normalized
	nfile = n_elements(filenames)
	flat = mrdfits(flatname, 0, hflat, /silent)

	med = median(flat)
	if med lt 0.5 or med gt 2.0 then begin
  		print, 'Warning: flat is not normalized to zero'
	endif

	;For each file, read in fits
	for f = 0, nfile-1 do begin
   		data = mrdfits(filenames[f],0,h, /silent)
   		if n_elements(data) ne n_elements(flat) then begin
      		print, 'WARNING - Data and flat field are different dimensions!!'
         		; runs, but doesn't work; check this later (is this still true????????)
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
      		w = where(flat gt flatmaxval, ct)
       		if ct gt 0 then flat[w] = !values.f_nan
   		endif
 
 		;Divides out flattened field and adds keywords to header to show change	
   		fdata = data / flat 

   		sxaddpar, h, 'FLATFLD', flatname
   		skycts = median(fdata[goodsignal])
   		sxaddpar, h, 'SKYCTS', skycts, 'Sky counts'
   		if sxpar(h,'EXPTIME') gt 0 then $
      		sxaddpar, h, 'CTRATE', skycts/sxpar(h,'EXPTIME'), after='SKYCTS', 'Sky counts per second'

   		get_date, now
   		hist = 'Processed by flatproc '+now
   		sxaddhist, hist, h
   	
   		;Saves flat fielded fits file with new header 
   		lastslashpos = strpos(filenames[f],'/',/reverse_search)
   		if lastslashpos gt 0 then $
     		outname = strmid(filenames[f],0,lastslashpos) + '/' + 'f' + strmid(filenames[f],lastslashpos+1) $
   		else  $
     		outname = 'f'+filenames[f]

   		mwrfits, fdata, outname, h, /create
   	
	endfor

end
