;+
; NAME:
;	skypipeproc
;
; PURPOSE:
;	Subtracts sky flat from data and then subtracts median of that from remaining data.  
;	Adds 1000 offset afterwards.  Then crops and saves new fits file. NOTE: crop is specific
;	to RATIR
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
; Modified by Vicki Toy 11/18/2013
;
; FUTURE IMPROVEMENTS:
;	automated cropping of image?
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
	;stop, must be same dimensions).  If there is a minimum or maximum flat value set, forces values
	;outside of that range to NaN.  Use finite values above 0.1 to determine skycounts, and subtract out
	;flat along with median of flattened data (add 1000 offset).  Crops data based on hand chosen datapoints
	;NEEDS TO BE CHANGED FOR RIMAS VLT
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
   		
   		;Subtract out skyflat, and subtract the median of the subsequent flat subtracted data
   		;Add 1000 so no negative data? (unclear from original program, but just an offset so not crucial)
   		;Then calculate skycounts from data (above a minimum, or by default above 0.1)
   		fdata  = data  - flat 
   		tmp    = median(fdata)
   		fdata  = fdata - tmp
   		fdata  = fdata + 1000.
   		skycts = median(fdata[goodsignal])
   
   
   		;Adds header keywords to denote new median counts and that we flatfielded with a sky flat
   		sxaddpar, h, 'SFLATFLD', flatname
   		skycts = median(fdata[goodsignal])
   		sxaddpar, h, 'COUNTS', skycts
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

		;Different cropping regions for RATIR's CCDs vs. H2RGs CHANGE FOR RIMAS VLT
		;Extracts subset of data and saves to outfile with same altered header
   		s=size(fdata)
   		camera = strcompress(sxpar(h, 'WAVELENG'), /REMOVE_ALL)

   		if camera eq 'OPT' then begin
   			hextract,fdata,h,fdata,h,50,975,50,975
   		endif else begin
   			hextract,fdata,h,fdata,h,75,899,0,s(2)-1
   		endelse

   		mwrfits, fdata, outname, h, /create
	endfor

end
