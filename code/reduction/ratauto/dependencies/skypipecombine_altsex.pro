;+
; NAME:
;	skypipecombine
;
; PURPOSE:
;	Create sigma clipped median sky flat.  Scales each file based on the overall median (sigma clipped),
;	then will remove objects in each file that are well beyond median sigma clipped value (after removing
;	those beyond saturation and count limits) for 'r' and 'i' bands this is done with sextractor selection 
;	instead (2x 100% flux fraction radius). Then takes sigma clipped median of each pixel and save as with outfile name.
;
; INPUT:
;	filelist - files to be processed
;	outfile  - name for output fits file
;	filt	 - filter of files
;	pipevar  - pipeline parameters in structure
;
; KEYWORDS:
;	removeobjects 	- specifies if you want objects remove, need to run if you want median or mean trimmed
;	mediani 		- sets algorithm to median (default)
;	meani			- sets algorithm to mean
;	objthresh		- sets sigma to determine if point source (default is 6) taken after sigma clipping
;	trimlo			- how much to trim off bottom of data in mean algorithm mode (default is 25%)
;	trimhi			- how much to trim off top of data in mean algorithm mode (default is 25%)
;	mincounts		- sets minimum counts allowed (default is 1)
;	maxcounts		- sets maximum counts allowed (default is 55000)
;	satlevel		- sets saturation level (default is 30000)
;	type			- sets 'SKYTYPE' keyword in header of outfile to this string
;
; EXAMPLE:
;	skypipecombine, list, 'outfile.fits', /removeobjects, type='sky'
;
; DEPENDENCIES:
;	sigclipstats_vt.pro, djs_iterstat.pro
;
; Written by Dan Perley 
; Modified by Vicki Toy 11/18/2013
;-

pro skypipecombine_altsex, filelist, outfile, filt, pipevar, removeobjects=removeobjects, objthresh=objthresh, $
	mediani=mediani, meani=meani, trimlo=trimlo, trimhi=trimhi, mincounts=mincounts, $
	maxcounts=maxcounts, type=type, satlevel=satlevel
	
	;Set algorithm
	algorithm = 'median' ; default
	if keyword_set(mediani) then algorithm = 'median'
	if keyword_set(meani)   then algorithm = 'mean' 
	
	;Sets defaults for trimming (25% of list), counts (1-55,000), and saturation level and radius
	if n_elements(trimlo)    eq 0 then trimlo    = (n_elements(filelist)+1)/4
	if n_elements(trimhi)    eq 0 then trimhi    = trimlo
	if n_elements(mincounts) eq 0 then mincounts = 1
	if n_elements(maxcounts) eq 0 then maxcounts = 55000
	if n_elements(satlevel)  eq 0 then satlevel	 = 30000.

	;If outfile isn't set, then just save it with same root as filelist with fits extension
	if n_elements(outfile) eq 0 then begin
  		filearr = strsplit(filelist,'.', /extract)
  		fileroot = filearr[0]
  		outfile = fileroot + '.fits'
	endif

	;If given list, then grab all filenames, saved to files
	if n_elements(filelist) eq 1 then begin
  		files = grabcolumn(filelist,0,/str)
  		files = files[where(files ne '')]
	endif else begin
  		files = filelist
	endelse

	nfiles = n_elements(files)
	nmid = nfiles/2
	
	;Read in middle file and initialize arrays
	refdata = mrdfits(files[nmid], 0, h, /silent, /fscale)
	s  = size(refdata)
	nx = s[1]
	ny = s[2]

	data     = fltarr(nx, ny, nfiles) + !values.f_nan
	inmeds   = fltarr(nfiles)
	usefiles = strarr(nfiles)
	skymed   = fltarr(nfiles)

	z = 0
	keptfiles = []
	;For each file and make sure size matches middle file, calculate sigma clipped median(3sig, 6 iter), then if within 
	;counts limit save data into 3d data cube and save clipped median into skymed, and mark file as usable
	;Increment z by one when this is true
	for f = 0, nfiles-1 do begin
   		indata = mrdfits(files[f], 0, h, /silent, /fscale)
   		ins = size(indata)
   		inx = ins[1]
   		iny = ins[2]
   
   		if inx ne nx or ny ne ny then begin
      		print, 'File ', files[f], ' has wrong dimensions ('+clip(inx)+'x'+clip(iny)+'; should have '+clip(nx)+'x'+clip(ny)+')'
      		continue
   		endif

		;Perform 3 sigma clipped median and save to inmeds
   		sigclipstats_vt, indata, median=inmed, sigmahi=3, sigmalo=-3
   		inmeds[f] = inmed
   		
   		;If median is within limits save data, otherwise exclude files
   		if inmed ge mincounts and inmed le maxcounts then begin
      		
      		print, '  ' + removepath(files[f]) + ' (' + strtrim(long(inmed),2) + ' counts/pix)'
      		skymed[z]   = inmed
      		data[*,*,z] = float(indata)
      		usefiles[z] = files[f]
      		z = z + 1

   		endif else begin
      		if inmed lt mincounts then $
      			print, '  ' + removepath(files[f]) + ' (' + strtrim(long(inmed),2) + ' counts/pix) - too few counts; excluding'
      		if inmed gt maxcounts then $
      			print, '  ' + removepath(files[f]) + ' (' + strtrim(long(inmed),2) + ' counts/pix) - too many counts; excluding'
   		endelse

		keptfiles = [keptfiles, files[f]]
	endfor
	
	;Need more than one sky file to make a skyflat
	if z lt 2 then begin
    	print, 'ERROR - Not enough counts to make a flat with these data!'
    	return
	endif
	
	;Median of sigma clipped medians
	medsky=median(skymed, /even)
	
	;Multiply each fits file data by the median of sigma clipped medians divided by the median of the data (NOT sigma clipped)
	;Simply scales each data set to the median of all datasets
	for f = 0, z-1 do begin
   		inmed = median(data[*,*,f], /even) ; if each flat has a changing sky background
   		factor=medsky/inmed
   		data[*,*,f]=data[*,*,f]*factor(0)
	endfor

	;Chops off empty data
	if z ne nfiles then data = data[*,*,0:z-1]

	;Removes objects from field if keyword set by calculating median sigma clipping (5 sigma, 6 iter) and using
	;the calculated stddev to remove 6sigma (or non-default object threshold) data from the median along with 
	;values above the saturation limit.
	;Then find 3 sigma clipped median of each pixel with remaining values or mean of middle 50% of data
	if keyword_set(removeobjects) then begin
   		print, '  Identifying objects...'

		;Sets object threshold (sigma limit) and buffer
   		if n_elements(objthresh) eq 0 then objthresh = 6
   
   		;For each file within count limits find 5 sigma clipped stats and set data that are
   		;6 sigma (objethresh) from median to NAN and data beyond saturation level to NAN
   		for f = 0, z-1 do begin

      		indata = data[*,*,f]
      		satmask = fltarr(nx, ny)
      		
      		;Find sources from sextractor for RATIR's CCDs, 
      		;otherwise use previous source extraction method using iterative sigma clipping
      		if (filt eq 'i') or (filt eq 'r') then begin
      		
      			;Runs sextractor with special commands for quick and dirty source extraction
      			print,pipevar.sexcommand + ' ' +files[f]+  ' -c sex_source.config'
      			spawn,pipevar.sexcommand + ' ' +files[f]+  ' -c sex_source.config'
      			readcol, 'skysource.cat', x,y,ra,dec,mag,magerr,rad, fwhms,ell, cstar, flag
      			
      			;Finds sources with reasonable FWHM (i.e. >1 pixel) and no raised flags, 
      			;except if source near edge, still want to sky subtract without bright source
      			realsources = where(fwhms gt 1 and (flag eq 0 or flag eq 16))

				rad = rad[realsources]
				x   = x[realsources]
				y   = y[realsources]

      			sdata = size(indata)
      			scale = 1.0
				
				;For each good source, set everything within given radius to NAN
      			for so=0,n_elements(rad)-1 do begin
      		
      				;Source removal using box with radius from sextractor
      				rd = abs(rad[so])*scale
      				minx = x[so]-rd
      				maxx = x[so]+rd
      				miny = y[so]-rd
      				maxy = y[so]+rd
 
 					;If at an edge just go to the edge of the image
 					if minx lt 0 then minx = 0
 					if miny lt 0 then miny = 0
 					if maxx gt sdata[1]-1 then maxx = sdata[1]-1
      				if maxy gt sdata[2]-1 then maxy = sdata[2]-1
      			
      				indata[ minx:maxx, miny:maxy ] = !Values.F_NAN
      				
      			endfor

      		endif else begin
      		
      			sigclipstats_vt, indata, sigmahi=5, sigmalo=-5, median=datamed, stdevi=datastdev
      			sourcepixels = where(indata ge datamed + objthresh * datastdev, ctsourcepix)
      
      			;Make a mask from objects that are 6 sigma (or non-default object threshold) 
      			;from the median (from sigma clipped stats) and set to NAN
      			if ctsourcepix gt 0 then begin
         			indata[sourcepixels] = !Values.F_NAN
      			endif
      			
      		endelse
      	
      		ctsatpix = 0
      	
      		;Find data that is beyond saturation limit and set to NAN
      		if n_elements(satlevel) gt 0 then begin
         		satpixels = where(indata ge satlevel, ctsatpix)
         	
         		if ctsatpix gt 0 then begin
            		satmask[satpixels] = 1
            		indata[nearsatpixels] = !Values.F_NAN
         		endif
      		endif

      		data[*,*,f] = indata
      	
		endfor

   		reflat = fltarr(nx, ny)

		;If algorithm set to median, find 3 sigma clipped median of each pixel 
		;(excluding saturated or above 6 sigma - NAN values which are eventually set to 1)
   		if algorithm eq 'median' then begin
     		print, '  Median-combining...'
     		i = 0l
     	
     		for y = 0, ny-1 do begin
       			for x = 0, nx-1 do begin
         			i = i + 1
         			vector=data[x,y,*]
         			tmp=where(finite(vector))
         			djs_iterstat, vector(tmp), sigrej=3, median=me
         			reflat[x,y] = me
       			endfor		
       		endfor
       	
       		;Replace bad pixels with 1
     		bad = where(finite(reflat) eq 0, ct)
     		if ct gt 0 then reflat[bad] = 1
     	
		endif

		;If algorithm set to mean, takes mean of sorted values that have been trimmed
		;default is to trim 25% off top and bottom, if not enough good data, set trimming to 0
   		if algorithm eq 'mean' then begin
     		print,  '  Combining via trimmed mean...'
     		i = 0l
     
     		for y = 0, ny-1 do begin
       			for x = 0, nx-1 do begin
         			slice = data[x,y,*]
         			good = where(finite(slice) eq 1, ctgood)
         
         			if ctgood eq 0 then begin
             			reflat[x,y] = 1
         			endif else begin
             			slice = slice[good]
             			itrimlo = trimlo
             			itrimhi = trimhi
             	
             			while ctgood-itrimlo-itrimhi lt 1 do begin
               				itrimlo = (itrimlo - 1) > 0
               				itrimhi = (itrimhi - 1) > 0
             			endwhile
             	
             			slice = slice[sort(slice)]
             			slice = slice[itrimlo:ctgood-1-itrimhi]
             			reflat[x,y] = mean(slice)
         			endelse
         	
       			endfor
     		endfor
		endif

   		flat = reflat

	endif

	;Adds header information to signify what files we used 
	for f = 0, z-1 do begin
  		sxaddpar, h, 'SKY'+strtrim(string(f),2), removepath(usefiles[f])
	endfor

	if keyword_set(type) then sxaddpar,  h, 'SKYTYPE', type

	get_date, now
	hist = 'Processed by flatcombine '+now
	sxaddhist, hist, h

	;Saves new fits files with new sky flat
	mwrfits, flat, outfile, h, /create

	print, '  Written to ', outfile
	
end
