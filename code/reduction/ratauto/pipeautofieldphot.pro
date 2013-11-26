;Finds the minimum (cartesian) distance from one object to another in the same
;array.  Returns an array of the minimum distance of that indices' object to any other object
function nndist, x, y
	allmindist = fltarr(n_elements(x))
   	
   	for i = 0L, n_elements(x)-1 do begin
    	mindist = 1e9
      
      	for j = 0, n_elements(x)-1 do begin
         	if i eq j then continue
         	dist = sqrt((x[i]-x[j])^2+(y[i]-y[j])^2)
         	if dist lt mindist then mindist = dist
      	endfor
      	
      	allmindist[i] = mindist
      	
   	endfor
   	return, allmindist
end


function findneighbors, x, y, i, mindist=mindist, count=count
	neighbors = [-1]
    for j = 0L, n_elements(x)-1 do begin
		if i eq j then continue
        dist = sqrt((x[i]-x[j])^2+(y[i]-y[j])^2)
        if dist lt mindist then neighbors = [neighbors,j]
    endfor
      
    count = n_elements(neighbors)-1
    if count gt 0 then neighbors = neighbors[1:*]
    return, neighbors
end

function nearest, x, y, xarr, yarr, mindist=mindist, count=count

    count = 0
    index = -1
    imindist = mindist

    for j = 0L, n_elements(xarr)-1 do begin
    	dist = sqrt((x-xarr[j])^2+(y-yarr[j])^2)
        if dist lt imindist then begin
        	imindist = dist
            index = j
            count = count + 1
        endif
    endfor

    return, index
end

;Find maximum pixel within +/- r range in both x and Y, for r=1 that is a 5x5 square
;Return maximum pixel value for each file within +/- r
function maxpixel, imagename, x, y, r=r
	data = mrdfits(imagename, /silent)
   	if n_elements(r) eq 0 then r = 2
   	nx = (size(data)) [1]
   	ny = (size(data)) [2]
   	m  = fltarr(n_elements(x))
   	
   	for i = 0, n_elements(x)-1 do begin
     	xi = fix(x[i])
     	yi = fix(y[i])
     	if xi-r ge nx or yi-r ge ny then m[i] = !values.f_nan else $
     		m[i] = max(data[(xi-r)>0:(xi+r)<(nx-1),  (yi-r)>0:(yi+r)<(ny-1)])
   	endfor

   	return, m
end



pro fieldphot, imagename, ra, dec, mag, magerr, zpt1=zpt1, airtermi=airtermi, rinner=rinner, router=router, seeing=seeing, relative=relative, getseeing=getseeing

; r1 is the optimum aperture.
; r2 is a large aperture that can be treated as infinite.

	common lrisconfig, lrisautopath, refdatapath, defaultspath

	;VLT USES DIFFERENT PARAM FILES THAN ASTROMETRY, CHANGE TO MATCH?  MUST ADD FLAG
	;Check if sextractor files are in current folder, if not save them
	if file_test('temp.param') eq 0 then spawn, 'cp ' + defaultspath+'/temp.param ./'
	if file_test('sex.config') eq 0 then spawn, 'cp ' + defaultspath+'/sex.config ./'
	if file_test('sex.conv')   eq 0 then spawn, 'cp ' + defaultspath+'/sex.conv ./'

	;Make sure air term is positive
	if airtermi lt 0 then airterm = abs(airtermi[0]) else airterm = airtermi[0]

	;Read in header and grab keywords
	data = mrdfits(imagename,0,h,/silent)
	filter = sxpar(h,'FILTER')
	exptime = sxpar(h,'EXPTIME')
	
	;CHANGE TO VARIABLE VALUES VLT
	saturate = 100000.
	airmass = 1

	sexpath = getsexpath()

	;Remove sextractor catalog if one already exists, then run sextractor
	if file_test('temp.cat') ge 1 then spawn, 'rm -f temp.cat'
	cmd = sexpath+'sex '+imagename+' -c sex.config'
	spawn, cmd

	;If sextractor catalog not made, stop program to alert user
	if file_test('temp.cat') eq 0 then begin
   		print, 'Sextractor failed on '+ imagename
   		print, 'Command was ', cmd
   		stop
   		return
	endif

	;Read in sextractor catalog
	readcol, 'temp.cat', x, y, ra, dec, sexmag, sexmagerr, ellip, fwhm, /silent

	;Find the minimum distance to other objects in same catalog (returns array with minimum
	;distance, where the index specifies which object it refers to) 
	mindist = nndist(x,y)

	skypix = sxpar(h,'COUNTS')
	
	;Finds maximum pixel value in image in 5x5 cube around x, y coordinates
	maxpix = maxpixel(imagename, x,y, r=2)

	;ALWAYS SET TO 48000 VLT, CHANGE
	satval = 48000 < saturate
	
	;Hard coded abovesky values for each filter NEEDS TO BE CHANGED FOR RIMAS VLT
	IF strcmp(filter,'H') THEN abovesky = 6000 ELSE IF strcmp(filter,'J') THEN abovesky = 4500 ELSE IF strcmp(filter,'Y') THEN abovesky = 800 ELSE abovesky = 100

	;VLT HARD CODED NUMBERS REMOVE FOR RIMAS?
	
	;
	ct = 0
	while ct lt 7 do begin
   		if abovesky lt 0 then break
   		pstars = where(fwhm gt seeing[0] and fwhm lt seeing[1]*airmass and maxpix lt satval and maxpix gt abovesky+skypix, ct)
   		if ct lt 7 and abovesky gt 1600 then abovesky = abovesky - 1000
   		if ct lt 7 and abovesky le 1600 then abovesky = abovesky - 100
	endwhile
	
	;If no stars found, stop program for user to retry other inputs
	if ct eq 0 then begin
   		print, "Couldn't identify any calibration stars in the image."
   		print, "Check input seeing and try again."
   		stop ; crash deliberately here
	endif
	
	;Find the median of the FWHM of the potential stars.  Choose the smaller of 1.2*med(FWHM)
	;or airmass*seeing upper limit for the FWHM threshold
	medstarfwhm = median(fwhm[pstars], /even)
	fwhmthresh  = medstarfwhm*1.2 < airmass*seeing[1]

	;Set get seeing to the median FWHM of the potential stars
	getseeing = medstarfwhm
	; could do the others too, but no need - watch out for this however

	;If inner radius is not set use 1.2*med(FWHM) otherwise use set radius
	if n_elements(rinner) eq 0 then r1 = medstarfwhm*1.2 else r1 = rinner

	;If relative keyword set then set the large aperture to the value of r1 (either rinner of 1.2*med(fwhm)
	if keyword_set(relative) then begin
  		r2 = r1
	endif else begin
  		r2 = router
	endelse

	;Using r1 as aperture radius
	print, 'Using aperture radius ', clip(r1)

	rcheck1 = medstarfwhm     
	rcheck2 = medstarfwhm*1.5 

	;Trust sextractor's centers and use IDL astro library function, aper, that takes a 
	;dataset with (x,y) coordinates and returns the magnitude, magnitude error, sky value,
	;and sky error.  This command sets the photons per analog digital units to 1 (aper assumes poisson statistics),
	;sets 4 photometry aperture radii <-r1,r2,rcheck1,rcheck2, inner and outer radii of sky annulus (40 and 50), 
	;and minimum and maximum values of a good pixel (-100 to saturation value)
	;Zeropoint magnitude set to 25 initially
	aper, data, x, y, apmagarr, apmagerrarr, sky, skyerr, 1.0, [r1,r2,rcheck1,rcheck2], [40,50], [-100,satval], /silent
	
	;Adjusts the magnitude for exposure time, why multiply???? VLT
	apmagarr = apmagarr + 2.5*alog10(exptime) ; adjust for exposure time!!!
	
	;VLT arbitrary magnitude difference?
	dmfloor = 0.03 ; some danger here if the fwhm-estimated aperture is way too large or small

	;Magnitude difference between two apertures (rcheck1 and rcheck2)
	dm = (apmagarr[2,*] - apmagarr[3,*]) [*]; flux difference between two apertures
	
	;Choose objects with FWHM less than threshold (1.2*med(fwhm) or airmass*maxseeing) 
	;and difference in mags larger than 0.03.  Separate bright stars and create dm threshold
	;using the median of bright stars VLT WHY THESE NUMBERS.  Make sure threshold is at least 0.04
	stars       = where(fwhm lt fwhmthresh and dm gt dmfloor)
	bstars      = where(fwhm lt fwhmthresh and dm gt dmfloor and maxpix ge median(maxpix[stars], /even) )
	dmthresh    = (median(dm[bstars],/even) * [0.8, 1.2]) + [-0.1, 0.1]
	dmthresh[0] = dmthresh[0] > 0.04 ; some danger here if the fwhm-estimated aperture is way too large or small

	;Separate magnitude results for each aperture
	apmag1     = (apmagarr[0,*]) [*]
	apmag2     = (apmagarr[1,*]) [*]
	apmagerr1  = (apmagerrarr[0,*]) [*]
	apmagerr2  = (apmagerrarr[1,*]) [*]

	;Make sure magnitudes are within reasonable bounds and that there is a difference between 
	;magnitudes using med(fwhm) and 1.5*med(fwhm).  If there's not, maybe too faint? VLT
	isgood1mag = apmag1 gt -50 and apmag1 lt 50 and dm gt 0
	isgood2mag = apmag2 gt -50 and apmag2 lt 50 and dm gt 0

	;If relative keyword not set,
	if keyword_set(relative) eq 0 then begin
  		brightthresh = 0.02-0.0025
  		minsep = r2+2*r1
  
  		;Error lower than brightness threshold, and 
  		ct = 0
  		while ct lt 5 do begin ; used to be 3...
     		brightthresh = brightthresh + 0.0025
     		brightstars = where(apmagerr2 lt brightthresh and mindist gt minsep and fwhm lt fwhmthresh and $
                         	dm gt dmthresh[0] and dm lt dmthresh[1] and isgood2mag, ct)
                         	
     		if brightthresh gt 0.09 then begin ; not willing to accept any more than a 10 percent error, at this point start degrading the search radius
        		
        		if minsep gt r2 then begin 
           			minsep = minsep - r1
           			brightthresh = 0.02
           			print, 'WARNING - not enough isolated bright stars ('+clip(ct)+') found; decreasing minsp to ', fpr(minsep,2.2), ' pixels.'
        		endif else begin
           			if ct ge 1 then begin
              			print, 'WARNING - very few isolated bright stars ('+clip(ct)+') found.'
              			break
           			endif else begin
              			print, 'ERROR - found no isolated stars.'
              			return
           			endelse
        		endelse
     		endif
  		endwhile
  
  		apcorr = median(apmag1[brightstars]-apmag2[brightstars], /even)
  		
  		if ct gt 1 then begin
    		apunc  = stdev(apmag1[brightstars]-apmag2[brightstars])
    		apdev  = median(abs(apmag1[brightstars]-apmag2[brightstars] - apcorr),/even)
  		endif else begin
    		apunc = 1.0
    		apdev = 1.0
  		endelse
  			
  		napstars = n_elements(brightstars)
		print, 'Aperture correction:',fpr(apcorr,2.3), ' +/-', fpr(apunc,2.3), '(', fpr(apdev,1.3), '):  ', clip(napstars), ' stars'
	
	endif else begin
  		apcorr = 0.
  		apunc = 0.
  		apdev = 0.
 	 	napstars = 0
	endelse

	apmag1corr = apmag1 - apcorr

	zpteff = zpt1[0] - airterm*(airmass-1)

	printf, 3, fpr(medstarfwhm,2.3), ' ', fpr(r1,2.3), ' ', fpr(dmthresh[0],2.3), ' ', fpr(dmthresh[1],1.3), ' | ', fpr(zpteff,2.3), ' ', fpr(apcorr,2.3), ' ', fpr(apunc,1.3), ' ', fpr(apdev,1.3), ' ', clip(napstars,3)

	calthresh = 0.1
	calstars = where(apmagerr1 lt calthresh and mindist gt 3*r1 and fwhm lt fwhmthresh and dm gt dmthresh[0] and dm lt dmthresh[1] and isgood2mag, ct)
	ra = ra[calstars]
	dec = dec[calstars]
	mag = apmag1corr[calstars] + zpteff
	magerr = apmagerr1[calstars]

end


pro multifieldphot, imdata, stardata, zptdata, seeing=seeing, router=router

	allra = [0.D]
  	alldec = [0.D]
  	allmag = [0.]
  	allmagerr = [0.]
  	allfilt = ['']
  	allchip = ['']
  	alli = [-1]

	
  	for i = 0, n_elements(imdata)-1 do begin
     	print, imdata[i].filename
     	filt = imdata[i].filt
     	
     	;For each image corresponding to the target, write the filename, exposure time, filter, air term
     	printf, 3, clip(imdata[i].filename,23)+' '+clip(fix(imdata[i].exp),3)+' '+clip(filt,2)+' '+fpr(imdata[i].air,1.2)+' | ', format='($,A)'
     
     	;Find if zeropoint data matches filter
     	fi = where(zptdata.filt eq filt, ct)
     
     	;If zeropoint filter doesn't match file filter run fieldphot
     	if ct gt 0 then begin
        	zpt1 = zptdata[fi].zpt1
        	airterm = zptdata[fi].airterm
        	fieldphot, imdata[i].filename, ra, dec, mag, magerr, zpt1=zpt1, airterm=airterm, seeing=seeing, router=router, getseeing=getseeing
     	endif else begin
        	fieldphot, imdata[i].filename, ra, dec, mag, magerr, zpt1=0, airterm=0, seeing=seeing, router=router, /relative, getseeing=getseeing
     	endelse
     
     	if n_elements(mag) eq 0 then begin
        	print, 'Unable to determine field photometry for ', imdata[i].filename
        	continue
     	endif

     	imdata[i].seeing = getseeing

     	good      = where(mag lt 50 and mag gt -50)
     	good      = good[sort(mag[good])]
     	allra     = [allra, ra[good]]   ; append good stars to megalist
     	alldec    = [alldec, dec[good]]
     	allmag    = [allmag, mag[good]]
     	allmagerr = [allmagerr, magerr[good]]
     	allfilt   = [allfilt, replicate(imdata[i].filt,n_elements(good))]
     	allchip   = [allchip, replicate(imdata[i].chip,n_elements(good))]
     	alli      = [alli, replicate(i,n_elements(good))]
  	endfor

  	nstardata = n_elements(allra)-1
  	stardata  = replicate({ra:0.D, dec:0.D, mag:0., magerr:0., filt:'', chip:'', image:0, starid:0}, nstardata)
 
  	stardata.ra     = allra[1:*]
  	stardata.dec    = alldec[1:*]
  	stardata.mag    = allmag[1:*]
  	stardata.magerr = allmagerr[1:*]
  	stardata.filt   = allfilt[1:*]
  	stardata.chip   = allchip[1:*]
  	stardata.image  = alli[1:*]
  	stardata.starid = intarr(nstardata)-1

  	cosdec = cos(median(alldec, /even)*!pi/180.)
  
  	s = 0
  	for a = 0, nstardata-1 do begin
     	
     	if stardata[a].starid eq -1 then begin
        	stardata[a].starid = s
        	s = s + 1
     	endif
     
     	neighbors = findneighbors(stardata.ra*cosdec, stardata.dec, a, mindist=1./3600., count=ct) ;1./12000. = 0.3 arcsec
     	if ct gt 0 then stardata[neighbors].starid = stardata[a].starid
  	endfor
  	
  	ns = s
	return

end

pro makecal, cat, stardata, imdata, calname, object, uncthresh
;imdata should be imdata[thisobj]
	
	common lrisauto, autoastrocommand, swarpcommand, datadir, imworkingdir, spworkingdir, imfinaldir, spfinaldir, lrisversion, wildchar, overwrite
    common lrisfail, flatfail, catastrofail, relastrofail, fullastrofail, extractfail, wavsolfail, wavsolwarn
	
	filters = unique(stardata.filt)
    nfilt = n_elements(filters)
    maxs = max(stardata.starid)
    nexp = n_elements(unique(stardata.image))
    cosdec = cos(median(stardata.dec, /even)*!pi/180.)

    ; DIRECT cal - just get stars straight from the catalog that are in all images / not saturated

    caloutfile = imworkingdir + strlowcase(object) + '.'+strlowcase(calname)+'.direct.cal'

    openw, 1, caloutfile
    tstr = '#'
    tstr = tstr +      clip('RA',10-1)
    tstr = tstr + ' '+ ' ' + clip('dec',10-1)
    
    for f = 0, n_elements(filters)-1 do begin
    	tstr = tstr + ' '+ clip(filters[f],6)
       	tstr = tstr + ' '+ 'unc  '
    endfor
    
    printf, 1, tstr

    nminf = intarr(n_elements(filters)) ; star must be catalogued in this many images to be included
    
    for f = 0, n_elements(filters)-1 do begin
    	maxoverlap = 0
       	
       	for s = 0, maxs do begin
        	w = where(stardata.starid eq s and stardata.filt eq filters[f], n)
          	if n gt maxoverlap then maxoverlap = n  ; new code block
       	endfor
       
       	nminf[f] = maxoverlap*0.5
    endfor

    nprinted = 0
    for s = 0, maxs do begin
       	w = where(stardata.starid eq s, n)
       	if n lt maxoverlap or n lt 1 then continue   ; n lt nexp-nfilt-1 
       
       	if n eq 1 then begin
         	avgra = stardata[w[0]].ra
         	avgdec = stardata[w[0]].dec
       	endif else begin
         	avgra = median(stardata[w].ra, /even)
         	avgdec = median(stardata[w].dec, /even)
      	endelse

       	smatch = nearest(avgra*cosdec, avgdec, cat.ra*cosdec, cat.dec, mindist=1./3600., count=ct)
       	print,'SMATCH = ',smatch,s

       	if smatch eq -1 then continue

       	sstr = ''
       	sstr = sstr +       fpr(cat[smatch].ra,3.6)
       	sstr = sstr + ' ' + fpr(cat[smatch].dec,3.6)
       	
       	usestar = 0
       	for f = 0, n_elements(filters)-1 do begin
          	unc = 99.
          	filt = filters[f]
          	mf = (where(filt eq cat[smatch].filts, ct)) [0]
          	if ct gt 0 then begin
            	mag = cat[smatch].mags[mf]
            	unc = cat[smatch].maguncs[mf]
          	endif else begin
            	mag = 999
            	unc = 999
          	endelse

          	nexpfilt = n_elements(unique(stardata[where(stardata.filt eq filt)].image)) ; could use imdata again
          	dum = where(stardata.starid eq s and stardata.filt eq filters[f], ngoodfilt)
          
			print,s,f,mag,unc,filters[f],ngoodfilt

                                               ; formerly nexpfilt-1.  need to come up with a better solution for RL chips.
          	if unc ge -1. and unc lt uncthresh and ngoodfilt ge nminf[f] and mag ge 0. and mag lt 99. then begin  ; print the mag
              	usestar = 1
              	sstr = sstr + ' ' + fpr(mag,2.3) + ' ' + fpr(unc,1.3)
              	nprinted = nprinted + 1
          	endif else begin
              	sstr = sstr + ' ' + '-     '      + ' ' + '-    '
           	endelse
		endfor
		
       	if usestar then printf, 1, sstr
	endfor
	
    close, 1
    if nprinted gt 0 then begin
    	print, 'Catalog printed to ', caloutfile, ' (',clip(nprinted),' stars)'
    endif else begin
       	print, 'WARNING - no stars overlap all fields!  Cannot make direct catalog.'
       ;return  - no, keep going, although this might be perilous
    endelse

    ; EXTENDED cal using the shortest exposure
    ; need to make this look at both left and right, as usual.  
    ; loop (right, left, concatenate)?

	if strlowcase(object) eq 'grb051008' or object eq '051008' then return ; skip this due to bright star
    if (strlowcase(object) eq 'grb070810b' or strlowcase(object) eq '070810b') and strpos(calname,'sdss') ge 0 then return ; sloan is apparently near the field but not near enough
    if (strlowcase(object) eq 'grb080210' or object eq '080210') and strpos(calname,'sdss') ge 0 then return ; sloan is apparently near the field but
    if (strlowcase(object) eq 'grb080603a' or strlowcase(object) eq '080603a') and strpos(calname,'sdss') ge 0 then return ; sloan is apparently near the field but

    noshort = 0
    fdzpt = fltarr(n_elements(filters))
    fdzptstdev = fltarr(n_elements(filters))
    
    for f = 0, n_elements(filters)-1 do begin
    	filt = filters[f]

       	checkcatfilt = where(filt eq cat[0].filts, ctf)
       	
       	if ctf eq 0 then begin
          	noshort = noshort+1
          	print, 'Filter ', filt, ' not in catalog.'
          	continue
       	endif

       	ff = where(imdata.filt eq filt)

       	shortestexp = ff[(where(imdata[ff].exp eq min(imdata[ff].exp))) [0]]  ; the i index of the shortest file, hopefully
       	if min(imdata[ff].exp) gt 1000 then begin
           	noshort = noshort + 1
           	continue
       	endif

       	print, 'Calibrating from ', imdata[shortestexp].filename

       	dzpt = [-1]
       	for s = 0, maxs do begin
        	w = where(stardata.starid eq s and stardata.image eq shortestexp, n)
           	if n eq 0 then continue
           	if n gt 1 then continue ; n > 1 means two matching sources - sketchy (shouldn't happen?)
           	avgra = stardata[w].ra
           	avgdec = stardata[w].dec
           	smatch = nearest(avgra*cosdec, avgdec, cat.ra*cosdec, cat.dec, mindist=1./3600, count=ct)

           	if ct eq 0 then continue

           	mag = -1
           	mf = where(filt eq cat[smatch].filts, ctf)
           	if ctf gt 0 then begin
             	mag = cat[smatch].mags[mf]
             	unc = cat[smatch].maguncs[mf]
           	endif

           	if ctf gt 0 and mag gt 0. and unc ge 0. and unc lt 99. then dzpt = [dzpt, mag-stardata[w].mag]
       	endfor
       	
       	if n_elements(dzpt) eq 1 then continue ; some sort of failure has occurred
       	dzpt = dzpt[1:*]
       	ntot = n_elements(dzpt)
       	meddzpt = median(dzpt, /even)
       	if ntot gt 1 then stddzpt = stdev(dzpt) else stddzpt = 2
       	clipdzpt = dzpt[where(abs(dzpt-meddzpt) lt 3*stddzpt, ngood)]
       	nout = ntot - ngood
       	fdzpt[f] = median(clipdzpt, /even)
       	if ngood gt 1 then fdzptstdev[f] = stdev(clipdzpt) else  fdzptstdev[f] = 2
       	absdev = median(abs(clipdzpt - fdzpt[f]), /even)

       	print, filt, ' catalog zeropoint adjustment:', fpr(fdzpt[f],3.3), ' +/- ', fpr(fdzptstdev[f],1.3), $
       	             ' (', clip(ntot), ' stars, ', clip(nout), ' outliers)'
       	printf, 3, '# dzpt_'+strmid(calname,5)+'_'+filt+' =', fpr(fdzpt[f],2.3), ' ',fpr(fdzptstdev[f],1.3), ' ', fpr(absdev,1.3), ' ', clip(ngood,3), clip(nout,3)

	endfor

	if noshort lt n_elements(filters) then begin
    	caloutfile = imworkingdir + strlowcase(object) + '.'+strlowcase(calname)+'.extend.cal'
      	openw, 1, caloutfile
      	tstr = '#'
      	tstr = tstr +      clip('RA',10-1)
      	tstr = tstr + ' '+ ' ' + clip('dec',10-1)
      
      	for f = 0, n_elements(filters)-1 do begin
         	tstr = tstr + ' '+ clip(filters[f],6)
         	tstr = tstr + ' '+ 'unc  '
      	endfor
      
      	printf, 1, tstr

      	for s = 0, maxs do begin
         	usestar = 0
         	w = where(stardata.starid eq s, n)
         	if n eq 0 then continue ; this shouldn't happen but it does?
         	avgra = median([stardata[w].ra], /even)
         	avgdec = median([stardata[w].dec], /even)
         	sstr = ''
         	sstr = sstr +       fpr(avgra,3.6)
         	sstr = sstr + ' ' + fpr(avgdec,3.6)
         	
         	for f = 0, n_elements(filters)-1 do begin
            
        		shortestexp = ff[(where(imdata[ff].exp eq min(imdata[ff].exp))) [0]]
           		w = where(stardata.starid eq s and stardata.image eq shortestexp, n)

            	nexpfilt = n_elements(unique(stardata[where(stardata.filt eq filt)].image)) ; could use imdata again
            	dum = where(stardata.starid eq s and stardata.filt eq filters[f], ngoodfilt)

            	if n eq 1 then begin
              		mag = stardata[w].mag + fdzpt[f]
              		unc = stardata[w].magerr
              	
              		if ngoodfilt ge nminf[f] and unc lt 0.05 then begin   ; again, formerly nexpfilt-1
                 		sstr = sstr + ' ' + fpr(mag,2.3) + ' ' + fpr(unc,1.3)
                 		usestar = usestar + 1
              		endif else begin
                 		sstr = sstr + ' ' + '-    ' + ' ' + '-    '
              		endelse
        		endif

    		endfor
        
        	if usestar then printf, 1, sstr
      	endfor
      	
      	close, 1
	endif
    print, 'Catalog printed to ', caloutfile

	return
end

;Seeing need to be set???

;	seeingarcsec - seeing in arcsecs, if not set looks for seeing.txt, if that is not found set to default [0.4, 3.0]
;	objcat - catalog of objects looking for, if not set, use target name

pro pipeautofieldphot, red=red, blue=blue, seeingarcsec=seeingarcsec, forcephotometric=forcephotometric, objcat=objcat, chip=chip, outpipevar=outpipevar, inpipevar=inpipevar

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

	;If imworkingdir is not set, use current working directory
	if n_elements(pipevar.imworkingdir) eq 0 then pipevar.imworkingdir = ''

	;VLT Change for RIMAS
	prefchar = '2'
	
	;Finds images that have been prepped, flatfielded, and astrometry corrected
	imagefiles = choosefiles(pipevar.imworkingdir+'a*f*p'+prefchar+'*_img_?.fits')

	;End program if there are no correctly processed images
	if n_elements(unique(imagefiles)) eq 1 and imagefiles[0] eq '' then return
	imagefiles = imagefiles[where(imagefiles ne '')]
	
	close, /all

	;If seeingarcsec keyword not set either use seeing.txt file for seeing information 
	;or use default of 0.4 and 3.0 if text file does not exist.  Uses specified seeing if keyword set
	if keyword_set(seeingarcsec) eq 0 then begin
  		if file_test('seeing.txt') then begin
    		rdfloat, 'seeing.txt', see1, see2
    		seeingarcsec = [see1[0], see2[0]]
  		endif else begin
    		seeingarcsec = [0.4, 3.0]
  		endelse
	endif
	
	seeing  = seeingarcsec/0.135
	zptfile = pipevar.imworkingdir+'zeropoints.dat'
	
	;If zeropoint file is not in imworkingdir then use the file in the defaultspath
	if file_test(zptfile) eq 0 then begin
   		zptfile = pipevar.defaultspath+'/'+'zeropoints.dat'
   		photometric = 0
	endif else begin
	   	photometric = 1
	endelse
	
	;Open the zeropoint file
	openr, 1, zptfile
	utdatestr = ''
	aperstr   = ''
	readf, 1, utdatestr ;date zeropoint.dat uses
	readf, 1, aperstr	;aperture zeropoints use
	close, 1
	
	;Full aperture radius
	fullaperrad = float((strsplit(aperstr,/extract)) [0])

	;Read in zeropoint file
	readcol, zptfile, zfilt, zzpt1, zzpt1unc, zairterm, zairtermunc, format='a,a,f,f,f', /silent

	;Create structure with zeropoint data (filter, zeropoint, airterm) for as many filters in file
	;Save filter, zeropoint, and airterm into structure
	zptdata = replicate({filt:'', zpt1:0., airterm:0.},n_elements(zfilt))
	zptdata.filt    = zfilt
	zptdata.zpt1    = zzpt1
	zptdata.airterm = zairterm

	;For each image file read the header, save the filename, filter, exposure time, and airmass into structure
	;If object catalog is set, read the files and find the 
	nim = n_elements(imagefiles)
	imdata = replicate({filename:'', isgrb:0, object:'', filt:'', exp:0., air:0., seeing:0., fluxratio:0.}, nim)
	for i = 0, nim-1 do begin
  		h   = headfits(imagefiles[i])

  		imdata[i].filename = imagefiles[i]
  		imdata[i].filt     = clip(sxpar(h,'FILTER'))
  		imdata[i].exp      = clip(sxpar(h,'ELAPTIME'))
  		imdata[i].air      = sxpar(h,'AIRMASS')

		;If objcat keyword set then read in file with target name, ra, and dec, if it isn't, 
		;assume object exists and use target name
  		if n_elements(objcat) gt 0 then begin
  		   	grbsposfile = objcat
   			readcol, grbsposfile, grbname, grbra, grbdec, format='a,f,f', /silent
   			ngrb = n_elements(grbra)
    		distance = fltarr(ngrb)
    		
    		ra  = sxpar(h,'CRVAL1')  ; use the WCS header keywords to get the real center
  			dec = sxpar(h,'CRVAL2')
    
    		;For each object, find the cartesian distance in arcseconds, if this is small
    		;then calculate the great circle distance.  Save the distance
    		for r = 0, ngrb-1 do begin
       			dist = 3600*(abs(dec-grbdec[r]) + abs(ra-grbra[r]))
       			if dist lt 300 then gcirc, 2, ra, dec, grbra[r], grbdec[r], dist
       			distance[r] = dist
    		endfor
    		
    		;Find objects that are the closest to the pointing center
    		mindist = min(distance)
    		minr = (where(distance eq mindist)) [0]
    
    		;If the minimum distance to the object is less than 300 arcseconds then set
    		;set object as GRB and save object information
    		if mindist lt 300. then begin
      			imdata[i].isgrb = 1
      			print, clip(sxpar(h,'OBJECT')), ' -> ',  grbname[minr]
      			imdata[i].object = grbname[minr]     ;clip(sxpar(h,'OBJECT'))
    		endif
    		
  		endif else begin
    		imdata[i].isgrb = 1
    		imdata[i].object = clip(sxpar(h,'TARGNAME'))
  		endelse

	endfor

	;Only keep files that have our desired object in them, if objcat not set then this is all of the files
	;Objects is each target (RATIR's case will typically only be one target)
	imdata = imdata[where(imdata.isgrb)]
	objects = unique(imdata.object)

	stop

	;Create photometry summary files
	openw, 3, imworkingdir+'autophotsummary.txt'
	openw, 4, imworkingdir+'autophotsummaryflux.txt'

	;For each target, find the corresponding image files.
	for b = 0, n_elements(objects)-1 do begin
  		object = objects[b]
  		if object eq '' then continue

  		imfiles = where(imdata.object eq object, nobjim)  
  		objimdata = imdata[imfiles]                           ; only look at images of the relevant field

  		print
  		print, object
  		printf, 3, '# ', object
  		printf, 4, '# ', object

  		multifieldphot, objimdata, stardata, zptdata, seeing=seeing, router=fullaperrad
     	; stardata is the only *output*: list of ra, dec, mag, magerr, filt:'', image (#), starid (#)

  		; Compare exposures to tell if this was a photometric observation

  		filters = unique(stardata.filt)
  		for f = 0, n_elements(filters)-1 do begin
     		print, filters[f], ' transmission'
     		thisfilt = where(stardata.filt eq filters[f], ct)

     		if ct eq 0 then continue
     		fstardata = stardata[thisfilt]
     
     		fimagelist = unique(fstardata.image)   ; fimagelist is f subscripts, fstardata.image is o subscripts
     		fimagedata = objimdata[fimagelist]     ; fimagedata is f subscripts
     		fstarlist = unique(fstardata.starid)

     		nfi = n_elements(fimagelist) ; number of filter/chip images
     		dmagarr = fltarr(nfi) - 999
     		dmagij = fltarr(nfi,nfi)
     		ncompij = intarr(nfi,nfi)

     		; probably should do something here when nstars gets extremely large (>1000) to boost speed
     		if nfi gt 1000 then begin
      			print, 'Crowded field, may run slowly...'
     		endif

     		for fi = 0, nfi-1 do begin ; primary image loop          ; fi and fj indices are for fimagelist
       			for fj = 0, nfi-1 do begin ; compare to this image     
         			if fi eq fj then continue
         			oi = fimagelist[fi]                           ; oi and oi indices are for the object/starlist (all filters, chips)
         			oj = fimagelist[fj]
         			dmagarr = [-1]
         
         			for fs = 0, n_elements(fstarlist)-1 do begin ; collect the stars
            			os = fstarlist[fs]
            			wi = where(fstardata.starid eq os and fstardata.image eq oi, cti)
            			wj = where(fstardata.starid eq os and fstardata.image eq oj, ctj)
            
            			if cti eq 1 and ctj eq 1 then begin
               				sdm = (fstardata[wi].mag - fstardata[wj].mag) [0]
               				dmagarr = [dmagarr, sdm]
            			endif
         			endfor
         		
         			ncompstars = n_elements(dmagarr)-1
         			ncompij[fi,fj] = ncompstars
         			if ncompstars gt 0 then dmagij[fi,fj] = median(dmagarr[1:*], /even) else dmagij[fi,fj] = !values.f_nan

       			endfor
			endfor

     		dmag = fltarr(nfi)
     		for fi = 0, nfi-1 do begin 
       			for fj = fi+1, nfi-1 do begin 
          			if dmagij[fi,fj] lt -10 or finite(dmagij[fi,fj]) eq 0 then continue
          			newdmagij = dmagij[fi,fj] - dmag[fi] + dmag[fj]
          			dmag[fj] = dmag[fj] - newdmagij
       			endfor
    		endfor
     
     		fluxr = 10.0^(-dmag/2.5)
     		for fi = 0, nfi-1 do begin
       			objimdata[fimagelist[fi]].fluxratio = fluxr[fi] ; hopefully assigns correctly
       			print, fimagedata[fi].filename, ' :  dmag = ', fpr(dmag[fi],2.3), '  relflux = ', fpr(fluxr[fi],2.3), '  seeing = ', fpr(fimagedata[fi].seeing,2.3)
       			printf, 4, clip(fimagedata[fi].filename,23)+' '+clip(fix(fimagedata[fi].exp),3)+' '+clip(filters[f],2)+' '+clip(fimagedata[fi].dich,3)+' '+fpr(fimagedata[fi].air,1.2)+' | ', fpr(dmag[fi],2.3), ' ', fpr(fluxr[fi],2.3), ' ', fpr(fimagedata[fi].seeing,2.3)
     		endfor

     		maxabsdmag = max(abs(dmag))
			print,stdev(dmag)

     		if n_elements(dmag) gt 1 then stdevdmag = stdev(dmag) else stdevdmag = 0
     		if stdevdmag gt 0.05 then begin
       			print, 'Not photometric!'  
       			photometricobs = 0
     		endif else begin
       			photometricobs = 1
     		endelse
		endfor

  		; Look up 2MASS for this field and find matches
  		cra = median(stardata.ra, /even)
  		cdec =  median(stardata.dec, /even)
  		craref=cra
  		cdecref=cdec
  		tmpsc = gettmpsc(cra, cdec, 360.)

  		mindist = 3600. * min(sqrt(((stardata.ra-craref)*cos(cdecref*!pi/180.))^2 + (stardata.dec-cdecref)^2))
  		uncthresh = 0.12
  		if n_elements(tmpsc) gt 10 and mindist lt 120 then begin
    		tmpsccat = replicate({ra:0.D, dec:0.D, filts:strarr(3), mags:fltarr(3), maguncs:fltarr(3)}, n_elements(tmpsc))
    		tmpsccat.ra = tmpsc.ra
    		tmpsccat.dec = tmpsc.dec
    		tmpsccat.filts = ['J','H','K'] ; kind of stupidly inefficient, but whatever
    		tmpsccat.mags[0] = tmpsc.j
    		tmpsccat.mags[1] = tmpsc.h
    		tmpsccat.mags[2] = tmpsc.k
    		tmpsccat.maguncs[0] = 0
    		tmpsccat.maguncs[1] = 0
    		tmpsccat.maguncs[2] = 0
			makecal, tmpsccat, stardata, objimdata, 'tmpsc', strlowcase(object), uncthresh
		endif else begin
    		print, 'Not a 2MASS field.'
  		endelse

  		; Look up the Nickel-based catalog and use that as a calibrator
  		; in the future, point to a general catalog and look at everything
  		gipos = strpos(object,'0') < strpos(object,'1')
  		grb = strlowcase(strmid(object,gipos))
  		grb = (strsplit(grb,'_',/extract)) [0]
	
  		grbcalfile = '~/research/redux/nickelcal/grb' + grb + '.cat' 
  		if strlen(grb) eq 6 and file_test(grbcalfile) eq 0 then $
    		grbcalfile = '~/research/redux/nickelcal/grb' + grb + 'a.cat' 

  		if file_test(grbcalfile) then begin
     		ns = countlines(grbcalfile)-1
     		nickelcat = replicate({ra:0.D, dec:0.D, filts:strarr(5), mags:(fltarr(5)-1), maguncs:(fltarr(5)-1)}, ns)

     		inline = ''
     		openr, 1, grbcalfile
     		readf, 1, inline
     		headings = strsplit(strmid(inline,1),/extract)
     		nh = n_elements(headings)
     
     		f = 0
     		for i = 5, nh-1, 4 do begin
        		nickelcat[*].filts[f] = headings[i]
        		f = f + 1
     		endfor
     	
     		nfilt = f
     		for s = 0, ns-1 do begin
        		readf, 1, inline
        		inarr = strsplit(inline,/extract)
        		nickelcat[s].ra = double(inarr[0])
        		nickelcat[s].dec = double(inarr[1])
        
        		for f = 0, nfilt-1 do begin
           			if clip(inarr[5+f*4]) ne '-' then nickelcat[s].mags[f] = float(inarr[5+f*4])
           			if clip(inarr[5+f*4+2]) ne '-' and inarr[5+f*4+1] ne '-' then nickelcat[s].maguncs[f] = float(inarr[5+f*4+1])/3. > float(inarr[5+f*4+2]) 
        		endfor
     		endfor
     	
     		close, 1

     		makecal, nickelcat, stardata, objimdata, camera + '.nickel', strlowcase(object), uncthresh

  		endif else begin
     		print, 'Not a Nickel field.'
  		endelse

  		; Do an absolute calibration
  		if (keyword_set(photometric) and photometricobs) or keyword_set(forcephotometric) then begin
			landoltoutfile = imworkingdir + strlowcase(object)+'.'+camera + '.landolt.cal'
    		openw, 1, landoltoutfile

    		ns = max(stardata.starid)+1
    		filters = unique(stardata.filt)
    		cosdec = cos(median(stardata.dec, /even)*!pi/180.)
    		tstr = '#'
    		tstr = tstr +      clip('RA',10-1)
    		tstr = tstr + ' '+ ' ' + clip('dec',10-1)
    		tstr = tstr + ' '+ clip('unc_RA',8)
    		tstr = tstr + ' '+ clip('unc_dec',8)
    		tstr = tstr + ' '+ clip('n',2)
    
    		for f = 0, n_elements(filters)-1 do begin
       			tstr = tstr + ' '+ clip(filters[f],6)
       			tstr = tstr + ' '+ clip('var',5)
       			tstr = tstr + ' '+ clip('unc',5)
       			tstr = tstr + ' '+ clip('n',2)
    		endfor
    	
    		printf, 1, tstr
    
    		for s = 0, ns-1 do begin
       			w = where(stardata.starid eq s, n)
       			if n eq 0 then continue ; should almost never happen, but occasionally a flpt error...?
       
       			if n gt 1 then begin
          			avgra = median(stardata[w].ra, /even)
          			avgdec = median(stardata[w].dec, /even)
          			uncra = stdev(1.0D*stardata[w].ra)*3600.*cosdec
          			uncdec = stdev(1.0D*stardata[w].dec)*3600.
       			endif else begin
          			avgra = stardata[w[0]].ra
          			avgdec = stardata[w[0]].dec
          			uncra = 1
          			uncdec = 1
       			endelse
   
       			sstr = ''
       			sstr = sstr +       fpr(avgra,3.6)
       			sstr = sstr + ' ' + fpr(avgdec,3.6)
       			sstr = sstr + ' ' + fpr(uncra,1.6)
       			sstr = sstr + ' ' + fpr(uncdec,1.6)
       			sstr = sstr + ' ' + clip(n,2)
       		
       			for f = 0, n_elements(filters)-1 do begin
          			thisfilt = where(stardata[w].filt eq filters[f], ct)
          
          			if ct eq 0 then begin
             			sstr = sstr + ' ' + clip('-',6)
             			sstr = sstr + ' ' + clip('-',5)
             			sstr = sstr + ' ' + clip('-',5)
             			sstr = sstr + ' ' + clip(0,2)
          			endif else begin
          		
             			if ct gt 1 then begin
                			avgmag = median(stardata[w[thisfilt]].mag, /even)
                			avgmagerr = median(stardata[w[thisfilt]].magerr, /even)
                			varmag = stdev(stardata[w[thisfilt]].mag)
             			endif else begin
                			avgmag = stardata[w[thisfilt[0]]].mag
                			avgmagerr = stardata[w[thisfilt[0]]].magerr
                			varmag = 1
             			endelse
             		
             			sstr = sstr + ' ' + fpr(avgmag,2.3)
             			sstr = sstr + ' ' + fpr(varmag,1.3)
             			sstr = sstr + ' ' + fpr(avgmagerr,1.3)
             			sstr = sstr + ' ' + clip(ct,2)
             		
          			endelse
       			endfor
       		 
       			printf, 1, sstr
    		endfor
    
    		close, 1
    		print, 'Catalog printed to ', landoltoutfile
  		endif else begin
    		print, 'Not photometric.'
  		endelse
  
  		print

	endfor
	close, 3
	close, 4

end



