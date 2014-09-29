;+
; NAME:
;	calc_zpt
;
; PURPOSE:
;	Finds sextractor objects with optional input parameters and default values.  Saves
;	object mask, starfile, and estimates the seeing from the stars found. 
;
; INPUTS:
;	trumag  - list of catalog magnitudes (array of stars by observations - [ [ob1st1, ob1st2, ob1st3], [ob2st1, ob2st2, ob2st3] ])
;	obsmag  - list of observed magnitudes (array of stars by observations - [ [ob1st1, ob1st2, ob1st3], [ob2st1, ob2st2, ob2st3] ])
;	wts     - list of weights for each object (array of stars by observations - [ [ob1st1, ob1st2, ob1st3], [ob2st1, ob2st2, ob2st3] ])
;
; OPTIONAL KEYWORDS:
;	sigma   - sigma value for how far values can be from robust scatter value
;   plotter - filename to save zeropoint plot
;
; OUTPUTS:
;	z2     - zeropoint values after 2 zeropoint corrections
;	scats  - robust scatter after 2 zeropoint corrections
;	rmss   - rms values after 2 zeropoint corrections
;
; EXAMPLE:
;	mix = calc_zpt(reflist, maglist, errlist, sigma=3.0, plotter='test.ps')
;
; DEPENDENCIES:
;	None
;
; Written by Brad Cenko 
; Abridged version written by Vicki Toy 7/7/2014
;
; FUTURE IMPROVEMENTS:
;	Less redundant code?
;-

function calc_zpt, trumag, obsmag, wts, sigma=sigma, plotter=plotter
	
	if (not KEYWORD_SET(sigma)) then sigma = 3.0
	
	;Find difference between catalog and observed magnitudes
	diff   = trumag - obsmag

	;Find number of observations and stars	
	sizes  = size(obsmag)
	nstars = sizes[1]
	
	if sizes[0] eq 1 then begin
		nobs = 1
	endif else begin
		nobs   = sizes[2]
	endelse
	
	;For each observation (i.e. frame) find the weighted difference and store zeropoint
	;and new magnitude with zeropoint correction
	z = []
	modmag = obsmag
	for i = 0, nobs-1 do begin
		indz = total(diff[*,i]*wts[*,i])/total(wts[*,i])
		z = [z, indz]
		modmag[*,i] = obsmag[*,i]+indz
	endfor
	
	;Find difference of catalog and zeropoint corrected values. Remove any values with weights set
	;to 0 or lower.  Calculate robust scatter on these values.  If difference with these weights
	;is not within sigma*robust scatter then set weight to 0
	adiff1 = trumag - modmag
	for i = 0, nobs-1 do begin
		new1 = where(wts[*,i] gt 0, n1ct)
		if n1ct eq 0 then continue
		newd1 = adiff1[new1,i]
		scat1 = 1.48 * median(abs(newd1-median(newd1)))
		for j = 0, nstars-1 do begin
			if abs(adiff1[j,i] - median(newd1)) gt (sigma * scat1) then begin
				wts[j,i] = 0
			endif
		endfor
	endfor
	
	;Recalculate zeropoint and corrected magnitudes
	z2 = []
	modmag2 = obsmag
	for i = 0, nobs-1 do begin
		indz = total(diff[*,i]*wts[*,i])/total(wts[*,i])
		z2 = [z2, indz]
		modmag2[*,i] = obsmag[*,i]+indz
	endfor
		
	adiff2 = trumag - modmag2
	scats  = fltarr(nobs)
	rmss   = fltarr(nobs)
	
	;Recalculate robust scatter as well as rms scatter value on twice zeropoint corrected mags
	for p = 0, nobs-1 do begin
		new2 = where(wts[*,p] gt 0, n2ct)
		if n2ct eq 0 then continue
		newd2 = adiff2[new2,p]
		scats[p] = 1.48 * median(abs(newd2-median(newd2)))
		rmss[p] = stddev(newd2)
	endfor
	

	if keyword_set(plotter) then begin
	    set_plot, 'ps'
	    device, filename=plotter
		keep = where(wts ne 0)
	    plot, trumag[keep], adiff2[keep], psym=2
	    errplot, trumag[keep],adiff2[keep]-1.0/sqrt(wts[keep]), adiff2[keep]+1.0/sqrt(wts[keep])
	    device, /close
	    set_plot, 'x'
	endif
	
	return, [ [z2], [scats], [rmss] ]
end
	