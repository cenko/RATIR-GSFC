function calc_zpt, trumag, obsmag, wts, sigma=sigma
	
	if (not KEYWORD_SET(sigma)) then sigma = 3.0
	
	diff   = trumag - obsmag
	
	sizes  = size(obsmag)
	nstars = sizes[1]
	
	if sizes[0] eq 1 then begin
		nobs = 1
	endif else begin
		nobs   = sizes[2]
	endelse
	
	z = []
	modmag = obsmag
	for i = 0, nobs-1 do begin
		indz = total(diff[*,i]*wts[*,i])/total(wts[*,i])
		z = [z, indz]
		modmag[*,i] = obsmag[*,i]+indz
	endfor
			
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
	
	for p = 0, nobs-1 do begin
		new2 = where(wts[*,p] gt 0, n2ct)
		if n2ct eq 0 then continue
		newd2 = adiff2[new2,p]
		scats[p] = 1.48 * median(abs(newd2-median(newd2)))
		rmss[p] = stddev(newd2)
	endfor
	
	return, [ [z2], [scats], [rmss] ]
end

function sMAD, vec
	return, 1.48 * median(abs(vec-median(vec)))
end
	