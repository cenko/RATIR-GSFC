function nearest, x, y, xarr, yarr, maxdist
	
	dist = sqrt( (x-xarr)^2 + (y-yarr)^2)
	val = min(dist, ind)
	if val lt maxdist then begin
		return, ind
	endif else begin
		return, -1
	endelse

end