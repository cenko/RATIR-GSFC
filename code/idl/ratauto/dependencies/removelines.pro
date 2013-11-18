pro removelines, wav, flux, newflux, lines=lines, vel=vel, interpolate=interpolate, nan=nan, linezone=linezone

 ; remove lines in a flux star.

 good = intarr(n_elements(wav))+1

 for l = 0, n_elements(lines)-1 do begin
   cwav = lines[l]
   if n_elements(vel) gt 1 then v = vel[l] else v = vel
   linez = where(wav gt cwav*(1-v) and wav lt cwav*(1+v), ct)
   if ct gt 0 then good[linez] = 0
 endfor

 linezone = abs(1 - good) ; not a where, a bool arr

 if keyword_set(interpolate) eq 0 and keyword_set(nan) eq 0 then begin
    wav = wav[where(good)]
    flux = flux[where(good)]
 endif else begin
    if keyword_set(nan) then begin
       flux[where(good eq 0)] = !values.f_nan
    endif
    if keyword_set(interpolate) then begin
       newflux = interpol2(flux[where(good)], wav[where(good)], wav)
    endif
 endelse

end

