function interpol2, y, x, xnew, _extra

  ; interpolation with a more reasonable treatment of values outside the range for noisy data.

   v = interpol(y, x, xnew, _extra=_extra)

   wbelow = where(xnew lt min(x), ctbelow)
   wabove = where(xnew gt max(x), ctabove)

   ; really crude right now, assumes positive order, etc.

   n = n_elements(y)
   if ctbelow gt 0 then v[wbelow] = median(y[0:2])
   if ctabove gt 0 then v[wabove] = median(y[n-2:n-1])
  
   return, v

end

