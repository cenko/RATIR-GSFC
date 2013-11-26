function sigclipmedian, datain, sigmahi=sigmahi, sigmalo=sigmalo
   data = datain ; make a copy!
   if n_elements(maxiter) eq 0 then maxiter = 6
   endclipfrac = 0.01
 
   ndata = n_elements(data)
   if ndata eq 1 then return, data[0]
   if ndata le 3 then return, median(data, /even)
   if n_elements(sigmahi) eq 0 then sigmahi = 3.5
   if n_elements(sigmalo) eq 0 then sigmalo = 3.5

   good = where(finite(data), ctgood)
   if ctgood eq 0 then return, !values.f_nan
   if ctgood eq 1 then return, data[good[0]]
   data = data[good]
   clipdata = data[sort(data)] ; the going theory is that this speeds things up

   nremoved = ndata
   niter = 0
   while niter lt maxiter and nremoved gt 0 do begin
      runningmedian = median(clipdata, /even)
      runningstdev = stdev(clipdata)
      lolimit = runningmedian - sigmalo * runningstdev
      hilimit = runningmedian + sigmahi * runningstdev
      runningn = n_elements(clipdata)
      pass = where(clipdata ge lolimit and clipdata le hilimit, ct)
      if ct gt 0 then clipdata = clipdata[pass] else break
      nremoved = runningn - ct
      niter = niter + 1
   endwhile

   return, median(clipdata, /even)
end

