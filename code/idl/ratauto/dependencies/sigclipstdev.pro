function sigclipstdev, datain, sigmahi=sigmahi, sigmalo=sigmalo, maxiter=maxiter
   data = datain ; make a copy!
   if n_elements(maxiter) eq 0 then maxiter = 6
   endclipfrac = 0.01
 
   ndata = n_elements(data)
   if ndata eq 1 then return, data
   if ndata le 3 then return, median(data, /even)
   if n_elements(sigmahi) eq 0 then sigmahi = 3.6
   if n_elements(sigmalo) eq 0 then sigmalo = 3.6

   data = data[where(finite(data))]
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

   return, stdev(clipdata)
end

