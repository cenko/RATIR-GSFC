function sigclipabsdev, datain, sigmahi=sigmahi, sigmalo=sigmalo
   data = datain ; make a copy!
   maxiter = 6
   endclipfrac = 0.01
 
   ndata = n_elements(data)
   if ndata eq 1 then return, data
   if ndata le 3 then return, median(data, /even)
   if n_elements(sigmahi) eq 0 then sigmahi = 3.5
   if n_elements(sigmalo) eq 0 then sigmalo = 3.5

   data = data[where(finite(data))]
   clipdata = data[sort(data)] ; the going theory is that this speeds things up

   nremoved = ndata
   niter = 0
   while nremoved gt (ndata*0.01 < 50) and niter lt maxiter do begin
      runningmedian = median(clipdata)
      runningabsdev = mean(abs(clipdata-runningmedian))
      lolimit = runningmedian - sigmalo * runningabsdev
      hilimit = runningmedian + sigmahi * runningabsdev
      runningn = n_elements(clipdata)
      pass = where(clipdata ge lolimit and clipdata le hilimit, ct)
      if ct gt 0 then clipdata = clipdata[pass]
      nremoved = runningn - ct
      niter = niter + 1
   endwhile

   return, mean(abs(clipdata-runningmedian))
end

